################################################################################
### DESCRIPTION
### This file defines the Log-Logistic-Marginal function. It contains the 
### expressions of the pdf, cdf, its derivatives wrt the two parameters, names,
### response functions, offset functions, check_y functions required to create an mboost Family object.
### This function is handed over to a copula function and the expressions and functions
### are merged with the respective copula via the rlang package, in order to create 
### appropriate Families objects from gamboostLSS.


### libraries 
library(mboost)
library(gamboostLSS)
library(rlang)


### Inverse Gaussian 
InverseGaussian_Mar <- function(loc = NULL, offset_mu = NULL, offset_sigma = NULL){
  
  
  # check for appropriate offset values for parameter mu, sigma and tau
  
  if ((!is.null(offset_mu) && offset_mu <= 0))
    stop(sQuote("mu"), paste(" must be greater than zero in Marginal", loc))
  
  if ((!is.null(offset_sigma) && offset_sigma <= 0))                                          
    stop(sQuote("sigma"), paste(" must be greater than zero in Marginal", loc))
  
  
  
  
  # creating the location specific expressions
  y <- parse_expr(paste("y[,", loc, "]", sep = ""))                      
  mu <- parse_expr(paste("mu", loc, sep = ""))                           
  sigma <- parse_expr(paste("sigma", loc, sep = ""))                     
  
  
  
  ### generic functions
  
  # pdf
  
  pdf_gen <- expr(  exp(-0.5 * log(2 * pi) - log( !!sigma ) - (3/2) * log( !!y ) - (( !!y - !!mu )^2)/(2 * (!!sigma)^2 * ((!!mu)^2) * (!!y))) ) 
  # cdf
  cdf_gen <- expr(  pnorm( ( ((!!y)/(!!mu)) - 1)/( (!!sigma) * sqrt((!!y))) ) + exp(2/( (!!mu)*(!!sigma)^2)) * pnorm(-( ((!!y)/(!!mu)) + 1)/((!!sigma) * sqrt((!!y)))) )
  
  
  l_generic <- list(pdf = pdf_gen,
                    cdf = cdf_gen)
  
  
  
  ### mu functions
  
  # pdf
  pdf_mu  <- expr( exp(-0.5 * log(2 * pi) - log( !!sigma ) - (3/2) * log( !!y ) - (( !!y - exp(f) )^2)/(2 * (!!sigma)^2 * (exp(f)^2) * (!!y))) )
  
  # derivativ logpdf
  derlpdf1.deretamu   <- expr( ( ( !!y - exp(f) ) / ( exp(f)^3 * (!!sigma)^2 ) ) * exp(f) )
  
  # cdf
  cdf_mu  <- expr( pnorm( ( ((!!y)/(exp(f))) - 1)/( (!!sigma) * sqrt((!!y))) ) + exp(2/( (exp(f))*(!!sigma)^2)) * pnorm(-( ((!!y)/(exp(f))) + 1)/((!!sigma) * sqrt((!!y)))) )
  
  # derivative cdf
  dercdf.deretamu  <- expr(-( exp(2/(exp(f) * (!!sigma)^2)) * ((!!y) * dnorm(-((1 + (!!y)/exp(f))/((!!sigma) * sqrt((!!y)))))/(exp(f)^2 * (!!sigma) * sqrt((!!y))) - 2 * ((!!sigma)^2 * pnorm(-((1 + (!!y)/exp(f))/((!!sigma) * sqrt((!!y)))))/(exp(f) * (!!sigma)^2)^2)) - (!!y) * dnorm(((!!y)/exp(f) -  1)/((!!sigma) * sqrt((!!y))))/(exp(f)^2 * (!!sigma) * sqrt((!!y))) ) * exp(f) )
  
  
  l_mu <- list(pdf = pdf_mu,
               derlpdf = derlpdf1.deretamu,
               cdf = cdf_mu,
               dercdf = dercdf.deretamu)
  
  
  
  ### sigma functions
  
  # pdf
  pdf_sig <- expr( exp(-0.5 * log(2 * pi) - log( exp(f) ) - (3/2) * log( !!y ) - (( !!y - !!mu )^2)/(2 * (exp(f))^2 * ((!!mu)^2) * (!!y))) )
  
  # derivative pdf
  derlpdf1.deretasig  <- expr( ( - (1 / exp(f)) + ( !!y - !!mu )^2 / ( (!!mu)^2 * exp(f)^3 * !!y ) ) * exp(f) )
  
  # cdf 
  cdf_sig <- expr( pnorm( ( ((!!y)/(!!mu)) - 1)/( (exp(f)) * sqrt((!!y))) ) + exp(2/( (!!mu)*(exp(f))^2)) * pnorm(-( ((!!y)/(!!mu)) + 1)/((exp(f)) * sqrt((!!y))))  )
  
  # derivative cdf 
  dercdf.deretasig <- expr( ( expr( ((1 + (!!y)/(!!mu)) * dnorm(-((1 + (!!y)/(!!mu))/(exp(f) * sqrt((!!y))))) * sqrt((!!y))/(exp(f) * sqrt((!!y)))^2 - 4 * ((!!mu) * exp(f) * pnorm(-((1 + (!!y)/(!!mu))/(exp(f) *  sqrt((!!y)))))/((!!mu) * exp(f)^2)^2)) * exp(2/((!!mu) * exp(f)^2)) - dnorm(((!!y)/(!!mu) - 1)/(exp(f) * sqrt((!!y)))) * sqrt((!!y)) * ((!!y)/(!!mu) - 1)/(exp(f) * sqrt((!!y)))^2        ) ) * exp(f) )
  
  
  l_sigma <- list(pdf = pdf_sig,
                  derlpdf = derlpdf1.deretasig,
                  cdf = cdf_sig,
                  dercdf = dercdf.deretasig)
  
  
  
  
  
  ### response functions 
  
  response_mu <- function(f) exp(f)
  response_sigma <- function(f) exp(f)
  
  
  l_response <- list(mu = response_mu,
                     sigma = response_sigma)
  
  ### offset functions     
  
  offset_mu <- expr({
    if (!is.null(!!mu)) {
      RET <- log(!!mu)
    }
    else {
      RET <- log(weighted.mean((!!y + weighted.mean(!!y, w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE)) # weighted version
    }
    return(RET)
  })
  
  offset_sigma <- expr({            # taken from gamlss
    if (!is.null(!!sigma)) {
      RET <- log(!!sigma)
    }
    else {
      
      #sigma <- rep(0.1, length(!!y))      
      sigma <- sqrt(weighted.var( !!y, w = w, na.rm = TRUE))/(weighted.mean(!!y, w = w, na.rm = TRUE) )^1.5 
      
      RET <- log(mean(sigma))
    }
    return(RET)
  })
  
  
  
  
  l_offset <- list(mu = offset_mu,
                   sigma = offset_sigma)
  
  ### names 
  
  name_mu <- "Inverse Gaussian distribution: mu(log link)"
  name_sigma <- "Inverse Gaussian distribution: sigma(log link)"
  
  l_names <- list(mu = name_mu,
                  sigma = name_sigma)
  
  ### check y function   
  
  
  check_y_mu <- expr({
    if (!is.numeric(!!y))
      stop("response is not numeric but ", sQuote("WeibullMar()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("WeibullMar()"))
    y
  })
  
  check_y_sigma <- expr({
    if (!is.numeric(!!y))
      stop("response is not numeric but ", sQuote("WeibullMar()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("WeibullMar()"))
    y
  })
  
  
  
  l_check_y <- list(mu = check_y_mu,
                    sigma = check_y_sigma)   
  
  
  
  
  # return list
  
  l = list(parameters = c(mu, sigma),
           parameter_names = c("mu", "sigma"),
           
           generic = l_generic,
           mu = l_mu,
           sigma = l_sigma,
           
           response = l_response,
           offset = l_offset,
           name = l_names,
           check_y = l_check_y,
           marg_name = "InverseGaussianMarg")
  
  
  return(l)
  
}


# mar1 <- LogLogistic_Mar(loc = 1)
