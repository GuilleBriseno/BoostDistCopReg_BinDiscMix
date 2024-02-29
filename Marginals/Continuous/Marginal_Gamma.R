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


### Gamma  Marginal function

Gamma_Mar <- function(loc = NULL, offset_mu = NULL, offset_sigma = NULL){
  
  
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
  
  pdf_gen <- expr( dgamma(!!y, shape = 1/(!!sigma)^2, scale = (!!mu) * (!!sigma)^2) )
  # cdf
  cdf_gen <- expr( pgamma(!!y, shape = 1/(!!sigma)^2, scale = (!!mu) * (!!sigma)^2) )
  
  
  l_generic <- list(pdf = pdf_gen,
                    cdf = cdf_gen)
  
  
  
  ### mu functions
  
  # pdf
  pdf_mu  <- expr( dgamma(!!y, shape = 1/(!!sigma)^2, scale = exp(f) * (!!sigma)^2) )
  
  # derivativ logpdf
  derlpdf1.deretamu   <- expr( (-((!!sigma)^2 * exp((log(!!y) - (2 * log((!!sigma)) + log((exp(f))) + (!!y)/(exp(f))))/(!!sigma)^2 - (lgamma(1/(!!sigma)^2) + log((!!y)))) * ((exp(f)) - (!!y))/((exp(f)) * (!!sigma)^2)^2) ) * exp(f) )
  
  # cdf
  cdf_mu  <- expr( pgamma(!!y, shape = 1/(!!sigma)^2, scale = exp(f) * (!!sigma)^2) )
  
  # derivative cdf
  dercdf.deretamu  <- expr( (-( (exp(-((!!y)/((exp(f)) * (!!sigma)^2))) * (!!y) *((!!y)/((exp(f)) * (!!sigma)^2))^(-1 + 1/(!!sigma)^2))/((exp(f))^2 *(!!sigma)^2*gamma(1/(!!sigma)^2))) ) * exp(f) ) 
  
  
  l_mu <- list(pdf = pdf_mu,
               derlpdf = derlpdf1.deretamu,
               cdf = cdf_mu,
               dercdf = dercdf.deretamu)
  
  
  
  ### sigma functions
  
  # pdf
  pdf_sig <- expr( dgamma(!!y, shape = 1/(exp(f))^2, scale = (!!mu) * (exp(f))^2) )
  
  # derivative pdf
  derlpdf1.deretasig  <- expr( (  -(((2 * (log((!!y)) - (2 * log((exp(f))) + log((!!mu)))) - 2 * digamma(1/(exp(f))^2))/(exp(f))^3 + (!!mu) * (exp(f)) * (2 * (!!mu) - 2 * (!!y))/((!!mu) * (exp(f))^2)^2) * exp((log((!!y)) - (2 * log((exp(f))) + log((!!mu)) + (!!y)/(!!mu)))/(exp(f))^2 - (lgamma(1/(exp(f))^2) +  log((!!y)))))  ) * exp(f) )
  
  # cdf 
  cdf_sig <- expr( pgamma(!!y, shape = 1/(exp(f))^2, scale = (!!mu) * (exp(f))^2) )
  # derivative cdf 
  dercdf.deretasig <- expr( ( derGAMMA_CDF.dersigma(resp = !!y, param1 = !!mu, param2 = exp(f)) ) * exp(f))
  
  
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
      sigma <- rep(1, length(!!y))             
      RET <- log(mean(sigma))
    }
    return(RET)
  })
  
  
  
  
  l_offset <- list(mu = offset_mu,
                   sigma = offset_sigma)
  
  ### names 
  
  name_mu <- "Gamma distribution: mu(log link)"
  name_sigma <- "Gamma distribution: sigma(log link)"
  
  l_names <- list(mu = name_mu,
                  sigma = name_sigma)
  
  ### check y function   
  
  
  check_y_mu <- expr({
    if (!is.numeric(!!y))
      stop("response is not numeric but ", sQuote("Gamma()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("Gamma()"))
    y
  })
  
  check_y_sigma <- expr({
    if (!is.numeric(!!y))
      stop("response is not numeric but ", sQuote("Gamma()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("Gamma()"))
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
           marg_name = "GammaMarg")
  
  
  return(l)
  
}


# mar1 <- LogLogistic_Mar(loc = 1)
