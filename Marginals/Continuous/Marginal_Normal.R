################################################################################
### DESCRIPTION
### This file defines the Log-Normal-Marginal function. It contains the 
### expressions of the pdf, cdf, its derivatives wrt the two parameters, names,
### response functions, offset functions, check_y functions required to create an mboost Family object.
### This function is handed over to a copula function and the expressions and functions
### are merged with the respective copula via the rlang package, in order to create 
### appropriate Families objects from gamboostLSS.


### libraries 
library(mboost)
library(gamboostLSS)
library(rlang)


Normal_Mar <- function(loc = NULL, offset_sigma = NULL){
  
  
  # check for appropriate offset values for parameter sigma
  if ((!is.null(offset_sigma) && offset_sigma <= 0))
    stop(sQuote("sigma"), paste(" must be greater than zero in Marginal", loc))
  
  
  
  # creating the location specific expressions
  y <- parse_expr(paste("y[,", loc, "]", sep = ""))                      
  mu <- parse_expr(paste("mu", loc, sep = ""))                           
  sigma <- parse_expr(paste("sigma", loc, sep = ""))                     
  
  
  
  ### generic functions
  
  # pdf
  
  pdf_gen <- expr(dnorm(x = !!y, mean = !!mu, sd = !!sigma))
  # cdf
  cdf_gen <- expr(pnorm(q = !!y, mean = !!mu, sd = !!sigma))
  
  
  l_generic <- list(pdf = pdf_gen,
                    cdf = cdf_gen)
  
  
  
  ### mu functions
  
  # pdf
  pdf_mu  <- expr(dnorm(x = !!y, mean = f, sd = !!sigma))
  # der pdf 
  derlpdf1.deretamu   <- expr((!!y - f)/(!!sigma)^2 * 1)
  # cdf
  cdf_mu  <- expr(pnorm(q = !!y, mean = f, sd = !!sigma))
  # derivative cdf
  dercdf.deretamu  <- expr(-dnorm(!!y, mean = f, sd = !!sigma) * 1 )   
  
  
  l_mu <- list(pdf = pdf_mu,
               #derpdf = derpdf.deretamu,
               derlpdf = derlpdf1.deretamu,
               cdf = cdf_mu,
               dercdf = dercdf.deretamu)
  
  
  
  ### sigma functions
  
  # pdf
  pdf_sig <- expr( dnorm(x = !!y, mean = !!mu, sd = exp(f)) )
  # derivative pdf
  derlpdf1.deretasig  <- expr( (-1/exp(f) + (!!y - !!mu)^2/exp(f)^3 ) * exp(f) )
  # cdf 
  cdf_sig <- expr( pnorm(q = !!y, mean = !!mu, sd = exp(f)) )
  # derivative cdf 
  dercdf.deretasig <- expr( -dnorm( (!!y - !!mu )/exp(f) ) * (!!y - !!mu)/exp(f)^2 * exp(f) )#expr(-dnorm(x = ( log(!!y) - !!mu )/exp(f)) * ( log(!!y) - !!mu )/exp(f)^2 * exp(f))
  
  
  
  l_sigma <- list(pdf = pdf_sig,
                  derlpdf = derlpdf1.deretasig,
                  cdf = cdf_sig,
                  dercdf = dercdf.deretasig)
  
  
  
  
  ### response functions 
  
  response_mu <- function(f) f
  response_sigma <- function(f) exp(f)
  
  l_response <- list(mu = response_mu,
                     sigma = response_sigma)
  
  ### offset functions     
  
  offset_mu <- expr({
    if (!is.null(!!mu)) {
      RET <- !!mu
    }
    else {
      RET <- weighted.mean((!!y + weighted.mean(!!y, w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
    }
    return(RET)
  })

  offset_sigma <- expr({
    if (!is.null(!!sigma)) {
      RET <- log(!!sigma)
    }
    else {
      sigma <- rep(weighted.sd( !!y , w = w, na.rm = TRUE), length(!!y)) # taken from gamlss
      RET <- log(mean(sigma))                                                 # but weighted version, taken from gamboostLSS
    }
    return(RET)
  })
  
  
  
  l_offset <- list(mu = offset_mu,
                   sigma = offset_sigma)
  
  ### names 
  
  name_mu <- "Normal distribution: mu(id link)"
  name_sigma <- "Normal distribution: sigma(log link)"
  
  l_names <- list(mu = name_mu,
                  sigma = name_sigma)
  
  ### check y function   
  
  #check_y <- function(y) y                                                      
  check_y_mu <- expr({
    if (!is.numeric(!!y))
      stop("response is not but ", sQuote("NormalMar()"))
    y
  })
  
  check_y_sigma <- expr({
    if (!is.numeric(!!y))
      stop("response is not but ", sQuote("NormalLSS()"))
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
           marg_name = "NormMarg")
  
  
  return(l)
  
}

# marg <- LogNormal_Mar(loc = 1)
