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


BernoulliLogit_Mar <- function(loc = NULL, offset_mu = NULL){
  
  
  # check for appropriate offset values for parameter mu
  if ((!is.null(offset_mu) && (offset_mu < 0 || offset_mu > 1) ))
    stop(sQuote("mu"), paste(" must be within the interval (0,1) in Marginal", loc))
  
  
  
  # creating the location specific expressions
  y <- parse_expr(paste("y[,", loc, "]", sep = ""))                      
  mu <- parse_expr(paste("mu", loc, sep = ""))                      
  
  
  ### generic functions
  
  # pdf
  # PDF is not required here
  
  # cdf
  cdf_gen <- expr( !!mu ) 
  
  # In case of mixed responses
  cdf_gen_mixresp <- expr( 1 - (!!mu) )
  
  
  l_generic <- list(cdf = cdf_gen,
                    cdf_mixresp = cdf_gen_mixresp)
  
  
  
  ### mu functions
  
  # pdf
  # PDF is not required here
  
  # der pdf 
  # PDF is not required here
  # cdf
  cdf_mu  <- expr( plogis(f) ) 
  # derivative cdf
  dercdf.deretamu  <- expr( dlogis(f) )   
  
  
  # mixed responses
  cdf_mu_mixresp  <- expr( 1 - plogis(f) ) 
  # derivative cdf 
  dercdf.deretamu_mixresp  <- expr( -( ( 1 - exp(-f)/(1 + exp(-f)) ) * exp(-f)/(1 + exp(-f)) )  )  
  
  
  l_mu <- list(cdf = cdf_mu,
               dercdf = dercdf.deretamu,
               cdf_mixresp = cdf_mu_mixresp,
               dercdf_mixresp = dercdf.deretamu_mixresp)
  
  
  ### response functions 
  
  response_mu <- function(f) plogis(f)

  l_response <- list(mu = response_mu)
  
  ### offset functions     
  
  offset_mu <- expr({
    if (!is.null(!!mu)) {
      RET <- qlogis(!!mu)
    }
    else {
      #RET <- weighted.mean((!!y + weighted.mean(!!y, w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
    RET <- 0
      }
    return(RET)
  })
  
  
  
  offset_mu_mixresp <- expr({
    if (!is.null(!!mu)) {
      RET <- qlogis(!!mu)
    }
    else {
      
      RET <- qlogis( weighted.mean( ( !!y + 0.5 )/2, w = w, na.rm = TRUE) ) # weighted version
      
    }
    return(RET)
  })
  
  
  
  l_offset <- list(mu = offset_mu)
  
  l_offset_mixresp <- list(mu = offset_mu_mixresp)
  
  ### names 
  
  name_mu <- "Logit link function for binary responses: mu(logit link)"
  
  
  l_names <- list(mu = name_mu)
  
  ### check y function   
  
  #check_y <- function(y) y                                                      
  check_y_mu <- expr({
    if (!is.numeric(!!y))
      #stop("response is not but ", sQuote("LogNormLSS()"))
      stop("response is not but ", sQuote("LogitMargin()"))
    if (any(!!y < 0))
      #stop("response is not positive but ", sQuote("LogNormLSS()"))
      stop("response is not positive but ", sQuote("LogitMargin()"))
    if (any(!!y > 1))
      stop("response is above 1 but ", sQuote("LogitMargin()"))
    if ( any(unique(!!y) %in% c(0,1) == FALSE) )
      stop("response is not binary but ", sQuote("LogitMargin()"))
    y
  })
  
  
  l_check_y <- list(mu = check_y_mu)
  
  
  # return list
  
  l = list(parameters = c(mu),
           parameter_names = c("mu"),
           
           generic = l_generic,
           mu = l_mu,
           #sigma = l_sigma,
           
           response = l_response,
           offset = l_offset,
           offset_mixresp = l_offset_mixresp,
           name = l_names,
           check_y = l_check_y,
           marg_name = "LogitMarg")
  
  
  return(l)
  
}

# marg <- LogNormal_Mar(loc = 1)
