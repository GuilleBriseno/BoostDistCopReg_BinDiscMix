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


### Log-Logistic-Marginal function

LogLogistic_Mar_GENCENS <- function(loc = NULL, offset_mu = NULL, offset_sigma = NULL){
  
  
  # check for appropriate offset values for parameter mu, sigma and tau
  
  if ((!is.null(offset_mu) && offset_mu <= 0))
    stop(sQuote("mu"), paste(" must be greater than zero in Marginal", loc))
  
  if ((!is.null(offset_sigma) && offset_sigma <= 0))                                          
    stop(sQuote("sigma"), paste(" must be greater than zero in Marginal", loc))
  
  
  
  
  # creating the location specific expressions
  #y <- parse_expr(paste("y[,", loc, "]", sep = ""))                      
  
  if(loc == 1){
    
    yleft <- parse_expr(paste("y[,", 1, "]", sep = ""))
    
    yright <- parse_expr(paste("y[,", 2, "]", sep = ""))
    
    y <- parse_expr(paste("y[,", 2, "]", sep = "")) # this y is used only for the offset. 
    
  }
  
  if(loc == 2){
    
    yleft <- parse_expr(paste("y[,", 3, "]", sep = ""))
    
    yright <- parse_expr(paste("y[,", 4, "]", sep = ""))
    
    y <- parse_expr(paste("y[,", 4, "]", sep = "")) # this y is used only for the offset. 
    
  }
  
  mu <- parse_expr(paste("mu", loc, sep = ""))                           
  sigma <- parse_expr(paste("sigma", loc, sep = ""))                     
  
  
  
  ### generic functions
  
  # pdf
  pdf_gen <- expr(dLogLogistic(x = !!y, mu = !!mu, sigma = !!sigma))
  
  # cdf
  cdf_gen_LEFT <- expr(pLogLogistic(x = !!yleft, mu = !!mu, sigma = !!sigma))
  
  cdf_gen_RIGHT <- expr(pLogLogistic(x = !!yright, mu = !!mu, sigma = !!sigma))
  
  l_generic <- list(pdf = pdf_gen,
                    cdfleft = cdf_gen_LEFT,
                    cdfright = cdf_gen_RIGHT)
  
  
  
  ### mu functions
  
  # pdf
  pdf_mu  <- expr(dLogLogistic(x = !!y, mu = exp(f), sigma = !!sigma))
  # derivativ logpdf
  derlpdf1.deretamu   <- expr(( -(!!sigma)/exp(f) + 2 * (!!sigma) * (!!y / exp(f))^(!!sigma) / (exp(f) * (1 + (!!y / exp(f))^(!!sigma))) ) *
                                exp(f))
  # cdf
  cdf_mu_LEFT  <- expr(pLogLogistic(x = !!yleft, mu = exp(f), sigma = !!sigma))
  
  cdf_mu_RIGHT  <- expr(pLogLogistic(x = !!yright, mu = exp(f), sigma = !!sigma))
  
  # derivative cdf
  dercdf.deretamu_LEFT  <- expr(- ( (!!sigma)/exp(f) * ((!!yleft)/exp(f))^(-(!!sigma)) / (1 + ((!!yleft)/exp(f))^(-(!!sigma)))^2 ) *
                             exp(f)) 
  
  dercdf.deretamu_RIGHT  <- expr(- ( (!!sigma)/exp(f) * ((!!yright)/exp(f))^(-(!!sigma)) / (1 + ((!!yright)/exp(f))^(-(!!sigma)))^2 ) *
                             exp(f))
  
  
  l_mu <- list(pdf = pdf_mu,
               derlpdf = derlpdf1.deretamu,
               cdfleft = cdf_mu_LEFT,
               cdfright = cdf_mu_RIGHT,
               dercdfleft = dercdf.deretamu_LEFT,
               dercdfright = dercdf.deretamu_RIGHT
               )
  
  
  
  ### sigma functions
  
  # pdf
  pdf_sig <- expr(dLogLogistic(x = !!y, mu = !!mu, sigma = exp(f)))
  # derivative pdf
  derlpdf1.deretasig  <- expr(( (1/exp(f)) + log(!!y) - log(!!mu) - 2/(1 + ((!!y) / (!!mu))^(exp(f))) * ((!!y) / (!!mu))^(exp(f)) * log((!!y)/(!!mu)) ) *
                                exp(f))
  
  # cdf 
  cdf_sig_LEFT <- expr(pLogLogistic( x = !!yleft, mu = !!mu, sigma = exp(f)))
  
  cdf_sig_RIGHT <- expr(pLogLogistic( x = !!yright, mu = !!mu, sigma = exp(f)))
  
  # derivative cdf 
  dercdf.deretasig_LEFT <- expr(( ((!!yleft) / (!!mu))^(-(exp(f))) * log((!!yleft)/(!!mu)) / (1 + ((!!yleft) / (!!mu))^(-(exp(f))))^2  ) *
                             exp(f))
  
  dercdf.deretasig_RIGHT <- expr(( ((!!yright) / (!!mu))^(-(exp(f))) * log((!!yright)/(!!mu)) / (1 + ((!!yright) / (!!mu))^(-(exp(f))))^2  ) *
                             exp(f))
  
  

    
  l_sigma <- list(pdf = pdf_sig,
                  derlpdf = derlpdf1.deretasig,
                  cdfleft = cdf_sig_LEFT,
                  cdfright = cdf_sig_RIGHT,
                  dercdfleft = dercdf.deretasig_LEFT,
                  dercdfright = dercdf.deretasig_RIGHT
                  )
  
  
  
  
  
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
      sigma <- rep(0.1, length(!!y))             
      RET <- log(mean(sigma))
    }
    return(RET)
  })
  
  
  
  
  l_offset <- list(mu = offset_mu,
                   sigma = offset_sigma)
  
  ### names 
  
  name_mu <- "Log-Logistic distribution: mu(log link)"
  name_sigma <- "Log-Logistic distribution: sigma(log link)"
  
  l_names <- list(mu = name_mu,
                  sigma = name_sigma)
  
  ### check y function   
  
  
  check_y_mu <- expr({
    if (!is.numeric(!!y))
      stop("response is not numeric but ", sQuote("LogLogisticLSS()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("LogLogisticLSS()"))
    y
  })
  
  check_y_sigma <- expr({
    if (!is.numeric(!!y))
      stop("response is not numeric but ", sQuote("LogLogisticLSS()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("LogLogisticLSS()"))
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
           marg_name = "LogLogisticMarg")
  
  
  return(l)
  
}


# mar1 <- LogLogistic_Mar(loc = 1)
