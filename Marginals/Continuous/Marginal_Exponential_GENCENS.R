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


Exponential_Mar_GENCENS <- function(loc = NULL, offset_mu = NULL){
  
  
  # check for appropriate offset values for parameter sigma
  if ((!is.null(offset_mu) && offset_mu <= 0))
    stop(sQuote("mu"), paste(" must be greater than zero in Marginal", loc))
  
  
  
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
  #sigma <- parse_expr(paste("sigma", loc, sep = ""))                     
  
  
  
  ### generic functions
  
  # pdf
  pdf_gen <- expr( dexp(x = !!y, rate = !!mu ) ) #expr(dlnorm(x = !!y, meanlog = !!mu, sdlog = !!sigma))
  
  # cdf
  cdf_gen_LEFT <- expr( pexp(q = !!yleft, rate = !!mu ) ) #expr(plnorm(q = !!y, meanlog = !!mu, sdlog = !!sigma))
  
  cdf_gen_RIGHT <- expr( pexp(q = !!yright, rate = !!mu ) ) #expr(plnorm(q = !!y, meanlog = !!mu, sdlog = !!sigma))
  
  l_generic <- list(pdf = pdf_gen,
                    cdfleft = cdf_gen_LEFT,
                    cdfright = cdf_gen_RIGHT
                    )
  
  
  
  ### mu functions
  
  # pdf
  pdf_mu  <- expr( dexp(x = !!y, rate = exp(f)) ) #expr(dlnorm(x = !!y, meanlog = f, sdlog = !!sigma))
  # der pdf 
  derlpdf1.deretamu   <- expr(  1 - !!y * exp(f) ) #expr((log(!!y) - f)/(!!sigma)^2 * 1)
  
  # cdf
  cdf_mu_LEFT  <- expr( pexp(q = !!yleft, rate = exp(f)) )  #expr(plnorm(q = !!y, meanlog = f, sdlog = !!sigma))
  
  cdf_mu_RIGHT  <- expr( pexp(q = !!yright, rate = exp(f)) )  #expr(plnorm(q = !!y, meanlog = f, sdlog = !!sigma))
  
  # derivative cdf
  dercdf.deretamu_LEFT  <- expr( exp( - exp(f) * !!yleft ) * !!yleft * exp(f) )   
  
  dercdf.deretamu_RIGHT <- expr( exp( - exp(f) * !!yright ) * !!yright * exp(f) )   
  
  
  l_mu <- list(pdf = pdf_mu,
               #derpdf = derpdf.deretamu,
               derlpdf = derlpdf1.deretamu,
               cdfleft = cdf_mu_LEFT,
               cdfright = cdf_mu_RIGHT,
               dercdfleft = dercdf.deretamu_LEFT,
               dercdfright = dercdf.deretamu_RIGHT
               )
  
  
  
  ### sigma functions
  
  # # pdf
  # pdf_sig <- expr(dlnorm(x = !!y, meanlog = !!mu, sdlog = exp(f)))
  # # derivative pdf
  # derlpdf1.deretasig  <- expr((1/(exp(f)^3)) * ((log(!!y) - !!mu)^2 - exp(f)^2) * exp(f))
  # # cdf 
  # cdf_sig <- expr(plnorm(q = !!y, meanlog = !!mu, sdlog = exp(f)))
  # # derivative cdf 
  # dercdf.deretasig <- expr(-dnorm(x = ( log(!!y) - !!mu )/exp(f))*
  #                            ( log(!!y) - !!mu )/exp(f)^2 * exp(f))
  # 
  # 
  # 
  # l_sigma <- list(pdf = pdf_sig,
  #                 derlpdf = derlpdf1.deretasig,
  #                 cdf = cdf_sig,
  #                 dercdf = dercdf.deretasig)
  # 
  
  
  
  ### response functions 
  
  response_mu <- function(f) exp(f)
  #response_sigma <- function(f) exp(f)
  
  l_response <- list(mu = response_mu)
  
  ### offset functions     
  
  offset_mu <- expr({
    if (!is.null(!!mu)) {
      RET <- !!mu
    }
    else {
      RET <- weighted.mean(((!!y) + weighted.mean((!!y), w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
    }
    return(RET)
  })
  
  # offset_sigma <- expr({
  #   if (!is.null(!!sigma)) {
  #     RET <- log(!!sigma)
  #   }
  #   else {
  #     sigma <- rep(weighted.sd( log(!!y) , w = w, na.rm = TRUE), length(!!y)) # taken from gamlss
  #     RET <- log(mean(sigma))                                                 # but weighted version, taken from gamboostLSS
  #   }
  #   return(RET)
  # })
  
  
  
  l_offset <- list(mu = offset_mu)
  
  ### names 
  
  name_mu <- "Exponential distribution: mu(log link)"
  #name_sigma <- "Log-Normal distribution: sigma(log link)"
  
  l_names <- list(mu = name_mu)
  
  ### check y function   
  
  #check_y <- function(y) y                                                      
  # check_y_mu <- expr({
  #   if (!is.numeric(!!y))
  #     stop("response is not but ", sQuote("LogNormLSS()"))
  #   if (any(!!y < 0))
  #     stop("response is not positive but ", sQuote("LogNormLSS()"))
  #   y
  # })
  
  check_y_mu <- expr({
    if (!is.numeric(!!y))
      stop("response is not but ", sQuote("LogNormLSS()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("LogNormLSS()"))
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
           name = l_names,
           check_y = l_check_y,
           marg_name = "ExponentialMarg")
  
  
  return(l)
  
}

# marg <- LogNormal_Mar(loc = 1)
