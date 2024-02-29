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


LogNormal_Mar_GENCENS <- function(loc = NULL, offset_sigma = NULL){
  
  
  # check for appropriate offset values for parameter sigma
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
  pdf_gen <- expr(dlnorm(x = !!y, meanlog = !!mu, sdlog = !!sigma))
  
  # cdf
  cdf_gen_LEFT <- expr(plnorm(q = !!yleft, meanlog = !!mu, sdlog = !!sigma))
  
  cdf_gen_RIGHT <- expr(plnorm(q = !!yright, meanlog = !!mu, sdlog = !!sigma))
  
  
  l_generic <- list(pdf = pdf_gen,
                    cdfleft = cdf_gen_LEFT,
                    cdfright = cdf_gen_RIGHT)
  
  
  
  ### mu functions
  
  # pdf
  pdf_mu  <- expr(dlnorm(x = !!y, meanlog = f, sdlog = !!sigma))
  # der pdf 
  derlpdf1.deretamu   <- expr((log(!!y) - f)/(!!sigma)^2 * 1)
  
  # cdf
  cdf_mu_LEFT  <- expr(plnorm(q = !!yleft, meanlog = f, sdlog = !!sigma))
  
  cdf_mu_RIGHT  <- expr(plnorm(q = !!yright, meanlog = f, sdlog = !!sigma))
  
  # derivative cdf
  dercdf.deretamu_LEFT  <- expr(-dnorm(log(!!yleft), mean = f, sd = !!sigma) * 1)   
  
  dercdf.deretamu_RIGHT  <- expr(-dnorm(log(!!yright), mean = f, sd = !!sigma) * 1)   
  
  l_mu <- list(pdf = pdf_mu,
               #derpdf = derpdf.deretamu,
               derlpdf = derlpdf1.deretamu,
               cdfleft = cdf_mu_LEFT,
               cdfright = cdf_mu_RIGHT,
               dercdfleft = dercdf.deretamu_LEFT,
               dercdfright = dercdf.deretamu_RIGHT
               )
  
  
  
  ### sigma functions
  
  # pdf
  pdf_sig <- expr(dlnorm(x = !!y, meanlog = !!mu, sdlog = exp(f)))
  # derivative pdf
  derlpdf1.deretasig  <- expr((1/(exp(f)^3)) * ((log(!!y) - !!mu)^2 - exp(f)^2) * exp(f))
  
  # cdf 
  cdf_sig_LEFT <- expr(plnorm(q = !!yleft, meanlog = !!mu, sdlog = exp(f)))
  
  cdf_sig_RIGHT <- expr(plnorm(q = !!yright, meanlog = !!mu, sdlog = exp(f)))
  
  # derivative cdf 
  dercdf.deretasig_LEFT <- expr(-dnorm(x = ( log(!!yleft) - !!mu )/exp(f))*
                             ( log(!!yleft) - !!mu )/exp(f)^2 * exp(f))
  
  dercdf.deretasig_RIGHT <- expr(-dnorm(x = ( log(!!yright) - !!mu )/exp(f))*
                             ( log(!!yright) - !!mu )/exp(f)^2 * exp(f))
  
  
  l_sigma <- list(pdf = pdf_sig,
                  derlpdf = derlpdf1.deretasig,
                  cdfleft = cdf_sig_LEFT,
                  cdfright = cdf_sig_RIGHT,
                  dercdfleft = dercdf.deretasig_LEFT,
                  dercdfright = dercdf.deretasig_RIGHT
                  )
  
  
  
  
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
      RET <- weighted.mean((log(!!y) + weighted.mean(log(!!y), w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
    }
    return(RET)
  })
  
  offset_sigma <- expr({
    if (!is.null(!!sigma)) {
      RET <- log(!!sigma)
    }
    else {
      sigma <- rep(weighted.sd( log(!!y) , w = w, na.rm = TRUE), length(!!y)) # taken from gamlss
      RET <- log(mean(sigma))                                                 # but weighted version, taken from gamboostLSS
    }
    return(RET)
  })
  
  
  
  l_offset <- list(mu = offset_mu,
                   sigma = offset_sigma)
  
  ### names 
  
  name_mu <- "Log-Normal distribution: mu(id link)"
  name_sigma <- "Log-Normal distribution: sigma(log link)"
  
  l_names <- list(mu = name_mu,
                  sigma = name_sigma)
  
  ### check y function   
  
  #check_y <- function(y) y                                                      
  check_y_mu <- expr({
    if (!is.numeric(!!y))
      stop("response is not but ", sQuote("LogNormLSS()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("LogNormLSS()"))
    y
  })
  
  check_y_sigma <- expr({
    if (!is.numeric(!!y))
      stop("response is not but ", sQuote("LogNormLSS()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("LogNormLSS()"))
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
           marg_name = "LogNormMarg")
  
  
  return(l)
  
}

# marg <- LogNormal_Mar(loc = 1)
