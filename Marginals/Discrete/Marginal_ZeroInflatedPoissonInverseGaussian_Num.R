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


ZIPIG_Mar <- function(loc = NULL, offset_mu = NULL, offset_sigma = NULL, offset_nu = NULL, aux_mat = NULL, dim_aux_mat = NULL){
  
  
  # check for appropriate offset values for parameter mu
  if ((!is.null(offset_mu) && offset_mu <= 0))
    stop(sQuote("mu"), paste(" must be greater than zero in Marginal", loc))
  
  # check for appropriate offset values for parameter sigma
  if ((!is.null(offset_sigma) && offset_sigma <= 0))
    stop(sQuote("sigma"), paste(" must be greater than zero in Marginal", loc))
  
  # check for appropriate offset values for parameter sigma
  if ((!is.null(offset_nu) && offset_nu <= 0))
    stop(sQuote("nu"), paste(" must be greater than zero in Marginal", loc))
  
  
  # creating the location specific expressions
  y <- parse_expr(paste("y[,", loc, "]", sep = ""))                      
  mu <- parse_expr(paste("mu", loc, sep = ""))  
  sigma <- parse_expr(paste("sigma", loc, sep = ""))  
  nu <- parse_expr(paste("nu", loc, sep = ""))
  
  
  ### generic functions
  # pdf
  pdf_gen <- expr( dZIPIG(x = !!y, mu = !!mu, sigma = !!sigma, nu = !!nu) )
  # cdf
  cdf_gen <- expr( pZIPIG(q = !!y, mu = !!mu, sigma = !!sigma, nu = !!nu) )  
  
  
  l_generic <- list(pdf = pdf_gen,
                    cdf = cdf_gen)
  
  
  ### mu functions
  # pdf
  pdf_mu  <- expr( dZIPIG(x = !!y, mu = exp(f), sigma = !!sigma, nu = !!nu) )
  # der pdf 
  derpdf1.deretamu   <- expr({  derZIPIG_PDF.dermu2(resp = !!y, param1 = exp(f), param2 = !!sigma, param3 = !!nu)   })
  
  # cdf
  cdf_mu  <- expr( pZIPIG(q = !!y, mu = exp(f), sigma = !!sigma, nu = !!nu) )  
  
  dercdf.deretamu <- expr( get_derCDF_3param_etaparam(func = derZIPIG_PDF.dermu2, resp = !!y, par1 = exp(f), par2 = !!sigma, par3 = !!nu) )
  
  
  l_mu <- list(pdf = pdf_mu,
               derpdf = derpdf1.deretamu,
               cdf = cdf_mu,
               dercdf = dercdf.deretamu)
  
  
  ### sigma functions
  # pdf
  pdf_sigma <- expr( dZIPIG(x = !!y, mu = !!mu, sigma = exp(f), nu = !!nu)  )
  
  # der pdf
  derpdf1.deretasigma <- expr({ derZIPIG_PDF.dersigma2(resp = !!y, param1 = !!mu, param2 = exp(f), param3 = !!nu)  })
  
  # cdf
  cdf_sigma <- expr( pZIPIG(q = !!y, mu = !!mu, sigma = exp(f), nu = !!nu) )  
  
  # der cdf
  dercdf.deretasigma <- expr( get_derCDF_3param_etaparam(func = derZIPIG_PDF.dersigma2, resp = !!y, par1 = !!mu, par2 = exp(f), par3 = !!nu) )
  
  l_sigma <- list(pdf = pdf_sigma,
                  derpdf = derpdf1.deretasigma,
                  cdf = cdf_sigma,
                  dercdf = dercdf.deretasigma)
  
  
  ### nu functions
  # pdf
  pdf_nu <- expr(  dZINBI(x = !!y, mu = !!mu, sigma = !!sigma, nu = plogis(f))   )
  
  # der pdf
  derpdf1.deretanu <- expr({ derZIPIG_PDF.dernu2(resp = !!y, param1 = !!mu, param2 = !!sigma, param3 = plogis(f)) })
  
  # cdf
  cdf_nu <- expr( pZINBI(q = !!y, mu = !!mu, sigma = !!sigma, nu = plogis(f)) )  
  
  # der cdf
  dercdf.deretanu <- expr( get_derCDF_3param_etaparam(func = derZIPIG_PDF.dernu2, resp = !!y, par1 = !!mu, par2 = !!sigma, par3 = plogis(f)) )
  
  l_nu <- list(pdf = pdf_nu,
               derpdf = derpdf1.deretanu,
               cdf = cdf_nu,
               dercdf = dercdf.deretanu)
  
  
  ### response functions 
  
  response_mu <- function(f) exp(f)
  response_sigma <- function(f) exp(f)
  response_nu <- function(f) plogis(f)
  
  l_response <- list(mu = response_mu, 
                     sigma = response_sigma,
                     nu = response_nu)
  
  ### offset functions     
  
  offset_mu <- expr({
    if (!is.null(!!mu)) {
      RET <- log( !!mu  )
    }
    else {
      
      #RET <- weighted.mean((log(!!y) + weighted.mean(log(!!y), w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
      #temp <-  weighted.mean((log(!!y) + weighted.mean(log(!!y), w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
      #RET <- temp
      temp <- weighted.mean((!!y + weighted.mean(!!y, w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
      
      RET <- log(temp)     
      
    }
    return(RET)
  })
  
  
  offset_sigma <- expr({
    if (!is.null(!!sigma)) {
      RET <- log( !!sigma )
    }
    else {
      
      #RET <- weighted.mean((log(!!y) + weighted.mean(log(!!y), w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
      #temp <-  weighted.mean((log(!!y) + weighted.mean(log(!!y), w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
      #RET <- temp
      #temp <- weighted.mean((!!y + weighted.mean(!!y, w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
      temp <- rep(max( ( (weighted.sd((!!y), w = w, na.rm = TRUE)^2 - weighted.mean((!!y), w = w, na.rm = TRUE)) / (weighted.mean((!!y), w = w, na.rm = TRUE)^2) ), 0.1), length(!!y) ) 
      
      RET <- mean( log(temp) )
      
    }
    return(RET)
  })
  
  
  offset_nu <- expr({
    if (!is.null(!!nu)) {
      RET <- log(!!nu / (1 - !!nu) )
    }
    else {
      
      #RET <- weighted.mean((log(!!y) + weighted.mean(log(!!y), w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
      #temp <-  weighted.mean((log(!!y) + weighted.mean(log(!!y), w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
      #RET <- temp
      
      #temp <- weighted.mean((!!y + weighted.mean(!!y, w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
      temp <- rep(((sum(!!y == 0)/length(!!y))+0.01)/2, length(!!y))
      
      RET <- mean(log(temp / (1 - temp)))
      
    }
    return(RET)
  })
  
  
  
  l_offset <- list(mu = offset_mu,
                   sigma = offset_sigma,
                   nu = offset_nu)
  
  ### names 
  
  name_mu <- "Zero Inflated Poisson Inverse Gaussian distribution: mu(log link)"
  name_sigma <- "Zero Inflated Poisson Inverse Gaussian distribution: sigma(log link)"
  name_nu <- "Zero Inflated Poisson Inverse Gaussian distribution: nu(logit link)"
  
  
  l_names <- list(mu = name_mu,
                  sigma = name_sigma,
                  nu = name_nu)
  
  ### check y function   
  
  #check_y <- function(y) y                                                      
  check_y_mu <- expr({
    if (!is.numeric(!!y))
      stop("response is not but ", sQuote("ZIPIG()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("ZIPIG()"))
    y
  })
  
  
  check_y_sigma <- expr({
    if (!is.numeric(!!y))
      stop("response is not but ", sQuote("ZIPIG()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("ZIPIG()"))
    y
  })
  
  
  check_y_nu <- expr({
    if (!is.numeric(!!y))
      stop("response is not but ", sQuote("ZIPIG()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("ZIPIG()"))
    y
  })
  
  l_check_y <- list(mu = check_y_mu,
                    sigma = check_y_sigma,
                    nu = check_y_nu)
  
  
  # return list
  
  l = list(parameters = c(mu, sigma, nu),
           parameter_names = c("mu", "sigma", "nu"),
           
           generic = l_generic,
           mu = l_mu,
           sigma = l_sigma,
           nu = l_nu,
           
           response = l_response,
           offset = l_offset,
           name = l_names,
           check_y = l_check_y,
           marg_name = "ZIPIGMarg")
  
  
  return(l)
  
}

# marg <- LogNormal_Mar(loc = 1)
