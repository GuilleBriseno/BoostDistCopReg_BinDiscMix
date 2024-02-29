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


PIG_Mar <- function(loc = NULL, offset_mu = NULL, offset_sigma = NULL, aux_mat = NULL, dim_aux_mat = NULL){
  
  
  # check for appropriate offset values for parameter mu
  if ((!is.null(offset_mu) && offset_mu <= 0))
    stop(sQuote("mu"), paste(" must be greater than zero in Marginal", loc))
  
  # check for appropriate offset values for parameter sigma
  if ((!is.null(offset_sigma) && offset_sigma <= 0))
    stop(sQuote("sigma"), paste(" must be greater than zero in Marginal", loc))
  
  
  # creating the location specific expressions
  y <- parse_expr(paste("y[,", loc, "]", sep = ""))                      
  mu <- parse_expr(paste("mu", loc, sep = ""))  
  sigma <- parse_expr(paste("sigma", loc, sep = ""))  
  
  
  
  ### generic functions
  # pdf
  pdf_gen <- expr( dPIG(x = !!y, mu = !!mu, sigma = !!sigma) )
  # cdf
  cdf_gen <- expr( pPIG(q = !!y, mu = !!mu, sigma = !!sigma) )  
  
  
  l_generic <- list(pdf = pdf_gen,
                    cdf = cdf_gen)
  
  
  ### mu functions
  # pdf
  pdf_mu  <- expr( dPIG(x = !!y, mu = exp(f), sigma = !!sigma ) )
  # der pdf 
  derpdf1.deretamu   <- expr({  derPIG_PDF.dermu2(resp = !!y, param1 = exp(f), param2 = !!sigma)   })
  
  # cdf
  cdf_mu  <- expr( pPIG(q = !!y, mu = exp(f), sigma = !!sigma ) )  
  
  dercdf.deretamu <-  expr( get_derCDF_etaparam(func = derPIG_PDF.dermu2, resp = !!y, par1 = exp(f), par2 = !!sigma) )
  
  
  l_mu <- list(pdf = pdf_mu,
               derpdf = derpdf1.deretamu,
               cdf = cdf_mu,
               dercdf = dercdf.deretamu)
  
  
  ### sigma functions
  # pdf
  pdf_sigma <- expr( dPIG(x = !!y, mu = !!mu, sigma = exp(f) )  )
  
  # der pdf
  derpdf1.deretasigma <- expr({ derPIG_PDF.dersigma2(resp = !!y, param1 = !!mu, param2 = exp(f))  })
  
  # cdf
  cdf_sigma <- expr( pPIG(q = !!y, mu = !!mu, sigma = exp(f) ) )  
  
  # der cdf
  dercdf.deretasigma <- expr( get_derCDF_etaparam(func = derPIG_PDF.dersigma2, resp = !!y, par1 = !!mu, par2 = exp(f)) )
  
  l_sigma <- list(pdf = pdf_sigma,
                  derpdf = derpdf1.deretasigma,
                  cdf = cdf_sigma,
                  dercdf = dercdf.deretasigma)
  
  
  ### response functions 
  
  response_mu <- function(f) exp(f)
  response_sigma <- function(f) exp(f)
  
  l_response <- list(mu = response_mu, 
                     sigma = response_sigma)
  
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
  
  
  
  l_offset <- list(mu = offset_mu,
                   sigma = offset_sigma)
  
  ### names 
  
  name_mu <- "Poisson Inverse Gaussian distribution: mu(log link)"
  name_sigma <- "Poisson Inverse Gaussian distribution: sigma(log link)"
  
  
  l_names <- list(mu = name_mu,
                  sigma = name_sigma)
  
  ### check y function   
  
  #check_y <- function(y) y                                                      
  check_y_mu <- expr({
    if (!is.numeric(!!y))
      stop("response is not but ", sQuote("PIG()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("PIG()"))
    y
  })
  
  
  check_y_sigma <- expr({
    if (!is.numeric(!!y))
      stop("response is not but ", sQuote("PIG()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("PIG()"))
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
           marg_name = "PIGMarg")
  
  
  return(l)
  
}

# marg <- LogNormal_Mar(loc = 1)
