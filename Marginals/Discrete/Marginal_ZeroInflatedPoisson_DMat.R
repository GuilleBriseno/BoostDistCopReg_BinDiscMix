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


ZIP_Mar <- function(loc = NULL, offset_mu = NULL, offset_sigma = NULL, aux_mat = NULL, dim_aux_mat = NULL){
  
  
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
  pdf_gen <- expr( dZIP(x = !!y, mu = !!mu, sigma = !!sigma) )
  # cdf
  cdf_gen <- expr( pZIP(q = !!y, mu = !!mu, sigma = !!sigma) )  
  
  
  l_generic <- list(pdf = pdf_gen,
                    cdf = cdf_gen)
  
  
  ### mu functions
  # pdf
  pdf_mu  <- expr( dZIP(x = !!y, mu = exp(f), sigma = !!sigma ) )
  # der pdf 
  derpdf1.deretamu   <- expr({  ifelse( !!y > 0, 
                                          (1 - (!!sigma)) * ((exp(f))^(((!!y)) - 1) * ((!!y)))/gamma((!!y) + 1) * exp(-(exp(f))) - 
                                           (1 - (!!sigma)) * (exp(f))^((!!y))/gamma((!!y) + 1) * exp(-(exp(f))), 
                                          -((1 - (!!sigma)) * exp(-(exp(f))))) * exp(f)  })
  
  # cdf
  cdf_mu  <- expr( pZIP(q = !!y, mu = exp(f), sigma = !!sigma ) )  
  
  dercdf.deretamu <-  expr( get_derCDF_etaparam(func = derZIP_PDF.dermu2, resp = !!y, par1 = exp(f), par2 = !!sigma) )
  
  
  l_mu <- list(pdf = pdf_mu,
               derpdf = derpdf1.deretamu,
               cdf = cdf_mu,
               dercdf = dercdf.deretamu)
  
  
  ### sigma functions
  # pdf
  pdf_sigma <- expr( dZIP(x = !!y, mu = !!mu, sigma = exp(f) )  )
  
  # der pdf
  derpdf1.deretasigma <- expr({  ifelse(!!y > 0, 
                                          -((!!mu)^((!!y))/gamma((!!y) + 1) * exp(-(!!mu))), 
                                          1 - exp(-(!!mu))) * ((plogis(f)/(1-plogis(f))*(1+plogis(f)/(1-plogis(f)))-(plogis(f)/(1-plogis(f)))^2)/(1+plogis(f)/(1-plogis(f)))^2)   })
  
  # cdf
  cdf_sigma <- expr( pZIP(q = !!y, mu = !!mu, sigma = plogis(f) ) )  
  
  # der cdf
  dercdf.deretasigma <- expr( get_derCDF_etaparam(func = derZIP_PDF.dersigma2, resp = !!y, par1 = !!mu, par2 = plogis(f)) )
  
  l_sigma <- list(pdf = pdf_sigma,
                  derpdf = derpdf1.deretasigma,
                  cdf = cdf_sigma,
                  dercdf = dercdf.deretasigma)
  
  
  ### response functions 
  
  response_mu <- function(f) exp(f)
  response_sigma <- function(f) plogis(f)
  
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
      RET <- log(!!sigma / (1 - !!sigma) )
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
                   sigma = offset_sigma)
  
  ### names 
  
  name_mu <- "Zero Inflated Poisson distribution: mu(log link)"
  name_sigma <- "Zero Inflated Poisson distribution: sigma(log link)"
  
  
  l_names <- list(mu = name_mu,
                  sigma = name_sigma)
  
  ### check y function   
  
  #check_y <- function(y) y                                                      
  check_y_mu <- expr({
    if (!is.numeric(!!y))
      stop("response is not but ", sQuote("ZIP()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("ZIP()"))
    y
  })
  
  
  check_y_sigma <- expr({
    if (!is.numeric(!!y))
      stop("response is not but ", sQuote("ZIP()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("ZIP()"))
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
           marg_name = "ZIPMarg")
  
  
  return(l)
  
}

# marg <- LogNormal_Mar(loc = 1)
