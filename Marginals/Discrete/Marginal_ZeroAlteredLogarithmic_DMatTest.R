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


ZALG_Mar <- function(loc = NULL, offset_mu = NULL, offset_sigma = NULL, aux_mat = NULL, dim_aux_mat = NULL){
  
  
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
  pdf_gen <- expr( dZALG(x = !!y, mu = !!mu, sigma = !!sigma) )
  # cdf
  cdf_gen <- expr( pZALG(q = !!y, mu = !!mu, sigma = !!sigma) )  
  
  
  l_generic <- list(pdf = pdf_gen,
                    cdf = cdf_gen)
  
  
  ### mu functions
  # pdf
  pdf_mu  <- expr( dZALG(x = !!y, mu = plogis(f), sigma = !!sigma ) )
  # der pdf 
  #derpdf1.deretamu   <- expr( ( ( exp( -exp(f) ) * (exp(f)^(!!y - 1 ) * (!!y) - exp(f)^(!!y) ) ) / factorial(!!y) ) * exp(f) ) 
  derpdf1.deretamu   <- expr({  ifelse( !!y > 0, 
                                        (plogis(f)^(-1 + (!!y)) * (-1 + (!!sigma)) * (-plogis(f) + (!!y) * (-1 + plogis(f)) * log(1 - plogis(f))))/((!!y) * (-1 + plogis(f)) * log(1 - plogis(f))^2), 
                                        0) * (1 - plogis(f)) * plogis(f) 
  })
  
  # cdf
  cdf_mu  <- expr( pZALG(q = !!y, mu = plogis(f), sigma = !!sigma ) )  
  
  #dercdf.deretamu <-  expr( get_derCDF_etaparam(func = derZALG_PDF.dermu2, resp = !!y, par1 = plogis(f), par2 = !!sigma) )
  dercdf.deretamu <-  expr( derDISCDIST_PDF.derparam2(func = derZALG_PDF.dermu2, auxmatrix = !!aux_mat, dim_aux_matrix = !!dim_aux_mat, resp = !!y, par1 = plogis(f), par2 = !!sigma) )
  
  l_mu <- list(pdf = pdf_mu,
               derpdf = derpdf1.deretamu,
               cdf = cdf_mu,
               dercdf = dercdf.deretamu)
  
  
  ### sigma functions
  # pdf
  pdf_sigma <- expr( dZALG(x = !!y, mu = !!mu, sigma = plogis(f) )  )
  
  # der pdf
  derpdf1.deretasigma <- expr({ ifelse( !!y > 0,
                                        -((-(log(1 - (!!mu)))^(-1)) * (!!mu)^((!!y))/(!!y)),
                                        1) * ( (plogis(f)/(1 - plogis(f))*(1 + plogis(f)/(1 - plogis(f)))-(plogis(f)/(1 - plogis(f)))^2) / (1 + plogis(f)/(1 - plogis(f)))^2 )  
  })
  
  # cdf
  cdf_sigma <- expr( pZALG(q = !!y, mu = !!mu, sigma = plogis(f) ) )  
  
  # der cdf
  #dercdf.deretasigma <- expr( get_derCDF_etaparam(func = derZALG_PDF.dersigma2, resp = !!y, par1 = !!mu, par2 = plogis(f)) )
  dercdf.deretasigma <- expr( derDISCDIST_PDF.derparam2(func = derZALG_PDF.dersigma2, auxmatrix = !!aux_mat, dim_aux_matrix = !!dim_aux_mat, resp = !!y, par1 = !!mu, par2 = plogis(f)) )
  
  l_sigma <- list(pdf = pdf_sigma,
                  derpdf = derpdf1.deretasigma,
                  cdf = cdf_sigma,
                  dercdf = dercdf.deretasigma)
  
  
  ### response functions 
  
  response_mu <- function(f) plogis(f)
  response_sigma <- function(f) plogis(f)
  
  l_response <- list(mu = response_mu, 
                     sigma = response_sigma)
  
  ### offset functions     
  
  offset_mu <- expr({
    if (!is.null(!!mu)) {
      RET <- log( !!mu / (1 - !!mu) )
    }
    else {
      
      #RET <- weighted.mean((log(!!y) + weighted.mean(log(!!y), w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
      #temp <-  weighted.mean((log(!!y) + weighted.mean(log(!!y), w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
      #RET <- temp
      temp <- weighted.mean((!!y + weighted.mean(!!y, w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
      
      temp <- temp / ( temp + 1)
      
      RET <- log(temp/(1 - temp))     
      
    }
    return(RET)
  })
  
  
  offset_sigma <- expr({
    if (!is.null(!!sigma)) {
      RET <- log( !!sigma / (1 - !!sigma) )
    }
    else {
      
      #RET <- weighted.mean((log(!!y) + weighted.mean(log(!!y), w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
      #temp <-  weighted.mean((log(!!y) + weighted.mean(log(!!y), w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
      #RET <- temp
      #temp <- weighted.mean((!!y + weighted.mean(!!y, w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
      temp <- rep(0.3, length(!!y))
      
      RET <- mean( log(temp/(1 - temp)) )
      
    }
    return(RET)
  })
  
  
  
  l_offset <- list(mu = offset_mu,
                   sigma = offset_sigma)
  
  ### names 
  
  name_mu <- "Zero-altered Logarithmic distribution: mu(logit link)"
  name_sigma <- "Zero-altered Logarithmic distribution: sigma(logit link)"
  
  
  l_names <- list(mu = name_mu,
                  sigma = name_sigma)
  
  ### check y function   
  
  #check_y <- function(y) y                                                      
  check_y_mu <- expr({
    if (!is.numeric(!!y))
      stop("response is not but ", sQuote("ZALG()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("ZALG()"))
    y
  })
  
  
  check_y_sigma <- expr({
    if (!is.numeric(!!y))
      stop("response is not but ", sQuote("ZALG()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("ZALG()"))
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
           marg_name = "ZALGMarg")
  
  
  return(l)
  
}

# marg <- LogNormal_Mar(loc = 1)
