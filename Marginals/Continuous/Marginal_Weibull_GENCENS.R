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


### Weibull function

Weibull_Mar_GENCENS <- function(loc = NULL, offset_mu = NULL, offset_sigma = NULL){
  
  
  # check for appropriate offset values for parameter mu, sigma and tau
  
  if ((!is.null(offset_mu) && offset_mu <= 0))
    stop(sQuote("mu"), paste(" must be greater than zero in Marginal", loc))
  
  if ((!is.null(offset_sigma) && offset_sigma <= 0))                                          
    stop(sQuote("sigma"), paste(" must be greater than zero in Marginal", loc))
  
  
  
  
  # creating the location specific expressions ----> the format of the response in R is: y = [ y1_left, y1_right, y2_left, y2_right, cens1, cens2 ]
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
  
  ## Parameters mu and sigma stay unmodified: 
  mu <- parse_expr(paste("mu", loc, sep = ""))                           
  sigma <- parse_expr(paste("sigma", loc, sep = ""))                     
  
  
  
  ### generic functions (all have a LEFT and RIGHT version)
  
  # pdf
  #pdf_gen_LEFT <- expr(dweibull(x = !!yleft, scale = !!mu, shape = !!sigma))
  
  pdf_gen <- expr(dweibull(x = !!yright, scale = !!mu, shape = !!sigma))
  
  
  # cdf
  cdf_gen_LEFT <- expr(pweibull(q = !!yleft, scale = !!mu, shape = !!sigma))
  
  cdf_gen_RIGHT <- expr(pweibull(q = !!yright, scale = !!mu, shape = !!sigma))
  
  
  l_generic <- list(#pdfleft = pdf_gen_LEFT,
                    pdf = pdf_gen,
                    cdfleft = cdf_gen_LEFT,
                    cdfright = cdf_gen_RIGHT)
  
  
  
  ### mu functions
  
  # pdf
  #pdf_mu_LEFT  <- expr(dweibull(x = !!yleft, scale = exp(f), shape = !!sigma))
  
  pdf_mu  <- expr(dweibull(x = !!yright, scale = exp(f), shape = !!sigma))
  
  
  # derivativ logpdf
  #derlpdf1.deretamu_LEFT   <- expr( (- 1 / (exp(f)) + (!!sigma - 1) * (- 1 / (exp(f))) + (!!sigma) * (!!yleft)^(!!sigma)*(exp(f))^(- (!!sigma) - 1) ) * exp(f) )
  
  derlpdf1.deretamu   <- expr( (- 1 / (exp(f)) + (!!sigma - 1) * (- 1 / (exp(f))) + (!!sigma) * (!!yright)^(!!sigma)*(exp(f))^(- (!!sigma) - 1) ) * exp(f) )
  
  # cdf
  cdf_mu_LEFT  <- expr(pweibull(q = !!yleft, scale = exp(f), shape = !!sigma))
  
  cdf_mu_RIGHT  <- expr(pweibull(q = !!yright, scale = exp(f), shape = !!sigma))
  
  # derivative cdf
  dercdf.deretamu_LEFT  <- expr(-(!!sigma * !!yleft * exp(- (!!yleft / (exp(f)) )^(!!sigma) ) * ( !!yleft/ (exp(f)) )^(!!sigma - 1) * (1 / (exp(f))^2 ) ) * exp(f) )
  
  dercdf.deretamu_RIGHT  <- expr(-(!!sigma * !!yright * exp(- (!!yright / (exp(f)) )^(!!sigma) ) * ( !!yright/ (exp(f)) )^(!!sigma - 1) * (1 / (exp(f))^2 ) ) * exp(f) )
  
  l_mu <- list(#pdfleft = pdf_mu_LEFT,
               pdf = pdf_mu,
               #derlpdfleft = derlpdf1.deretamu_LEFT,
               derlpdf = derlpdf1.deretamu,
               cdfleft = cdf_mu_LEFT,
               cdfright = cdf_mu_RIGHT,
               dercdfleft = dercdf.deretamu_RIGHT,
               dercdfright = dercdf.deretamu_LEFT)
  
  
  
  ### sigma functions
  
  # pdf
  #pdf_sig_LEFT <- expr(dweibull(x = !!y, scale = !!mu, shape = exp(f)))
  
  pdf_sig <- expr(dweibull(x = !!yright, scale = !!mu, shape = exp(f)))
  
  # derivative pdf
  #derlpdf1.deretasig_LEFT  <- expr( ( (1/(exp(f))) + log(!!y) - log(!!mu) - (!!y/ (!!mu))^(exp(f)) * log(!!y / (!!mu) )  ) * exp(f) )
  
  derlpdf1.deretasig  <- expr( ( (1/(exp(f))) + log(!!yright) - log(!!mu) - (!!yright / (!!mu))^(exp(f)) * log(!!yright / (!!mu) )  ) * exp(f) )
  
  # cdf 
  cdf_sig_LEFT <- expr(pweibull(q = !!yleft, scale = !!mu, shape = exp(f)))
  
  cdf_sig_RIGHT <- expr(pweibull(q = !!yright, scale = !!mu, shape = exp(f)))
  
  # derivative cdf 
  dercdf.deretasig_LEFT <- expr( ( exp(- (!!yleft / (!!mu) )^(exp(f)) ) * ( (!!yleft / (!!mu) )^(exp(f)) * log(!!yleft / (!!mu) ) ) ) * exp(f) )
  
  dercdf.deretasig_RIGHT <- expr( ( exp(- (!!yright / (!!mu) )^(exp(f)) ) * ( (!!yright / (!!mu) )^(exp(f)) * log(!!yright / (!!mu) ) ) ) * exp(f) )
  
  
  l_sigma <- list(#pdfleft = pdf_sig_LEFT,
                  pdf = pdf_sig,
                  #derlpdfleft = derlpdf1.deretasig_LEFT,
                  derlpdf = derlpdf1.deretasig,
                  cdfleft = cdf_sig_LEFT,
                  cdfright = cdf_sig_RIGHT,
                  dercdfleft = dercdf.deretasig_LEFT,
                  dercdfright = dercdf.deretasig_RIGHT)
  
  
  
  
  
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
  
  name_mu <- "Weibull distribution: mu(log link)"
  name_sigma <- "Weibull distribution: sigma(log link)"
  
  l_names <- list(mu = name_mu,
                  sigma = name_sigma)
  
  ### check y function   
  
  
  check_y_mu <- expr({
    if (!is.numeric(!!y))
      stop("response is not numeric but ", sQuote("WeibullMar()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("WeibullMar()"))
    y
  })
  
  check_y_sigma <- expr({
    if (!is.numeric(!!y))
      stop("response is not numeric but ", sQuote("WeibullMar()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("WeibullMar()"))
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
           marg_name = "WeibullMarg")
  
  
  return(l)
  
}


# mar1 <- LogLogistic_Mar(loc = 1)
