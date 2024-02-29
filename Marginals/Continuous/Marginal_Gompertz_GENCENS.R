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


### Gompertz marginal

Gompertz_Mar_GENCENS <- function(loc = NULL, offset_mu = NULL, offset_sigma = NULL){
  
  
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
  
  pdf_gen <- expr( (!!sigma) * (!!mu) * exp( (!!mu) + (!!sigma) * (!!y) - (!!mu) * exp( (!!sigma) * (!!y) ) ) )
  
  # cdf
  cdf_gen_LEFT <- expr( 1 - exp(- (!!mu) * ( exp( (!!sigma) * (!!yleft) ) - 1 ) ) )
  
  cdf_gen_RIGHT <- expr( 1 - exp(- (!!mu) * ( exp( (!!sigma) * (!!yright) ) - 1 ) ) )
  
  
  l_generic <- list(pdf = pdf_gen,
                    cdfleft = cdf_gen_LEFT,
                    cdfright = cdf_gen_RIGHT)
  
  
  
  ### mu functions
  
  # pdf
  pdf_mu  <- expr( (!!sigma) * (exp(f)) * exp( (exp(f)) + (!!sigma) * (!!y) - (exp(f)) * exp( (!!sigma) * (!!y) ) ) )
  
  # derivativ logpdf
  derlpdf1.deretamu   <- expr( ( (1/exp(f)) + 1 - exp((!!sigma) * (!!y) ) ) * exp(f) )
  
  # cdf
  cdf_mu_LEFT  <- expr( 1 - exp(- exp(f) * ( exp( (!!sigma) * (!!yleft) ) - 1 ) ) )
  
  cdf_mu_RIGHT  <- expr( 1 - exp(- exp(f) * ( exp( (!!sigma) * (!!yright) ) - 1 ) ) )
  
  # derivative cdf
  dercdf.deretamu_LEFT  <- expr(-( ( exp((!!sigma) * (!!yleft)) - 1) * exp( exp(f) - exp(f) * exp( (!!sigma) * (!!yleft)) ) ) * exp(f) )
  
  dercdf.deretamu_RIGHT <- expr(-( ( exp((!!sigma) * (!!yright)) - 1) * exp( exp(f) - exp(f) * exp( (!!sigma) * (!!yright)) ) ) * exp(f) )
  
  
  l_mu <- list(pdf = pdf_mu,
               derlpdf = derlpdf1.deretamu,
               cdfleft = cdf_mu_LEFT,
               cdfright = cdf_mu_RIGHT,
               dercdfleft = dercdf.deretamu_LEFT,
               dercdfright = dercdf.deretamu_RIGHT
               )
  
  
  
  ### sigma functions
  
  # pdf
  pdf_sig <- expr( (exp(f)) * (!!mu) * exp( (!!mu) + (exp(f)) * (!!y) - (!!mu) * exp( (exp(f)) * (!!y) ) ) )
  
  # derivative pdf
  derlpdf1.deretasig  <- expr( ( (1/exp(f)) + (!!y) - (!!y) * (!!mu) * exp(exp(f) * (!!y) )  ) * exp(f) )
  
  # cdf 
  cdf_sig_LEFT <- expr( 1 - exp(- (!!mu) * ( exp( (exp(f)) * (!!yleft) ) - 1 ) )  ) 
  
  cdf_sig_RIGHT <- expr( 1 - exp(- (!!mu) * ( exp( (exp(f)) * (!!yright) ) - 1 ) )  ) 
  
  # derivative cdf 
  dercdf.deretasig_LEFT <- expr( ( (!!mu) * (!!yleft) * exp( (!!mu) * (- exp(exp(f) * (!!yleft))) + (!!mu) + (!!sigma) * (!!yleft) ) ) * exp(f) )
  
  dercdf.deretasig_RIGHT <- expr( ( (!!mu) * (!!yright) * exp( (!!mu) * (- exp(exp(f) * (!!yright))) + (!!mu) + (!!sigma) * (!!yright) ) ) * exp(f) )
  
  
  l_sigma <- list(pdf = pdf_sig,
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
  
  name_mu <- "Gompertz distribution: mu(log link)"
  name_sigma <- "Gompertz distribution: sigma(log link)"
  
  l_names <- list(mu = name_mu,
                  sigma = name_sigma)
  
  ### check y function   
  
  
  check_y_mu <- expr({
    if (!is.numeric(!!y))
      stop("response is not numeric but ", sQuote("GompertzMar()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("GompertzMar()"))
    y
  })
  
  check_y_sigma <- expr({
    if (!is.numeric(!!y))
      stop("response is not numeric but ", sQuote("GompertzMar()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("GompertzMar()"))
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
           marg_name = "GompertzMarg")
  
  
  return(l)
  
}


# mar1 <- LogLogistic_Mar(loc = 1)
