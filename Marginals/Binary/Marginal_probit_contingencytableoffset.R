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



BernoulliProbit_Mar <- function(loc = NULL, offset_mu = NULL){
  
  
  # check for appropriate offset values for parameter mu
  if ((!is.null(offset_mu) && (offset_mu < 0 || offset_mu > 1) ))
    stop(sQuote("mu"), paste(" must be within the interval (0,1) in Marginal", loc))
  
  
  
  # creating the location specific expressions
  y <- parse_expr(paste("y[,", loc, "]", sep = ""))                      
  mu <- parse_expr(paste("mu", loc, sep = ""))                      
  
  
  #### ADD the odds ratio as an offset: 
  if(loc == 1){
    
    y_other <- parse_expr(paste("y[,", 2, "]", sep = ""))    
    
  }else{
    
    y_other <- parse_expr(paste("y[,", 1, "]", sep = ""))    
    
  }
  
  
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
  cdf_mu  <- expr( pnorm(f) ) 
  # derivative cdf
  dercdf.deretamu  <- expr( dnorm(f) )   
  
  # Mixed responses 
  cdf_mu_mixresp  <- expr( 1 - pnorm(f) ) 
  # derivative cdf
  dercdf.deretamu_mixresp  <- expr( -dnorm(-f) )
  
  
  l_mu <- list(cdf = cdf_mu,
               dercdf = dercdf.deretamu,
               cdf_mixresp = cdf_mu_mixresp,
               dercdf_mixresp = dercdf.deretamu_mixresp)
  
  
  ### response functions 
  
  response_mu <- function(f) pnorm(f)
  
  l_response <- list(mu = response_mu)
  
  ### offset functions     
  
  # offset_mu <- expr({
  #   if (!is.null(!!mu)) {
  #     RET <- qnorm(!!mu)
  #   }else{
  #     #RET <- weighted.mean((!!y + weighted.mean(!!y, w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
  #     #RET <- weighted.mean(!!y, w = w, na.rm = TRUE)
  #     RET <- 0 
  #   }
  #   return(RET)
  # })
  # ADDED A SELECTOR OF OFFSET FUNCTIONS: this will overwrite the previous offset. More thoughts for BIVARIATE MIXED BINARY CONTINUOUS!
  if(loc == 1){
    
    
    offset_mu <- expr({
      if (!is.null(!!mu)) {
        RET <- qnorm(!!mu)
      }
      else {
        
        # Weighted mean, based on gamlss.distr
        #RET <- qlogis( weighted.mean( ( !!y + 0.5 )/2, w = w, na.rm = TRUE) ) # weighted version
        
        # CONTINGENCY TABLE OFFSET
        a <- sum(rep(1, length(which(y == 0)[which(y==0) %in% which(!!y_other==0)])) * w[which(y == 0)[which(y ==0) %in% which(!!y_other==0)]])
        b <- sum(rep(1, length(which(y == 0)[which(y==0) %in% which(!!y_other!=0)])) * w[which(y == 0)[which(y ==0) %in% which(!!y_other!=0)]])
        c <- sum(rep(1, length(which(y != 0)[which(y!=0) %in% which(!!y_other==0)])) * w[which(y != 0)[which(y !=0) %in% which(!!y_other==0)]])
        d <- sum(rep(1, length(which(y != 0)[which(y!=0) %in% which(!!y_other!=0)])) * w[which(y != 0)[which(y !=0) %in% which(!!y_other!=0)]])
        
        temp1 <- (c+d)/sum(a+b+c+d)
        RET <-  qnorm(temp1)    # log(temp1/(1-temp1))
        # Warning! Offset at zero is detrimental to the oobag approach to optimise mstop!
        #RET <- 0
      }
      return(RET)
    })
    
  }else{
    
    offset_mu <- expr({
      if (!is.null(!!mu)) {
        RET <- qnorm(!!mu)
      }
      else {
        
        # Weighted mean, based on gamlss.distr
        #RET <- qlogis( weighted.mean( ( !!y + 0.5 )/2, w = w, na.rm = TRUE) ) # weighted version
        
        # CONTINGENCY TABLE OFFSET
        a <- sum(rep(1, length(which(!!y_other == 0)[which(!!y_other==0) %in% which(y==0)])) * w[which(!!y_other == 0)[which(!!y_other==0) %in% which(y==0)]])
        b <- sum(rep(1, length(which(!!y_other == 0)[which(!!y_other==0) %in% which(y!=0)])) * w[which(!!y_other == 0)[which(!!y_other==0) %in% which(y!=0)]])
        c <- sum(rep(1, length(which(!!y_other != 0)[which(!!y_other!=0) %in% which(y==0)])) * w[which(!!y_other != 0)[which(!!y_other!=0) %in% which(y==0)]])
        d <- sum(rep(1, length(which(!!y_other != 0)[which(!!y_other!=0) %in% which(y!=0)])) * w[which(!!y_other != 0)[which(!!y_other!=0) %in% which(y!=0)]])
        
        temp2 <- (b+d)/sum(a+b+c+d)
        RET <- qnorm(temp2)      # log(temp2/(1-temp2))
        # Warning! Offset at zero is detrimental to the oobag approach to optimise mstop!
        #RET <- 0
      }
      return(RET)
    })
    
    
  }
  
  
  offset_mu_mixresp <-  expr({
    
    if(!is.null(!!mu)) {
      
      RET <- qlogis(!!mu)
      
    }
    else {
      
      # Weighted mean, based on gamlss.distr
      RET <- qlogis( weighted.mean( ( !!y + 0.5 )/2, w = w, na.rm = TRUE) ) # weighted version
      
    }
    
    return(RET)
    
  })
  
  
  l_offset <- list(mu = offset_mu)
  
  l_offset_mixresp <- list(mu = offset_mu_mixresp)
  
  ### names 
  
  name_mu <- "Probit link function for binary responses: mu(probit link)"
  
  
  l_names <- list(mu = name_mu)
  
  ### check y function   
  
  #check_y <- function(y) y                                                      
  check_y_mu <- expr({
    if (!is.numeric(!!y))
      #stop("response is not but ", sQuote("LogNormLSS()"))
      stop("response is not but ", sQuote("ProbitMargin()"))
    if (any(!!y < 0))
      #stop("response is not positive but ", sQuote("LogNormLSS()"))
      stop("response is not positive but ", sQuote("ProbitMargin()"))
    if (any(!!y > 1))
      stop("response is above 1 but ", sQuote("ProbitMargin()"))
    if ( any(unique(!!y) %in% c(0,1) == FALSE) )
      stop("response is not binary but ", sQuote("ProbitMargin()"))
    y
  })
  
  
  l_check_y <- list(mu = check_y_mu)
  
  
  # return list
  
  l = list(parameters = c(mu),
           parameter_names = c("mu"),
           generic = l_generic,
           mu = l_mu,
           response = l_response,
           offset = l_offset,
           offset_mixresp = l_offset_mixresp,
           name = l_names,
           check_y = l_check_y,
           marg_name = "ProbitMarg")
  
  
  return(l)
  
}

# marg <- LogNormal_Mar(loc = 1)
