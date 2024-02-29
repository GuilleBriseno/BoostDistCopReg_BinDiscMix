################################################################################
### DESCRIPTION
### This file defines a gamboostLSS Families object for the Log-Logistic distribution
### Ready for use in univariate gamboostLSS fitting.


library(mboost)
library(gamboostLSS)
library(gamlss)
library(gamlss.dist)
library(VGAM)



### helpers

check_stabilization <- function(stabilization = c("none", "MAD", "L2")) {
  stabilization <- match.arg(stabilization)
  ## check if old stabilization interface is used and issue a warning
  if (getOption("gamboostLSS_stab_ngrad")) {
    warning("Usage of ", sQuote("options(gamboostLSS_stab_ngrad = TRUE)"),
            " is deprecated.\n", "Use argument ", sQuote("stabilization"),
            " in the fitting family. See ?Families for details.")
    if (stabilization == "none")
      warning(sQuote("stabilization"), " is set to ", dQuote("MAD"))
  }
  stabilization
}

stabilize_ngradient <- function(ngr, w = 1, stabilization) {
  ## set which to MAD if gamboostLSS_stab_ngrad = TRUE and which == "none"
  if (stabilization == "none" && getOption("gamboostLSS_stab_ngrad"))
    stabilization <- "MAD"
  ## stabilization using the mean absolute deviation (MAD)
  if (stabilization == "MAD") {
    div <- weighted.median(abs(ngr - weighted.median(ngr, w = w, na.rm = TRUE)),
                           w = w, na.rm = TRUE)
    div <- ifelse(div < 0.0001, 0.0001, div)
    ngr <- ngr / div
  }
  if (stabilization == "L2") {
    div <- sqrt(weighted.mean(ngr^2, w =w,  na.rm = TRUE))
    div <- ifelse(div < 1e-04, 1e-04, div)
    div <- ifelse(div > 1e+04, 1e+04, div)
    ngr <- ngr / div
  }
  ngr
}

# weighted.median <- gamboostLSS::weighted.median

## weighted median
weighted.median <- function (x, w = 1, na.rm = FALSE) {
  if (length(w) == 1)
    w <- rep(w, length(x))
  
  ## remove observations with zero weights
  x <- x[w != 0]
  w <- w[w != 0]
  
  ## remove NAs if na.rm = TRUE
  if (na.rm) {
    keep <- !is.na(x) & !is.na(w)
    x <- x[keep]
    w <- w[keep]
  } else {
    if (any(is.na(x)) | any(is.na(w)))
      return(NA)
  }
  
  ## sort data and weights
  ind <- order(x)
  x <- x[ind]
  w <- w[ind]
  
  ## first time that fraction of weights is above 0.5
  ind1 <- min(which(cumsum(w)/sum(w) > 0.5))
  
  ## first time that fraction of weights is below 0.5
  ind2 <- ifelse(ind1 == 1, 1, max(which(cumsum(w)/sum(w) <= 0.5)))
  
  ## if sum of weights is an even integer
  if(sum(w) %% 1 == 0 && sum(w) %% 2 == 0)
    return(mean(c(x[ind1], x[ind2])))
  
  ## else return
  return(max(c(x[ind1], x[ind2])))
}







######### mu-model ##########

LogLogisticMu <- function(mu = NULL, sigma = NULL, stabilization){
  
  loss <- function(y, f, sigma){
    -(log(sigma) + sigma * (log(y) - log(exp(f))) - log(y) - 2*(log(1 + (y/exp(f))^sigma)))
  }
  
  
  risk <- function(y, f, w = 1){
    sum(w * loss(y = y, f = f, sigma = sigma))
  }
  
  
  ngradient <- function(y, f, w = 1){
    
    ngr <- (-sigma/exp(f) + (2 * sigma * (y/exp(f))^(sigma))/(exp(f) * (1+(y/exp(f))^sigma))) * exp(f)
    ngr <- stabilize_ngradient(ngr, w = w, stabilization)
    ngr
    
  }
  
  
  offset <- function(y, w = 1){
    if (!is.null(mu)) {
      RET <- log(mu)
    }
    else {
      mu <- (y + mean(y))/2
      RET <- log(mean(mu))
      # mu <-  rep(1, length(y))
      # RET <- log(mean(mu))
    }
    return(RET)
  }
  
  
  check_y <- function(y) {
    if (!is.numeric(y) || !is.null(dim(y)))
      stop("response is not a numeric vector but ", sQuote("LogLogisticLSS()"))
    if (any(y < 0))
      stop("response is not positive but ", sQuote("LogLogisticLSS()"))
    y
  }
  
  
  response <- function(f){
    exp(f)
  }
  
  
  mboost::Family(ngradient = ngradient, 
                 loss = loss, 
                 risk = risk,
                 offset = offset,
                 response = response,
                 check_y = check_y,
                 name = "LogLogistic: mu(log link)")
  
}


######## Sigma-model ########## 


LogLogisticSigma <- function(mu = NULL, sigma = NULL, stabilization){
  
  
  
  loss <- function(y, mu, f){
    -(log(exp(f)) + exp(f) * (log(y) - log(mu)) - log(y) - 2*(log(1 + (y/mu)^exp(f))))
  }
  
  
  
  risk <- function(y, f, w = 1){
    sum(w * loss(y = y, mu = mu, f = f))
  }
  
  
  
  
  ngradient <- function(y, f, w = 1){
    ngr <- (1/exp(f) + log(y) - log(mu) - 2/(1 + (y/mu)^exp(f)) * log(y/mu) * (y/mu)^exp(f)) * exp(f)
    ngr <- stabilize_ngradient(ngr, w = w, stabilization)
    ngr
    
  }
  
  
  
  
  offset <- function (y, w = 1) {
    if (!is.null(sigma)) {
      RET <- log(sigma)
    }
    else {
      sigma <- rep(0.1, length(y))
      RET <- log(mean(sigma))
    }
    return(RET)
  }
  
  
  
  check_y <- function(y) {
    if (!is.numeric(y) || !is.null(dim(y)))
      stop("response is not a numeric vector but ", sQuote("LOgLogisticLSS()"))
    if (any(y < 0))
      stop("response is not positive but ", sQuote("LogLogisticLSS()"))
    y
  }
  
  
  response <- function(f){
    exp(f)
  }
  
  
  mboost::Family(ngradient = ngradient,
                 loss = loss,
                 risk = risk,
                 offset = offset,
                 check_y = check_y,
                 response = response,
                 name = "LogLogistic distribution: sigma(log link)")
  
}



######## creation of the Log-Norm Families object #########


LogLogisticLSS <- function(mu = NULL, sigma = NULL, 
                           stabilization = c("none", "MAD", "L2")){
  
  if ((!is.null(mu) && mu < 0))
    stop(sQuote("mu"), " must be greater than zero")
  
  if ((!is.null(sigma) && sigma < 0))
    stop(sQuote("sigma"), " must be greater than zero")
  
  
  stabilization <- check_stabilization(stabilization)
  
  
  gamboostLSS::Families(mu = LogLogisticMu(mu = mu, sigma = sigma, stabilization = stabilization),
                        sigma = LogLogisticSigma(mu = mu, sigma = sigma, stabilization = stabilization),
                        name = "LogLogistic")
  
  
}

