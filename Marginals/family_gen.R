################################################################################
### DESCRIPTION
### This is a helper function that creates for the different
### parameter sub-models, i.e. the parameters of the marginals and the 
### copula dependence parameter, the respective loss, risk and gradient functions.
### Finally, it creates and returns a parameter specific mboost Family objects .


### libraries 
library(numDeriv)



family_gen <- function(mu1 = NULL, sigma1 = NULL, nu1 = NULL, tau1 = NULL,
                       mu2 = NULL, sigma2 = NULL, nu2 = NULL, tau2 = NULL, 
                       rho = NULL,
                       stabilization,
                       loss_body, loss_args, 
                       risk_body, risk_args,
                       grad_body, grad_args,
                       offset_body, offset_args,
                       check_y_body, check_y_args,
                       response, name){
  
  
  
################################################################################
######################### Helpers from gamboostLSS #############################
################### To be deleted when integration into package ################
  
  ## function for weighted sd
  weighted.sd <- function(x, w, ...) {
    if (missing(w))
      w <- rep(1, length(x))
    m <- weighted.mean(x, w, ...)
    var <- weighted.mean((x - m)^2, w, ...) * sum(w) / (sum(w) - 1)
    return(sqrt(var))
  }
  
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
  
  ## helper function that stabilizes the negative gradient if requested by the user
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
  
  
  ### pdf Singh-Maddala
  
  dSinghMaddala <- function (x, mu = 1, sigma = 1, tau = 1, log = FALSE){
    
    if (any(mu < 0)) 
      stop(paste("mu must be positive", "\n", ""))
    
    if (any(tau < 0)) 
      stop(paste("tau must be positive", "\n", 
                 ""))
    
    z <- (x/mu)^sigma
    
    loglik <- log(z) + log(abs(sigma)) - log(x) + log(tau) - (1 + tau) * log(1 + z)
    
    if (log == FALSE) 
      ft <- exp(loglik)
    else ft <- loglik
    ft
  }
  
  ### cdf Singh-Maddala
  
  pSinghMaddala <- function (x, mu = 1, sigma = 1, tau = 1, log = FALSE){
    
    if (any(mu < 0)) 
      stop(paste("mu must be positive", "\n", ""))
    
    if (any(tau < 0)) 
      stop(paste("tau must be positive", "\n", 
                 ""))
    
    z <- (x/mu)^sigma
    
    log_cdf <- log(1 - (1 + z)^(-tau))
    
    if (log == FALSE) 
      ft <- exp(log_cdf)
    else ft <- log_cdf
    ft
  }
  
  # pdf Dagum
  
  dDagum <- function (x, mu = 1, sigma = 1, nu = 1, log = FALSE){
    
    if (any(mu < 0)) 
      stop(paste("mu must be positive", "\n", ""))
    
    if (any(nu < 0)) 
      stop(paste("nu must be positive", "\n", ""))
    
    z <- (x/mu)^sigma
    loglik <- nu * log(z) + log(abs(sigma)) - log(x) - lgamma(nu) - 
      lgamma(1) + lgamma(nu + 1) - (nu + 1) * log(1 + z)
    
    if (log == FALSE) 
      ft <- exp(loglik)
    else ft <- loglik
    
    ft
  }
  
  ## cdf Dagum
  
  pDagum <- function (x, mu = 1, sigma = 1, nu = 1, log = FALSE){
    
    if (any(mu < 0))
      stop(paste("mu must be positive", "\n", ""))
    
    if (any(nu < 0))
      stop(paste("nu must be positive", "\n",
                 ""))
    
    log_cdf <- log((1 + (x/mu)^-sigma)^-nu)
    
    if (log == FALSE)
      ft <- exp(log_cdf)
    else ft <- log_cdf
    ft
  }
  
  
  ### pdf LogLogistic
  
  dLogLogistic <- function (x, mu = 1, sigma = 1, log = FALSE){
    
    if (any(mu < 0)) 
      stop(paste("mu must be positive", "\n", ""))
    
    if (any(sigma < 0)) 
      stop(paste("sigma must be positive", "\n", ""))
    
    
    z <- (x/mu)^(sigma)
    
    loglik <- log(sigma) + (sigma) * (log(x) - log(mu)) - log(x) - 2 * log(1 + z)
    
    if (log == FALSE) 
      ft <- exp(loglik)
    else ft <- loglik
    ft
  }
  
  ### cdf LogLogistic
  
  pLogLogistic <- function (x, mu = 1, sigma = 1, log = FALSE){
    
    if (any(mu < 0)) 
      stop(paste("mu must be positive", "\n", ""))
    
    if (any(sigma < 0)) 
      stop(paste("sigma must be positive", "\n", ""))
    
    
    log_cdf <- log(1/(1 + (x/mu)^(-sigma)))
    
    if (log == FALSE) 
      ft <- exp(log_cdf)
    else ft <- log_cdf
    ft
  }
  
  weighted.sd <- function(x, w, ...) {
    if (missing(w))
      w <- rep(1, length(x))
    m <- weighted.mean(x, w, ...)
    var <- weighted.mean((x - m)^2, w, ...) * sum(w) / (sum(w) - 1)
    return(sqrt(var))
  }
  
  
  
################################################################################
######## Helper for numerical approximation of the gradient from gjrm ##########
################################################################################
  
  # helper for numerical approximation from gjrm
  
  num_grad <- function(funcD, para){
    
    para <- c(para)
    
    # if(length(para) == 1){
    #   
    #   fi <- jacobian(function(para) funcD(para), para)
    #   
    # }
    # 
    # if(length(para) > 1){
      
      fi <- grad(function(para) funcD(para), para)
      
    # }
    
    return(fi)
    
    
  }
  
  
  derFunc.derparameter1 <- function(func, y2, param1, paramconstant){ num_grad(func, param1) }
  derFunc.derparameter2 <- function(func, y2, paramconstat, param2){ num_grad(func, param2) }  
  
  
  GetNumDeriv_DiscCDFOLD_DONOTUSE <- function(resp, param1, param2, funcName, which_param){
    
    if(which_param == 1){
      
      numder <- derFunc.derparameter1(func = function(param1) funcName(q = resp, mu = param1, sigma = param2), 
                                      y2 = resp, 
                                      param1 = param1, 
                                      paramconstant = param2)
      
    }
    
    if(which_param == 2){
      
      numder <- derFunc.derparameter2(func = function(param2) funcName(q = resp, mu = param1, sigma = param2),
                                      y2 = resp,
                                      paramconstat = param1, 
                                      param2 = param2)
      
    }
    
    return(numder)
    
  }

  GetNumDeriv_DiscCDFP1 <- function(resp, quantity1, quantity2, funcName, derParamderEta){
    
    par1 <- quantity1
    
    numgd <- num_grad(function(par1) funcName(resp, mu = par1, sigma = quantity2), para = par1)
    
    numgdeta <- numgd * derParamderEta
    
    return(numgdeta)
  }
  
  GetNumDeriv_DiscCDFP2 <- function(resp, quantity1, quantity2, funcName, derParamderEta){
    
    par2 <- quantity2
    
    numgd <- num_grad(function(par2) funcName(resp, mu = quantity1, sigma = par2), para = par2)
    
    numgdeta <- numgd * derParamderEta
    
    return(numgdeta)
  }
  
    
################################################################################
######## Helper for computation of gradient of discrete and other ugly CDFs w.r.t parameter ###
################################################################################
  
  
  # Should be helpful with some tight differences of probabilities:
  # min.pr = 1e-16, max.pr = 0.999999
  
  pdffz <- function(input){
    
    pdiff <- ifelse(input > 0.999999, 0.999999, input) 
    pdiff <- ifelse(pdiff < 1e-16, 1e-16, pdiff) 
    
    return(pdiff)
    
  }
  
  
  dvffz <- function(derivative){
    
    deriv <- ifelse(is.na(derivative), .Machine$double.eps, derivative)
    deriv <- ifelse(deriv == Inf, 8.218407e+20, deriv)
    deriv <- ifelse(deriv == -Inf, -8.218407e+20, deriv)
    
    return(deriv)
  }
  
  
  check_eta_theta <- function(additivepredictor, coptype){
    
    if(coptype == "GAUSS"){
    
      checked_predictor <- ifelse(abs(additivepredictor) > 8.75, sign(additivepredictor) * 8.75, additivepredictor)
    
    }
    
    if(coptype == "FRANK"){
      
      epsilon <- c(sqrt(.Machine$double.eps))
      
      checked_predictor <- ifelse(abs(additivepredictor) < epsilon, epsilon, additivepredictor)
      
    }
    
    if(coptype %in% c("GUMBEL", "CLAYTON", "JOE")){
      
      checked_predictor <- ifelse(additivepredictor > 20, 20, additivepredictor)
      
      checked_predictor <- ifelse(additivepredictor < -17, -17, additivepredictor)
    }
    
    return(checked_predictor)
  }
  
  
  # GBS: Credit where it is due: Another incredibly useful piece of software from
  # gjrm:
  # 
  aux_matrix_disc <- function(response){

    ly1 <- length(response)
    y1m <- list()
    my1 <- max(response)

    for(i in 1:ly1){

      y1m[[i]] <- seq(0, response[i]);

      length(y1m[[i]]) <- my1 + 1

    }


    y1m <- do.call(rbind, y1m)


    return(y1m)
  }
  
  
  
  ###### 
  derDISCDIST_PDF.derparam2 <- function(auxmatrix, dim_aux_matrix, func, resp, par1, par2){
    
    ### Check first if we are at iteration 1 (i.e. there is only one value for mu and sigma due to them being offsets...)
    if(length(par1) == 1){
      
      param1 <- rep(par1, length(resp))
      
    }else{ param1 <- par1 }
    
    if(length(par2) == 1){
      
      param2 <- rep(par2, length(resp))
      
    }else{ param2 <- par2 }
    
    # Change the stuff below to use the matrix, so hopefully it will be faster now:
    #CDF_der <- rep(0, length(resp))
    # for(i in 1:length(resp)){
    #   allval <- seq(0, resp[i])
    #   CDF_der[i] <- sum(func(allval, param1 = param1[i], param2 = param2[i]))
    # }
    
    CDF_der <- rowSums(t(sapply(1:(dim_aux_matrix[1]), function(i) func(auxmatrix[i,], param1[i], param2[i]))), na.rm = TRUE) 
    
    return(CDF_der)
  }
  
  derDISCDIST_PDF.derparam2_3PARAMS <- function(auxmatrix, dim_aux_matrix, func, resp, par1, par2, par3){
    
    ### Check first if we are at iteration 1 (i.e. there is only one value for mu and sigma due to them being offsets...)
    if(length(par1) == 1){
      
      param1 <- rep(par1, length(resp))
      
    }else{ param1 <- par1 }
    
    if(length(par2) == 1){
      
      param2 <- rep(par2, length(resp))
      
    }else{ param2 <- par2 }

    if(length(par3) == 1){
      
      param3 <- rep(par3, length(resp))
      
    }else{ param3 <- par3 }    
    # Change the stuff below to use the matrix, so hopefully it will be faster now:
    #CDF_der <- rep(0, length(resp))
    # for(i in 1:length(resp)){
    #   allval <- seq(0, resp[i])
    #   CDF_der[i] <- sum(func(allval, param1 = param1[i], param2 = param2[i]))
    # }
    
    CDF_der <- rowSums(t(sapply(1:(dim_aux_matrix[1]), function(i) func(auxmatrix[i,], param1[i], param2[i], param3[i]))), na.rm = TRUE) 
    
    return(CDF_der)
  }
  
  # # Discrete PDFs as functions (useful for discrete marginal distributions!)
  # # Poisson PDF
  derPOISSON_PDF.dermu2 <- function(resp, param1){ ( exp(-param1) * (param1^(resp - 1) * resp - param1^resp)/factorial(resp) ) * param1 }
  #derPOISSON_PDF.dermu2 <- function(resp, param1){  ( ( exp( -param1 ) * ( param1^(resp - 1 ) * (resp) - param1^(resp) ) ) / factorial(resp) )  }
  # derPOISSON_PDF.dermu2 <- function(resp, param1){
  #   
  #   gamma <- exp( -param1 ) * param1^(resp)/( dpois(resp, param1))
  #   
  #   deriv <- as.vector(1/gamma*(-exp(-param1)*param1^(resp) + exp(-param1)*resp*param1^(resp-1))*param1) 
  #   
  #   return(deriv)
  # }
  
  
  # Geometric PDF
  derGEOMETRIC_PDF.dermu2 <- function(resp, param){ (resp - param)* param^(resp-1) * (param + 1)^(-resp -2) }
  
  # ZERO INFLATED POISSON PDF
  derZIP_PDF.dermu2_OLD <- function(resp, param1, param2){
    
    
    derivpdf2 <- ifelse( resp > 0, 
                         (1 - param2) * (param1^((resp) - 1) * (resp))/gamma(resp + 1) * exp(-param1) - (1 - param2) * param1^(resp)/gamma(resp + 1) * exp(-param1), 
                         -( (1 - param2) * exp(-param1) )
                         )
    
    
    derivpdf2 <- as.vector(derivpdf2) * param1
    
    
    return(derivpdf2)
  }
  
  derZIP_PDF.dersigma2_OLD <- function(resp, param1, param2){
    
    derivpdf2 <- ifelse(resp > 0, 
                        -(param1^(resp)/gamma(resp + 1) * exp(-param1)), 
                        1 - exp(-param1)) 
    
    
    derivpdf2 <- as.vector(derivpdf2) * ( (param2/(1-param2)*(1+param2/(1-param2))-(param2/(1-param2))^2)/(1+param2/(1-param2))^2 )
    
    
    return(derivpdf2)
  }
  
  
  
  #### Zero inflated poisson 
  derZIP_PDF.derparams <- function(resp, param1, param2, which_param, zeroIndices){
    
    derZIP <- vector("numeric", length = length(resp))
    
    if(which_param == 1){
      

      derZIP[zeroIndices] <- -( (1 - param2[zeroIndices]) * exp(-param1[zeroIndices]) )
      derZIP[-zeroIndices] <-  (1 - param2[-zeroIndices]) * (param1[-zeroIndices]^((resp[-zeroIndices]) - 1) * (resp[-zeroIndices]))/gamma(resp[-zeroIndices] + 1) * exp(-param1[-zeroIndices]) - (1 - param2[-zeroIndices]) * param1[-zeroIndices]^(resp[-zeroIndices])/gamma(resp[-zeroIndices] + 1) * exp(-param1[-zeroIndices]) 
      #(1 - param2[-zeroIndices]) * ( exp(-param1[-zeroIndices]) * (param1[-zeroIndices]^(resp[-zeroIndices] - 1) * resp[-zeroIndices] - param1[-zeroIndices]^resp[-zeroIndices])/factorial(resp[-zeroIndices]) ) 
      
      derParamderEta <- param1
      
      derZIP <- derZIP * derParamderEta
      
    }
    
    if(which_param == 2){
      
      derZIP[zeroIndices] <- ( 1 - exp(-param1[zeroIndices]) )
      derZIP[-zeroIndices] <-  -(param1[-zeroIndices]^(resp[-zeroIndices])/gamma(resp[-zeroIndices] + 1) * exp(-param1[-zeroIndices])) #- dpois(x = resp[-zeroIndices], lambda = param1[-zeroIndices])
      
      derParamderEta <- ( (param2/(1-param2)*(1+param2/(1-param2))-(param2/(1-param2))^2)/(1+param2/(1-param2))^2 )
      
      derZIP <- derZIP * derParamderEta
      
    }
    
    return(derZIP)
    
  }
  getDiscCDFDerivs <- function(resp, param, aux_dimensions){
    
    derivs <- rowSums(matrix(as.numeric(derPOISSON_PDF.dermu2(resp, param) ), 
                             aux_dimensions[1],
                             aux_dimensions[2]), na.rm = TRUE)
    
    return(derivs)
  }
  
  
  
  
  ##### NBI
  derNBI_PDF.dermu2 <- function(resp, param1, param2){
    
    ( gamma( resp + (1/(param2)) )/ ( (gamma( 1/(param2) ) )*( gamma( resp + 1) )) ) * ( ((param2) * (param1)) / (1 + (param2)*(param1)) )^(resp) * (1/(1 + (param2)*(param1) ))^(1/(param2)) * 
      ( resp*(param1*param2)^(-1)*param2*param1-(resp+1/param2)*(param1*param2+1)^(-1)*param2*param1 ) 
    
  }
  
  derNBI_PDF.dersigma2 <- function(resp, param1, param2){
    
    ( gamma( resp + (1/(param2)) )/ ( (gamma( 1/(param2) ) )*( gamma( resp + 1) )) ) * ( ((param2) * (param1)) / (1 + (param2)*(param1)) )^(resp) * (1/(1 + (param2)*(param1) ))^(1/(param2)) * 
      ( digamma(resp+1/param2)*(-1/param2^2)+resp*(param1*param2)^(-1)*param1-(digamma(1/param2)*(-1/param2^2)+(-1/param2^2)*log(param1*param2+1)+(resp+1/param2)*(1/(param1*param2+1))*param1) ) * 
      param2 
    
  }
  
  
  ##### ZALG
  derZALG_PDF.dermu2 <- function(resp, param1, param2){
    
    
    derivpdf2 <- ifelse( resp > 0, 
                         (param1^(-1 + resp) * (-1 + param2) * (-param1 + resp * (-1 + param1) * log(1 - param1)))/(resp * (-1 + param1) * log(1 - param1)^2), 
                         0)
    
    
    derivpdf2 <- as.vector(derivpdf2) * (1 - param1) * param1
    
    
    return(derivpdf2)
  }
  
  derZALG_PDF.dersigma2 <- function(resp, param1, param2){
    
    
    derivpdf2 <- ifelse(resp > 0,
                        -((-(log(1 - param1))^(-1)) * param1^(resp)/resp),
                        1)
    
    derivpdf2 <- as.vector(derivpdf2) * ( (param2/(1 - param2)*(1 + param2/(1 - param2))-(param2/(1 - param2))^2) / (1 + param2/(1 - param2))^2 )  
    
    return(derivpdf2)
  }
  
  ##### ZINBI
  derZINBI_PDF.dermu2 <- function(resp, param1, param2, param3){
    
    
    derivpdf2 <- ifelse( resp > 0, 
                         #
                         -(((-1 + param3) * (1/(1 + param2 * param1))^(1/param2) * ((param1 * param2/(1 + param1 * param2))^(resp - 1) * (resp * (param2/(1 + param1 * param2) - 
                                                                                                                                                    param1 * param2 * param2/(1 + param1 * param2)^2))) - (-1 + param3) * ((1/(1 + param2 * param1))^((1/param2) - 1) * ((1/param2) * (param2/(1 + 
                                                                                                                                                                                                                                                                                                 param2 * param1)^2))) * (param1 * param2/(1 + param1 * param2))^resp) * gamma(1/param2 + resp)/(gamma(1/param2) * gamma(1 + resp))), 
                         #
                         (-1 + param3) * ((1/(1 + param1 * param2))^((1/param2) - 1) * ((1/param2) * (param2/(1 + param1 * param2)^2))) )
    
    
    derivpdf2 <- as.vector(derivpdf2) #* param1
    
    return(derivpdf2)
    
  }
  
  derZINBI_PDF.dersigma2 <- function(resp, param1, param2, param3){
    
    derivpdf2 <- ifelse( resp > 0, 
                         #
                         ((-1 + param3) * (1/(1 + param1 * param2))^(1 + 1/param2) * ((param1 * param2)/(1 + param1 * param2))^
                            resp * gamma(resp + 1/param2) * (-resp * param2 + param1 * param2 + log(1/(1 + param1 * param2)) + 
                                                               param1 * param2 * log(1/(1 + param1 * param2)) + (1 + param1 * param2) * digamma(
                                                                 resp + 1/param2) - (1 + param1 * param2) * digamma(1/param2)))/(param2^2 * gamma(1 + resp) * gamma(1/param2)), 
                         #
                         ((-1 + param3) * (1/(1 + param1 * param2))^(1 + 1/param2) * (param1 * param2 + (1 + param1 * param2) * log(1/(1 + param1 * param2))))/param2^2)
    
    derivpdf2 <- as.vector(derivpdf2) #* param2
    
    return(derivpdf2)
  }
  
  derZINBI_PDF.dernu2 <- function(resp, param1, param2, param3){
    
    derivpdf2  <- ifelse(resp > 0, 
                         #
                         -((1/(1 + param2 * param1))^(1/param2) * (param1 * param2/(1 + param1 * param2))^resp * gamma(1/param2 + resp)/(gamma(1/param2) * gamma(1 + resp))), 
                         #
                         1 - (1/(1 + param1 * param2))^(1/param2)) 
    
    derivpdf2 <- as.vector(derivpdf2) # * ( (param3/(1-param3)*(1+param3/(1-param3))-(param3/(1-param3))^2)/(1+param3/(1-param3))^2 )   
    
    return(derivpdf2)
  }
  
  
  ##### ZANBI
  derZANBI_PDF.dermu2 <- function(resp, param1, param2, param3){
    
    derivpdf2 <- ifelse(resp > 0, 
                        ((-1 + param3) * (1/(1 + param1 * param2))^(1 + 1/param2) * ((param1 * param2)/(1 + param1 * param2))^resp * (param1 + 
                                        resp * (-1 + (1/(1 + param1 * param2))^(1/param2))) * gamma(resp + 1/param2))/(param1 * (-1 + (1/(1 + param1 * param2))^(1/param2))^2 * gamma(1 + resp) * gamma(1/param2)),
                        0)
    
    
    derivpdf2 <- as.vector(derivpdf2) * param1
    
    return(derivpdf2)
  }
  
  derZANBI_PDF.dersigma2 <- function(resp, param1, param2, param3){
    
    
    derivpdf2 <-  ifelse(resp > 0, 
                         -(((-1 + param3) * (1/(1 + param1 * param2))^(1 + 1/param2) * ((param1 * param2)/(1 + param1 * param2))^resp * 
                                    gamma(resp + 1/param2) * (resp * param2 - param1 * param2 - 
                                                             resp * param2 * (1/(1 + param1 * param2))^(1/param2) - log(1/(1 + param1 * param2)) -
                                                             param1 * param2 * log(1/(1 + param1 * param2)) + (1 + param1 * param2) * (-1 + (1/(1 + param1 * param2))^(1/
                                                            param2)) * digamma(resp + 1/param2) - (1 + param1 * param2) * (-1 + (1/(1 + param1 * param2))^(1/param2)) * digamma(1/param2)))/(param2^2 * (-1 + (1/(1 + param1 * param2))^(1/param2))^2 * gamma(1 + resp) * 
                                                                                                                                                                                                                                                                        gamma(1/param2))), 
                         0)
    
    
    derivpdf2 <- as.vector(derivpdf2) * param2 
    
    return(derivpdf2)
  }
  
  derZANBI_PDF.dernu2 <- function(resp, param1, param2, param3){
    
    derivpdf2 <- ifelse(resp > 0, 
                        (1/(1 + param1 * param2))^(1/param2) * (param1 * param2/(1 + param1 * param2))^resp * 
                          gamma(1/param2 + resp)/((-1 + (1/(1 + param2 * param1))^(1/param2)) * 
                                                 gamma(1/param2) * gamma(1 + resp)),
                        1)
    
    derivpdf2 <- as.vector(derivpdf2) *  ( (param3/(1-param3)*(1+param3/(1-param3))-(param3/(1-param3))^2)/(1+param3/(1-param3))^2 )
    
    return(derivpdf2)
  }
  
  
  ############ Poisson Inverse Gaussian
  derPIG_PDF.dermu2 <- function(resp, param1, param2){
    
    fd.prec <- 10^(-7)
    
    f2 <- dPIG(resp, mu = param1, sigma = param2)
    
    f2.fd.mu <- dPIG(resp, mu = (param1 + fd.prec), sigma = param2)
    
    f2.fd.mu <- ifelse(f2.fd.mu >  10^(-8), f2.fd.mu,  10^(-8))
    
    derivpdf2 <- (f2.fd.mu - f2) / (fd.prec)
    
    derivpdf2 <- as.vector(derivpdf2) * param1 
    
    return(derivpdf2)
    
  }
  
  derPIG_PDF.dersigma2 <- function(resp, param1, param2){
    
    fd.prec <- 10^(-7)
    
    f2 <- dPIG(resp, mu = param1, sigma = param2)
    
    f2.fd.sigma <- dPIG(resp, mu = param1, sigma=(param2 + fd.prec))
    
    f2.fd.sigma <- ifelse(f2.fd.sigma > 10^(-8), f2.fd.sigma, 10^(-8))
    
    derivpdf2 <- (f2.fd.sigma - f2) / (fd.prec)
    
    derivpdf2 <- as.vector(derivpdf2) * param2 
    
    return(derivpdf2)
    
  }
  
  
  ############ ZERO INFLATED Poisson Inverse Gaussian
  derZIPIG_PDF.dermu2 <- function(resp, param1, param2, param3){
    
    fd.prec <- 10^(-7)
    
    f2 <- dZIPIG(resp, mu = param1, sigma = param2, nu = param3)
    
    f2.fd.mu <- dZIPIG(resp, mu = (param1 + fd.prec), sigma = param2, nu = param3)
    
    f2.fd.mu <- ifelse(f2.fd.mu >  10^(-8), f2.fd.mu,  10^(-8))
    
    derivpdf2 <- (f2.fd.mu - f2) / (fd.prec)
    
    derivpdf2 <- as.vector(derivpdf2) * param1 
    
    return(derivpdf2)
    
  }
  
  derZIPIG_PDF.dersigma2 <- function(resp, param1, param2, param3){
    
    fd.prec <- 10^(-7)
    
    f2 <- dZIPIG(resp, mu = param1, sigma = param2, nu = param3)
    
    f2.fd.sigma <- dZIPIG(resp, mu = param1, sigma = (param2 + fd.prec), nu = param3)
    
    f2.fd.sigma <- ifelse(f2.fd.sigma > 10^(-8), f2.fd.sigma, 10^(-8))
    
    derivpdf2 <- (f2.fd.sigma - f2) / (fd.prec)
    
    derivpdf2 <- as.vector(derivpdf2) * param2 
    
    return(derivpdf2)
    
  }
  
  
  derZIPIG_PDF.dernu2 <- function(resp, param1, param2, param3){
    
    fd.prec <- 10^(-7)
    
    f2 <- dZIPIG(resp, mu = param1, sigma = param2, nu = param3)
    
    f2.fd.nu <- dZIPIG(resp, mu = param1, sigma = param2, nu = (param3 + fd.prec))
    
    f2.fd.nu <- ifelse(f2.fd.nu > 10^(-8), f2.fd.nu, 10^(-8))
    
    derivpdf2 <- (f2.fd.nu - f2) / (fd.prec)
    
    derivpdf2 <- as.vector(derivpdf2) * ((param3/(1-param3)*(1+param3/(1-param3))-(param3/(1-param3))^2)/(1+param3/(1-param3))^2) 
    
    return(derivpdf2)
    
  }
  
  
  ###### Zero Inflated Poisson
  derZIP_PDF.dermu2 <- function(resp, param1, param2){
    
    derivpdf2 <-  ifelse(resp > 0, 
                         (1 - param2) * (param1^((resp) - 1) * (resp))/gamma(resp + 1) * exp(-param1) - 
                           (1 - param2) * param1^(resp)/gamma(resp + 1) * exp(-param1), 
                         -((1 - param2) * exp(-param1)))
    
    derivpdf2 <- as.vector(derivpdf2) * param1
    
  }
  
  derZIP_PDF.dersigma2 <- function(resp, param1, param2){
    
    derivpdf2 <- ifelse(resp > 0,
                        -(param1^(resp)/gamma(resp + 1) * exp(-param1)), 
                        1 - exp(-param1)) 
    
    derivpdf2 <- as.vector(derivpdf2) * (1 - param2) * param2
  }
  
  
  ### New approach: 
  # funcD should be from gamlss.dist
  getDiscFuncDerivs_P1 <- function(inputf, funcD, resp, other){
    
    input <- inputf
    
    grdnt <- grad(func  = function(input) funcD(resp, mu = input, sigma = other), x = input)
    
    return(grdnt)
  }
  
  getDiscFuncDerivs_P2 <- function(inputf, funcD, resp, other){
  
    input <- inputf
    
    grdnt <- grad(func  = function(input) funcD(resp, mu = other, sigma = input), x = input)
  
    return(grdnt)
    }
  # Something like this should work
  # rowSums(t(sapply(1:dim_aux_mat[1], function(i) derPOISSON_PDF.dermu2(aux_mat[i,], exp(f)[i]))), na.rm = TRUE)
  
  
  ##### Wrappers
  getDerivs_TWOParams <- function(func, p1, p2, resp){
    
    CDF_der <- rep(0, length(resp))
    
    #print(length(resp))
    
    for(i in 1:length(resp)){
      #y.y <- y[i]
      #mm <- nmu[i]
      #ms <- nsigma[i]
      allval <- seq(0, resp[i])
      CDF_der[i] <- sum(func(allval, p1[i], p2[i]))
      #CDF_der[i] <- sum(pdfall)
    }
    #CDF_der <- sapply(1:length(p1), function(i)  sum(func(seq(0, y[i]), p1[i], p2[i])) )
    
    
    return(CDF_der)
  }
  
  getDerivs_THREEParams <- function(func, p1, p2, p3, resp, dimension){
    
    #CDF_der <- rep(0, length(p1))
    CDF_der <- as.numeric(sapply(1:dimension, function(i)  sum(func(seq(0, resp[i]), p1[i], p2[i])), simplify = TRUE ))
    
    return(CDF_der)
  }
  
  get_derCDF_1param_etaparam <- function(func, resp, par1){
    
    # if(length(par1) == 1){
    #   
    #   param1 <- rep(par1, length(resp))
    #   
    # }else{ param1 <- par1 }
    param1 <- par1
    
    #derZALG_PDF.dermu2(seq(0, resp_i), param1, param2)
    CDF_der <- rep(0, length(resp))
    
    
    for(i in 1:length(resp)){
      #y.y <- y[i]
      #mm <- nmu[i]
      #ms <- nsigma[i]
      allval <- seq(0, resp[i])
      CDF_der[i] <- sum(func(allval, param1 = param1[i]))
      #CDF_der[i] <- sum(pdfall)
    }
    
    return(CDF_der)
    
  }
  
  get_derCDF_etaparam <- function(func, resp, par1, par2){
    
    if(length(par1) == 1){
      
      param1 <- rep(par1, length(resp))
      
    }else{ param1 <- par1 }
    
    if(length(par2) == 1){
      
      param2 <- rep(par2, length(resp))
      
    }else{ param2 <- par2 }
    
    #derZALG_PDF.dermu2(seq(0, resp_i), param1, param2)
    CDF_der <- rep(0, length(resp))
    
    
    for(i in 1:length(resp)){
      #y.y <- y[i]
      #mm <- nmu[i]
      #ms <- nsigma[i]
      allval <- seq(0, resp[i])
      CDF_der[i] <- sum(func(allval, param1 = param1[i], param2 = param2[i]))
      #CDF_der[i] <- sum(pdfall)
    }
    
    return(CDF_der)
    
  }
  
  get_derCDF_3param_etaparam <- function(func, resp, par1, par2, par3){
    
    if(length(par1) == 1){
      
      param1 <- rep(par1, length(resp))
      
    }else{ param1 <- par1 }
    
    if(length(par2) == 1){
      
      param2 <- rep(par2, length(resp))
      
    }else{ param2 <- par2 }
    
    if(length(par3) == 1){
      
      param3 <- rep(par3, length(resp))
      
    }else{ param3 <- par3 }
    
    #derZALG_PDF.dermu2(seq(0, resp_i), param1, param2)
    CDF_der <- rep(0, length(resp))
    
    
    for(i in 1:length(resp)){
      #y.y <- y[i]
      #mm <- nmu[i]
      #ms <- nsigma[i]
      allval <- seq(0, resp[i])
      CDF_der[i] <- sum(func(allval, param1 = param1[i], param2 = param2[i], param3 = param3[i]))
      #CDF_der[i] <- sum(pdfall)
    }
    
    return(CDF_der)
    
  }
  
  
  ###### GAMMA DISTRIBUTION: Derivative of the CDF w.r.t. sigma
  derGAMMA_CDF.dersigma <- function(resp, param1, param2){
    
    fd.prec <- 10^(-7)
    
    f2 <- pgamma(resp, shape = 1/param2^2, scale = param1 * param2^2) 
    
    f2.fd.sig <- pgamma(resp, shape = 1/(param2 + fd.prec)^2, scale = param1 * (param2 + fd.prec)^2) 
    
    f2.fd.sig <- ifelse(f2.fd.sig >  10^(-8), f2.fd.sig,  10^(-8))
    
    derivcdf <- (f2.fd.sig - f2) / (fd.prec)
    
    derivcdf <- as.vector(derivcdf)  
  
    return(derivcdf)
    
  }
  
  
  
  
################################################################################
######################### Function definition ##################################
################################################################################
  
  
  # defining the loss function
  loss <- function(){}
  body(loss) <- loss_body
  formals(loss) <- loss_args
  
  
  # defining the risk function
  risk <- function(){}   
  body(risk) <- risk_body
  formals(risk) <- risk_args
  
  
  # defining the gradient function
  ngradient <- function(){}
  body(ngradient) <- as.call(c(as.name("{"), grad_body)) # creating a multi-expression body (see: https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/body)
  formals(ngradient) <- grad_args
  
  # defining offset function
  offset <- function(){}
  body(offset) <- offset_body
  formals(offset) <- offset_args
  
  # defining the check_y function
  
  check_y <- function(){}
  body(check_y) <- check_y_body
  formals(check_y) <- check_y_args
  
  # setting the appropriate environment for offset, y_check and response function
  #environment(offset) <- current_env()
  environment(response) <- current_env()
  #environment(check_y) <- current_env()
  
  
  # creating the Family object
  mboost::Family(ngradient = ngradient,
                 risk = risk,
                 loss = loss,
                 response = response,
                 offset = offset,
                 name = name,
                 check_y = check_y)
  
  
}