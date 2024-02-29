#########################################
library("VGAM")
library("gamboostLSS")
library("mvtnorm")
library("pROC")
library("scoringRules")
#library("GJRM")
#setwd("BOOSTCOPFILES/")
#########################################
### load Copula


# Binary continuous case!!! (get copulas of 90 and 270 degrees rotation to work)
#source("Copulas/BinaryContinuous/Copula_Clayton_BinCont.R")
#source("Copulas/BinaryContinuous/Copula_Clayton270_BinCont.R")
#source("Copulas/BinaryContinuous/Copula_Clayton270_NEW_BinCont.R")

source("Copulas/BinaryContinuous/Copula_Clayton180_BinCont.R")
source("Copulas/BinaryContinuous/Copula_Clayton270_NEW2_BinCont.R")
source("Copulas/BinaryContinuous/Copula_Frank2_BinCont.R")


### local directories
# source("Copulas/Mixed_BinaryContinuous/Copula_Clayton270_NEW2_BinCont.R")
# source("Copulas/Mixed_BinaryContinuous/Copula_Frank2_BinCont.R")
# source("Copulas/Mixed_BinaryContinuous/Copula_Clayton180_BinCont.R")
# source("Copulas/Mixed_BinaryContinuous/Copula_Clayton90_BinCont.R")
# source("Copulas/Mixed_BinaryContinuous/Copula_Clayton_BinCont.R")

### Gumbel 
# source("Copulas/Mixed_BinaryContinuous/Copula_Gumbel_BinCont.R")
# source("Copulas/Mixed_BinaryContinuous/Copula_Gumbel90_BinCont.R")
# source("Copulas/Mixed_BinaryContinuous/Copula_Gumbel180_BinCont.R")
# source("Copulas/Mixed_BinaryContinuous/Copula_Gumbel270_BinCont.R")


# Useful trycatch variation:
myTryCatch <- function(expr){
  
  warn <- err <- NULL
  
  value <- withCallingHandlers(
    
    tryCatch(expr, error=function(e) {
      
      err <<- e
      
      NULL
    }), warning=function(w) {
      
      warn <<- w
      
      invokeRestart("muffleWarning")
    })
  
  list(value=value, warning=warn, error=err)
}

pdffz <- function(input){
  
  pdiff <- ifelse(input > 0.999999, 0.999999, input) 
  pdiff <- ifelse(pdiff < 1e-16, 1e-16, pdiff) 
  
  return(pdiff)
  
}

# some checks for the copula parameter / kendalls tau because it may be the case that a perfect 1 (or some other value at the edge of the parameter space)
# gets predicted. This funciton is harmless in 99.9999999% of the cases. 
check_extremes <- function(input, copula = "GAUSS"){
  
  if(copula == "GAUSS"){
    inp <- ifelse(input == 1, 0.99999999, input)
    inp <- ifelse(inp == -1, -0.99999999, inp)
  }
  
  if(copula == "FRANK"){
    
    epsilon <- c(sqrt(.Machine$double.eps))
    
    inp <- ifelse(abs(input) < epsilon, epsilon, input)
  }
  
  if(copula == "CLAYTON"){
    
    inp <- ifelse(input > 20, 20, input)
    inp <- ifelse(inp < -17, -17, inp)
    
  }
  
  return(inp)
  
}




### Rotated Clayton by 270 degrees
sim <- function(seed, n.train, n.mstop, p, n.test, corr, boost.nu = 0.1){
  
  
  ## Likelihood of the bivariate bernoulli distribution
  loss_bivbern <- function(mu1, mu2, or, y){
    y1 <- y[,1]
    y2 <- y[,2]
    
    
    loss_vec <- c()
    for(i in 1:length(y1)){
      
      psiminone = or[i] - 1
      
      hilfs1 = 1 + (mu1[i] + mu2[i])*psiminone
      hilfs2 = -4*or[i]*psiminone*mu1[i]*mu2[i]
      
      if(or[i] == 1) {
        p11 <- mu1[i] * mu2[i]
      }else{
        p11 = 0.5*(psiminone^(-1))*( hilfs1 - sqrt((hilfs1^2 + hilfs2)))
      }
      
      if(y1[i] == 0 && y2[i] == 0){
        loss_vec[i] <- -log(1 + p11 - mu2[i] - mu1[i])
      } else if(y1[i] == 0 && y2[i] == 1){
        loss_vec[i] <- -log(mu2[i] - p11)
      } else if(y1[i] == 1 && y2[i] == 0){
        loss_vec[i] <- -log(mu1[i] - p11)
      } else{
        loss_vec[i] <- -log(p11)
      }
      
      
    }
    loss_vec
  }
  
  loss <- function(mu1, mu2, sigma2, rho, y){
    
    F1 <- 1 - mu1 # F(0)
    
    F2 <- 1 - pnorm(y[,2], mean = mu2, sd = sigma2)
    
    dens2 <- dnorm(y[,2], mean = mu2, sd = sigma2)
    
    thet <- rho
    
    HFunc <-  pdffz( (F1^(-thet) + F2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (F2^((-thet) - 1) * (-thet))) )
    
    return( -( (1 - y[,1]) * log( HFunc ) + y[,1] * log( pdffz(1 - HFunc) ) + log(dens2) ) )
    
  }
  
  loss_VC <- function(mu1, mu2, sigma2, rho, y, FAM = 33){
    
    F1 <- 1 - mu1 # F(0)
    
    F2 <- pnorm(y[,2], mean = mu2, sd = sigma2)
    
    dens2 <- dnorm(y[,2], mean = mu2, sd = sigma2)
    
    thet <- rho
    
    HFunc <-  pdffz( VineCopula:::BiCopHfunc2(u1 = F1, u2 = F2, family = FAM, par = thet) ) #pdffz( (F1^(-thet) + F2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (F2^((-thet) - 1) * (-thet))) )
    
    return( -( (1 - y[,1]) * log( HFunc ) + y[,1] * log( pdffz(1 - HFunc) ) + log(dens2) ) )
    
  }
  
  
  # SAMPLING OF BIVARIATE RESPONSE FROM COPULA
  # simulate data from bivariate binary DGP:
  data.gen.bivbin <- function(FAM, mu1, mu2, sigma2, theta){
    
    
    u1u2 <- VineCopula::BiCopSim(length(mu1), family = FAM, par = theta)
    
    y1 <- qbinom(p = u1u2[,1], size = 1, prob = mu1)
    
    y2 <- qnorm(p = u1u2[,2], mean = mu2, sd = sigma2)
    
    dat <- data.frame(y1, y2)
    
    return(dat)
  }
  
  data.gen.bivbin_energyscore <- function(n, FAM, mu1, mu2, sigma2, theta, from_univariate = FALSE){
    
    #if(!from_univariate){
    
    u1u2 <- VineCopula::BiCopSim(n, family = FAM, par = theta)
    
    y1 <- qbinom(p = u1u2[,1], size = 1, prob = mu1)
    
    y2 <- qnorm(p = u1u2[,2],  mean = mu2, sd = sigma2)
    
    #}else{
    
    #y1 <- rbinom(n, size = 1, prob = mu1)
    
    #y2 <- rnorm(n,  mean = mu2, sd = sigma2)
    
    #}
    
    
    dat <- data.frame(y1, y2)
    
    return(dat)
  }
  
  
  # 1.5x2 − 1x3 + 1.5x4 − √x1ix1i,
  #mu1 <- c(1, 1.5, -1, 1.5, 0, 0, rep(0,p-6))
  #
  # x1, x2 ,x3 x4, x5, x6, x7, x8, x9, x10
  #mu1 <- c( 0, 1.5, -1, 1.5, 0, 0, 0, 0, 0, 0)
  mu1 <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  # 0.5x2 + 1.5x3
  #mu2 <- c(2, -1, 1, 0, 0, rep(0,p-5)) 
  #
  mu2 <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  
  # −0.5 + cos(2 * x2)
  #
  #sigma2 <- c(2, -1, 1, 0, 0, rep(0,p-5)) 
  sigma2 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  
  # −1.5 + 1x5 − 1.5x6 + sin(x3)
  #
  #or  <- c(0, 0, 1, 1.5, 0, 0, 0, 0,rep(0,p-8))
  or <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  
  nl_function_mu1 <- function(x){
    
    
    nl <- sqrt(x)*x^2 - 2*cos(3*x)
    
    return(nl*0.5)
    
  }
  
  nl_function_mu2 <- function(x){
    
    nl <- -0.7*exp(x^2) + exp(x^0.4) - 0*log(x+0.2) 
    
    return(nl)
  }
  
  nl_function_sigma2 <- function(x){
    
    return(cos(2 * x))
  }
  
  nl_function_rho <- function(x){
    
    return(sin(x*4))
  }
  
  TrueBeta <-  vector('list')
  TrueBeta$mu1 <- mu1
  TrueBeta$mu2 <- mu2
  TrueBeta$sigma2 <- sigma2
  TrueBeta$or <- or
  
  set.seed(seed)
  
  n.mstop = 1500
  n <- n.train + n.mstop
  weight.mstop <- c(rep(1, times = n.train),rep(0, times = n.mstop)) 
  
  # - training data
  #x.train <- rmvnorm(n = n, mean = rep(0,p), sigma =  toeplitz(sapply(seq(0,p-1), function(x) corr^x)))
  x1  <- runif(n, 0, 1)
  x2 <- runif(n, 0, 1)
  x3 <- runif(n, 0, 1)
  x4 <- runif(n, 0, 1)
  x5 <- runif(n, 0, 1)
  x6 <- runif(n, 0, 1)
  x7 <- runif(n, 0, 1)
  x8 <- runif(n, 0, 1)
  x9 <- runif(n, 0, 1) 
  x10 <- runif(n, 0, 1)
  
  x.train <- cbind(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
  
  
  ##########################################################################################################################################
  # GENERATION OF DISTRIBUTION PARAMETERS AND STUFF 
  # predictor
  train.eta.mu1_center <-  (x.train %*% mu1)
  
  ### Add non-linear effects (additive): − √x1 x1,
  train.eta.mu1_center <- train.eta.mu1_center + -1*0 + nl_function_mu1(x3)
  
  # apply response function: PROBIT
  train.eta.mu1 <-  pnorm(train.eta.mu1_center) 
  
  ##########################################################################################################################################
  # predictor and apply response function: IDENTITY
  train.eta.mu2 <-  x.train %*% mu2 + nl_function_mu2(x1)
  
  # no non-linear effect here
  train.eta.mu2 <- train.eta.mu2 + 0
  train.eta.mu2 <-  train.eta.mu2 # Identity link
  
  ##########################################################################################################################################
  #### SIGMA FOR CONTINOUS RESPONSE
  train.eta.sigma2 <-  x.train %*% sigma2
  
  train.eta.sigma2 <- train.eta.sigma2 + nl_function_sigma2(x2) - 0.5
  train.eta.sigma2 <- exp(train.eta.sigma2)
  
  ##########################################################################################################################################
  ### COPULA PARAMETER FOR ROTATED CLAYTON COPULA (270 DEGREES)
  train.eta.or_center <- (x.train %*% or)
  
  train.eta.or_center <- train.eta.or_center + - 1 + 3*nl_function_rho(x3)
  train.eta.or_center <- -exp(train.eta.or_center)
  train.copula.parameter <- train.eta.or_center
  
  ##########################################################################################################################################
  ##########################################################################################################################################
  
  y.train <- data.gen.bivbin(FAM = 33, mu1 = train.eta.mu1, mu2 = train.eta.mu2, sigma2 = train.eta.sigma2, theta = train.copula.parameter)
  
  colnames(y.train)<- c('y1','y2')
  dat.train <- data.frame(y.train,x.train)
  
  ############################################################################################################# - test data
  #x.test <- rmvnorm(n = n.test, mean = rep(0,p), sigma = toeplitz(sapply(seq(0,p-1), function(x) corr^x)))
  x1.test <- runif(n.test, 0, 1)
  x2.test <- runif(n.test, 0, 1)
  x3.test <- runif(n.test, 0, 1)
  x4.test <- runif(n.test, 0, 1)
  x5.test <- runif(n.test, 0, 1)
  x6.test <- runif(n.test, 0, 1)
  x7.test <- runif(n.test, 0, 1)
  x8.test <- runif(n.test, 0, 1)
  x9.test <- runif(n.test, 0, 1) 
  x10.test <- runif(n.test, 0, 1)
  
  x.test <- cbind(x1.test, x2.test, x3.test, x4.test, x5.test, x6.test, x7.test, x8.test, x9.test, x10.test)
  
  colnames(x.test) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
  
  ##########################################################################################################################################
  test.eta.mu1 <-  (x.test %*% mu1)
  test.eta.mu1_center <- test.eta.mu1 + -1*0 + nl_function_mu1(x3.test)
  test.eta.mu1 <-  pnorm(test.eta.mu1_center) 
  
  ##########################################################################################################################################
  test.eta.mu2 <-  x.test %*% mu2
  test.eta.mu2 <- test.eta.mu2 + 0 + nl_function_mu2(x1.test)
  test.eta.mu2 <-  test.eta.mu2 
  
  ##########################################################################################################################################
  test.eta.sigma2 <-  x.test %*% sigma2
  test.eta.sigma2 <- test.eta.sigma2 + nl_function_sigma2(x2.test) - 0.5
  test.eta.sigma2 <- exp(test.eta.sigma2)
  
  ##########################################################################################################################################
  test.eta.or_center <- (x.test %*% or)
  test.eta.or_center <- test.eta.or_center + - 1 + 3*nl_function_rho(x3.test)
  test.eta.or_center <- -exp(test.eta.or_center)
  test.copula.parameter <- test.eta.or_center
  
  
  ##########################################################################################################################################
  ##########################################################################################################################################
  # SAMPLE TEST OBSERVATIONS
  y.test <- data.gen.bivbin(FAM = 33, mu1 = test.eta.mu1, mu2 = test.eta.mu2, sigma2 = test.eta.sigma2, theta = test.copula.parameter)
  
  colnames(y.test)<- c('y1','y2')
  dat.test <- data.frame(y.test,x.test)
  
  ##########################################################################################################################################
  ##########################################################################################################################################
  
  mstop.bivBern <-  vector('list')
  
  
  AUCbivBern <-  vector('list')
  BrierbivBern <-  vector('list')
  lik <- vector('list')
  
  
  # # Predict distributional quantities
  model_equation <- formula(cbind(y1,y2) ~ bbs(x1, knots = 20, degree = 3, difference = 2) + 
                              bbs(x2, knots = 20, degree = 3, difference = 2) + 
                              bbs(x3, knots = 20, degree = 3, difference = 2) + 
                              bbs(x4, knots = 20, degree = 3, difference = 2) + 
                              bbs(x5, knots = 20, degree = 3, difference = 2) + 
                              bbs(x6, knots = 20, degree = 3, difference = 2) +
                              bbs(x7, knots = 20, degree = 3, difference = 2) + 
                              bbs(x8, knots = 20, degree = 3, difference = 2) + 
                              bbs(x9, knots = 20, degree = 3, difference = 2) + 
                              bbs(x10, knots = 20, degree = 3, difference = 2))
  
  model_equation <- formula(cbind(y1,y2) ~ bbs(x1) + 
                              bbs(x2) + 
                              bbs(x3) + 
                              bbs(x4) + 
                              bbs(x5) + 
                              bbs(x6) +
                              bbs(x7) + 
                              bbs(x8) + 
                              bbs(x9) + 
                              bbs(x10))
  
  bivModel_Formula <- list(mu1 = model_equation,
                           mu2 = model_equation,
                           sigma2 = model_equation,
                           rho = model_equation
  )
  
  #############################################
  ################# - COPULA - ################
  #############################################
  # Bivariate copula model
  bivBinCopula <- gamboostLSS(bivModel_Formula, data = dat.train, 
                              families = Clayton270_Cop_BinCont(marg1 = "PROBIT", 
                                                                marg2 = "NORM", 
                                                                stabilization = "L2"), 
                              control = boost_control(mstop = 2000, risk = 'oobag', nu = boost.nu, trace = TRUE), 
                              method = 'noncyclic', weights = weight.mstop)
  
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 1990){
    mstop(bivBinCopula) <- 4000
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 3990){
    bivBinCopula[8000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 7990){
    bivBinCopula[10000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 9990){
    bivBinCopula[15000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 14990){
    bivBinCopula[20000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula, merge = T))
  
  if(MSTOP_COP >= 19990){
    bivBinCopula[25000]
  }
  
  
  
  MSTOP_COP <- which.min(risk(bivBinCopula, merge = T))
  oobag.risk.cop <- risk(bivBinCopula,merge = T)
  
  #rm(bivBinCopula)
  #dat.train_biv <- dat.train[weight.mstop == 1, ]
  
  # RE-FIT until OPTIMAL MSTOP
  #bivBinCopula <- gamboostLSS(bivModel_Formula, data = dat.train_biv, 
  #                           families = Clayton270_Cop_BinCont(marg1 = "PROBIT", marg2 = "NORM"), 
  #                           control = boost_control(mstop = MSTOP_COP, nu  = boost.nu, trace = TRUE), method = 'noncyclic')
  bivBinCopula <- bivBinCopula[MSTOP_COP]
  
  
  mstop.bivBinCopula <-  vector('list')
  mstop.bivBinCopula$mstop <- MSTOP_COP
  mstop.bivBinCopula$mu1 <- bivBinCopula$mu1$mstop()
  mstop.bivBinCopula$mu2 <- bivBinCopula$mu2$mstop()
  mstop.bivBinCopula$sigma2 <- bivBinCopula$sigma2$mstop()
  mstop.bivBinCopula$rho  <- bivBinCopula$rho$mstop()
  
  # EXTRACT ALL COEFFICIENTS!
  coef.bivBinCopula <- coef(bivBinCopula, which = "")
  
  binaryMarginMetrics <-  vector('list')
  continuousMarginMetrics <-  vector('list')
  likBinCopula <- vector('list')
  
  
  ################ Get predictions from each covariate?
  newX <- matrix(rep(seq(from = 0, to = 1,length.out = 500),times = 10), nrow = 500)
  colnames(newX) <-  c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
  
  
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ### Copula predictions / selection
  Copula_MU1_predictions <- list()
  Copula_MU2_predictions <- list()
  Copula_SIGMA2_predictions <- list()
  Copula_RHO_predictions <- list()
  
  ### MU 1
  Copula_MU1_predictions$x1 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = 1, type = "link" )
  Copula_MU1_predictions$x2 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x2", type = "link" )
  Copula_MU1_predictions$x3 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x3", type = "link" )
  Copula_MU1_predictions$x4 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x4", type = "link" )
  Copula_MU1_predictions$x5 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x5", type = "link" )
  Copula_MU1_predictions$x6 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x6", type = "link" )
  Copula_MU1_predictions$x7 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x7", type = "link" )
  Copula_MU1_predictions$x8 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x8", type = "link" )
  Copula_MU1_predictions$x9 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x9", type = "link" )
  Copula_MU1_predictions$x10 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x10", type = "link" )
  
  
  #### MU 2
  Copula_MU2_predictions$x1 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = 1, type = "link" )
  Copula_MU2_predictions$x2 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x2", type = "link" )
  Copula_MU2_predictions$x3 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x3", type = "link" )
  Copula_MU2_predictions$x4 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x4", type = "link" )
  Copula_MU2_predictions$x5 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x5", type = "link" )
  Copula_MU2_predictions$x6 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x6", type = "link" )
  Copula_MU2_predictions$x7 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x7", type = "link" )
  Copula_MU2_predictions$x8 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x8", type = "link" )
  Copula_MU2_predictions$x9 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x9", type = "link" )
  Copula_MU2_predictions$x10 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x10", type = "link" )
  
  #### SIGMA 2
  Copula_SIGMA2_predictions$x1 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = 1, type = "link" )
  Copula_SIGMA2_predictions$x2 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x2", type = "link" )
  Copula_SIGMA2_predictions$x3 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x3", type = "link" )
  Copula_SIGMA2_predictions$x4 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x4", type = "link" )
  Copula_SIGMA2_predictions$x5 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x5", type = "link" )
  Copula_SIGMA2_predictions$x6 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x6", type = "link" )
  Copula_SIGMA2_predictions$x7 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x7", type = "link" )
  Copula_SIGMA2_predictions$x8 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x8", type = "link" )
  Copula_SIGMA2_predictions$x9 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x9", type = "link" )
  Copula_SIGMA2_predictions$x10 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x10", type = "link" )
  
  
  #### RHO
  Copula_RHO_predictions$x1 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = 1, type = "link" )
  Copula_RHO_predictions$x2 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x2", type = "link" )
  Copula_RHO_predictions$x3 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x3", type = "link" )
  Copula_RHO_predictions$x4 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x4", type = "link" )
  Copula_RHO_predictions$x5 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x5", type = "link" )
  Copula_RHO_predictions$x6 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x6", type = "link" )
  Copula_RHO_predictions$x7 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x7", type = "link" )
  Copula_RHO_predictions$x8 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x8", type = "link" )
  Copula_RHO_predictions$x9 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x9", type = "link" )
  Copula_RHO_predictions$x10 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x10", type = "link" )
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  
  # Predict distributional quantities
  predCopula.mu1 <- predict(bivBinCopula$mu1, newdata = dat.test, type = 'response')
  predCopula.mu2 <- predict(bivBinCopula$mu2, newdata = dat.test, type = 'response')
  predCopula.sigma2 <- predict(bivBinCopula$sigma2, newdata = dat.test, type = 'response')
  predCopula.rho <- predict(bivBinCopula$rho, newdata = dat.test, type = 'response')
  
  ######### COMPUTE PERFORMANCE: 
  # BINARY MARGIN
  # AUC (margin 1)
  binaryMarginMetrics$AUC <- roc(dat.test$y1, as.numeric(predCopula.mu1))$auc
  # Brier score (margin 1, then margin 2)
  binaryMarginMetrics$BrierScore <- mean((as.numeric(predCopula.mu1)-(as.numeric(dat.test$y1)))^2)
  
  ### CONTINUOUS MARGIN: 
  continuousMarginMetrics$MSE <- mean((as.numeric(predCopula.mu2)-as.numeric(dat.test$y2))^2)
  
  
  # LOG-SCORE (LOG-LIKELIHOOD) (change to the copula likelihood)
  lik$biv <- sum(loss(mu1 = predCopula.mu1, mu2 = predCopula.mu2, sigma2 = predCopula.sigma2, rho = predCopula.rho, y = y.test))
  
  
  ###########################################################################################################################
  ################ - univariate fits - ######################################################################################
  ###########################################################################################################################
  
  univariateBinaryEquation <- formula(as.factor(y1) ~  bbs(x1) + 
                                        bbs(x2) + 
                                        bbs(x3) + 
                                        bbs(x4) + 
                                        bbs(x5) + 
                                        bbs(x6) +
                                        bbs(x7) + 
                                        bbs(x8) + 
                                        bbs(x9) + 
                                        bbs(x10))
  
  univariateContinuousEquation <- list(mu = formula(y2 ~ bbs(x1) + 
                                                      bbs(x2) + 
                                                      bbs(x3) + 
                                                      bbs(x4) + 
                                                      bbs(x5) + 
                                                      bbs(x6) +
                                                      bbs(x7) + 
                                                      bbs(x8) + 
                                                      bbs(x9) + 
                                                      bbs(x10)),
                                       sigma = formula(y2 ~ bbs(x1) + 
                                                         bbs(x2) + 
                                                         bbs(x3) + 
                                                         bbs(x4) + 
                                                         bbs(x5) + 
                                                         bbs(x6) +
                                                         bbs(x7) + 
                                                         bbs(x8) + 
                                                         bbs(x9) + 
                                                         bbs(x10))
  )
  
  
  
  UNIVARIATE_MU1_predictions <- list()
  UNIVARIATE_MU2_predictions <- list()
  UNIVARIATE_SIGMA2_predictions <- list()
  
  mstop.uni <- vector('list')
  coef.uni  <- vector('list')
  coef.uni_1  <- vector('list')
  AUC.uni   <- vector('list')
  Brier.uni <- vector('list')
  
  # - mu1
  dat.train.mu1 = dat.train[, !(names(dat.train) %in% "y2")]
  dat.test.mu1  = dat.test[ , !(names(dat.test) %in% "y2")]
  
  glm.uni.mu1 <- gamboost(as.factor(y1) ~. , data = dat.train.mu1, 
                          family = Binomial(link = "probit", type = "glm"), 
                          control = boost_control(risk = 'oobag', nu = boost.nu, mstop = 1000, trace = T), 
                          weights = weight.mstop)
  
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  if(mstop.uni$mu1 >= 90){
    glm.uni.mu1[1000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 990){
    glm.uni.mu1[5000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 4990){
    glm.uni.mu1[10000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  #dat.train_mu1 <- dat.train.mu1[weight.mstop == 1, ]
  
  #glm.uni.mu1 <- gamboost(as.factor(y1) ~. , data = dat.train_mu1, family = Binomial(link = "probit", type = "glm"), 
  #                       control = boost_control(mstop = mstop.uni$mu1, nu  = 0.01))
  glm.uni.mu1 <- glm.uni.mu1[mstop.uni$mu1]
  
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 1  
  coef.uni$mu1 <- coef(glm.uni.mu1, which = '')
  
  binaryMarginMetrics$Univariate_AUC <- roc(dat.test.mu1$y1, as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response")))$auc
  binaryMarginMetrics$Univariate_BrierScore <- mean((as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1, type = "response"))-(as.numeric(dat.test.mu1$y1)))^2)
  
  ### MU 1
  UNIVARIATE_MU1_predictions$x1 <- predict(glm.uni.mu1, data.frame(newX), which = 1, type = "link" )
  UNIVARIATE_MU1_predictions$x2 <- predict(glm.uni.mu1, data.frame(newX), which = "x2", type = "link" )
  UNIVARIATE_MU1_predictions$x3 <- predict(glm.uni.mu1, data.frame(newX), which = "x3", type = "link" )
  UNIVARIATE_MU1_predictions$x4 <- predict(glm.uni.mu1, data.frame(newX), which = "x4", type = "link" )
  UNIVARIATE_MU1_predictions$x5 <- predict(glm.uni.mu1, data.frame(newX), which = "x5", type = "link" )
  UNIVARIATE_MU1_predictions$x6 <- predict(glm.uni.mu1, data.frame(newX), which = "x6", type = "link" )
  UNIVARIATE_MU1_predictions$x7 <- predict(glm.uni.mu1, data.frame(newX), which = "x7", type = "link" )
  UNIVARIATE_MU1_predictions$x8 <- predict(glm.uni.mu1, data.frame(newX), which = "x8", type = "link" )
  UNIVARIATE_MU1_predictions$x9 <- predict(glm.uni.mu1, data.frame(newX), which = "x9", type = "link" )
  UNIVARIATE_MU1_predictions$x10 <- predict(glm.uni.mu1, data.frame(newX), which = "x10", type = "link" )
  
  
  ### CONTINUOUS MARGIN
  # - mu2
  dat.train.mu2 =  dat.train[, !(names(dat.train) %in% "y1")]
  dat.test.mu2 = dat.test[, !(names(dat.test) %in% "y1")]
  
  glm.uni.mu2 <- gamboostLSS(y2 ~., data = dat.train.mu2, 
                             families = GaussianLSS(stabilization = "L2"), 
                             method = "noncyclic", 
                             control = boost_control(risk = 'oobag', nu = boost.nu, mstop = 1000, trace = TRUE), 
                             weights = weight.mstop)
  
  
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  if(mstop.uni$mu2 >= 90){
    glm.uni.mu2[1000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  if(mstop.uni$mu2 >= 990){
    glm.uni.mu2[10000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  if(mstop.uni$mu2 >= 9990){
    glm.uni.mu2[15000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  #dat.train_mu2 <- dat.train.mu2[weight.mstop == 1, ]
  
  #glm.uni.mu2 <- gamboostLSS(y2 ~., data = dat.train_mu2, families = GaussianLSS(), method = "noncyclic", 
  #                           control = boost_control(mstop = mstop.uni$mu2, nu  = boost.nu, trace = TRUE))
  glm.uni.mu2 <- glm.uni.mu2[mstop.uni$mu2]
  
  
  mstop.uni$sigma2 <- mstop(glm.uni.mu2$sigma)
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 2  
  coef.uni$mu2 <- coef(glm.uni.mu2$mu, which = '')
  coef.uni$sigma2 <- coef(glm.uni.mu2$sigma, which = '')
  
  continuousMarginMetrics$Univariate_MSE <- mean(((predict(glm.uni.mu2$mu, newdata = dat.test.mu2, type = "response"))-as.numeric(dat.test$y2))^2)
  
  
  
  # Likelihood
  pred.mu1.uni <- as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1, type = "response"))
  pred.mu2.uni <- as.numeric(predict(glm.uni.mu2$mu, newdata = dat.test.mu2, type = "response"))
  pred.sigma2.uni <- as.numeric(predict(glm.uni.mu2$sigma, newdata = dat.test.mu2, type = "response"))
  
  mu1.uni.loglik <- sum(-dbinom(x = dat.test.mu1$y1, size = 1, prob = pred.mu1.uni, log = T))
  mu2.uni.loglik <- sum(-dnorm(x = dat.test.mu2$y2, mean = pred.mu2.uni, sd = pred.sigma2.uni, log = T))
  
  lik$uni_usingCopula <- sum(loss(mu1 = pred.mu1.uni, mu2 = pred.mu2.uni, sigma2 = pred.sigma2.uni, rho = rep((.Machine$double.eps)^2, length(nrow(y.test))), y = y.test))
  
  lik$uni<- mu1.uni.loglik + mu2.uni.loglik
  
  
  n.ges = c(n.train, n.mstop, n.test)
  pred.ges <- NULL#list(pred.mu1, pred.mu2, pred.or)
  pred.uni <- list(mu1 = pred.mu1.uni, mu2 = pred.mu2.uni, sigma2 = pred.sigma2.uni)
  pred.cop <- list(mu1 = predCopula.mu1, mu2 = predCopula.mu2, sigma2 = predCopula.sigma2, rho = predCopula.rho)
  
  
  
  #### MU 2
  UNIVARIATE_MU2_predictions$x1 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = 1, type = "link" )
  UNIVARIATE_MU2_predictions$x2 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x2", type = "link" )
  UNIVARIATE_MU2_predictions$x3 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x3", type = "link" )
  UNIVARIATE_MU2_predictions$x4 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x4", type = "link" )
  UNIVARIATE_MU2_predictions$x5 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x5", type = "link" )
  UNIVARIATE_MU2_predictions$x6 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x6", type = "link" )
  UNIVARIATE_MU2_predictions$x7 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x7", type = "link" )
  UNIVARIATE_MU2_predictions$x8 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x8", type = "link" )
  UNIVARIATE_MU2_predictions$x9 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x9", type = "link" )
  UNIVARIATE_MU2_predictions$x10 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x10", type = "link" )
  
  #### SIGMA 2
  UNIVARIATE_SIGMA2_predictions$x1 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = 1, type = "link" )
  UNIVARIATE_SIGMA2_predictions$x2 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x2", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x3 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x3", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x4 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x4", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x5 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x5", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x6 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x6", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x7 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x7", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x8 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x8", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x9 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x9", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x10 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x10", type = "link" )
  
  
  
  ###############################################################################################################################
  ######## Energy Score ########################################################################################################
  ##############################################################################################################################
  es_biv <- vector()
  es_uni<- vector()
  es_cop <- vector()
  es_gjrm <- vector()
  
  for(i in 1:length(pred.mu1.uni)){
    
    pred_sample_uni <- matrix(NA, nrow = 2, ncol = 1000)
    pred_sample_biv <- matrix(NA, nrow = 2, ncol = 1000)
    pred_sample_cop <- matrix(NA, nrow = 2, ncol = 1000)
    
    
    
    # univariate
    #sample_uni <- rbinom2.or(10000, mu1 =  pred.mu1.uni[i], mu2 = pred.mu2.uni[i], oratio = 1)
    sample_uni <- data.gen.bivbin_energyscore(n = 1000, 
                                              mu1 = pred.mu1.uni[i], 
                                              mu2 = pred.mu2.uni[i], 
                                              sigma2 = pred.sigma2.uni[i], 
                                              FAM = 33, theta = -(.Machine$double.eps)^2,
                                              from_univariate = TRUE)
    
    pred_sample_uni[1,] <- sample_uni[,1]
    pred_sample_uni[2,] <- sample_uni[,2]
    
    es_uni[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_uni) 
    
    
    
    
    # copula
    sample_cop <- data.gen.bivbin_energyscore(n = 1000, 
                                              mu1 = predCopula.mu1[i],  
                                              mu2 = predCopula.mu2[i], 
                                              sigma2 = predCopula.sigma2[i], 
                                              FAM = 33, theta = -predCopula.rho[i])
    
    pred_sample_cop[1, ] <- sample_cop[,1]
    pred_sample_cop[2, ] <- sample_cop[,2]
    
    es_cop[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_cop) 
    
    
  }
  
  energy_score <- list()
  #energy_score$biv <- mean(es_biv)
  energy_score$uni <- mean(es_uni)
  energy_score$cop <- mean(es_cop)
  
  
  ##################################### Collect: mstop, predictions, 
  #                                             AUC for binary margin, Brier score for binary margin, 
  #                                             MSE for continuous margin
  #                                             energy score and negative log likelihood 
  output <- list(TrueBeta = TrueBeta, 
                 n = n.ges, 
                 p  = p, 
                 Likelihood = lik,
                 predict = pred.ges, 
                 predict.uni = pred.uni, 
                 predict.cop = pred.cop,
                 energy_score = energy_score,  
                 ####
                 Coefficients.uni = coef.uni, 
                 mstop.uni = mstop.uni,
                 ################################################
                 BinaryMarginMetrics = binaryMarginMetrics,
                 ContinuousMarginMetrics = continuousMarginMetrics,
                 ################################################
                 CoefficientsCOPULA = coef.bivBinCopula, 
                 oobag.riskCOPULA = oobag.risk.cop, 
                 mstopCOPULA = mstop.bivBinCopula, 
                 ###
                 Copula_Predictions = list(MU1 = Copula_MU1_predictions, 
                                           MU2 = Copula_MU2_predictions,
                                           SIGMA2 = Copula_SIGMA2_predictions,
                                           RHO = Copula_RHO_predictions),
                 Univariate_Predictions = list(MU1 = UNIVARIATE_MU1_predictions, 
                                               MU2 = UNIVARIATE_MU2_predictions, 
                                               SIGMA2 = UNIVARIATE_SIGMA2_predictions)
  )
  
  
  
  return(output)
}

### FRANK COPULA
sim_FRANK <- function(seed, n.train, n.mstop, p, n.test, corr, boost.nu = 0.1){
  
  
  ## Likelihood of the bivariate bernoulli distribution
  loss_bivbern <- function(mu1, mu2, or, y){
    y1 <- y[,1]
    y2 <- y[,2]
    
    
    loss_vec <- c()
    for(i in 1:length(y1)){
      
      psiminone = or[i] - 1
      
      hilfs1 = 1 + (mu1[i] + mu2[i])*psiminone
      hilfs2 = -4*or[i]*psiminone*mu1[i]*mu2[i]
      
      if(or[i] == 1) {
        p11 <- mu1[i] * mu2[i]
      }else{
        p11 = 0.5*(psiminone^(-1))*( hilfs1 - sqrt((hilfs1^2 + hilfs2)))
      }
      
      if(y1[i] == 0 && y2[i] == 0){
        loss_vec[i] <- -log(1 + p11 - mu2[i] - mu1[i])
      } else if(y1[i] == 0 && y2[i] == 1){
        loss_vec[i] <- -log(mu2[i] - p11)
      } else if(y1[i] == 1 && y2[i] == 0){
        loss_vec[i] <- -log(mu1[i] - p11)
      } else{
        loss_vec[i] <- -log(p11)
      }
      
      
    }
    loss_vec
  }
  
  loss <- function(mu1, mu2, sigma2, rho, y){
    
    F1 <- 1 - mu1 # F(0)
    
    F2 <- pnorm(y[,2], mean = mu2, sd = sigma2)
    
    dens2 <- dnorm(y[,2], mean = mu2, sd = sigma2)
    
    thet <- rho
    
    HFunc <- pdffz( (exp(thet)* (-1 + exp(F1* thet)))/(-exp((F1 + F2)* thet) + exp(thet)* (-1 + exp(F2* thet) + exp(F1* thet))) )
    
    return( -( (1 - y[,1]) * log( HFunc ) + y[,1] * log( pdffz(1 - HFunc) ) + log(dens2) ) )
    
  }
  
  loss_VC <- function(mu1, mu2, sigma2, rho, y, FAM = 5){
    
    F1 <- 1 - mu1 # F(0)
    
    F2 <- pnorm(y[,2], mean = mu2, sd = sigma2)
    
    dens2 <- dnorm(y[,2], mean = mu2, sd = sigma2)
    
    thet <- rho
    
    HFunc <-  pdffz( VineCopula:::BiCopHfunc2(u1 = F1, u2 = F2, family = FAM, par = thet) ) #pdffz( (F1^(-thet) + F2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (F2^((-thet) - 1) * (-thet))) )
    
    return( -( (1 - y[,1]) * log( HFunc ) + y[,1] * log( pdffz(1 - HFunc) ) + log(dens2) ) )
    
  }
  
  
  # SAMPLING OF BIVARIATE RESPONSE FROM COPULA
  # simulate data from bivariate binary DGP:
  data.gen.bivbin <- function(FAM, mu1, mu2, sigma2, theta){
    
    
    u1u2 <- VineCopula::BiCopSim(length(mu1), family = FAM, par = theta)
    
    y1 <- qbinom(p = u1u2[,1], size = 1, prob = mu1)
    
    y2 <- qnorm(p = u1u2[,2], mean = mu2, sd = sigma2)
    
    dat <- data.frame(y1, y2)
    
    return(dat)
  }
  
  data.gen.bivbin_energyscore <- function(n, FAM, mu1, mu2, sigma2, theta, from_univariate = FALSE){
    
    #if(!from_univariate){
    
    u1u2 <- VineCopula::BiCopSim(n, family = FAM, par = theta)
    
    y1 <- qbinom(p = u1u2[,1], size = 1, prob = mu1)
    
    y2 <- qnorm(p = u1u2[,2],  mean = mu2, sd = sigma2)
    
    #}else{
    
    #y1 <- rbinom(n, size = 1, prob = mu1)
    
    #y2 <- rnorm(n,  mean = mu2, sd = sigma2)
    
    #}
    
    
    dat <- data.frame(y1, y2)
    
    return(dat)
  }
  
  
  # 1.5x2 − 1x3 + 1.5x4 − √x1ix1i,
  #mu1 <- c(1, 1.5, -1, 1.5, 0, 0, rep(0,p-6))
  #
  # x1, x2 ,x3 x4, x5, x6, x7, x8, x9, x10
  #mu1 <- c( 0, 1.5, -1, 1.5, 0, 0, 0, 0, 0, 0)
  mu1 <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  # 0.5x2 + 1.5x3
  #mu2 <- c(2, -1, 1, 0, 0, rep(0,p-5)) 
  #
  mu2 <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  
  # −0.5 + cos(2 * x2)
  #
  #sigma2 <- c(2, -1, 1, 0, 0, rep(0,p-5)) 
  sigma2 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  
  # −1.5 + 1x5 − 1.5x6 + sin(x3)
  #
  #or  <- c(0, 0, 1, 1.5, 0, 0, 0, 0,rep(0,p-8))
  or <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  
  nl_function_mu1 <- function(x){
    
    
    nl <- sqrt(x)*x^2 - 2*cos(3*x)
    
    return(nl*0.5)
    
  }
  
  nl_function_mu2 <- function(x){
    
    nl <- -0.7*exp(x^2) + exp(x^0.4) - 0*log(x+0.2) 
    
    return(nl)
  }
  
  nl_function_sigma2 <- function(x){
    
    return(cos(2 * x))
  }
  
  nl_function_rho <- function(x){
    
    return(10* sin(x*4) )
  }
  
  TrueBeta <-  vector('list')
  TrueBeta$mu1 <- mu1
  TrueBeta$mu2 <- mu2
  TrueBeta$sigma2 <- sigma2
  TrueBeta$or <- or
  
  set.seed(seed)
  
  n.mstop = 1500
  n <- n.train + n.mstop
  weight.mstop <- c(rep(1, times = n.train),rep(0, times = n.mstop)) 
  
  # - training data
  #x.train <- rmvnorm(n = n, mean = rep(0,p), sigma =  toeplitz(sapply(seq(0,p-1), function(x) corr^x)))
  x1  <- runif(n, 0, 1)
  x2 <- runif(n, 0, 1)
  x3 <- runif(n, 0, 1)
  x4 <- runif(n, 0, 1)
  x5 <- runif(n, 0, 1)
  x6 <- runif(n, 0, 1)
  x7 <- runif(n, 0, 1)
  x8 <- runif(n, 0, 1)
  x9 <- runif(n, 0, 1) 
  x10 <- runif(n, 0, 1)
  
  x.train <- cbind(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
  
  
  ##########################################################################################################################################
  # GENERATION OF DISTRIBUTION PARAMETERS AND STUFF 
  # predictor
  train.eta.mu1_center <-  (x.train %*% mu1)
  
  ### Add non-linear effects (additive): − √x1 x1,
  train.eta.mu1_center <- train.eta.mu1_center + -1*0 + nl_function_mu1(x3)
  
  # apply response function: PROBIT
  train.eta.mu1 <-  pnorm(train.eta.mu1_center) 
  
  ##########################################################################################################################################
  # predictor and apply response function: IDENTITY
  train.eta.mu2 <-  x.train %*% mu2 + nl_function_mu2(x1)
  
  # no non-linear effect here
  train.eta.mu2 <- train.eta.mu2 + 0
  train.eta.mu2 <-  train.eta.mu2 # Identity link
  
  ##########################################################################################################################################
  #### SIGMA FOR CONTINOUS RESPONSE
  train.eta.sigma2 <-  x.train %*% sigma2
  
  train.eta.sigma2 <- train.eta.sigma2 + nl_function_sigma2(x2) - 0.5
  train.eta.sigma2 <- exp(train.eta.sigma2)
  
  ##########################################################################################################################################
  ### COPULA PARAMETER FOR ROTATED CLAYTON COPULA (270 DEGREES)
  train.eta.or_center <- (x.train %*% or)
  
  train.eta.or_center <- train.eta.or_center + - 1 + 3*nl_function_rho(x3)
  #train.eta.or_center <- -exp(train.eta.or_center)
  train.copula.parameter <- train.eta.or_center
  
  ##########################################################################################################################################
  ##########################################################################################################################################
  
  y.train <- data.gen.bivbin(FAM = 5, mu1 = train.eta.mu1, mu2 = train.eta.mu2, sigma2 = train.eta.sigma2, theta = train.copula.parameter)
  
  colnames(y.train)<- c('y1','y2')
  dat.train <- data.frame(y.train,x.train)
  
  ############################################################################################################# - test data
  #x.test <- rmvnorm(n = n.test, mean = rep(0,p), sigma = toeplitz(sapply(seq(0,p-1), function(x) corr^x)))
  x1.test <- runif(n.test, 0, 1)
  x2.test <- runif(n.test, 0, 1)
  x3.test <- runif(n.test, 0, 1)
  x4.test <- runif(n.test, 0, 1)
  x5.test <- runif(n.test, 0, 1)
  x6.test <- runif(n.test, 0, 1)
  x7.test <- runif(n.test, 0, 1)
  x8.test <- runif(n.test, 0, 1)
  x9.test <- runif(n.test, 0, 1) 
  x10.test <- runif(n.test, 0, 1)
  
  x.test <- cbind(x1.test, x2.test, x3.test, x4.test, x5.test, x6.test, x7.test, x8.test, x9.test, x10.test)
  
  colnames(x.test) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
  
  ##########################################################################################################################################
  test.eta.mu1 <-  (x.test %*% mu1)
  test.eta.mu1_center <- test.eta.mu1 + -1*0 + nl_function_mu1(x3.test)
  test.eta.mu1 <-  pnorm(test.eta.mu1_center) 
  
  ##########################################################################################################################################
  test.eta.mu2 <-  x.test %*% mu2
  test.eta.mu2 <- test.eta.mu2 + 0 + nl_function_mu2(x1.test)
  test.eta.mu2 <-  test.eta.mu2 
  
  ##########################################################################################################################################
  test.eta.sigma2 <-  x.test %*% sigma2
  test.eta.sigma2 <- test.eta.sigma2 + nl_function_sigma2(x2.test) - 0.5
  test.eta.sigma2 <- exp(test.eta.sigma2)
  
  ##########################################################################################################################################
  test.eta.or_center <- (x.test %*% or)
  test.eta.or_center <- test.eta.or_center + - 1 + 3*nl_function_rho(x3.test)
  #test.eta.or_center <- -exp(test.eta.or_center)
  test.copula.parameter <- test.eta.or_center
  
  
  
  
  TrueKendallRange <- VineCopula::BiCopPar2Tau(family = 5, par = range(train.copula.parameter))
  
  
  ##########################################################################################################################################
  ##########################################################################################################################################
  # SAMPLE TEST OBSERVATIONS
  y.test <- data.gen.bivbin(FAM = 5, mu1 = test.eta.mu1, mu2 = test.eta.mu2, sigma2 = test.eta.sigma2, theta = test.copula.parameter)
  
  colnames(y.test)<- c('y1','y2')
  dat.test <- data.frame(y.test,x.test)
  
  ##########################################################################################################################################
  ##########################################################################################################################################
  
  mstop.bivBern <-  vector('list')
  
  
  AUCbivBern <-  vector('list')
  BrierbivBern <-  vector('list')
  lik <- vector('list')
  
  
  # # Predict distributional quantities
  model_equation <- formula(cbind(y1,y2) ~ bbs(x1, knots = 20, degree = 3, difference = 2) + 
                              bbs(x2, knots = 20, degree = 3, difference = 2) + 
                              bbs(x3, knots = 20, degree = 3, difference = 2) + 
                              bbs(x4, knots = 20, degree = 3, difference = 2) + 
                              bbs(x5, knots = 20, degree = 3, difference = 2) + 
                              bbs(x6, knots = 20, degree = 3, difference = 2) +
                              bbs(x7, knots = 20, degree = 3, difference = 2) + 
                              bbs(x8, knots = 20, degree = 3, difference = 2) + 
                              bbs(x9, knots = 20, degree = 3, difference = 2) + 
                              bbs(x10, knots = 20, degree = 3, difference = 2))
  
  model_equation <- formula(cbind(y1,y2) ~ bbs(x1) + 
                              bbs(x2) + 
                              bbs(x3) + 
                              bbs(x4) + 
                              bbs(x5) + 
                              bbs(x6) +
                              bbs(x7) + 
                              bbs(x8) + 
                              bbs(x9) + 
                              bbs(x10))
  
  bivModel_Formula <- list(mu1 = model_equation,
                           mu2 = model_equation,
                           sigma2 = model_equation,
                           rho = model_equation
  )
  
  #############################################
  ################# - COPULA - ################
  #############################################
  # Bivariate copula model
  bivBinCopula <- gamboostLSS(bivModel_Formula, data = dat.train, 
                              families = Frank_Cop_BinCont(marg1 = "PROBIT", 
                                                                marg2 = "NORM", 
                                                                stabilization = "L2"), 
                              control = boost_control(mstop = 2000, risk = 'oobag', nu = boost.nu, trace = TRUE), 
                              method = 'noncyclic', weights = weight.mstop)
  
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 1990){
    mstop(bivBinCopula) <- 4000
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 3990){
    bivBinCopula[8000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 7990){
    bivBinCopula[10000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 9990){
    bivBinCopula[15000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 14990){
    bivBinCopula[20000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula, merge = T))
  
  if(MSTOP_COP >= 19990){
    bivBinCopula[25000]
  }
  
  
  
  MSTOP_COP <- which.min(risk(bivBinCopula, merge = T))
  oobag.risk.cop <- risk(bivBinCopula,merge = T)
  
  #rm(bivBinCopula)
  #dat.train_biv <- dat.train[weight.mstop == 1, ]
  
  # RE-FIT until OPTIMAL MSTOP
  #bivBinCopula <- gamboostLSS(bivModel_Formula, data = dat.train_biv, 
  #                           families = Clayton270_Cop_BinCont(marg1 = "PROBIT", marg2 = "NORM"), 
  #                           control = boost_control(mstop = MSTOP_COP, nu  = boost.nu, trace = TRUE), method = 'noncyclic')
  bivBinCopula <- bivBinCopula[MSTOP_COP]
  
  
  mstop.bivBinCopula <-  vector('list')
  mstop.bivBinCopula$mstop <- MSTOP_COP
  mstop.bivBinCopula$mu1 <- bivBinCopula$mu1$mstop()
  mstop.bivBinCopula$mu2 <- bivBinCopula$mu2$mstop()
  mstop.bivBinCopula$sigma2 <- bivBinCopula$sigma2$mstop()
  mstop.bivBinCopula$rho  <- bivBinCopula$rho$mstop()
  
  # EXTRACT ALL COEFFICIENTS!
  coef.bivBinCopula <- coef(bivBinCopula, which = "")
  
  binaryMarginMetrics <-  vector('list')
  continuousMarginMetrics <-  vector('list')
  likBinCopula <- vector('list')
  
  
  ################ Get predictions from each covariate?
  newX <- matrix(rep(seq(from = 0, to = 1,length.out = 500),times = 10), nrow = 500)
  colnames(newX) <-  c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
  
  
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ### Copula predictions / selection
  Copula_MU1_predictions <- list()
  Copula_MU2_predictions <- list()
  Copula_SIGMA2_predictions <- list()
  Copula_RHO_predictions <- list()
  
  ### MU 1
  Copula_MU1_predictions$x1 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = 1, type = "link" )
  Copula_MU1_predictions$x2 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x2", type = "link" )
  Copula_MU1_predictions$x3 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x3", type = "link" )
  Copula_MU1_predictions$x4 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x4", type = "link" )
  Copula_MU1_predictions$x5 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x5", type = "link" )
  Copula_MU1_predictions$x6 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x6", type = "link" )
  Copula_MU1_predictions$x7 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x7", type = "link" )
  Copula_MU1_predictions$x8 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x8", type = "link" )
  Copula_MU1_predictions$x9 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x9", type = "link" )
  Copula_MU1_predictions$x10 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x10", type = "link" )
  
  
  #### MU 2
  Copula_MU2_predictions$x1 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = 1, type = "link" )
  Copula_MU2_predictions$x2 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x2", type = "link" )
  Copula_MU2_predictions$x3 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x3", type = "link" )
  Copula_MU2_predictions$x4 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x4", type = "link" )
  Copula_MU2_predictions$x5 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x5", type = "link" )
  Copula_MU2_predictions$x6 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x6", type = "link" )
  Copula_MU2_predictions$x7 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x7", type = "link" )
  Copula_MU2_predictions$x8 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x8", type = "link" )
  Copula_MU2_predictions$x9 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x9", type = "link" )
  Copula_MU2_predictions$x10 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x10", type = "link" )
  
  #### SIGMA 2
  Copula_SIGMA2_predictions$x1 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = 1, type = "link" )
  Copula_SIGMA2_predictions$x2 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x2", type = "link" )
  Copula_SIGMA2_predictions$x3 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x3", type = "link" )
  Copula_SIGMA2_predictions$x4 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x4", type = "link" )
  Copula_SIGMA2_predictions$x5 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x5", type = "link" )
  Copula_SIGMA2_predictions$x6 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x6", type = "link" )
  Copula_SIGMA2_predictions$x7 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x7", type = "link" )
  Copula_SIGMA2_predictions$x8 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x8", type = "link" )
  Copula_SIGMA2_predictions$x9 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x9", type = "link" )
  Copula_SIGMA2_predictions$x10 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x10", type = "link" )
  
  
  #### RHO
  Copula_RHO_predictions$x1 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = 1, type = "link" )
  Copula_RHO_predictions$x2 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x2", type = "link" )
  Copula_RHO_predictions$x3 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x3", type = "link" )
  Copula_RHO_predictions$x4 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x4", type = "link" )
  Copula_RHO_predictions$x5 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x5", type = "link" )
  Copula_RHO_predictions$x6 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x6", type = "link" )
  Copula_RHO_predictions$x7 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x7", type = "link" )
  Copula_RHO_predictions$x8 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x8", type = "link" )
  Copula_RHO_predictions$x9 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x9", type = "link" )
  Copula_RHO_predictions$x10 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x10", type = "link" )
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  
  # Predict distributional quantities
  predCopula.mu1 <- predict(bivBinCopula$mu1, newdata = dat.test, type = 'response')
  predCopula.mu2 <- predict(bivBinCopula$mu2, newdata = dat.test, type = 'response')
  predCopula.sigma2 <- predict(bivBinCopula$sigma2, newdata = dat.test, type = 'response')
  predCopula.rho <- predict(bivBinCopula$rho, newdata = dat.test, type = 'response')
  
  ######### COMPUTE PERFORMANCE: 
  # BINARY MARGIN
  # AUC (margin 1)
  binaryMarginMetrics$AUC <- roc(dat.test$y1, as.numeric(predCopula.mu1))$auc
  # Brier score (margin 1, then margin 2)
  binaryMarginMetrics$BrierScore <- mean((as.numeric(predCopula.mu1)-(as.numeric(dat.test$y1)))^2)
  
  ### CONTINUOUS MARGIN: 
  continuousMarginMetrics$MSE <- mean((as.numeric(predCopula.mu2)-as.numeric(dat.test$y2))^2)
  
  
  # LOG-SCORE (LOG-LIKELIHOOD) (change to the copula likelihood)
  lik$biv <- sum(loss(mu1 = predCopula.mu1, mu2 = predCopula.mu2, sigma2 = predCopula.sigma2, rho = predCopula.rho, y = y.test))
  
  
  ###########################################################################################################################
  ################ - univariate fits - ######################################################################################
  ###########################################################################################################################
  
  univariateBinaryEquation <- formula(as.factor(y1) ~  bbs(x1) + 
                                        bbs(x2) + 
                                        bbs(x3) + 
                                        bbs(x4) + 
                                        bbs(x5) + 
                                        bbs(x6) +
                                        bbs(x7) + 
                                        bbs(x8) + 
                                        bbs(x9) + 
                                        bbs(x10))
  
  univariateContinuousEquation <- list(mu = formula(y2 ~ bbs(x1) + 
                                                      bbs(x2) + 
                                                      bbs(x3) + 
                                                      bbs(x4) + 
                                                      bbs(x5) + 
                                                      bbs(x6) +
                                                      bbs(x7) + 
                                                      bbs(x8) + 
                                                      bbs(x9) + 
                                                      bbs(x10)),
                                       sigma = formula(y2 ~ bbs(x1) + 
                                                         bbs(x2) + 
                                                         bbs(x3) + 
                                                         bbs(x4) + 
                                                         bbs(x5) + 
                                                         bbs(x6) +
                                                         bbs(x7) + 
                                                         bbs(x8) + 
                                                         bbs(x9) + 
                                                         bbs(x10))
  )
  
  
  
  UNIVARIATE_MU1_predictions <- list()
  UNIVARIATE_MU2_predictions <- list()
  UNIVARIATE_SIGMA2_predictions <- list()
  
  mstop.uni <- vector('list')
  coef.uni  <- vector('list')
  coef.uni_1  <- vector('list')
  AUC.uni   <- vector('list')
  Brier.uni <- vector('list')
  
  # - mu1
  dat.train.mu1 = dat.train[, !(names(dat.train) %in% "y2")]
  dat.test.mu1  = dat.test[ , !(names(dat.test) %in% "y2")]
  
  glm.uni.mu1 <- gamboost(as.factor(y1) ~. , data = dat.train.mu1, 
                          family = Binomial(link = "probit", type = "glm"), 
                          control = boost_control(risk = 'oobag', nu = boost.nu, mstop = 1000, trace = T), 
                          weights = weight.mstop)
  
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  if(mstop.uni$mu1 >= 90){
    glm.uni.mu1[1000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 990){
    glm.uni.mu1[5000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 4990){
    glm.uni.mu1[10000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  #dat.train_mu1 <- dat.train.mu1[weight.mstop == 1, ]
  
  #glm.uni.mu1 <- gamboost(as.factor(y1) ~. , data = dat.train_mu1, family = Binomial(link = "probit", type = "glm"), 
  #                       control = boost_control(mstop = mstop.uni$mu1, nu  = 0.01))
  glm.uni.mu1 <- glm.uni.mu1[mstop.uni$mu1]
  
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 1  
  coef.uni$mu1 <- coef(glm.uni.mu1, which = '')
  
  binaryMarginMetrics$Univariate_AUC <- roc(dat.test.mu1$y1, as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response")))$auc
  binaryMarginMetrics$Univariate_BrierScore <- mean((as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1, type = "response"))-(as.numeric(dat.test.mu1$y1)))^2)
  
  ### MU 1
  UNIVARIATE_MU1_predictions$x1 <- predict(glm.uni.mu1, data.frame(newX), which = 1, type = "link" )
  UNIVARIATE_MU1_predictions$x2 <- predict(glm.uni.mu1, data.frame(newX), which = "x2", type = "link" )
  UNIVARIATE_MU1_predictions$x3 <- predict(glm.uni.mu1, data.frame(newX), which = "x3", type = "link" )
  UNIVARIATE_MU1_predictions$x4 <- predict(glm.uni.mu1, data.frame(newX), which = "x4", type = "link" )
  UNIVARIATE_MU1_predictions$x5 <- predict(glm.uni.mu1, data.frame(newX), which = "x5", type = "link" )
  UNIVARIATE_MU1_predictions$x6 <- predict(glm.uni.mu1, data.frame(newX), which = "x6", type = "link" )
  UNIVARIATE_MU1_predictions$x7 <- predict(glm.uni.mu1, data.frame(newX), which = "x7", type = "link" )
  UNIVARIATE_MU1_predictions$x8 <- predict(glm.uni.mu1, data.frame(newX), which = "x8", type = "link" )
  UNIVARIATE_MU1_predictions$x9 <- predict(glm.uni.mu1, data.frame(newX), which = "x9", type = "link" )
  UNIVARIATE_MU1_predictions$x10 <- predict(glm.uni.mu1, data.frame(newX), which = "x10", type = "link" )
  
  
  ### CONTINUOUS MARGIN
  # - mu2
  dat.train.mu2 =  dat.train[, !(names(dat.train) %in% "y1")]
  dat.test.mu2 = dat.test[, !(names(dat.test) %in% "y1")]
  
  glm.uni.mu2 <- gamboostLSS(y2 ~., data = dat.train.mu2, 
                             families = GaussianLSS(stabilization = "L2"), 
                             method = "noncyclic", 
                             control = boost_control(risk = 'oobag', nu = boost.nu, mstop = 1000, trace = TRUE), 
                             weights = weight.mstop)
  
  
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  if(mstop.uni$mu2 >= 90){
    glm.uni.mu2[1000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  if(mstop.uni$mu2 >= 990){
    glm.uni.mu2[10000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  if(mstop.uni$mu2 >= 9990){
    glm.uni.mu2[15000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  #dat.train_mu2 <- dat.train.mu2[weight.mstop == 1, ]
  
  #glm.uni.mu2 <- gamboostLSS(y2 ~., data = dat.train_mu2, families = GaussianLSS(), method = "noncyclic", 
  #                           control = boost_control(mstop = mstop.uni$mu2, nu  = boost.nu, trace = TRUE))
  glm.uni.mu2 <- glm.uni.mu2[mstop.uni$mu2]
  
  
  mstop.uni$sigma2 <- mstop(glm.uni.mu2$sigma)
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 2  
  coef.uni$mu2 <- coef(glm.uni.mu2$mu, which = '')
  coef.uni$sigma2 <- coef(glm.uni.mu2$sigma, which = '')
  
  continuousMarginMetrics$Univariate_MSE <- mean(((predict(glm.uni.mu2$mu, newdata = dat.test.mu2, type = "response"))-as.numeric(dat.test$y2))^2)
  
  
  
  # Likelihood
  pred.mu1.uni <- as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1, type = "response"))
  pred.mu2.uni <- as.numeric(predict(glm.uni.mu2$mu, newdata = dat.test.mu2, type = "response"))
  pred.sigma2.uni <- as.numeric(predict(glm.uni.mu2$sigma, newdata = dat.test.mu2, type = "response"))
  
  mu1.uni.loglik <- sum(-dbinom(x = dat.test.mu1$y1, size = 1, prob = pred.mu1.uni, log = T))
  mu2.uni.loglik <- sum(-dnorm(x = dat.test.mu2$y2, mean = pred.mu2.uni, sd = pred.sigma2.uni, log = T))
  
  lik$uni_usingCopula <- sum(loss(mu1 = pred.mu1.uni, mu2 = pred.mu2.uni, sigma2 = pred.sigma2.uni, rho = rep((.Machine$double.eps)^2, length(nrow(y.test))), y = y.test))
  
  lik$uni<- mu1.uni.loglik + mu2.uni.loglik
  
  
  n.ges = c(n.train, n.mstop, n.test)
  pred.ges <- NULL#list(pred.mu1, pred.mu2, pred.or)
  pred.uni <- list(mu1 = pred.mu1.uni, mu2 = pred.mu2.uni, sigma2 = pred.sigma2.uni)
  pred.cop <- list(mu1 = predCopula.mu1, mu2 = predCopula.mu2, sigma2 = predCopula.sigma2, rho = predCopula.rho)
  
  
  
  #### MU 2
  UNIVARIATE_MU2_predictions$x1 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = 1, type = "link" )
  UNIVARIATE_MU2_predictions$x2 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x2", type = "link" )
  UNIVARIATE_MU2_predictions$x3 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x3", type = "link" )
  UNIVARIATE_MU2_predictions$x4 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x4", type = "link" )
  UNIVARIATE_MU2_predictions$x5 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x5", type = "link" )
  UNIVARIATE_MU2_predictions$x6 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x6", type = "link" )
  UNIVARIATE_MU2_predictions$x7 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x7", type = "link" )
  UNIVARIATE_MU2_predictions$x8 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x8", type = "link" )
  UNIVARIATE_MU2_predictions$x9 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x9", type = "link" )
  UNIVARIATE_MU2_predictions$x10 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x10", type = "link" )
  
  #### SIGMA 2
  UNIVARIATE_SIGMA2_predictions$x1 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = 1, type = "link" )
  UNIVARIATE_SIGMA2_predictions$x2 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x2", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x3 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x3", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x4 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x4", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x5 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x5", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x6 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x6", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x7 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x7", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x8 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x8", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x9 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x9", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x10 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x10", type = "link" )
  
  
  
  ###############################################################################################################################
  ######## Energy Score ########################################################################################################
  ##############################################################################################################################
  es_biv <- vector()
  es_uni<- vector()
  es_cop <- vector()
 
  
  for(i in 1:length(pred.mu1.uni)){
    
    pred_sample_uni <- matrix(NA, nrow = 2, ncol = 1000)
    pred_sample_biv <- matrix(NA, nrow = 2, ncol = 1000)
    pred_sample_cop <- matrix(NA, nrow = 2, ncol = 1000)
    
    
    
    # univariate
    #sample_uni <- rbinom2.or(10000, mu1 =  pred.mu1.uni[i], mu2 = pred.mu2.uni[i], oratio = 1)
    sample_uni <- data.gen.bivbin_energyscore(n = 1000, 
                                              mu1 = pred.mu1.uni[i], 
                                              mu2 = pred.mu2.uni[i], 
                                              sigma2 = pred.sigma2.uni[i], 
                                              FAM = 5, 
                                              theta = sqrt(.Machine$double.eps),
                                              from_univariate = TRUE)
    
    pred_sample_uni[1,] <- sample_uni[,1]
    pred_sample_uni[2,] <- sample_uni[,2]
    
    es_uni[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_uni) 
    
    
    
    
    # copula
    sample_cop <- data.gen.bivbin_energyscore(n = 1000, 
                                              mu1 = predCopula.mu1[i],  
                                              mu2 = predCopula.mu2[i], 
                                              sigma2 = predCopula.sigma2[i], 
                                              FAM = 5, 
                                              theta = predCopula.rho[i])
    
    pred_sample_cop[1, ] <- sample_cop[,1]
    pred_sample_cop[2, ] <- sample_cop[,2]
    
    es_cop[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_cop) 
    
    
  }
  
  energy_score <- list()
  #energy_score$biv <- mean(es_biv)
  energy_score$uni <- mean(es_uni)
  energy_score$cop <- mean(es_cop)
  
  
  ##################################### Collect: mstop, predictions, 
  #                                             AUC for binary margin, Brier score for binary margin, 
  #                                             MSE for continuous margin
  #                                             energy score and negative log likelihood 
  output <- list(TrueBeta = TrueBeta, 
                 n = n.ges, 
                 p  = p, 
                 Likelihood = lik,
                 predict = pred.ges, 
                 predict.uni = pred.uni, 
                 predict.cop = pred.cop,
                 energy_score = energy_score,  
                 ####
                 Coefficients.uni = coef.uni, 
                 mstop.uni = mstop.uni,
                 ################################################
                 BinaryMarginMetrics = binaryMarginMetrics,
                 ContinuousMarginMetrics = continuousMarginMetrics,
                 ################################################
                 CoefficientsCOPULA = coef.bivBinCopula, 
                 oobag.riskCOPULA = oobag.risk.cop, 
                 mstopCOPULA = mstop.bivBinCopula, 
                 ###
                 Copula_Predictions = list(MU1 = Copula_MU1_predictions, 
                                           MU2 = Copula_MU2_predictions,
                                           SIGMA2 = Copula_SIGMA2_predictions,
                                           RHO = Copula_RHO_predictions),
                 Univariate_Predictions = list(MU1 = UNIVARIATE_MU1_predictions, 
                                               MU2 = UNIVARIATE_MU2_predictions, 
                                               SIGMA2 = UNIVARIATE_SIGMA2_predictions),
                 KendallRange = TrueKendallRange
  )
  
  
  
  return(output)
}

#### ARCHIMEDEAN COPULAE WITH ONE-DIRECTIONAL DEPENDENCE.
sim_CLAYTON180 <- function(seed, n.train, n.mstop, p, n.test, corr, boost.nu = 0.1){
  
  
  
  loss <- function(mu1, mu2, sigma2, rho, y){
    
    F1 <- 1 - mu1 # F(0)
    
    F2 <- 1 - pnorm(y[,2], mean = mu2, sd = sigma2)
    
    dens2 <- dnorm(y[,2], mean = mu2, sd = sigma2)
    
    thet <- rho
    
    HFunc <-  pdffz( (F1^(-thet) + F2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (F2^((-thet) - 1) * (-thet))) )
    
    return( -( (1 - y[,1]) * log( HFunc ) + y[,1] * log( pdffz(1 - HFunc) ) + log(dens2) ) )
    
  }
  
  loss_VC <- function(mu1, mu2, sigma2, rho, y, FAM = 13){
    
    F1 <- 1 - mu1 # F(0)
    
    F2 <- pnorm(y[,2], mean = mu2, sd = sigma2)
    
    dens2 <- dnorm(y[,2], mean = mu2, sd = sigma2)
    
    thet <- rho
    
    HFunc <-  pdffz( VineCopula:::BiCopHfunc2(u1 = F1, u2 = F2, family = FAM, par = thet) ) #pdffz( (F1^(-thet) + F2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (F2^((-thet) - 1) * (-thet))) )
    
    return( -( (1 - y[,1]) * log( HFunc ) + y[,1] * log( pdffz(1 - HFunc) ) + log(dens2) ) )
    
  }
  
  
  # SAMPLING OF BIVARIATE RESPONSE FROM COPULA
  # simulate data from bivariate binary DGP:
  data.gen.bivbin <- function(FAM, mu1, mu2, sigma2, theta){
    
    
    u1u2 <- VineCopula::BiCopSim(length(mu1), family = FAM, par = theta)
    
    y1 <- qbinom(p = u1u2[,1], size = 1, prob = mu1)
    
    y2 <- qnorm(p = u1u2[,2], mean = mu2, sd = sigma2)
    
    dat <- data.frame(y1, y2)
    
    return(dat)
  }
  
  data.gen.bivbin_energyscore <- function(n, FAM, mu1, mu2, sigma2, theta, from_univariate = FALSE){
    
    #if(!from_univariate){
    
    u1u2 <- VineCopula::BiCopSim(n, family = FAM, par = theta)
    
    y1 <- qbinom(p = u1u2[,1], size = 1, prob = mu1)
    
    y2 <- qnorm(p = u1u2[,2],  mean = mu2, sd = sigma2)
    
    #}else{
    
    #y1 <- rbinom(n, size = 1, prob = mu1)
    
    #y2 <- rnorm(n,  mean = mu2, sd = sigma2)
    
    #}
    
    
    dat <- data.frame(y1, y2)
    
    return(dat)
  }
  
  
  # 1.5x2 − 1x3 + 1.5x4 − √x1ix1i,
  #mu1 <- c(1, 1.5, -1, 1.5, 0, 0, rep(0,p-6))
  #
  # x1, x2 ,x3 x4, x5, x6, x7, x8, x9, x10
  #mu1 <- c( 0, 1.5, -1, 1.5, 0, 0, 0, 0, 0, 0)
  mu1 <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  # 0.5x2 + 1.5x3
  #mu2 <- c(2, -1, 1, 0, 0, rep(0,p-5)) 
  #
  mu2 <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  
  # −0.5 + cos(2 * x2)
  #
  #sigma2 <- c(2, -1, 1, 0, 0, rep(0,p-5)) 
  sigma2 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  
  # −1.5 + 1x5 − 1.5x6 + sin(x3)
  #
  #or  <- c(0, 0, 1, 1.5, 0, 0, 0, 0,rep(0,p-8))
  or <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  
  nl_function_mu1 <- function(x){
    
    
    nl <- sqrt(x)*x^2 - 2*cos(3*x)
    
    return(nl*0.5)
    
  }
  
  nl_function_mu2 <- function(x){
    
    nl <- -0.7*exp(x^2) + exp(x^0.4) - 0*log(x+0.2) 
    
    return(nl)
  }
  
  nl_function_sigma2 <- function(x){
    
    return(cos(2 * x))
  }
  
  nl_function_rho <- function(x){
    
    return(sin(x*4))
  }
  
  TrueBeta <-  vector('list')
  TrueBeta$mu1 <- mu1
  TrueBeta$mu2 <- mu2
  TrueBeta$sigma2 <- sigma2
  TrueBeta$or <- or
  
  set.seed(seed)
  
  #n.mstop = 1500
  n <- n.train + n.mstop
  weight.mstop <- c(rep(1, times = n.train),rep(0, times = n.mstop)) 
  
  # - training data
  #x.train <- rmvnorm(n = n, mean = rep(0,p), sigma =  toeplitz(sapply(seq(0,p-1), function(x) corr^x)))
  x1  <- runif(n, 0, 1)
  x2 <- runif(n, 0, 1)
  x3 <- runif(n, 0, 1)
  x4 <- runif(n, 0, 1)
  x5 <- runif(n, 0, 1)
  x6 <- runif(n, 0, 1)
  x7 <- runif(n, 0, 1)
  x8 <- runif(n, 0, 1)
  x9 <- runif(n, 0, 1) 
  x10 <- runif(n, 0, 1)
  
  x.train <- cbind(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
  
  #x.test <- rmvnorm(n = n.test, mean = rep(0,p), sigma = toeplitz(sapply(seq(0,p-1), function(x) corr^x)))
  x1.test <- runif(n.test, 0, 1)
  x2.test <- runif(n.test, 0, 1)
  x3.test <- runif(n.test, 0, 1)
  x4.test <- runif(n.test, 0, 1)
  x5.test <- runif(n.test, 0, 1)
  x6.test <- runif(n.test, 0, 1)
  x7.test <- runif(n.test, 0, 1)
  x8.test <- runif(n.test, 0, 1)
  x9.test <- runif(n.test, 0, 1) 
  x10.test <- runif(n.test, 0, 1)
  
  x.test <- cbind(x1.test, x2.test, x3.test, x4.test, x5.test, x6.test, x7.test, x8.test, x9.test, x10.test)
  
  colnames(x.test) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
  
  
  ##########################################################################################################################################
  # GENERATION OF DISTRIBUTION PARAMETERS AND STUFF 
  # predictor
  train.eta.mu1_center <-  (x.train %*% mu1)
  test.eta.mu1 <-  (x.test %*% mu1)
  
  ### Add non-linear effects (additive): − √x1 x1,
  train.eta.mu1_center  <- train.eta.mu1_center + -1*0 + nl_function_mu1(x3)
  test.eta.mu1_center   <- test.eta.mu1 + -1*0 + nl_function_mu1(x3.test)
  
  # apply response function: PROBIT
  train.eta.mu1 <-  pnorm(train.eta.mu1_center) 
  test.eta.mu1  <-  pnorm(test.eta.mu1_center) 
  
  ##########################################################################################################################################
  # predictor and apply response function: IDENTITY
  train.eta.mu2 <-  x.train %*% mu2 + nl_function_mu2(x1)
  test.eta.mu2  <-  x.test %*% mu2
  test.eta.mu2  <- test.eta.mu2 + 0 + nl_function_mu2(x1.test)
  
  # no non-linear effect here
  train.eta.mu2 <- train.eta.mu2 + 0
  train.eta.mu2 <-  train.eta.mu2 # Identity link
  test.eta.mu2  <-  test.eta.mu2 
  
  
  ##########################################################################################################################################
  #### SIGMA FOR CONTINOUS RESPONSE
  train.eta.sigma2 <-  x.train %*% sigma2
  test.eta.sigma2 <-  x.test %*% sigma2
  
  train.eta.sigma2  <- train.eta.sigma2 + nl_function_sigma2(x2) - 0.5
  test.eta.sigma2   <- test.eta.sigma2 + nl_function_sigma2(x2.test) - 0.5
  
  train.eta.sigma2  <- exp(train.eta.sigma2)
  test.eta.sigma2   <- exp(test.eta.sigma2)
  
  
  ##########################################################################################################################################
  ### COPULA PARAMETER FOR ROTATED CLAYTON COPULA (180 DEGREES)
  train.eta.or_center <- (x.train %*% or)
  test.eta.or_center <- (x.test %*% or)
  
  train.eta.or_center <- train.eta.or_center + - 1 + 3*nl_function_rho(x3)
  test.eta.or_center <- test.eta.or_center + - 1 + 3*nl_function_rho(x3.test)
  
  train.eta.or_center <- exp(train.eta.or_center) 
  test.eta.or_center  <- exp(test.eta.or_center) 
  
  train.copula.parameter <- train.eta.or_center
  test.copula.parameter <- test.eta.or_center
  
  
  TrueKendallTauRange <- VineCopula::BiCopPar2Tau(family = 13, par = range(train.copula.parameter))
  
  ##########################################################################################################################################
  ##########################################################################################################################################
  
  y.train <- data.gen.bivbin(FAM = 13, mu1 = train.eta.mu1, mu2 = train.eta.mu2, sigma2 = train.eta.sigma2, theta = train.copula.parameter)
  # SAMPLE TEST OBSERVATIONS
  y.test  <- data.gen.bivbin(FAM = 13, mu1 = test.eta.mu1, mu2 = test.eta.mu2, sigma2 = test.eta.sigma2, theta = test.copula.parameter)
  
  
  colnames(y.train)<- c('y1','y2')
  dat.train <- data.frame(y.train,x.train)
  
  colnames(y.test)<- c('y1','y2')
  dat.test <- data.frame(y.test,x.test)
  
  ##########################################################################################################################################
  ##########################################################################################################################################
  
  mstop.bivBern <-  vector('list')
  
  
  AUCbivBern <-  vector('list')
  BrierbivBern <-  vector('list')
  lik <- vector('list')
  
  
  # # Predict distributional quantities
  model_equation <- formula(cbind(y1,y2) ~ bbs(x1, knots = 20, degree = 3, difference = 2) + 
                              bbs(x2, knots = 20, degree = 3, difference = 2) + 
                              bbs(x3, knots = 20, degree = 3, difference = 2) + 
                              bbs(x4, knots = 20, degree = 3, difference = 2) + 
                              bbs(x5, knots = 20, degree = 3, difference = 2) + 
                              bbs(x6, knots = 20, degree = 3, difference = 2) +
                              bbs(x7, knots = 20, degree = 3, difference = 2) + 
                              bbs(x8, knots = 20, degree = 3, difference = 2) + 
                              bbs(x9, knots = 20, degree = 3, difference = 2) + 
                              bbs(x10, knots = 20, degree = 3, difference = 2))
  
  model_equation <- formula(cbind(y1,y2) ~ bbs(x1) + 
                              bbs(x2) + 
                              bbs(x3) + 
                              bbs(x4) + 
                              bbs(x5) + 
                              bbs(x6) +
                              bbs(x7) + 
                              bbs(x8) + 
                              bbs(x9) + 
                              bbs(x10))
  
  bivModel_Formula <- list(mu1 = model_equation,
                           mu2 = model_equation,
                           sigma2 = model_equation,
                           rho = model_equation
  )
  
  #############################################
  ################# - COPULA - ################
  #############################################
  # Bivariate copula model
  bivBinCopula <- gamboostLSS(bivModel_Formula, data = dat.train, 
                              families = Clayton180_Cop_BinCont(marg1 = "PROBIT", 
                                                                marg2 = "NORM", 
                                                                stabilization = "L2"), 
                              control = boost_control(mstop = 2000, risk = 'oobag', nu = boost.nu, trace = TRUE), 
                              method = 'noncyclic', weights = weight.mstop)
  
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 1990){
    mstop(bivBinCopula) <- 4000
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 3990){
    bivBinCopula[8000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 7990){
    bivBinCopula[10000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 9990){
    bivBinCopula[15000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 14990){
    bivBinCopula[20000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula, merge = T))
  
  if(MSTOP_COP >= 19990){
    bivBinCopula[25000]
  }
  
  
  
  MSTOP_COP <- which.min(risk(bivBinCopula, merge = T))
  oobag.risk.cop <- risk(bivBinCopula,merge = T)
  
  #rm(bivBinCopula)
  #dat.train_biv <- dat.train[weight.mstop == 1, ]
  
  # RE-FIT until OPTIMAL MSTOP
  #bivBinCopula <- gamboostLSS(bivModel_Formula, data = dat.train_biv, 
  #                           families = Clayton270_Cop_BinCont(marg1 = "PROBIT", marg2 = "NORM"), 
  #                           control = boost_control(mstop = MSTOP_COP, nu  = boost.nu, trace = TRUE), method = 'noncyclic')
  bivBinCopula <- bivBinCopula[MSTOP_COP]
  
  
  mstop.bivBinCopula <-  vector('list')
  mstop.bivBinCopula$mstop <- MSTOP_COP
  mstop.bivBinCopula$mu1 <- bivBinCopula$mu1$mstop()
  mstop.bivBinCopula$mu2 <- bivBinCopula$mu2$mstop()
  mstop.bivBinCopula$sigma2 <- bivBinCopula$sigma2$mstop()
  mstop.bivBinCopula$rho  <- bivBinCopula$rho$mstop()
  
  # EXTRACT ALL COEFFICIENTS!
  coef.bivBinCopula <- coef(bivBinCopula, which = "")
  
  binaryMarginMetrics <-  vector('list')
  continuousMarginMetrics <-  vector('list')
  likBinCopula <- vector('list')
  
  
  ################ Get predictions from each covariate?
  newX <- matrix(rep(seq(from = 0, to = 1,length.out = 500),times = 10), nrow = 500)
  colnames(newX) <-  c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
  
  
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ### Copula predictions / selection
  Copula_MU1_predictions <- list()
  Copula_MU2_predictions <- list()
  Copula_SIGMA2_predictions <- list()
  Copula_RHO_predictions <- list()
  
  ### MU 1
  Copula_MU1_predictions$x1 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = 1, type = "link" )
  Copula_MU1_predictions$x2 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x2", type = "link" )
  Copula_MU1_predictions$x3 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x3", type = "link" )
  Copula_MU1_predictions$x4 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x4", type = "link" )
  Copula_MU1_predictions$x5 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x5", type = "link" )
  Copula_MU1_predictions$x6 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x6", type = "link" )
  Copula_MU1_predictions$x7 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x7", type = "link" )
  Copula_MU1_predictions$x8 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x8", type = "link" )
  Copula_MU1_predictions$x9 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x9", type = "link" )
  Copula_MU1_predictions$x10 <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = "x10", type = "link" )
  
  
  #### MU 2
  Copula_MU2_predictions$x1 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = 1, type = "link" )
  Copula_MU2_predictions$x2 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x2", type = "link" )
  Copula_MU2_predictions$x3 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x3", type = "link" )
  Copula_MU2_predictions$x4 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x4", type = "link" )
  Copula_MU2_predictions$x5 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x5", type = "link" )
  Copula_MU2_predictions$x6 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x6", type = "link" )
  Copula_MU2_predictions$x7 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x7", type = "link" )
  Copula_MU2_predictions$x8 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x8", type = "link" )
  Copula_MU2_predictions$x9 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x9", type = "link" )
  Copula_MU2_predictions$x10 <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = "x10", type = "link" )
  
  #### SIGMA 2
  Copula_SIGMA2_predictions$x1 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = 1, type = "link" )
  Copula_SIGMA2_predictions$x2 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x2", type = "link" )
  Copula_SIGMA2_predictions$x3 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x3", type = "link" )
  Copula_SIGMA2_predictions$x4 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x4", type = "link" )
  Copula_SIGMA2_predictions$x5 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x5", type = "link" )
  Copula_SIGMA2_predictions$x6 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x6", type = "link" )
  Copula_SIGMA2_predictions$x7 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x7", type = "link" )
  Copula_SIGMA2_predictions$x8 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x8", type = "link" )
  Copula_SIGMA2_predictions$x9 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x9", type = "link" )
  Copula_SIGMA2_predictions$x10 <- predict(bivBinCopula, data.frame(newX), parameter = "sigma2", which = "x10", type = "link" )
  
  
  #### RHO
  Copula_RHO_predictions$x1 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = 1, type = "link" )
  Copula_RHO_predictions$x2 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x2", type = "link" )
  Copula_RHO_predictions$x3 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x3", type = "link" )
  Copula_RHO_predictions$x4 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x4", type = "link" )
  Copula_RHO_predictions$x5 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x5", type = "link" )
  Copula_RHO_predictions$x6 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x6", type = "link" )
  Copula_RHO_predictions$x7 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x7", type = "link" )
  Copula_RHO_predictions$x8 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x8", type = "link" )
  Copula_RHO_predictions$x9 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x9", type = "link" )
  Copula_RHO_predictions$x10 <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = "x10", type = "link" )
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  
  # Predict distributional quantities
  predCopula.mu1 <- predict(bivBinCopula$mu1, newdata = dat.test, type = 'response')
  predCopula.mu2 <- predict(bivBinCopula$mu2, newdata = dat.test, type = 'response')
  predCopula.sigma2 <- predict(bivBinCopula$sigma2, newdata = dat.test, type = 'response')
  predCopula.rho <- predict(bivBinCopula$rho, newdata = dat.test, type = 'response')
  
  ######### COMPUTE PERFORMANCE: 
  # BINARY MARGIN
  # AUC (margin 1)
  binaryMarginMetrics$AUC <- roc(dat.test$y1, as.numeric(predCopula.mu1))$auc
  # Brier score (margin 1, then margin 2)
  binaryMarginMetrics$BrierScore <- mean((as.numeric(predCopula.mu1)-(as.numeric(dat.test$y1)))^2)
  
  ### CONTINUOUS MARGIN: 
  continuousMarginMetrics$MSE <- mean((as.numeric(predCopula.mu2)-as.numeric(dat.test$y2))^2)
  
  
  # LOG-SCORE (LOG-LIKELIHOOD) (change to the copula likelihood)
  lik$biv <- sum(loss_VC(mu1 = predCopula.mu1, mu2 = predCopula.mu2, sigma2 = predCopula.sigma2, rho = predCopula.rho, y = y.test))
  
  
  ###########################################################################################################################
  ################ - univariate fits - ######################################################################################
  ###########################################################################################################################
  
  univariateBinaryEquation <- formula(as.factor(y1) ~  bbs(x1) + 
                                        bbs(x2) + 
                                        bbs(x3) + 
                                        bbs(x4) + 
                                        bbs(x5) + 
                                        bbs(x6) +
                                        bbs(x7) + 
                                        bbs(x8) + 
                                        bbs(x9) + 
                                        bbs(x10))
  
  univariateContinuousEquation <- list(mu = formula(y2 ~ bbs(x1) + 
                                                      bbs(x2) + 
                                                      bbs(x3) + 
                                                      bbs(x4) + 
                                                      bbs(x5) + 
                                                      bbs(x6) +
                                                      bbs(x7) + 
                                                      bbs(x8) + 
                                                      bbs(x9) + 
                                                      bbs(x10)),
                                       sigma = formula(y2 ~ bbs(x1) + 
                                                         bbs(x2) + 
                                                         bbs(x3) + 
                                                         bbs(x4) + 
                                                         bbs(x5) + 
                                                         bbs(x6) +
                                                         bbs(x7) + 
                                                         bbs(x8) + 
                                                         bbs(x9) + 
                                                         bbs(x10))
  )
  
  
  
  UNIVARIATE_MU1_predictions <- list()
  UNIVARIATE_MU2_predictions <- list()
  UNIVARIATE_SIGMA2_predictions <- list()
  
  mstop.uni <- vector('list')
  coef.uni  <- vector('list')
  coef.uni_1  <- vector('list')
  AUC.uni   <- vector('list')
  Brier.uni <- vector('list')
  
  # - mu1
  dat.train.mu1 = dat.train[, !(names(dat.train) %in% "y2")]
  dat.test.mu1  = dat.test[ , !(names(dat.test) %in% "y2")]
  
  glm.uni.mu1 <- gamboost(as.factor(y1) ~. , data = dat.train.mu1, 
                          family = Binomial(link = "probit", type = "glm"), 
                          control = boost_control(risk = 'oobag', nu = boost.nu, mstop = 1000, trace = T), 
                          weights = weight.mstop)
  
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  if(mstop.uni$mu1 >= 90){
    glm.uni.mu1[1000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 990){
    glm.uni.mu1[5000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 4990){
    glm.uni.mu1[10000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  #dat.train_mu1 <- dat.train.mu1[weight.mstop == 1, ]
  
  #glm.uni.mu1 <- gamboost(as.factor(y1) ~. , data = dat.train_mu1, family = Binomial(link = "probit", type = "glm"), 
  #                       control = boost_control(mstop = mstop.uni$mu1, nu  = 0.01))
  glm.uni.mu1 <- glm.uni.mu1[mstop.uni$mu1]
  
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 1  
  coef.uni$mu1 <- coef(glm.uni.mu1, which = '')
  
  binaryMarginMetrics$Univariate_AUC <- roc(dat.test.mu1$y1, as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response")))$auc
  binaryMarginMetrics$Univariate_BrierScore <- mean((as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1, type = "response"))-(as.numeric(dat.test.mu1$y1)))^2)
  
  ### MU 1
  UNIVARIATE_MU1_predictions$x1 <- predict(glm.uni.mu1, data.frame(newX), which = 1, type = "link" )
  UNIVARIATE_MU1_predictions$x2 <- predict(glm.uni.mu1, data.frame(newX), which = "x2", type = "link" )
  UNIVARIATE_MU1_predictions$x3 <- predict(glm.uni.mu1, data.frame(newX), which = "x3", type = "link" )
  UNIVARIATE_MU1_predictions$x4 <- predict(glm.uni.mu1, data.frame(newX), which = "x4", type = "link" )
  UNIVARIATE_MU1_predictions$x5 <- predict(glm.uni.mu1, data.frame(newX), which = "x5", type = "link" )
  UNIVARIATE_MU1_predictions$x6 <- predict(glm.uni.mu1, data.frame(newX), which = "x6", type = "link" )
  UNIVARIATE_MU1_predictions$x7 <- predict(glm.uni.mu1, data.frame(newX), which = "x7", type = "link" )
  UNIVARIATE_MU1_predictions$x8 <- predict(glm.uni.mu1, data.frame(newX), which = "x8", type = "link" )
  UNIVARIATE_MU1_predictions$x9 <- predict(glm.uni.mu1, data.frame(newX), which = "x9", type = "link" )
  UNIVARIATE_MU1_predictions$x10 <- predict(glm.uni.mu1, data.frame(newX), which = "x10", type = "link" )
  
  
  ### CONTINUOUS MARGIN
  # - mu2
  dat.train.mu2 =  dat.train[, !(names(dat.train) %in% "y1")]
  dat.test.mu2 = dat.test[, !(names(dat.test) %in% "y1")]
  
  glm.uni.mu2 <- gamboostLSS(y2 ~., data = dat.train.mu2, 
                             families = GaussianLSS(stabilization = "L2"), 
                             method = "noncyclic", 
                             control = boost_control(risk = 'oobag', nu = boost.nu, mstop = 1000, trace = TRUE), 
                             weights = weight.mstop)
  
  
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  if(mstop.uni$mu2 >= 90){
    glm.uni.mu2[1000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  if(mstop.uni$mu2 >= 990){
    glm.uni.mu2[10000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  if(mstop.uni$mu2 >= 9990){
    glm.uni.mu2[15000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  #dat.train_mu2 <- dat.train.mu2[weight.mstop == 1, ]
  
  #glm.uni.mu2 <- gamboostLSS(y2 ~., data = dat.train_mu2, families = GaussianLSS(), method = "noncyclic", 
  #                           control = boost_control(mstop = mstop.uni$mu2, nu  = boost.nu, trace = TRUE))
  glm.uni.mu2 <- glm.uni.mu2[mstop.uni$mu2]
  
  
  mstop.uni$sigma2 <- mstop(glm.uni.mu2$sigma)
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 2  
  coef.uni$mu2 <- coef(glm.uni.mu2$mu, which = '')
  coef.uni$sigma2 <- coef(glm.uni.mu2$sigma, which = '')
  
  continuousMarginMetrics$Univariate_MSE <- mean(((predict(glm.uni.mu2$mu, newdata = dat.test.mu2, type = "response"))-as.numeric(dat.test$y2))^2)
  
  
  
  # Likelihood
  pred.mu1.uni <- as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1, type = "response"))
  pred.mu2.uni <- as.numeric(predict(glm.uni.mu2$mu, newdata = dat.test.mu2, type = "response"))
  pred.sigma2.uni <- as.numeric(predict(glm.uni.mu2$sigma, newdata = dat.test.mu2, type = "response"))
  
  mu1.uni.loglik <- sum(-dbinom(x = dat.test.mu1$y1, size = 1, prob = pred.mu1.uni, log = T))
  mu2.uni.loglik <- sum(-dnorm(x = dat.test.mu2$y2, mean = pred.mu2.uni, sd = pred.sigma2.uni, log = T))
  
  lik$uni_usingCopula <- sum(loss_VC(mu1 = pred.mu1.uni, 
                                     mu2 = pred.mu2.uni, 
                                     sigma2 = pred.sigma2.uni, 
                                     rho = rep(sqrt(.Machine$double.eps), length(nrow(y.test))), y = y.test))
  
  lik$uni<- mu1.uni.loglik + mu2.uni.loglik
  
  
  n.ges = c(n.train, n.mstop, n.test)
  pred.ges <- NULL#list(pred.mu1, pred.mu2, pred.or)
  pred.uni <- list(mu1 = pred.mu1.uni, mu2 = pred.mu2.uni, sigma2 = pred.sigma2.uni)
  pred.cop <- list(mu1 = predCopula.mu1, mu2 = predCopula.mu2, sigma2 = predCopula.sigma2, rho = predCopula.rho)
  
  
  
  #### MU 2
  UNIVARIATE_MU2_predictions$x1 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = 1, type = "link" )
  UNIVARIATE_MU2_predictions$x2 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x2", type = "link" )
  UNIVARIATE_MU2_predictions$x3 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x3", type = "link" )
  UNIVARIATE_MU2_predictions$x4 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x4", type = "link" )
  UNIVARIATE_MU2_predictions$x5 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x5", type = "link" )
  UNIVARIATE_MU2_predictions$x6 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x6", type = "link" )
  UNIVARIATE_MU2_predictions$x7 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x7", type = "link" )
  UNIVARIATE_MU2_predictions$x8 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x8", type = "link" )
  UNIVARIATE_MU2_predictions$x9 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x9", type = "link" )
  UNIVARIATE_MU2_predictions$x10 <- predict(glm.uni.mu2, data.frame(newX), parameter = "mu", which = "x10", type = "link" )
  
  #### SIGMA 2
  UNIVARIATE_SIGMA2_predictions$x1 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = 1, type = "link" )
  UNIVARIATE_SIGMA2_predictions$x2 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x2", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x3 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x3", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x4 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x4", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x5 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x5", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x6 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x6", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x7 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x7", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x8 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x8", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x9 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x9", type = "link" )
  UNIVARIATE_SIGMA2_predictions$x10 <- predict(glm.uni.mu2, data.frame(newX), parameter = "sigma", which = "x10", type = "link" )
  
  
  
  ###############################################################################################################################
  ######## Energy Score ########################################################################################################
  ##############################################################################################################################
  es_biv <- vector()
  es_uni<- vector()
  es_cop <- vector()
  
  for(i in 1:length(pred.mu1.uni)){
    
    pred_sample_uni <- matrix(NA, nrow = 2, ncol = 1000)
    pred_sample_biv <- matrix(NA, nrow = 2, ncol = 1000)
    pred_sample_cop <- matrix(NA, nrow = 2, ncol = 1000)
    
    
    
    # univariate
    sample_uni <- data.gen.bivbin_energyscore(n = 1000, 
                                              mu1 = pred.mu1.uni[i], 
                                              mu2 = pred.mu2.uni[i], 
                                              sigma2 = pred.sigma2.uni[i], 
                                              FAM = 13, 
                                              theta = sqrt(.Machine$double.eps),
                                              from_univariate = TRUE)
    
    pred_sample_uni[1,] <- sample_uni[,1]
    pred_sample_uni[2,] <- sample_uni[,2]
    
    es_uni[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_uni) 
    
    
    
    
    # copula
    sample_cop <- data.gen.bivbin_energyscore(n = 1000, 
                                              mu1 = predCopula.mu1[i],  
                                              mu2 = predCopula.mu2[i], 
                                              sigma2 = predCopula.sigma2[i], 
                                              FAM = 13, 
                                              theta = predCopula.rho[i])
    
    pred_sample_cop[1, ] <- sample_cop[,1]
    pred_sample_cop[2, ] <- sample_cop[,2]
    
    es_cop[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_cop) 
    
    
  }
  
  energy_score <- list()
  energy_score$uni <- mean(es_uni)
  energy_score$cop <- mean(es_cop)
  
  
  ##################################### Collect: mstop, predictions, 
  #                                             AUC for binary margin, Brier score for binary margin, 
  #                                             MSE for continuous margin
  #                                             energy score and negative log likelihood 
  output <- list(TrueBeta = TrueBeta, 
                 n = n.ges, 
                 p  = p, 
                 Likelihood = lik,
                 predict = pred.ges, 
                 predict.uni = pred.uni, 
                 predict.cop = pred.cop,
                 energy_score = energy_score,  
                 ####
                 Coefficients.uni = coef.uni, 
                 mstop.uni = mstop.uni,
                 ################################################
                 BinaryMarginMetrics = binaryMarginMetrics,
                 ContinuousMarginMetrics = continuousMarginMetrics,
                 ################################################
                 CoefficientsCOPULA = coef.bivBinCopula, 
                 oobag.riskCOPULA = oobag.risk.cop, 
                 mstopCOPULA = mstop.bivBinCopula, 
                 ###
                 Copula_Predictions = list(MU1 = Copula_MU1_predictions, 
                                           MU2 = Copula_MU2_predictions,
                                           SIGMA2 = Copula_SIGMA2_predictions,
                                           RHO = Copula_RHO_predictions),
                 Univariate_Predictions = list(MU1 = UNIVARIATE_MU1_predictions, 
                                               MU2 = UNIVARIATE_MU2_predictions, 
                                               SIGMA2 = UNIVARIATE_SIGMA2_predictions),
                 KendallRange = TrueKendallTauRange
  )
  
  
  
  return(output)
}