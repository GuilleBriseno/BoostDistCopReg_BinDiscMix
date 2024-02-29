#########################################
library("VGAM")
library("gamboostLSS")
library("mvtnorm")
library("pROC")
library("scoringRules")
library("GJRM")
library("copula")
#setwd("BOOSTCOPFILES/")
#########################################
### load Copulas
#source("Copulas/BivariateDiscrete/Copula_Gaussian_BivDisc_p1m1.R")
#source("Copulas/BivariateDiscrete/Copula_FGM_BivDisc_p1m1.R")
source("Copulas/BivariateDiscrete/Copula_FGM2_BivDisc_p1m1.R")
source("Copulas/BivariateDiscrete/Copula_Joe_BivDisc_p1m1.R")
source("Copulas/BivariateDiscrete/Copula_Joe180_BivDisc_p1m1.R")
source("Copulas/BivariateDiscrete/Copula_Frank_BivDisc_p1m1.R")

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

# some checks for the copula parameter / kendalls tau
check_extremes <- function(input){
  
  inp <- ifelse(input == 1, 0.99999999, input)
  inp <- ifelse(inp == -1, -0.99999999, inp)
  
  return(inp)
  
}

check_extremes_nuZINBI <- function(input){
  
  inp <- ifelse(input == 1, 0.99999999, input)
  inp <- ifelse(inp == 0, 0.00000001, inp)
  
  return(inp)
  
}

#### USING FGM COPULA
sim_FGM <- function(seed, n.train, n.mstop, p, n.test, corr, boost.nu.steplength = 0.1){
  
  
  loss <- function(mu1, mu2, sigma1, sigma2, rho, y, nu2 = NULL){
    
    
    # Margin 1 is ZALG:
    F1 <- pdffz(pZALG(y[,1], mu = mu1, sigma = sigma1))
    
    pdf1 <- pdffz(dZALG(y[,1], mu = mu1, sigma = sigma1))
    
    
    # Margin 2 is NBI:
    # F2 <- pdffz(pNBI(y[,2], mu = mu2, sigma = sigma2))
    # 
    # pdf2 <- pdffz(dNBI(y[,2], mu = mu2, sigma = sigma2))
    F2 <- pdffz(pZINBI(y[,2], mu = mu2, sigma = sigma2, nu = nu2))
    
    pdf2 <- pdffz(dZINBI(y[,2], mu = mu2, sigma = sigma2, nu = nu2))
    
    # Copula parameter (same range as Gauss copula)
    thet <- check_extremes(rho) # FGM copula
    
    
    # Minus 1
    CDF1m1 <- pdffz(F1 - pdf1)
    
    CDF2m1 <- pdffz(F2 - pdf2)
    
    ### Copula Terms
    T1 <- pdffz(F1 * F2 * (1 + thet * (1 - F1) * (1 - F2)))
    
    T2 <- pdffz(CDF1m1 * F2 * (1 + thet * (1 - CDF1m1) * (1 - F2)))
    
    T3 <- pdffz(F1 * CDF2m1 * (1 + thet * (1 - F1) * (1 - CDF2m1)))
    
    T4 <- pdffz(CDF1m1 * CDF2m1 * (1 + thet * (1 - CDF1m1) * (1 - CDF2m1)))
    
    # Negative log-likelihood
    return(- log( T1 - T2 - T3 + T4) )
    
  }
  
  
  
  ## Coefficients:
  # x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
  mu1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  mu2 = c(0, 0, 0, 0, 0,  0, 0, 0, 0, 0) 
  
  # x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
  sigma1 <- c(0, 0, 0, 0, 0,  0, 0, 0, 0, 0) 
  sigma2 <- c(0, 0, 0, 0, 0,  0, 0, 0, 0, 0) 
  
  nu2 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0) 
  
  # x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
  or  = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  
  ### function for mu1
  nl_function_mu1 <- function(x){
    
    
    nl <- sqrt(x)*x^2 - 2*cos(3*x)
    
    return(nl*0.5)
    
  }
  
  ### function for sigma1
  nl_function_sigma1 <- function(x){
    
    func <- -40*(x^(3/2) - x^(4/3))
    
    func <- 2 * func
    
    return(func)
  }
  
  nl_func_sigma1 <- function(x){
    
    func <- -(log(x)) + cos(x)
    
    func <- 2 * func
    
    return(func)
  }
  
  ## function for mu2
  nl_function_mu2 <- function(x){
    
    nl <- -0.7*exp(x^2) + exp(x^0.4) - 0*log(x+0.2) 
    
    return(nl)
  }
  
  nl_function_sigma2 <- function(x){
    
    return(cos(2 * x))
  }
  
  nl_function_nu2 <- function(x){
    
    func <- -(sin(x) - exp(x)^2)
    
    func <- 0.7*func
    
    return(func)
  }
  
  nl_function_rho <- function(x){
    
    return(2*sin(x*4))
  }
  
  TrueBeta <-  vector('list')
  TrueBeta$mu1 <- mu1
  TrueBeta$mu2 <- mu2
  
  TrueBeta$sigma1 <- sigma1
  TrueBeta$sigma2 <- sigma2
  
  TrueBeta$nu2 <- nu2
  
  TrueBeta$or <- or
  
  set.seed(seed)
  
  #n.mstop <- 1500
  #n.mstop <- 500
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
  
  ############################################################################################################################################
  # GENERATION OF DISTRIBUTION PARAMETERS AND STUFF 
  # predictor
  train.eta.mu1_center <-  (x.train %*% mu1)
  train.eta.mu1_center <-  train.eta.mu1_center + nl_function_mu1(x1)
  
  # apply response function:
  train.eta.mu1<-  plogis(train.eta.mu1_center) 
  range(train.eta.mu1)
  
  #### SIGMA FOR ZALG RESPONSE
  train.eta.sigma1 <-  x.train %*% sigma1
  
  train.eta.sigma1 <- train.eta.sigma1 - 2 + nl_function_sigma1(x3)
  train.eta.sigma1 <- plogis(train.eta.sigma1)
  range(train.eta.sigma1)
  
  # predictor and apply response function: 
  train.eta.mu2 <-  x.train %*% mu2
  train.eta.mu2 <-  -2*0 + train.eta.mu2 + nl_function_mu2(x2)
  train.eta.mu2 <-  exp(train.eta.mu2) 
  range(train.eta.mu2)
  
  # predictor and apply response function: 
  train.eta.sigma2 <-  x.train %*% sigma2
  train.eta.sigma2 <-  train.eta.sigma2 + nl_function_sigma2(x4)
  train.eta.sigma2 <-  exp(train.eta.sigma2) 
  range(train.eta.sigma2)
  
  # predictor and apply response function: 
  train.eta.nu2 <-  x.train %*% nu2
  train.eta.nu2 <-  train.eta.nu2 + nl_function_nu2(x1) - 3
  train.eta.nu2 <-  plogis(train.eta.nu2) 
  range(train.eta.nu2)
  
  ### Copula parameter
  train.eta.or_center <- (x.train %*% or)
  train.eta.or_center <- train.eta.or_center - nl_function_rho(x3) 
  
  #train.eta.or <-  exp( train.eta.or_center)
  train.copula.parameter <- tanh(train.eta.or_center)
  range(train.copula.parameter)
  ############################################################################################################################################
  
  # SAMPLING OF BIVARIATE RESPONSE FROM COPULA
  # simulate data from bivariate binary DGP:
  data.gen.bivdisc <- function(mu1, mu2, sigma1, sigma2, theta, nu2 = NULL){
    
    y1y2 <- matrix(0, ncol = 2, nrow = length(mu1))
    
    for(i in 1:length(mu1)){
      
      paramlist1 <- list(mu = mu1[i], 
                         sigma = sigma1[i])
      
      
      paramlist2 <- list(mu = mu2[i], 
                         sigma = sigma2[i],
                         nu = nu2[i])
      
      copObj <- copula::fgmCopula(param = theta[i])
      
      copThing <- mvdc(copObj, margins = c("ZALG", "ZINBI"), paramMargins = list(paramlist1, paramlist2))
      
      y1y2[i,] <- rMvdc(copThing, n = 1)
    }
    
    
    dat <- data.frame(y1y2[,1], y1y2[,2])
    
    return(dat)
  }
  
  data.gen.bivdisc_energyscore <- function(n, mu1, mu2, sigma1, sigma2, theta, nu = NULL){
    
    
    paramlist1 <- list(mu = mu1, sigma = sigma1)
    
    paramlist2 <- list(mu = mu2, sigma = sigma2, nu = nu) 
    
    copObj <- copula::fgmCopula(param = theta)
    
    copThing <- mvdc(copObj, margins = c("ZALG", "ZINBI"), paramMargins = list(paramlist1, paramlist2))
    
    y1y2 <- rMvdc(copThing, n = n)
    
    dat <- data.frame(y1y2[,1], y1y2[,2])
    
    
    return(dat)
  }
  
  
  y.train <- data.gen.bivdisc(mu1 = train.eta.mu1, mu2 = train.eta.mu2, sigma1 = train.eta.sigma1, sigma2 = train.eta.sigma2, theta = train.copula.parameter, nu2 = train.eta.nu2)
  
  colnames(y.train)<- c('y1','y2')
  dat.train <- data.frame(y.train,x.train)
  
  # - test data
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
  
  ############################################################################################################################################
  # GENERATION OF DISTRIBUTION PARAMETERS AND STUFF 
  # predictor
  test.eta.mu1_center <- (x.test %*% mu1)
  test.eta.mu1_center <- test.eta.mu1_center + nl_function_mu1(x1.test)
  
  # apply response function:
  test.eta.mu1<-  plogis(test.eta.mu1_center) 
  
  #### SIGMA FOR ZALG RESPONSE
  test.eta.sigma1 <-  x.test %*% sigma1
  
  test.eta.sigma1 <- test.eta.sigma1 - 2 + nl_function_sigma1(x3.test)
  test.eta.sigma1 <- plogis(test.eta.sigma1)
  
  
  # predictor and apply response function: 
  test.eta.mu2 <-  x.test %*% mu2
  test.eta.mu2 <-  -2*0 + test.eta.mu2 + nl_function_mu2(x2.test)
  test.eta.mu2 <-  exp(test.eta.mu2) 
  
  # predictor and apply response function: 
  test.eta.sigma2 <-  x.test %*% sigma2
  test.eta.sigma2 <-  test.eta.sigma2 +  nl_function_sigma2(x4.test)
  test.eta.sigma2 <-  exp(test.eta.sigma2) 
  
  # predictor and apply response function: 
  test.eta.nu2 <-  x.test %*% nu2
  test.eta.nu2 <-  test.eta.nu2 + nl_function_nu2(x1.test) - 3
  test.eta.nu2 <-  plogis(test.eta.nu2) 
  
  
  ### Copula parameter
  test.eta.or_center <- (x.test %*% or)
  test.eta.or_center <- test.eta.or_center - nl_function_rho(x3.test)
  
  #test.eta.or <-  exp( test.eta.or_center)
  test.copula.parameter <- tanh(test.eta.or_center)
  ############################################################################################################################################
  
  y.test <- data.gen.bivdisc(mu1 = test.eta.mu1, mu2 = test.eta.mu2, sigma1 = test.eta.sigma1, sigma2 = test.eta.sigma2, theta = test.copula.parameter, nu2 = test.eta.nu2)
  
  colnames(y.test)<- c('y1','y2')
  dat.test <- data.frame(y.test,x.test)
  
  
  # FITTIIIINGG
  mstop.bivBern <-  vector('list')
  
  AUCbivBern <-  vector('list')
  BrierbivBern <-  vector('list')
  lik <- vector('list')
  MSEbivDisc <- vector("list")
  
  
  
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
                           sigma1 = model_equation,
                           sigma2 = model_equation,
                           rho = model_equation
  )
  
  
  bivModel_Formula <- list(mu1 = model_equation,
                           mu2 = model_equation,
                           sigma1 = model_equation,
                           sigma2 = model_equation,
                           nu2 = model_equation,
                           rho = model_equation
  )
  
  #############################################
  ################# - COPULA - ################
  #############################################
  # Bivariate copula model
  bivDiscCopula <- gamboostLSS(bivModel_Formula, data = dat.train, 
                               families = FGM_Cop_BivDiscrete(marg1 = "ZALG", 
                                                              marg2 = "ZINBI",
                                                              stabilization = "L2"), 
                               control = boost_control(mstop = 2000, risk = 'oobag',
                                                       nu  = boost.nu.steplength, trace = TRUE), 
                               method = 'noncyclic', weights = weight.mstop)
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 1900){
    mstop(bivDiscCopula) <- 3000
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 2900){
    mstop(bivDiscCopula) <- 4000
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 3990){
    bivDiscCopula[8000]
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 7990){
    bivDiscCopula[10000]
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 9990){
    bivDiscCopula[15000]
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 14990){
    bivDiscCopula[20000]
  }
  
  
  
  
  MSTOP_COP <- which.min(risk(bivDiscCopula, merge = T))
  oobag.risk.cop <- risk(bivDiscCopula,merge = T)
  
  # rm(bivDiscCopula)
  # dat.train_biv <- dat.train[weight.mstop == 1, ]
  # 
  # # RE-FIT until OPTIMAL MSTOP
  # bivDiscCopula <- gamboostLSS(bivModel_Formula, data = dat.train_biv, 
  #             families = FGM_Cop_BivDiscrete(marg1 = "ZALG", 
  #                                            marg2 = "ZINBI"),
  #             control = boost_control(mstop = MSTOP_COP, nu  = boost.nu.steplength, trace = TRUE), 
  #             method = 'noncyclic')
  bivDiscCopula <- bivDiscCopula[MSTOP_COP]
  
  
  mstop.bivCopula <-  vector('list')
  mstop.bivCopula$mstop <- MSTOP_COP
  mstop.bivCopula$mu1 <- bivDiscCopula$mu1$mstop()
  mstop.bivCopula$sigma1 <- bivDiscCopula$sigma1$mstop()
  
  mstop.bivCopula$mu2 <- bivDiscCopula$mu2$mstop()
  mstop.bivCopula$sigma2 <- bivDiscCopula$sigma2$mstop()
  mstop.bivCopula$nu2 <- bivDiscCopula$nu2$mstop()
  
  mstop.bivCopula$rho  <- bivDiscCopula$rho$mstop()
  
  # EXTRACT ALL COEFFICIENTS!
  coef.bivDiscCopula <- coef(bivDiscCopula, which = "")
  
  DiscreteMargins_Metrics <- vector("list")
  likBinCopula <- vector('list')
  
  
  ################ Get predictions from each covariate?
  newX <- matrix(rep(seq(from = 0, to = 1, length.out = 500),times = 10), nrow = 500)
  colnames(newX) <-  c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
  
  
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ### Copula predictions / selection
  Copula_MU1_predictions <- list()
  Copula_SIGMA1_predictions <- list()
  
  Copula_MU2_predictions <- list()
  Copula_SIGMA2_predictions <- list()
  Copula_NU2_predictions <- list()
  
  Copula_RHO_predictions <- list()
  
  ### MU 1
  Copula_MU1_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = 1, type = "link" )
  Copula_MU1_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x2", type = "link" )
  Copula_MU1_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x3", type = "link" )
  Copula_MU1_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x4", type = "link" )
  Copula_MU1_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x5", type = "link" )
  Copula_MU1_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x6", type = "link" )
  Copula_MU1_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x7", type = "link" )
  Copula_MU1_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x8", type = "link" )
  Copula_MU1_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x9", type = "link" )
  Copula_MU1_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x10", type = "link" )
  
  #### SIGMA 1
  Copula_SIGMA1_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = 1, type = "link" )
  Copula_SIGMA1_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x2", type = "link" )
  Copula_SIGMA1_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x3", type = "link" )
  Copula_SIGMA1_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x4", type = "link" )
  Copula_SIGMA1_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x5", type = "link" )
  Copula_SIGMA1_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x6", type = "link" )
  Copula_SIGMA1_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x7", type = "link" )
  Copula_SIGMA1_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x8", type = "link" )
  Copula_SIGMA1_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x9", type = "link" )
  Copula_SIGMA1_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x10", type = "link" )
  
  
  #### MU 2
  Copula_MU2_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = 1, type = "link" )
  Copula_MU2_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x2", type = "link" )
  Copula_MU2_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x3", type = "link" )
  Copula_MU2_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x4", type = "link" )
  Copula_MU2_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x5", type = "link" )
  Copula_MU2_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x6", type = "link" )
  Copula_MU2_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x7", type = "link" )
  Copula_MU2_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x8", type = "link" )
  Copula_MU2_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x9", type = "link" )
  Copula_MU2_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x10", type = "link" )
  
  #### SIGMA 2
  Copula_SIGMA2_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = 1, type = "link" )
  Copula_SIGMA2_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x2", type = "link" )
  Copula_SIGMA2_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x3", type = "link" )
  Copula_SIGMA2_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x4", type = "link" )
  Copula_SIGMA2_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x5", type = "link" )
  Copula_SIGMA2_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x6", type = "link" )
  Copula_SIGMA2_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x7", type = "link" )
  Copula_SIGMA2_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x8", type = "link" )
  Copula_SIGMA2_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x9", type = "link" )
  Copula_SIGMA2_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x10", type = "link" )
  
  
  #### NU 2
  Copula_NU2_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = 1, type = "link" )
  Copula_NU2_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x2", type = "link" )
  Copula_NU2_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x3", type = "link" )
  Copula_NU2_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x4", type = "link" )
  Copula_NU2_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x5", type = "link" )
  Copula_NU2_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x6", type = "link" )
  Copula_NU2_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x7", type = "link" )
  Copula_NU2_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x8", type = "link" )
  Copula_NU2_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x9", type = "link" )
  Copula_NU2_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x10", type = "link" )
  
  
  #### RHO
  Copula_RHO_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = 1, type = "link" )
  Copula_RHO_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x2", type = "link" )
  Copula_RHO_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x3", type = "link" )
  Copula_RHO_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x4", type = "link" )
  Copula_RHO_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x5", type = "link" )
  Copula_RHO_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x6", type = "link" )
  Copula_RHO_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x7", type = "link" )
  Copula_RHO_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x8", type = "link" )
  Copula_RHO_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x9", type = "link" )
  Copula_RHO_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x10", type = "link" )
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  
  
  
  # Predict distributional quantities
  predCopula.mu1 <- predict(bivDiscCopula$mu1, newdata = dat.test, type = 'response')
  predCopula.sigma1 <- predict(bivDiscCopula$sigma1, newdata = dat.test, type = 'response')
  
  predCopula.mu2 <- predict(bivDiscCopula$mu2, newdata = dat.test, type = 'response')
  predCopula.sigma2 <- predict(bivDiscCopula$sigma2, newdata = dat.test, type = 'response')
  predCopula.nu2 <- predict(bivDiscCopula$nu2, newdata = dat.test, type = 'response')
  
  
  predCopula.rho <- predict(bivDiscCopula$rho, newdata = dat.test, type = 'response')
  
  ######### COMPUTE PERFORMANCE: 
  
  # MSE (margin 1, then margin 2)
  DiscreteMargins_Metrics$mu1 <- mean((as.numeric(predCopula.mu1)-(as.numeric(dat.test$y1)))^2)
  DiscreteMargins_Metrics$mu2 <- mean((as.numeric(predCopula.mu2)-(as.numeric(dat.test$y2)))^2)
  
  # LOG-SCORE (LOG-LIKELIHOOD) (change to the copula likelihood)
  lik$biv <- sum(loss(mu1 = predCopula.mu1, 
                      mu2 = predCopula.mu2, 
                      sigma1 = predCopula.sigma1, 
                      sigma2 = predCopula.sigma2, 
                      nu2 = predCopula.nu2,
                      rho = predCopula.rho, 
                      y = y.test))
  
  
  ###########################################################################################################################
  ################ - univariat - ############################################################################################
  ###########################################################################################################################
  
  mstop.uni <- vector('list')
  coef.uni  <- vector('list')
  coef.uni_1  <- vector('list')
  
  univariateDiscreteEquation_1 <- list(mu = formula(y1 ~  bbs(x1) + 
                                                      bbs(x2) + 
                                                      bbs(x3) + 
                                                      bbs(x4) + 
                                                      bbs(x5) + 
                                                      bbs(x6) +
                                                      bbs(x7) + 
                                                      bbs(x8) + 
                                                      bbs(x9) + 
                                                      bbs(x10)),
                                       sigma = formula(y1 ~  bbs(x1) + 
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
  
  
  univariateDiscreteEquation_2 <- list(mu = formula(y2 ~  bbs(x1) + 
                                                      bbs(x2) + 
                                                      bbs(x3) + 
                                                      bbs(x4) + 
                                                      bbs(x5) + 
                                                      bbs(x6) +
                                                      bbs(x7) + 
                                                      bbs(x8) + 
                                                      bbs(x9) + 
                                                      bbs(x10)),
                                       sigma = formula(y2 ~  bbs(x1) + 
                                                         bbs(x2) + 
                                                         bbs(x3) + 
                                                         bbs(x4) + 
                                                         bbs(x5) + 
                                                         bbs(x6) +
                                                         bbs(x7) + 
                                                         bbs(x8) + 
                                                         bbs(x9) + 
                                                         bbs(x10),
                                       ),
                                       nu = formula(y2 ~  bbs(x1) + 
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
  
  
  
  
  # univariateDiscreteEquation_1 <- list(mu = formula(y1 ~  bols(x1) + 
  #                                                     bols(x2) + 
  #                                                     bols(x3) + 
  #                                                     bols(x4) + 
  #                                                     bols(x5) + 
  #                                                     bols(x6) +
  #                                                     bols(x7) + 
  #                                                     bols(x8) + 
  #                                                     bols(x9) + 
  #                                                     bols(x10)),
  #                                      sigma = formula(y1 ~  bols(x1) + 
  #                                                        bols(x2) + 
  #                                                        bols(x3) + 
  #                                                        bols(x4) + 
  #                                                        bols(x5) + 
  #                                                        bols(x6) +
  #                                                        bols(x7) + 
  #                                                        bols(x8) + 
  #                                                        bols(x9) + 
  #                                                        bols(x10))
  # )
  # 
  # 
  # univariateDiscreteEquation_2 <- list(mu = formula(y2 ~  bols(x1) + 
  #                                                     bols(x2) + 
  #                                                     bols(x3) + 
  #                                                     bols(x4) + 
  #                                                     bols(x5) + 
  #                                                     bols(x6) +
  #                                                     bols(x7) + 
  #                                                     bols(x8) + 
  #                                                     bols(x9) + 
  #                                                     bols(x10)),
  #                                      sigma = formula(y2 ~  bols(x1) + 
  #                                                        bols(x2) + 
  #                                                        bols(x3) + 
  #                                                        bols(x4) + 
  #                                                        bols(x5) + 
  #                                                        bols(x6) +
  #                                                        bols(x7) + 
  #                                                        bols(x8) + 
  #                                                        bols(x9) + 
  #                                                        bols(x10)),
  #                                      nu = formula(y2 ~  bols(x1) + 
  #                                                     bols(x2) + 
  #                                                     bols(x3) + 
  #                                                     bols(x4) + 
  #                                                     bols(x5) + 
  #                                                     bols(x6) +
  #                                                     bols(x7) + 
  #                                                     bols(x8) + 
  #                                                     bols(x9) + 
  #                                                     bols(x10))
  # )
  # 
  
  # - margin 1
  dat.train.mu1 = dat.train[, !(names(dat.train) %in% "y2")]
  dat.test.mu1  = dat.test[ , !(names(dat.test) %in% "y2")]
  
  glm.uni.mu1 <- gamboostLSS(univariateDiscreteEquation_1, 
                             data = dat.train.mu1, 
                             families = as.families(fname = "ZALG", stabilization = "L2"), 
                             control = boost_control(risk = 'oobag', nu = boost.nu.steplength, mstop = 2000, trace = T), 
                             method = "noncyclic",
                             weights = weight.mstop)
  
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  if(mstop.uni$mu1 >= 1990){
    glm.uni.mu1[5000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  if(mstop.uni$mu1 >= 4990){
    glm.uni.mu1[10000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  if(mstop.uni$mu1 >= 9990){
    glm.uni.mu1[20000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  # dat.train_mu1 <- dat.train.mu1[weight.mstop == 1, ]
  # 
  # glm.uni.mu1 <- gamboostLSS(univariateDiscreteEquation_1, 
  #                            data = dat.train_mu1, 
  #                            families = as.families(fname = "ZALG"), 
  #                         control = boost_control(mstop = mstop.uni$mu1, nu  = boost.nu.steplength, trace = T), 
  #                         method = "noncyclic")
  glm.uni.mu1 <- glm.uni.mu1[mstop.uni$mu1]
  
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 1  
  coef.uni$mu1 <- coef(glm.uni.mu1, which = '')
  
  DiscreteMargins_Metrics$Univariate_mu1 <- mean((as.numeric(predict(glm.uni.mu1$mu, newdata = dat.test.mu1, type = "response")) - (as.numeric(dat.test.mu1$y1)))^2)
  
  # - margin 2
  dat.train.mu2 =  dat.train[, !(names(dat.train) %in% "y1")]
  dat.test.mu2 = dat.test[, !(names(dat.test) %in% "y1")]
  
  glm.uni.mu2 <- gamboostLSS(univariateDiscreteEquation_2, 
                             data = dat.train.mu2, 
                             families = as.families(fname = "ZINBI", stabilization = "L2"), 
                             control = boost_control(risk = 'oobag', nu = boost.nu.steplength, trace = TRUE, mstop = 2000), 
                             method = "noncyclic",
                             weights = weight.mstop)
  
  
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  if(mstop.uni$mu2 >= 1990){
    glm.uni.mu2[5000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  if(mstop.uni$mu2 >= 4990){
    glm.uni.mu2[10000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  if(mstop.uni$mu2 >= 9990){
    glm.uni.mu2[20000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  # dat.train_mu2 <- dat.train.mu2[weight.mstop == 1, ]
  # 
  # glm.uni.mu2 <- gamboostLSS(univariateDiscreteEquation_2, 
  #                            data = dat.train_mu2, 
  #                            families = as.families(fname = "ZINBI"), 
  #                         control = boost_control(mstop = mstop.uni$mu2, nu  = boost.nu.steplength, trace=TRUE), method = "noncyclic")
  glm.uni.mu2 <- glm.uni.mu2[mstop.uni$mu2]
  
  mstop.uni$Margin1 <- mstop(glm.uni.mu1)
  
  mstop.uni$Margin2 <- mstop(glm.uni.mu2)
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 2  
  coef.uni$mu2 <- coef(glm.uni.mu2,which = '')
  
  DiscreteMargins_Metrics$Univariate_mu2 <- mean((as.numeric(predict(glm.uni.mu2$mu, newdata = dat.test.mu2, type = "response"))-(as.numeric(dat.test.mu2$y2)))^2)
  
  # Likelihood
  pred.mu1.uni <- as.numeric(predict(glm.uni.mu1$mu, newdata = dat.test.mu1, type = "response"))
  pred.sigma1.uni <- as.numeric(predict(glm.uni.mu1$sigma, newdata = dat.test.mu1, type = "response"))
  
  pred.mu2.uni <- as.numeric(predict(glm.uni.mu2$mu, newdata = dat.test.mu2, type = "response"))
  pred.sigma2.uni <- as.numeric(predict(glm.uni.mu2$sigma, newdata = dat.test.mu2, type = "response"))
  pred.nu2.uni <- as.numeric(predict(glm.uni.mu2$nu, newdata = dat.test.mu2, type = "response"))
  
  
  mu1.uni.loglik <- sum(-dZALG(x = dat.test.mu1$y1, mu = pred.mu1.uni, sigma = pred.sigma1.uni, log = T))
  mu2.uni.loglik <- sum(-dZINBI(x = dat.test.mu2$y2, mu = pred.mu2.uni, sigma = pred.sigma2.uni, nu = pred.nu2.uni, log = T))
  
  lik$uni<- mu1.uni.loglik + mu2.uni.loglik
  
  
  n.ges = c(n.train, n.mstop, n.test)
  pred.ges <- NULL#list(pred.mu1, pred.mu2, pred.or)
  pred.uni <- list(mu1 = pred.mu1.uni, 
                   sigma1 = pred.sigma1.uni,
                   mu2 = pred.mu2.uni,
                   sigma2 = pred.sigma2.uni,
                   nu2 = pred.nu2.uni
  )
  
  pred.cop <- list(mu1 = predCopula.mu1, 
                   sigma1 = predCopula.sigma1,
                   mu2 = predCopula.mu2,
                   sigma2 = predCopula.sigma2,
                   nu2 = predCopula.nu2,
                   rho = predCopula.rho)
  
  
  ##### Predictions:
  UNIVARIATE_MU1_predictions <- vector("list")
  UNIVARIATE_SIGMA1_predictions <- vector("list")
  
  UNIVARIATE_MU2_predictions <- vector("list")
  UNIVARIATE_SIGMA2_predictions <- vector("list")
  UNIVARIATE_NU2_predictions <- vector("list")
  
  
  #### MU 1
  UNIVARIATE_MU1_predictions$x1 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = 1, type = "link" )
  UNIVARIATE_MU1_predictions$x2 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x2", type = "link" )
  UNIVARIATE_MU1_predictions$x3 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x3", type = "link" )
  UNIVARIATE_MU1_predictions$x4 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x4", type = "link" )
  UNIVARIATE_MU1_predictions$x5 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x5", type = "link" )
  UNIVARIATE_MU1_predictions$x6 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x6", type = "link" )
  UNIVARIATE_MU1_predictions$x7 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x7", type = "link" )
  UNIVARIATE_MU1_predictions$x8 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x8", type = "link" )
  UNIVARIATE_MU1_predictions$x9 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x9", type = "link" )
  UNIVARIATE_MU1_predictions$x10 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x10", type = "link" )
  
  #### SIGMA 1
  UNIVARIATE_SIGMA1_predictions$x1 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = 1, type = "link" )
  UNIVARIATE_SIGMA1_predictions$x2 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x2", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x3 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x3", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x4 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x4", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x5 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x5", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x6 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x6", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x7 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x7", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x8 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x8", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x9 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x9", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x10 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x10", type = "link" )
  
  
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
  
  
  #### NU 2
  UNIVARIATE_NU2_predictions$x1 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = 1, type = "link" )
  UNIVARIATE_NU2_predictions$x2 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x2", type = "link" )
  UNIVARIATE_NU2_predictions$x3 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x3", type = "link" )
  UNIVARIATE_NU2_predictions$x4 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x4", type = "link" )
  UNIVARIATE_NU2_predictions$x5 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x5", type = "link" )
  UNIVARIATE_NU2_predictions$x6 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x6", type = "link" )
  UNIVARIATE_NU2_predictions$x7 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x7", type = "link" )
  UNIVARIATE_NU2_predictions$x8 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x8", type = "link" )
  UNIVARIATE_NU2_predictions$x9 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x9", type = "link" )
  UNIVARIATE_NU2_predictions$x10 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x10", type = "link" )
  
  
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
    sample_uni <- data.gen.bivdisc_energyscore(n = 1000, 
                                               mu1 = pred.mu1.uni[i], 
                                               mu2 = pred.mu2.uni[i], 
                                               sigma1 = pred.sigma1.uni[i], 
                                               sigma2 = pred.sigma2.uni[i],
                                               nu = check_extremes_nuZINBI(pred.nu2.uni[i]),
                                               theta = 0)
    
    pred_sample_uni[1,] <- sample_uni[,1]
    pred_sample_uni[2,] <- sample_uni[,2]
    
    es_uni[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_uni) 
    
    
    # copula
    sample_cop <- data.gen.bivdisc_energyscore(n = 1000, 
                                               mu1 = predCopula.mu1[i],  
                                               mu2 = predCopula.mu2[i], 
                                               sigma1 = predCopula.sigma1[i], 
                                               sigma2 = predCopula.sigma2[i],
                                               nu = check_extremes_nuZINBI(predCopula.nu2[i]),
                                               theta = check_extremes(predCopula.rho[i]))
    
    pred_sample_cop[1, ] <- sample_cop[,1]
    pred_sample_cop[2, ] <- sample_cop[,2]
    
    es_cop[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_cop) 
    
    
  }
  
  energy_score <- list()
  #energy_score$biv <- mean(es_biv)
  energy_score$uni <- mean(es_uni)
  energy_score$cop <- mean(es_cop)
  
  
  
  ########################################### We gather: mstop, predictions, MSE margin 1, margin 2, loglikelihood, energy score.
  output <- list(TrueBeta = TrueBeta, 
                 n = n.ges, 
                 p  = p, 
                 Likelihood = lik,
                 predict.uni = pred.uni, 
                 predict.cop = pred.cop,
                 energy_score = energy_score,  
                 ####
                 Coefficients.uni = coef.uni, 
                 mstop.uni = mstop.uni,
                 ###
                 DiscreteMarginMetrics = DiscreteMargins_Metrics,
                 ###
                 Copula_Predictions = list(MU1 = Copula_MU1_predictions, 
                                           SIGMA1 = Copula_SIGMA1_predictions,
                                           MU2 = Copula_MU2_predictions,
                                           SIGMA2 = Copula_SIGMA2_predictions,
                                           NU2 = Copula_NU2_predictions,
                                           RHO = Copula_RHO_predictions),
                 ###
                 Univariate_Predictions = list(MU1 = UNIVARIATE_MU1_predictions, 
                                               SIGMA1 = UNIVARIATE_SIGMA1_predictions,
                                               MU2 = UNIVARIATE_MU2_predictions, 
                                               SIGMA2 = UNIVARIATE_SIGMA2_predictions,
                                               NU2 = UNIVARIATE_NU2_predictions),
                 ###
                 CoefficientsCOPULA = coef.bivDiscCopula, 
                 oobag.riskCOPULA = oobag.risk.cop, 
                 energy_scoreCOPULA = energy_score,  
                 mstopCOPULA = mstop.bivCopula
  )
  
  
  return(output)
  
}

#### USING JOE COPULA
sim_JOE <- function(seed, n.train, n.mstop, p, n.test, corr, boost.nu.steplength = 0.1){
  
  
  loss <- function(mu1, mu2, sigma1, sigma2, rho, y, nu2 = NULL){
    
    
    # Margin 1 is ZALG:
    F1 <- pdffz(pZALG(y[,1], mu = mu1, sigma = sigma1))
    pdf1 <- pdffz(dZALG(y[,1], mu = mu1, sigma = sigma1))
    
    # Margin 2 is NBI:
    F2 <- pdffz(pZINBI(y[,2], mu = mu2, sigma = sigma2, nu = nu2))
    pdf2 <- pdffz(dZINBI(y[,2], mu = mu2, sigma = sigma2, nu = nu2))
    
    # Copula parameter 
    thet <- rho
    
    # Minus 1
    CDF1m1 <- pdffz(F1 - pdf1)
    CDF2m1 <- pdffz(F2 - pdf2)
    
    ### Copula Terms
    bit1 <- (1 - F1)^thet
    bit1m1 <- (1 - CDF1m1)^thet
    
    bit2 <- (1 - F2)^thet
    bit2m1 <- (1 - CDF2m1)^thet
    
    
    T1 <- pdffz(  1 - (bit1 + bit2 - bit1*bit2)^(1/thet)  )
    
    T2 <- pdffz(  1 - (bit1m1 + bit2 - bit1m1*bit2)^(1/thet)  ) 
    
    T3 <- pdffz(  1 - (bit1 + bit2m1 - bit1*bit2m1)^(1/thet)   )
    
    T4 <- pdffz(  1 - (bit1m1 + bit2m1 - bit1m1*bit2m1)^(1/thet)  ) 
    
    L_Total <- pdffz( T1 - T2 - T3 + T4 )
    
    # Negative log-likelihood
    return(- log( L_Total ) )
    
  }
  
  data.gen.bivdisc <- function(mu1, mu2, sigma1, sigma2, theta, nu2 = NULL){
    
    y1y2 <- matrix(0, ncol = 2, nrow = length(mu1))
    
    for(i in 1:length(mu1)){
      
      paramlist1 <- list(mu = mu1[i], 
                         sigma = sigma1[i])
      
      
      paramlist2 <- list(mu = mu2[i], 
                         sigma = sigma2[i],
                         nu = nu2[i])
      
      copObj <- copula::joeCopula(param = theta[i])
      
      copThing <- mvdc(copObj, margins = c("ZALG", "ZINBI"), paramMargins = list(paramlist1, paramlist2))
      
      y1y2[i,] <- rMvdc(copThing, n = 1)
    }
    
    
    dat <- data.frame(y1y2[,1], y1y2[,2])
    
    return(dat)
  }
  
  data.gen.bivdisc_energyscore <- function(n, mu1, mu2, sigma1, sigma2, theta, nu = NULL){
    
    
    paramlist1 <- list(mu = mu1, sigma = sigma1)
    
    paramlist2 <- list(mu = mu2, sigma = sigma2, nu = nu) 
    
    copObj <- copula::joeCopula(param = theta)
    
    copThing <- mvdc(copObj, margins = c("ZALG", "ZINBI"), paramMargins = list(paramlist1, paramlist2))
    
    y1y2 <- rMvdc(copThing, n = n)
    
    dat <- data.frame(y1y2[,1], y1y2[,2])
    
    
    return(dat)
  }
  
  
  
  
  ## Coefficients:
  # x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
  mu1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  mu2 = c(0, 0, 0, 0, 0,  0, 0, 0, 0, 0) 
  
  # x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
  sigma1 <- c(0, 0, 0, 0, 0,  0, 0, 0, 0, 0) 
  sigma2 <- c(0, 0, 0, 0, 0,  0, 0, 0, 0, 0) 
  
  nu2 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0) 
  
  # x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
  or  = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  
  ### function for mu1
  nl_function_mu1 <- function(x){
    
    
    nl <- sqrt(x)*x^2 - 2*cos(3*x)
    
    return(nl*0.5)
    
  }
  
  ### function for sigma1
  nl_function_sigma1 <- function(x){
    
    func <- -40*(x^(3/2) - x^(4/3))
    
    func <- 2 * func
    
    return(func)
  }
  
  nl_func_sigma1 <- function(x){
    
    func <- -(log(x)) + cos(x)
    
    func <- 2 * func
    
    return(func)
  }
  
  ## function for mu2
  nl_function_mu2 <- function(x){
    
    nl <- -0.7*exp(x^2) + exp(x^0.4) - 0*log(x+0.2) 
    
    return(nl)
  }
  
  nl_function_sigma2 <- function(x){
    
    return( -1.5*(1.5*cos(2 * x) + 0.5*log(x)) )
  }
  
  nl_function_nu2 <- function(x){
    
    func <- -(sin(x) - exp(x)^2)
    
    func <- 0.7*func
    
    return(func)
  }
  
  nl_function_rho <- function(x){
    
    #return(2*sin(x*4))
    return(2*sin(x*4))
  }
  
  TrueBeta <-  vector('list')
  TrueBeta$mu1 <- mu1
  TrueBeta$mu2 <- mu2
  
  TrueBeta$sigma1 <- sigma1
  TrueBeta$sigma2 <- sigma2
  
  TrueBeta$nu2 <- nu2
  
  TrueBeta$or <- or
  
  set.seed(seed)
  
 
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
  
  ############################################################################################################################################
  ############################################################################################################################################
  # - test data
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
  
  ############################################################################################################################################
  # GENERATION OF DISTRIBUTION PARAMETERS AND STUFF 
  # predictor
  
  ##### TRAIN DATA
  train.eta.mu1_center <- (x.train %*% mu1)
  train.eta.mu1_center <- train.eta.mu1_center + nl_function_mu1(x1)
  
  ##### TEST DATA
  test.eta.mu1_center <- (x.test %*% mu1)
  test.eta.mu1_center <- test.eta.mu1_center + nl_function_mu1(x1.test)

  # apply response function:
  train.eta.mu1<-  plogis(train.eta.mu1_center) 
  range(train.eta.mu1)
  
  test.eta.mu1<-  plogis(test.eta.mu1_center) 
  ########################################################################################################################
  #### SIGMA FOR ZALG RESPONSE
  train.eta.sigma1 <-  x.train %*% sigma1
  
  train.eta.sigma1 <- train.eta.sigma1 - 2 + nl_function_sigma1(x3) 
  train.eta.sigma1 <- plogis(train.eta.sigma1)
  range(train.eta.sigma1)
  
  ##### TEST DATA
  test.eta.sigma1 <-  x.test %*% sigma1
  test.eta.sigma1 <- test.eta.sigma1 - 2 + nl_function_sigma1(x3.test)
  test.eta.sigma1 <- plogis(test.eta.sigma1)
  ########################################################################################################################
  # predictor and apply response function: 
  train.eta.mu2 <-  x.train %*% mu2
  train.eta.mu2 <-  train.eta.mu2 + nl_function_mu2(x2)
  train.eta.mu2 <-  exp(train.eta.mu2) 
  range(train.eta.mu2)
  
  
  test.eta.mu2 <-  x.test %*% mu2
  test.eta.mu2 <-  test.eta.mu2 + nl_function_mu2(x2.test)
  test.eta.mu2 <-  exp(test.eta.mu2) 
  ########################################################################################################################
  # predictor and apply response function: 
  train.eta.sigma2 <-  x.train %*% sigma2
  train.eta.sigma2 <-  train.eta.sigma2 + nl_function_sigma2(x4)
  train.eta.sigma2 <-  exp(train.eta.sigma2) 
  range(train.eta.sigma2)
  
  test.eta.sigma2 <-  x.test %*% sigma2
  test.eta.sigma2 <-  test.eta.sigma2 +  nl_function_sigma2(x4.test)
  test.eta.sigma2 <-  exp(test.eta.sigma2) 
  ########################################################################################################################
  # predictor and apply response function: 
  train.eta.nu2 <-  x.train %*% nu2
  train.eta.nu2 <-  train.eta.nu2 + nl_function_nu2(x1) - 3
  train.eta.nu2 <-  plogis(train.eta.nu2) 
  range(train.eta.nu2)
  
  test.eta.nu2 <-  x.test %*% nu2
  test.eta.nu2 <-  test.eta.nu2 + nl_function_nu2(x1.test) - 3
  test.eta.nu2 <-  plogis(test.eta.nu2) 
  
  ########################################################################################################################
  ### Copula parameter
  train.eta.or_center <- (x.train %*% or)
  train.eta.or_center <- train.eta.or_center + nl_function_rho(x3) 
  
  train.copula.parameter <- (exp(train.eta.or_center) + 1)
  range(train.copula.parameter)
  
  TrueKendallTauRange <- VineCopula::BiCopPar2Tau(6, par = range(train.copula.parameter))  
  
  test.eta.or_center <- (x.test %*% or)
  test.eta.or_center <- test.eta.or_center + nl_function_rho(x3.test)
  
  test.copula.parameter <- (exp(test.eta.or_center) + 1)
  ############################################################################################################################################
  # SAMPLING OF BIVARIATE RESPONSE FROM COPULA
  y.train <- data.gen.bivdisc(mu1 = train.eta.mu1, mu2 = train.eta.mu2, sigma1 = train.eta.sigma1, sigma2 = train.eta.sigma2, theta = train.copula.parameter, nu2 = train.eta.nu2)
  
  colnames(y.train)<- c('y1','y2')
  dat.train <- data.frame(y.train,x.train)
  
  
  y.test <- data.gen.bivdisc(mu1 = test.eta.mu1, mu2 = test.eta.mu2, sigma1 = test.eta.sigma1, sigma2 = test.eta.sigma2, theta = test.copula.parameter, nu2 = test.eta.nu2)
  
  colnames(y.test)<- c('y1','y2')
  dat.test <- data.frame(y.test,x.test)
  
  ############################################################################################################################################
  # FITTIIIINGG
  mstop.bivBern <-  vector('list')
  
  AUCbivBern <-  vector('list')
  BrierbivBern <-  vector('list')
  lik <- vector('list')
  MSEbivDisc <- vector("list")
  
  
  
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
                           sigma1 = model_equation,
                           sigma2 = model_equation,
                           rho = model_equation
  )
  
  
  bivModel_Formula <- list(mu1 = model_equation,
                           mu2 = model_equation,
                           sigma1 = model_equation,
                           sigma2 = model_equation,
                           nu2 = model_equation,
                           rho = model_equation
  )
  
  #############################################
  ################# - COPULA - ################
  #############################################
  # Bivariate copula model
  bivDiscCopula <- gamboostLSS(bivModel_Formula, data = dat.train, 
                               families = Joe_Cop_BivDiscrete(marg1 = "ZALG", 
                                                              marg2 = "ZINBI",
                                                              stabilization = "L2"), 
                               control = boost_control(mstop = 2000, risk = 'oobag',
                                                       nu  = boost.nu.steplength, trace = TRUE), 
                               method = 'noncyclic', weights = weight.mstop)
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 1900){
    mstop(bivDiscCopula) <- 3000
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 2900){
    mstop(bivDiscCopula) <- 4000
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 3990){
    bivDiscCopula[8000]
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 7990){
    bivDiscCopula[10000]
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 9990){
    bivDiscCopula[15000]
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 14990){
    bivDiscCopula[20000]
  }
  
  
  
  
  MSTOP_COP <- which.min(risk(bivDiscCopula, merge = T))
  oobag.risk.cop <- risk(bivDiscCopula,merge = T)
  
  # rm(bivDiscCopula)
  # dat.train_biv <- dat.train[weight.mstop == 1, ]
  # 
  # # RE-FIT until OPTIMAL MSTOP
  # bivDiscCopula <- gamboostLSS(bivModel_Formula, data = dat.train_biv, 
  #             families = FGM_Cop_BivDiscrete(marg1 = "ZALG", 
  #                                            marg2 = "ZINBI"),
  #             control = boost_control(mstop = MSTOP_COP, nu  = boost.nu.steplength, trace = TRUE), 
  #             method = 'noncyclic')
  bivDiscCopula <- bivDiscCopula[MSTOP_COP]
  
  
  mstop.bivCopula <-  vector('list')
  mstop.bivCopula$mstop <- MSTOP_COP
  mstop.bivCopula$mu1 <- bivDiscCopula$mu1$mstop()
  mstop.bivCopula$sigma1 <- bivDiscCopula$sigma1$mstop()
  
  mstop.bivCopula$mu2 <- bivDiscCopula$mu2$mstop()
  mstop.bivCopula$sigma2 <- bivDiscCopula$sigma2$mstop()
  mstop.bivCopula$nu2 <- bivDiscCopula$nu2$mstop()
  
  mstop.bivCopula$rho  <- bivDiscCopula$rho$mstop()
  
  # EXTRACT ALL COEFFICIENTS!
  coef.bivDiscCopula <- coef(bivDiscCopula, which = "")
  
  DiscreteMargins_Metrics <- vector("list")
  likBinCopula <- vector('list')
  
  
  ################ Get predictions from each covariate?
  newX <- matrix(rep(seq(from = 0, to = 1, length.out = 500),times = 10), nrow = 500)
  colnames(newX) <-  c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
  
  
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ### Copula predictions / selection
  Copula_MU1_predictions <- list()
  Copula_SIGMA1_predictions <- list()
  
  Copula_MU2_predictions <- list()
  Copula_SIGMA2_predictions <- list()
  Copula_NU2_predictions <- list()
  
  Copula_RHO_predictions <- list()
  
  ### MU 1
  Copula_MU1_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = 1, type = "link" )
  Copula_MU1_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x2", type = "link" )
  Copula_MU1_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x3", type = "link" )
  Copula_MU1_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x4", type = "link" )
  Copula_MU1_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x5", type = "link" )
  Copula_MU1_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x6", type = "link" )
  Copula_MU1_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x7", type = "link" )
  Copula_MU1_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x8", type = "link" )
  Copula_MU1_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x9", type = "link" )
  Copula_MU1_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x10", type = "link" )
  
  #### SIGMA 1
  Copula_SIGMA1_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = 1, type = "link" )
  Copula_SIGMA1_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x2", type = "link" )
  Copula_SIGMA1_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x3", type = "link" )
  Copula_SIGMA1_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x4", type = "link" )
  Copula_SIGMA1_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x5", type = "link" )
  Copula_SIGMA1_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x6", type = "link" )
  Copula_SIGMA1_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x7", type = "link" )
  Copula_SIGMA1_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x8", type = "link" )
  Copula_SIGMA1_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x9", type = "link" )
  Copula_SIGMA1_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x10", type = "link" )
  
  
  #### MU 2
  Copula_MU2_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = 1, type = "link" )
  Copula_MU2_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x2", type = "link" )
  Copula_MU2_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x3", type = "link" )
  Copula_MU2_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x4", type = "link" )
  Copula_MU2_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x5", type = "link" )
  Copula_MU2_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x6", type = "link" )
  Copula_MU2_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x7", type = "link" )
  Copula_MU2_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x8", type = "link" )
  Copula_MU2_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x9", type = "link" )
  Copula_MU2_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x10", type = "link" )
  
  #### SIGMA 2
  Copula_SIGMA2_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = 1, type = "link" )
  Copula_SIGMA2_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x2", type = "link" )
  Copula_SIGMA2_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x3", type = "link" )
  Copula_SIGMA2_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x4", type = "link" )
  Copula_SIGMA2_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x5", type = "link" )
  Copula_SIGMA2_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x6", type = "link" )
  Copula_SIGMA2_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x7", type = "link" )
  Copula_SIGMA2_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x8", type = "link" )
  Copula_SIGMA2_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x9", type = "link" )
  Copula_SIGMA2_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x10", type = "link" )
  
  
  #### NU 2
  Copula_NU2_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = 1, type = "link" )
  Copula_NU2_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x2", type = "link" )
  Copula_NU2_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x3", type = "link" )
  Copula_NU2_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x4", type = "link" )
  Copula_NU2_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x5", type = "link" )
  Copula_NU2_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x6", type = "link" )
  Copula_NU2_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x7", type = "link" )
  Copula_NU2_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x8", type = "link" )
  Copula_NU2_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x9", type = "link" )
  Copula_NU2_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x10", type = "link" )
  
  
  #### RHO
  Copula_RHO_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = 1, type = "link" )
  Copula_RHO_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x2", type = "link" )
  Copula_RHO_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x3", type = "link" )
  Copula_RHO_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x4", type = "link" )
  Copula_RHO_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x5", type = "link" )
  Copula_RHO_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x6", type = "link" )
  Copula_RHO_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x7", type = "link" )
  Copula_RHO_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x8", type = "link" )
  Copula_RHO_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x9", type = "link" )
  Copula_RHO_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x10", type = "link" )
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  
  
  
  # Predict distributional quantities
  predCopula.mu1 <- predict(bivDiscCopula$mu1, newdata = dat.test, type = 'response')
  predCopula.sigma1 <- predict(bivDiscCopula$sigma1, newdata = dat.test, type = 'response')
  
  predCopula.mu2 <- predict(bivDiscCopula$mu2, newdata = dat.test, type = 'response')
  predCopula.sigma2 <- predict(bivDiscCopula$sigma2, newdata = dat.test, type = 'response')
  predCopula.nu2 <- predict(bivDiscCopula$nu2, newdata = dat.test, type = 'response')
  
  
  predCopula.rho <- predict(bivDiscCopula$rho, newdata = dat.test, type = 'response')
  
  ######### COMPUTE PERFORMANCE: 
  
  # MSE (margin 1, then margin 2)
  DiscreteMargins_Metrics$mu1 <- mean((as.numeric(predCopula.mu1)-(as.numeric(dat.test$y1)))^2)
  DiscreteMargins_Metrics$mu2 <- mean((as.numeric(predCopula.mu2)-(as.numeric(dat.test$y2)))^2)
  
  # LOG-SCORE (LOG-LIKELIHOOD) (change to the copula likelihood)
  lik$biv <- sum(loss(mu1 = predCopula.mu1, 
                      mu2 = predCopula.mu2, 
                      sigma1 = predCopula.sigma1, 
                      sigma2 = predCopula.sigma2, 
                      nu2 = predCopula.nu2,
                      rho = predCopula.rho, 
                      y = y.test))
  
  
  ###########################################################################################################################
  ################ - univariat - ############################################################################################
  ###########################################################################################################################
  
  mstop.uni <- vector('list')
  coef.uni  <- vector('list')
  coef.uni_1  <- vector('list')
  
  univariateDiscreteEquation_1 <- list(mu = formula(y1 ~  bbs(x1) + 
                                                      bbs(x2) + 
                                                      bbs(x3) + 
                                                      bbs(x4) + 
                                                      bbs(x5) + 
                                                      bbs(x6) +
                                                      bbs(x7) + 
                                                      bbs(x8) + 
                                                      bbs(x9) + 
                                                      bbs(x10)),
                                       sigma = formula(y1 ~  bbs(x1) + 
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
  
  
  univariateDiscreteEquation_2 <- list(mu = formula(y2 ~  bbs(x1) + 
                                                      bbs(x2) + 
                                                      bbs(x3) + 
                                                      bbs(x4) + 
                                                      bbs(x5) + 
                                                      bbs(x6) +
                                                      bbs(x7) + 
                                                      bbs(x8) + 
                                                      bbs(x9) + 
                                                      bbs(x10)),
                                       sigma = formula(y2 ~  bbs(x1) + 
                                                         bbs(x2) + 
                                                         bbs(x3) + 
                                                         bbs(x4) + 
                                                         bbs(x5) + 
                                                         bbs(x6) +
                                                         bbs(x7) + 
                                                         bbs(x8) + 
                                                         bbs(x9) + 
                                                         bbs(x10),
                                       ),
                                       nu = formula(y2 ~  bbs(x1) + 
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
  
  
  
  
  # univariateDiscreteEquation_1 <- list(mu = formula(y1 ~  bols(x1) + 
  #                                                     bols(x2) + 
  #                                                     bols(x3) + 
  #                                                     bols(x4) + 
  #                                                     bols(x5) + 
  #                                                     bols(x6) +
  #                                                     bols(x7) + 
  #                                                     bols(x8) + 
  #                                                     bols(x9) + 
  #                                                     bols(x10)),
  #                                      sigma = formula(y1 ~  bols(x1) + 
  #                                                        bols(x2) + 
  #                                                        bols(x3) + 
  #                                                        bols(x4) + 
  #                                                        bols(x5) + 
  #                                                        bols(x6) +
  #                                                        bols(x7) + 
  #                                                        bols(x8) + 
  #                                                        bols(x9) + 
  #                                                        bols(x10))
  # )
  # 
  # 
  # univariateDiscreteEquation_2 <- list(mu = formula(y2 ~  bols(x1) + 
  #                                                     bols(x2) + 
  #                                                     bols(x3) + 
  #                                                     bols(x4) + 
  #                                                     bols(x5) + 
  #                                                     bols(x6) +
  #                                                     bols(x7) + 
  #                                                     bols(x8) + 
  #                                                     bols(x9) + 
  #                                                     bols(x10)),
  #                                      sigma = formula(y2 ~  bols(x1) + 
  #                                                        bols(x2) + 
  #                                                        bols(x3) + 
  #                                                        bols(x4) + 
  #                                                        bols(x5) + 
  #                                                        bols(x6) +
  #                                                        bols(x7) + 
  #                                                        bols(x8) + 
  #                                                        bols(x9) + 
  #                                                        bols(x10)),
  #                                      nu = formula(y2 ~  bols(x1) + 
  #                                                     bols(x2) + 
  #                                                     bols(x3) + 
  #                                                     bols(x4) + 
  #                                                     bols(x5) + 
  #                                                     bols(x6) +
  #                                                     bols(x7) + 
  #                                                     bols(x8) + 
  #                                                     bols(x9) + 
  #                                                     bols(x10))
  # )
  # 
  
  # - margin 1
  dat.train.mu1 = dat.train[, !(names(dat.train) %in% "y2")]
  dat.test.mu1  = dat.test[ , !(names(dat.test) %in% "y2")]
  
  glm.uni.mu1 <- gamboostLSS(univariateDiscreteEquation_1, 
                             data = dat.train.mu1, 
                             families = as.families(fname = "ZALG", stabilization = "L2"), 
                             control = boost_control(risk = 'oobag', nu = boost.nu.steplength, mstop = 2000, trace = T), 
                             method = "noncyclic",
                             weights = weight.mstop)
  
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  if(mstop.uni$mu1 >= 1990){
    glm.uni.mu1[5000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  if(mstop.uni$mu1 >= 4990){
    glm.uni.mu1[10000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  if(mstop.uni$mu1 >= 9990){
    glm.uni.mu1[20000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  # dat.train_mu1 <- dat.train.mu1[weight.mstop == 1, ]
  # 
  # glm.uni.mu1 <- gamboostLSS(univariateDiscreteEquation_1, 
  #                            data = dat.train_mu1, 
  #                            families = as.families(fname = "ZALG"), 
  #                         control = boost_control(mstop = mstop.uni$mu1, nu  = boost.nu.steplength, trace = T), 
  #                         method = "noncyclic")
  glm.uni.mu1 <- glm.uni.mu1[mstop.uni$mu1]
  
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 1  
  coef.uni$mu1 <- coef(glm.uni.mu1, which = '')
  
  DiscreteMargins_Metrics$Univariate_mu1 <- mean((as.numeric(predict(glm.uni.mu1$mu, newdata = dat.test.mu1, type = "response")) - (as.numeric(dat.test.mu1$y1)))^2)
  
  # - margin 2
  dat.train.mu2 =  dat.train[, !(names(dat.train) %in% "y1")]
  dat.test.mu2 = dat.test[, !(names(dat.test) %in% "y1")]
  
  glm.uni.mu2 <- gamboostLSS(univariateDiscreteEquation_2, 
                             data = dat.train.mu2, 
                             families = as.families(fname = "ZINBI", stabilization = "L2"), 
                             control = boost_control(risk = 'oobag', nu = boost.nu.steplength, trace = TRUE, mstop = 2000), 
                             method = "noncyclic",
                             weights = weight.mstop)
  
  
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  if(mstop.uni$mu2 >= 1990){
    glm.uni.mu2[5000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  if(mstop.uni$mu2 >= 4990){
    glm.uni.mu2[10000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  if(mstop.uni$mu2 >= 9990){
    glm.uni.mu2[20000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  # dat.train_mu2 <- dat.train.mu2[weight.mstop == 1, ]
  # 
  # glm.uni.mu2 <- gamboostLSS(univariateDiscreteEquation_2, 
  #                            data = dat.train_mu2, 
  #                            families = as.families(fname = "ZINBI"), 
  #                         control = boost_control(mstop = mstop.uni$mu2, nu  = boost.nu.steplength, trace=TRUE), method = "noncyclic")
  glm.uni.mu2 <- glm.uni.mu2[mstop.uni$mu2]
  
  mstop.uni$Margin1 <- mstop(glm.uni.mu1)
  
  mstop.uni$Margin2 <- mstop(glm.uni.mu2)
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 2  
  coef.uni$mu2 <- coef(glm.uni.mu2,which = '')
  
  DiscreteMargins_Metrics$Univariate_mu2 <- mean((as.numeric(predict(glm.uni.mu2$mu, newdata = dat.test.mu2, type = "response"))-(as.numeric(dat.test.mu2$y2)))^2)
  
  # Likelihood
  pred.mu1.uni <- as.numeric(predict(glm.uni.mu1$mu, newdata = dat.test.mu1, type = "response"))
  pred.sigma1.uni <- as.numeric(predict(glm.uni.mu1$sigma, newdata = dat.test.mu1, type = "response"))
  
  pred.mu2.uni <- as.numeric(predict(glm.uni.mu2$mu, newdata = dat.test.mu2, type = "response"))
  pred.sigma2.uni <- as.numeric(predict(glm.uni.mu2$sigma, newdata = dat.test.mu2, type = "response"))
  pred.nu2.uni <- as.numeric(predict(glm.uni.mu2$nu, newdata = dat.test.mu2, type = "response"))
  
  
  mu1.uni.loglik <- sum(-dZALG(x = dat.test.mu1$y1, mu = pred.mu1.uni, sigma = pred.sigma1.uni, log = T))
  mu2.uni.loglik <- sum(-dZINBI(x = dat.test.mu2$y2, mu = pred.mu2.uni, sigma = pred.sigma2.uni, nu = pred.nu2.uni, log = T))
  
  lik$uni<- mu1.uni.loglik + mu2.uni.loglik
  
  
  n.ges = c(n.train, n.mstop, n.test)
  pred.ges <- NULL#list(pred.mu1, pred.mu2, pred.or)
  pred.uni <- list(mu1 = pred.mu1.uni, 
                   sigma1 = pred.sigma1.uni,
                   mu2 = pred.mu2.uni,
                   sigma2 = pred.sigma2.uni,
                   nu2 = pred.nu2.uni
  )
  
  pred.cop <- list(mu1 = predCopula.mu1, 
                   sigma1 = predCopula.sigma1,
                   mu2 = predCopula.mu2,
                   sigma2 = predCopula.sigma2,
                   nu2 = predCopula.nu2,
                   rho = predCopula.rho)
  
  
  ##### Predictions:
  UNIVARIATE_MU1_predictions <- vector("list")
  UNIVARIATE_SIGMA1_predictions <- vector("list")
  
  UNIVARIATE_MU2_predictions <- vector("list")
  UNIVARIATE_SIGMA2_predictions <- vector("list")
  UNIVARIATE_NU2_predictions <- vector("list")
  
  
  #### MU 1
  UNIVARIATE_MU1_predictions$x1 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = 1, type = "link" )
  UNIVARIATE_MU1_predictions$x2 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x2", type = "link" )
  UNIVARIATE_MU1_predictions$x3 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x3", type = "link" )
  UNIVARIATE_MU1_predictions$x4 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x4", type = "link" )
  UNIVARIATE_MU1_predictions$x5 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x5", type = "link" )
  UNIVARIATE_MU1_predictions$x6 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x6", type = "link" )
  UNIVARIATE_MU1_predictions$x7 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x7", type = "link" )
  UNIVARIATE_MU1_predictions$x8 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x8", type = "link" )
  UNIVARIATE_MU1_predictions$x9 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x9", type = "link" )
  UNIVARIATE_MU1_predictions$x10 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x10", type = "link" )
  
  #### SIGMA 1
  UNIVARIATE_SIGMA1_predictions$x1 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = 1, type = "link" )
  UNIVARIATE_SIGMA1_predictions$x2 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x2", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x3 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x3", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x4 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x4", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x5 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x5", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x6 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x6", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x7 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x7", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x8 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x8", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x9 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x9", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x10 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x10", type = "link" )
  
  
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
  
  
  #### NU 2
  UNIVARIATE_NU2_predictions$x1 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = 1, type = "link" )
  UNIVARIATE_NU2_predictions$x2 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x2", type = "link" )
  UNIVARIATE_NU2_predictions$x3 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x3", type = "link" )
  UNIVARIATE_NU2_predictions$x4 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x4", type = "link" )
  UNIVARIATE_NU2_predictions$x5 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x5", type = "link" )
  UNIVARIATE_NU2_predictions$x6 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x6", type = "link" )
  UNIVARIATE_NU2_predictions$x7 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x7", type = "link" )
  UNIVARIATE_NU2_predictions$x8 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x8", type = "link" )
  UNIVARIATE_NU2_predictions$x9 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x9", type = "link" )
  UNIVARIATE_NU2_predictions$x10 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x10", type = "link" )
  
  
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
    sample_uni <- data.gen.bivdisc_energyscore(n = 1000, 
                                               mu1 = pred.mu1.uni[i], 
                                               mu2 = pred.mu2.uni[i], 
                                               sigma1 = pred.sigma1.uni[i], 
                                               sigma2 = pred.sigma2.uni[i],
                                               nu = check_extremes_nuZINBI(pred.nu2.uni[i]),
                                               theta = 1)
    
    pred_sample_uni[1,] <- sample_uni[,1]
    pred_sample_uni[2,] <- sample_uni[,2]
    
    es_uni[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_uni) 
    
    
    # copula
    sample_cop <- data.gen.bivdisc_energyscore(n = 1000, 
                                               mu1 = predCopula.mu1[i],  
                                               mu2 = predCopula.mu2[i], 
                                               sigma1 = predCopula.sigma1[i], 
                                               sigma2 = predCopula.sigma2[i],
                                               nu = check_extremes_nuZINBI(predCopula.nu2[i]),
                                               theta = predCopula.rho[i])
    
    pred_sample_cop[1, ] <- sample_cop[,1]
    pred_sample_cop[2, ] <- sample_cop[,2]
    
    es_cop[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_cop) 
    
    
  }
  
  energy_score <- list()
  #energy_score$biv <- mean(es_biv)
  energy_score$uni <- mean(es_uni)
  energy_score$cop <- mean(es_cop)
  
  
  
  ########################################### We gather: mstop, predictions, MSE margin 1, margin 2, loglikelihood, energy score.
  output <- list(TrueBeta = TrueBeta, 
                 n = n.ges, 
                 p  = p, 
                 Likelihood = lik,
                 predict.uni = pred.uni, 
                 predict.cop = pred.cop,
                 energy_score = energy_score,  
                 ####
                 Coefficients.uni = coef.uni, 
                 mstop.uni = mstop.uni,
                 ###
                 DiscreteMarginMetrics = DiscreteMargins_Metrics,
                 ###
                 Copula_Predictions = list(MU1 = Copula_MU1_predictions, 
                                           SIGMA1 = Copula_SIGMA1_predictions,
                                           MU2 = Copula_MU2_predictions,
                                           SIGMA2 = Copula_SIGMA2_predictions,
                                           NU2 = Copula_NU2_predictions,
                                           RHO = Copula_RHO_predictions),
                 ###
                 Univariate_Predictions = list(MU1 = UNIVARIATE_MU1_predictions, 
                                               SIGMA1 = UNIVARIATE_SIGMA1_predictions,
                                               MU2 = UNIVARIATE_MU2_predictions, 
                                               SIGMA2 = UNIVARIATE_SIGMA2_predictions,
                                               NU2 = UNIVARIATE_NU2_predictions),
                 ###
                 CoefficientsCOPULA = coef.bivDiscCopula, 
                 oobag.riskCOPULA = oobag.risk.cop, 
                 energy_scoreCOPULA = energy_score,  
                 mstopCOPULA = mstop.bivCopula,
                 TrueKendallRange = TrueKendallTauRange
  )
  
  
  return(output)
  
}
sim_JOE_LINEAR <- function(seed, n.train, n.mstop, p, n.test, corr, boost.nu.steplength = 0.1){
  
  
  loss <- function(mu1, mu2, sigma1, sigma2, rho, y, nu2 = NULL){
    
    
    # Margin 1 is ZALG:
    F1 <- pdffz(pZALG(y[,1], mu = mu1, sigma = sigma1))
    pdf1 <- pdffz(dZALG(y[,1], mu = mu1, sigma = sigma1))
    
    # Margin 2 is NBI:
    F2 <- pdffz(pZINBI(y[,2], mu = mu2, sigma = sigma2, nu = nu2))
    pdf2 <- pdffz(dZINBI(y[,2], mu = mu2, sigma = sigma2, nu = nu2))
    
    # Copula parameter 
    thet <- rho
    
    # Minus 1
    CDF1m1 <- pdffz(F1 - pdf1)
    CDF2m1 <- pdffz(F2 - pdf2)
    
    ### Copula Terms
    # bit1 <- (1 - F1)^thet
    # bit1m1 <- (1 - CDF1m1)^thet
    # 
    # bit2 <- (1 - F2)^thet
    # bit2m1 <- (1 - CDF2m1)^thet
    # 
    # 
    # T1 <- pdffz(  1 - (bit1 + bit2 - bit1*bit2)^(1/thet)  )
    # 
    # T2 <- pdffz(  1 - (bit1m1 + bit2 - bit1m1*bit2)^(1/thet)  ) 
    # 
    # T3 <- pdffz(  1 - (bit1 + bit2m1 - bit1*bit2m1)^(1/thet)   )
    # 
    # T4 <- pdffz(  1 - (bit1m1 + bit2m1 - bit1m1*bit2m1)^(1/thet)  ) 
    T1 <- VineCopula::BiCopCDF(u1 = F1, u2 = F2, family = 6, par = thet)
    T2 <- VineCopula::BiCopCDF(u1 = CDF1m1, u2 = F2, family = 6, par = thet)
    T3 <- VineCopula::BiCopCDF(u1 = F1, u2 = CDF2m1, family = 6, par = thet)
    T4 <- VineCopula::BiCopCDF(u1 = CDF1m1, u2 = CDF2m1, family = 6, par = thet)
    
    L_Total <- pdffz( T1 - T2 - T3 + T4 )
    
    # Negative log-likelihood
    return(- log( L_Total ) )
    
  }
  
  data.gen.bivdisc <- function(mu1, mu2, sigma1, sigma2, theta, nu2 = NULL){
    
    y1y2 <- matrix(0, ncol = 2, nrow = length(mu1))
    
    for(i in 1:length(mu1)){
      
      paramlist1 <- list(mu = mu1[i], 
                         sigma = sigma1[i])
      
      
      paramlist2 <- list(mu = mu2[i], 
                         sigma = sigma2[i],
                         nu = nu2[i])
      
      copObj <- copula::joeCopula(param = theta[i])
      
      copThing <- mvdc(copObj, margins = c("ZALG", "ZINBI"), paramMargins = list(paramlist1, paramlist2))
      
      y1y2[i,] <- rMvdc(copThing, n = 1)
    }
    
    
    dat <- data.frame(y1y2[,1], y1y2[,2])
    
    return(dat)
  }
  
  data.gen.bivdisc_energyscore <- function(n, mu1, mu2, sigma1, sigma2, theta, nu = NULL){
    
    
    paramlist1 <- list(mu = mu1, sigma = sigma1)
    
    paramlist2 <- list(mu = mu2, sigma = sigma2, nu = nu) 
    
    copObj <- copula::joeCopula(param = theta)
    
    copThing <- mvdc(copObj, margins = c("ZALG", "ZINBI"), paramMargins = list(paramlist1, paramlist2))
    
    y1y2 <- rMvdc(copThing, n = n)
    
    dat <- data.frame(y1y2[,1], y1y2[,2])
    
    
    return(dat)
  }
  
  
  ### Only the first 5 have an effect:
  # x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
  mu1 <- c(-1, 0, +1, 0, 0, 0, 0, 0, 0, 0)
  mu2 <- c(+1.5, -1.5, 0, 0, 0, 0, 0, 0, 0, 0)
  
  # x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
  sigma1  <- c(0, 0, 0, +1, +1, 0, 0, -2, 0, 0)
  sigma2  <- c(0, -0.75, 0, +1, 0, 0, 0, 0, 0, 0)
  nu2     <- c(0, -0.75, +1, 0, 0, 0, 0, 0, 0, 0) 
  
  # x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
  or  <- c(0, -0.5, +1.5, 0, 1.5, 0, 0, 0, 0, 0)
  
  
  ### function for mu1
  nl_function_mu1 <- function(x){
    
    
    nl <- sqrt(x)*x^2 - 2*cos(3*x)
    
    return(nl*0.5)
    
  }
  
  ### function for sigma1
  nl_function_sigma1 <- function(x){
    
    func <- -40*(x^(3/2) - x^(4/3))
    
    func <- 2 * func
    
    return(func)
  }
  
  nl_func_sigma1 <- function(x){
    
    func <- -(log(x)) + cos(x)
    
    func <- 2 * func
    
    return(func)
  }
  
  ## function for mu2
  nl_function_mu2 <- function(x){
    
    nl <- -0.7*exp(x^2) + exp(x^0.4) - 0*log(x+0.2) 
    
    return(nl)
  }
  
  nl_function_sigma2 <- function(x){
    
    return( -1.5*(1.5*cos(2 * x) + 0.5*log(x)) )
  }
  
  nl_function_nu2 <- function(x){
    
    func <- -(sin(x) - exp(x)^2)
    
    func <- 0.7*func
    
    return(func)
  }
  
  nl_function_rho <- function(x){
    
    #return(2*sin(x*4))
    return(2*sin(x*4))
  }
  
  TrueBeta <-  vector('list')
  TrueBeta$mu1 <- mu1
  TrueBeta$mu2 <- mu2
  
  TrueBeta$sigma1 <- sigma1
  TrueBeta$sigma2 <- sigma2
  
  TrueBeta$nu2 <- nu2
  
  TrueBeta$or <- or
  
  set.seed(seed)
  
  
  n <- n.train + n.mstop
  weight.mstop <- c(rep(1, times = n.train), rep(0, times = n.mstop)) 
  
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
  
  ############################################################################################################################################
  ############################################################################################################################################
  # - test data
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
  
  ############################################################################################################################################
  # GENERATION OF DISTRIBUTION PARAMETERS AND STUFF 
  # predictor
  
  ##### TRAIN DATA
  train.eta.mu1_center <- (x.train %*% mu1)
  train.eta.mu1_center <- train.eta.mu1_center 
  
  ##### TEST DATA
  test.eta.mu1_center <- (x.test %*% mu1)
  test.eta.mu1_center <- test.eta.mu1_center 
  
  # apply response function:
  train.eta.mu1<-  plogis(train.eta.mu1_center) 
  range(train.eta.mu1)
  
  test.eta.mu1<-  plogis(test.eta.mu1_center) 
  ########################################################################################################################
  #### SIGMA FOR ZALG RESPONSE
  train.eta.sigma1 <-  x.train %*% sigma1
  
  train.eta.sigma1 <- train.eta.sigma1 #- 2
  train.eta.sigma1 <- plogis(train.eta.sigma1)
  range(train.eta.sigma1)
  
  ##### TEST DATA
  test.eta.sigma1 <-  x.test %*% sigma1
  test.eta.sigma1 <- test.eta.sigma1 #- 2
  test.eta.sigma1 <- plogis(test.eta.sigma1)
  ########################################################################################################################
  # predictor and apply response function: 
  train.eta.mu2 <-  x.train %*% mu2
  train.eta.mu2 <-  train.eta.mu2
  train.eta.mu2 <-  exp(train.eta.mu2) 
  range(train.eta.mu2)
  
  
  test.eta.mu2 <-  x.test %*% mu2
  test.eta.mu2 <-  test.eta.mu2 
  test.eta.mu2 <-  exp(test.eta.mu2) 
  ########################################################################################################################
  # predictor and apply response function: 
  train.eta.sigma2 <-  x.train %*% sigma2
  train.eta.sigma2 <-  train.eta.sigma2 
  train.eta.sigma2 <-  exp(train.eta.sigma2) 
  range(train.eta.sigma2)
  
  test.eta.sigma2 <-  x.test %*% sigma2
  test.eta.sigma2 <-  test.eta.sigma2
  test.eta.sigma2 <-  exp(test.eta.sigma2) 
  ########################################################################################################################
  # predictor and apply response function: 
  train.eta.nu2 <-  x.train %*% nu2
  train.eta.nu2 <-  train.eta.nu2 
  train.eta.nu2 <-  plogis(train.eta.nu2) 
  range(train.eta.nu2)
  
  test.eta.nu2 <-  x.test %*% nu2
  test.eta.nu2 <-  test.eta.nu2  
  test.eta.nu2 <-  plogis(test.eta.nu2) 
  
  ########################################################################################################################
  ### Copula parameter
  train.eta.or_center <- (x.train %*% or)
  train.eta.or_center <- train.eta.or_center 
  
  train.copula.parameter <- (exp(train.eta.or_center) + 1)
  range(train.copula.parameter)
  
  TrueKendallTauRange <- VineCopula::BiCopPar2Tau(6, par = range(train.copula.parameter))  
  
  test.eta.or_center <- (x.test %*% or)
  test.eta.or_center <- test.eta.or_center 
  
  test.copula.parameter <- (exp(test.eta.or_center) + 1)
  ############################################################################################################################################
  # SAMPLING OF BIVARIATE RESPONSE FROM COPULA
  y.train <- data.gen.bivdisc(mu1 = train.eta.mu1, mu2 = train.eta.mu2, sigma1 = train.eta.sigma1, sigma2 = train.eta.sigma2, theta = train.copula.parameter, nu2 = train.eta.nu2)
  
  colnames(y.train)<- c('y1','y2')
  dat.train <- data.frame(y.train,x.train)
  
  
  y.test <- data.gen.bivdisc(mu1 = test.eta.mu1, mu2 = test.eta.mu2, sigma1 = test.eta.sigma1, sigma2 = test.eta.sigma2, theta = test.copula.parameter, nu2 = test.eta.nu2)
  
  colnames(y.test)<- c('y1','y2')
  dat.test <- data.frame(y.test,x.test)
  
  ############################################################################################################################################
  # FITTIIIINGG
  mstop.bivBern <-  vector('list')
  
  AUCbivBern <-  vector('list')
  BrierbivBern <-  vector('list')
  lik <- vector('list')
  MSEbivDisc <- vector("list")
  
  
  
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
  
  model_equation <- formula(cbind(y1,y2) ~ bols(x1) + 
                              bols(x2) + 
                              bols(x3) + 
                              bols(x4) + 
                              bols(x5) + 
                              bols(x6) +
                              bols(x7) + 
                              bols(x8) + 
                              bols(x9) + 
                              bols(x10))
  
  bivModel_Formula <- list(mu1 = model_equation,
                           mu2 = model_equation,
                           sigma1 = model_equation,
                           sigma2 = model_equation,
                           rho = model_equation
  )
  
  
  bivModel_Formula <- list(mu1 = model_equation,
                           mu2 = model_equation,
                           sigma1 = model_equation,
                           sigma2 = model_equation,
                           nu2 = model_equation,
                           rho = model_equation
  )
  
  #############################################
  ################# - COPULA - ################
  #############################################
  # Bivariate copula model
  bivDiscCopula <- gamboostLSS(bivModel_Formula, data = dat.train, 
                               families = Joe_Cop_BivDiscrete(marg1 = "ZALG", 
                                                              marg2 = "ZINBI",
                                                              stabilization = "L2"), 
                               control = boost_control(mstop = 2000, risk = 'oobag',
                                                       nu  = boost.nu.steplength, trace = TRUE), 
                               method = 'noncyclic', weights = weight.mstop)
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 1900){
    mstop(bivDiscCopula) <- 3000
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 2900){
    mstop(bivDiscCopula) <- 4000
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 3990){
    bivDiscCopula[8000]
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 7990){
    bivDiscCopula[10000]
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 9990){
    bivDiscCopula[15000]
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 14990){
    bivDiscCopula[20000]
  }
  
  
  
  
  MSTOP_COP <- which.min(risk(bivDiscCopula, merge = T))
  oobag.risk.cop <- risk(bivDiscCopula,merge = T)
  
  # rm(bivDiscCopula)
  # dat.train_biv <- dat.train[weight.mstop == 1, ]
  # 
  # # RE-FIT until OPTIMAL MSTOP
  # bivDiscCopula <- gamboostLSS(bivModel_Formula, data = dat.train_biv, 
  #             families = FGM_Cop_BivDiscrete(marg1 = "ZALG", 
  #                                            marg2 = "ZINBI"),
  #             control = boost_control(mstop = MSTOP_COP, nu  = boost.nu.steplength, trace = TRUE), 
  #             method = 'noncyclic')
  bivDiscCopula <- bivDiscCopula[MSTOP_COP]
  
  
  mstop.bivCopula <-  vector('list')
  mstop.bivCopula$mstop <- MSTOP_COP
  mstop.bivCopula$mu1 <- bivDiscCopula$mu1$mstop()
  mstop.bivCopula$sigma1 <- bivDiscCopula$sigma1$mstop()
  
  mstop.bivCopula$mu2 <- bivDiscCopula$mu2$mstop()
  mstop.bivCopula$sigma2 <- bivDiscCopula$sigma2$mstop()
  mstop.bivCopula$nu2 <- bivDiscCopula$nu2$mstop()
  
  mstop.bivCopula$rho  <- bivDiscCopula$rho$mstop()
  
  # EXTRACT ALL COEFFICIENTS!
  coef.bivDiscCopula <- coef(bivDiscCopula, which = "")
  
  DiscreteMargins_Metrics <- vector("list")
  likBinCopula <- vector('list')
  
  
  ################ Get predictions from each covariate?
  newX <- matrix(rep(seq(from = 0, to = 1, length.out = 500),times = 10), nrow = 500)
  colnames(newX) <-  c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
  
  
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ### Copula predictions / selection
  Copula_MU1_predictions <- list()
  Copula_SIGMA1_predictions <- list()
  
  Copula_MU2_predictions <- list()
  Copula_SIGMA2_predictions <- list()
  Copula_NU2_predictions <- list()
  
  Copula_RHO_predictions <- list()
  
  ### MU 1
  Copula_MU1_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = 1, type = "link" )
  Copula_MU1_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x2", type = "link" )
  Copula_MU1_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x3", type = "link" )
  Copula_MU1_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x4", type = "link" )
  Copula_MU1_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x5", type = "link" )
  Copula_MU1_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x6", type = "link" )
  Copula_MU1_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x7", type = "link" )
  Copula_MU1_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x8", type = "link" )
  Copula_MU1_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x9", type = "link" )
  Copula_MU1_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x10", type = "link" )
  
  #### SIGMA 1
  Copula_SIGMA1_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = 1, type = "link" )
  Copula_SIGMA1_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x2", type = "link" )
  Copula_SIGMA1_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x3", type = "link" )
  Copula_SIGMA1_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x4", type = "link" )
  Copula_SIGMA1_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x5", type = "link" )
  Copula_SIGMA1_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x6", type = "link" )
  Copula_SIGMA1_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x7", type = "link" )
  Copula_SIGMA1_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x8", type = "link" )
  Copula_SIGMA1_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x9", type = "link" )
  Copula_SIGMA1_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x10", type = "link" )
  
  
  #### MU 2
  Copula_MU2_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = 1, type = "link" )
  Copula_MU2_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x2", type = "link" )
  Copula_MU2_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x3", type = "link" )
  Copula_MU2_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x4", type = "link" )
  Copula_MU2_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x5", type = "link" )
  Copula_MU2_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x6", type = "link" )
  Copula_MU2_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x7", type = "link" )
  Copula_MU2_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x8", type = "link" )
  Copula_MU2_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x9", type = "link" )
  Copula_MU2_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x10", type = "link" )
  
  #### SIGMA 2
  Copula_SIGMA2_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = 1, type = "link" )
  Copula_SIGMA2_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x2", type = "link" )
  Copula_SIGMA2_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x3", type = "link" )
  Copula_SIGMA2_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x4", type = "link" )
  Copula_SIGMA2_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x5", type = "link" )
  Copula_SIGMA2_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x6", type = "link" )
  Copula_SIGMA2_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x7", type = "link" )
  Copula_SIGMA2_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x8", type = "link" )
  Copula_SIGMA2_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x9", type = "link" )
  Copula_SIGMA2_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x10", type = "link" )
  
  
  #### NU 2
  Copula_NU2_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = 1, type = "link" )
  Copula_NU2_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x2", type = "link" )
  Copula_NU2_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x3", type = "link" )
  Copula_NU2_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x4", type = "link" )
  Copula_NU2_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x5", type = "link" )
  Copula_NU2_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x6", type = "link" )
  Copula_NU2_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x7", type = "link" )
  Copula_NU2_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x8", type = "link" )
  Copula_NU2_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x9", type = "link" )
  Copula_NU2_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x10", type = "link" )
  
  
  #### RHO
  Copula_RHO_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = 1, type = "link" )
  Copula_RHO_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x2", type = "link" )
  Copula_RHO_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x3", type = "link" )
  Copula_RHO_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x4", type = "link" )
  Copula_RHO_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x5", type = "link" )
  Copula_RHO_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x6", type = "link" )
  Copula_RHO_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x7", type = "link" )
  Copula_RHO_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x8", type = "link" )
  Copula_RHO_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x9", type = "link" )
  Copula_RHO_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x10", type = "link" )
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  
  
  
  # Predict distributional quantities
  predCopula.mu1 <- predict(bivDiscCopula$mu1, newdata = dat.test, type = 'response')
  predCopula.sigma1 <- predict(bivDiscCopula$sigma1, newdata = dat.test, type = 'response')
  
  predCopula.mu2 <- predict(bivDiscCopula$mu2, newdata = dat.test, type = 'response')
  predCopula.sigma2 <- predict(bivDiscCopula$sigma2, newdata = dat.test, type = 'response')
  predCopula.nu2 <- predict(bivDiscCopula$nu2, newdata = dat.test, type = 'response')
  
  
  predCopula.rho <- predict(bivDiscCopula$rho, newdata = dat.test, type = 'response')
  
  ######### COMPUTE PERFORMANCE: 
  
  # MSE (margin 1, then margin 2)
  DiscreteMargins_Metrics$mu1 <- mean((as.numeric(predCopula.mu1)-(as.numeric(dat.test$y1)))^2)
  DiscreteMargins_Metrics$mu2 <- mean((as.numeric(predCopula.mu2)-(as.numeric(dat.test$y2)))^2)
  
  # LOG-SCORE (LOG-LIKELIHOOD) (change to the copula likelihood)
  lik$biv <- sum(loss(mu1 = predCopula.mu1, 
                      mu2 = predCopula.mu2, 
                      sigma1 = predCopula.sigma1, 
                      sigma2 = predCopula.sigma2, 
                      nu2 = predCopula.nu2,
                      rho = predCopula.rho, 
                      y = y.test))
  
  
  ###########################################################################################################################
  ################ - univariat - ############################################################################################
  ###########################################################################################################################
  
  mstop.uni <- vector('list')
  coef.uni  <- vector('list')
  coef.uni_1  <- vector('list')
  
  # univariateDiscreteEquation_1 <- list(mu = formula(y1 ~  bbs(x1) + 
  #                                                     bbs(x2) + 
  #                                                     bbs(x3) + 
  #                                                     bbs(x4) + 
  #                                                     bbs(x5) + 
  #                                                     bbs(x6) +
  #                                                     bbs(x7) + 
  #                                                     bbs(x8) + 
  #                                                     bbs(x9) + 
  #                                                     bbs(x10)),
  #                                      sigma = formula(y1 ~  bbs(x1) + 
  #                                                        bbs(x2) + 
  #                                                        bbs(x3) + 
  #                                                        bbs(x4) + 
  #                                                        bbs(x5) + 
  #                                                        bbs(x6) +
  #                                                        bbs(x7) + 
  #                                                        bbs(x8) + 
  #                                                        bbs(x9) + 
  #                                                        bbs(x10))
  # )
  # 
  # 
  # univariateDiscreteEquation_2 <- list(mu = formula(y2 ~  bbs(x1) + 
  #                                                     bbs(x2) + 
  #                                                     bbs(x3) + 
  #                                                     bbs(x4) + 
  #                                                     bbs(x5) + 
  #                                                     bbs(x6) +
  #                                                     bbs(x7) + 
  #                                                     bbs(x8) + 
  #                                                     bbs(x9) + 
  #                                                     bbs(x10)),
  #                                      sigma = formula(y2 ~  bbs(x1) + 
  #                                                        bbs(x2) + 
  #                                                        bbs(x3) + 
  #                                                        bbs(x4) + 
  #                                                        bbs(x5) + 
  #                                                        bbs(x6) +
  #                                                        bbs(x7) + 
  #                                                        bbs(x8) + 
  #                                                        bbs(x9) + 
  #                                                        bbs(x10),
  #                                      ),
  #                                      nu = formula(y2 ~  bbs(x1) + 
  #                                                     bbs(x2) + 
  #                                                     bbs(x3) + 
  #                                                     bbs(x4) + 
  #                                                     bbs(x5) + 
  #                                                     bbs(x6) +
  #                                                     bbs(x7) + 
  #                                                     bbs(x8) + 
  #                                                     bbs(x9) + 
  #                                                     bbs(x10))
  # )
  # 
  
  
  
  univariateDiscreteEquation_1 <- list(mu = formula(y1 ~  bols(x1) +
                                                      bols(x2) +
                                                      bols(x3) +
                                                      bols(x4) +
                                                      bols(x5) +
                                                      bols(x6) +
                                                      bols(x7) +
                                                      bols(x8) +
                                                      bols(x9) +
                                                      bols(x10)),
                                       sigma = formula(y1 ~  bols(x1) +
                                                         bols(x2) +
                                                         bols(x3) +
                                                         bols(x4) +
                                                         bols(x5) +
                                                         bols(x6) +
                                                         bols(x7) +
                                                         bols(x8) +
                                                         bols(x9) +
                                                         bols(x10))
  )


  univariateDiscreteEquation_2 <- list(mu = formula(y2 ~  bols(x1) +
                                                      bols(x2) +
                                                      bols(x3) +
                                                      bols(x4) +
                                                      bols(x5) +
                                                      bols(x6) +
                                                      bols(x7) +
                                                      bols(x8) +
                                                      bols(x9) +
                                                      bols(x10)),
                                       sigma = formula(y2 ~  bols(x1) +
                                                         bols(x2) +
                                                         bols(x3) +
                                                         bols(x4) +
                                                         bols(x5) +
                                                         bols(x6) +
                                                         bols(x7) +
                                                         bols(x8) +
                                                         bols(x9) +
                                                         bols(x10)),
                                       nu = formula(y2 ~  bols(x1) +
                                                      bols(x2) +
                                                      bols(x3) +
                                                      bols(x4) +
                                                      bols(x5) +
                                                      bols(x6) +
                                                      bols(x7) +
                                                      bols(x8) +
                                                      bols(x9) +
                                                      bols(x10))
  )

  
  # - margin 1
  dat.train.mu1 = dat.train[, !(names(dat.train) %in% "y2")]
  dat.test.mu1  = dat.test[ , !(names(dat.test) %in% "y2")]
  
  glm.uni.mu1 <- gamboostLSS(univariateDiscreteEquation_1, 
                             data = dat.train.mu1, 
                             families = as.families(fname = "ZALG", stabilization = "L2"), 
                             control = boost_control(risk = 'oobag', nu = boost.nu.steplength, mstop = 2000, trace = T), 
                             method = "noncyclic",
                             weights = weight.mstop)
  
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  if(mstop.uni$mu1 >= 1990){
    glm.uni.mu1[5000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  if(mstop.uni$mu1 >= 4990){
    glm.uni.mu1[10000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  if(mstop.uni$mu1 >= 9990){
    glm.uni.mu1[20000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  # dat.train_mu1 <- dat.train.mu1[weight.mstop == 1, ]
  # 
  # glm.uni.mu1 <- gamboostLSS(univariateDiscreteEquation_1, 
  #                            data = dat.train_mu1, 
  #                            families = as.families(fname = "ZALG"), 
  #                         control = boost_control(mstop = mstop.uni$mu1, nu  = boost.nu.steplength, trace = T), 
  #                         method = "noncyclic")
  glm.uni.mu1 <- glm.uni.mu1[mstop.uni$mu1]
  
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 1  
  coef.uni$mu1 <- coef(glm.uni.mu1, which = '')
  
  DiscreteMargins_Metrics$Univariate_mu1 <- mean((as.numeric(predict(glm.uni.mu1$mu, newdata = dat.test.mu1, type = "response")) - (as.numeric(dat.test.mu1$y1)))^2)
  
  # - margin 2
  dat.train.mu2 =  dat.train[, !(names(dat.train) %in% "y1")]
  dat.test.mu2 = dat.test[, !(names(dat.test) %in% "y1")]
  
  glm.uni.mu2 <- gamboostLSS(univariateDiscreteEquation_2, 
                             data = dat.train.mu2, 
                             families = as.families(fname = "ZINBI", stabilization = "L2"), 
                             control = boost_control(risk = 'oobag', nu = boost.nu.steplength, trace = TRUE, mstop = 2000), 
                             method = "noncyclic",
                             weights = weight.mstop)
  
  
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  if(mstop.uni$mu2 >= 1990){
    glm.uni.mu2[5000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  if(mstop.uni$mu2 >= 4990){
    glm.uni.mu2[10000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  if(mstop.uni$mu2 >= 9990){
    glm.uni.mu2[20000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  # dat.train_mu2 <- dat.train.mu2[weight.mstop == 1, ]
  # 
  # glm.uni.mu2 <- gamboostLSS(univariateDiscreteEquation_2, 
  #                            data = dat.train_mu2, 
  #                            families = as.families(fname = "ZINBI"), 
  #                         control = boost_control(mstop = mstop.uni$mu2, nu  = boost.nu.steplength, trace=TRUE), method = "noncyclic")
  glm.uni.mu2 <- glm.uni.mu2[mstop.uni$mu2]
  
  mstop.uni$Margin1 <- mstop(glm.uni.mu1)
  
  mstop.uni$Margin2 <- mstop(glm.uni.mu2)
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 2  
  coef.uni$mu2 <- coef(glm.uni.mu2,which = '')
  
  DiscreteMargins_Metrics$Univariate_mu2 <- mean((as.numeric(predict(glm.uni.mu2$mu, newdata = dat.test.mu2, type = "response"))-(as.numeric(dat.test.mu2$y2)))^2)
  
  # Likelihood
  pred.mu1.uni <- as.numeric(predict(glm.uni.mu1$mu, newdata = dat.test.mu1, type = "response"))
  pred.sigma1.uni <- as.numeric(predict(glm.uni.mu1$sigma, newdata = dat.test.mu1, type = "response"))
  
  pred.mu2.uni <- as.numeric(predict(glm.uni.mu2$mu, newdata = dat.test.mu2, type = "response"))
  pred.sigma2.uni <- as.numeric(predict(glm.uni.mu2$sigma, newdata = dat.test.mu2, type = "response"))
  pred.nu2.uni <- as.numeric(predict(glm.uni.mu2$nu, newdata = dat.test.mu2, type = "response"))
  
  
  mu1.uni.loglik <- sum(-dZALG(x = dat.test.mu1$y1, mu = pred.mu1.uni, sigma = pred.sigma1.uni, log = T))
  mu2.uni.loglik <- sum(-dZINBI(x = dat.test.mu2$y2, mu = pred.mu2.uni, sigma = pred.sigma2.uni, nu = pred.nu2.uni, log = T))
  
  lik$uni<- mu1.uni.loglik + mu2.uni.loglik
  
  
  n.ges = c(n.train, n.mstop, n.test)
  pred.ges <- NULL#list(pred.mu1, pred.mu2, pred.or)
  pred.uni <- list(mu1 = pred.mu1.uni, 
                   sigma1 = pred.sigma1.uni,
                   mu2 = pred.mu2.uni,
                   sigma2 = pred.sigma2.uni,
                   nu2 = pred.nu2.uni
  )
  
  pred.cop <- list(mu1 = predCopula.mu1, 
                   sigma1 = predCopula.sigma1,
                   mu2 = predCopula.mu2,
                   sigma2 = predCopula.sigma2,
                   nu2 = predCopula.nu2,
                   rho = predCopula.rho)
  
  
  ##### Predictions:
  UNIVARIATE_MU1_predictions <- vector("list")
  UNIVARIATE_SIGMA1_predictions <- vector("list")
  
  UNIVARIATE_MU2_predictions <- vector("list")
  UNIVARIATE_SIGMA2_predictions <- vector("list")
  UNIVARIATE_NU2_predictions <- vector("list")
  
  
  #### MU 1
  UNIVARIATE_MU1_predictions$x1 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = 1, type = "link" )
  UNIVARIATE_MU1_predictions$x2 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x2", type = "link" )
  UNIVARIATE_MU1_predictions$x3 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x3", type = "link" )
  UNIVARIATE_MU1_predictions$x4 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x4", type = "link" )
  UNIVARIATE_MU1_predictions$x5 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x5", type = "link" )
  UNIVARIATE_MU1_predictions$x6 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x6", type = "link" )
  UNIVARIATE_MU1_predictions$x7 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x7", type = "link" )
  UNIVARIATE_MU1_predictions$x8 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x8", type = "link" )
  UNIVARIATE_MU1_predictions$x9 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x9", type = "link" )
  UNIVARIATE_MU1_predictions$x10 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x10", type = "link" )
  
  #### SIGMA 1
  UNIVARIATE_SIGMA1_predictions$x1 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = 1, type = "link" )
  UNIVARIATE_SIGMA1_predictions$x2 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x2", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x3 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x3", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x4 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x4", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x5 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x5", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x6 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x6", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x7 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x7", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x8 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x8", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x9 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x9", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x10 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x10", type = "link" )
  
  
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
  
  
  #### NU 2
  UNIVARIATE_NU2_predictions$x1 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = 1, type = "link" )
  UNIVARIATE_NU2_predictions$x2 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x2", type = "link" )
  UNIVARIATE_NU2_predictions$x3 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x3", type = "link" )
  UNIVARIATE_NU2_predictions$x4 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x4", type = "link" )
  UNIVARIATE_NU2_predictions$x5 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x5", type = "link" )
  UNIVARIATE_NU2_predictions$x6 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x6", type = "link" )
  UNIVARIATE_NU2_predictions$x7 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x7", type = "link" )
  UNIVARIATE_NU2_predictions$x8 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x8", type = "link" )
  UNIVARIATE_NU2_predictions$x9 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x9", type = "link" )
  UNIVARIATE_NU2_predictions$x10 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x10", type = "link" )
  
  
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
    sample_uni <- data.gen.bivdisc_energyscore(n = 1000, 
                                               mu1 = pred.mu1.uni[i], 
                                               mu2 = pred.mu2.uni[i], 
                                               sigma1 = pred.sigma1.uni[i], 
                                               sigma2 = pred.sigma2.uni[i],
                                               nu = check_extremes_nuZINBI(pred.nu2.uni[i]),
                                               theta = 1)
    
    pred_sample_uni[1,] <- sample_uni[,1]
    pred_sample_uni[2,] <- sample_uni[,2]
    
    es_uni[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_uni) 
    
    
    # copula
    sample_cop <- data.gen.bivdisc_energyscore(n = 1000, 
                                               mu1 = predCopula.mu1[i],  
                                               mu2 = predCopula.mu2[i], 
                                               sigma1 = predCopula.sigma1[i], 
                                               sigma2 = predCopula.sigma2[i],
                                               nu = check_extremes_nuZINBI(predCopula.nu2[i]),
                                               theta = predCopula.rho[i])
    
    pred_sample_cop[1, ] <- sample_cop[,1]
    pred_sample_cop[2, ] <- sample_cop[,2]
    
    es_cop[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_cop) 
    
    
  }
  
  energy_score <- list()
  #energy_score$biv <- mean(es_biv)
  energy_score$uni <- mean(es_uni)
  energy_score$cop <- mean(es_cop)
  
  
  
  ########################################### We gather: mstop, predictions, MSE margin 1, margin 2, loglikelihood, energy score.
  output <- list(TrueBeta = TrueBeta, 
                 n = n.ges, 
                 p  = p, 
                 Likelihood = lik,
                 predict.uni = pred.uni, 
                 predict.cop = pred.cop,
                 energy_score = energy_score,  
                 ####
                 Coefficients.uni = coef.uni, 
                 mstop.uni = mstop.uni,
                 ###
                 DiscreteMarginMetrics = DiscreteMargins_Metrics,
                 ###
                 Copula_Predictions = list(MU1 = Copula_MU1_predictions, 
                                           SIGMA1 = Copula_SIGMA1_predictions,
                                           MU2 = Copula_MU2_predictions,
                                           SIGMA2 = Copula_SIGMA2_predictions,
                                           NU2 = Copula_NU2_predictions,
                                           RHO = Copula_RHO_predictions),
                 ###
                 Univariate_Predictions = list(MU1 = UNIVARIATE_MU1_predictions, 
                                               SIGMA1 = UNIVARIATE_SIGMA1_predictions,
                                               MU2 = UNIVARIATE_MU2_predictions, 
                                               SIGMA2 = UNIVARIATE_SIGMA2_predictions,
                                               NU2 = UNIVARIATE_NU2_predictions),
                 ###
                 CoefficientsCOPULA = coef.bivDiscCopula, 
                 oobag.riskCOPULA = oobag.risk.cop, 
                 energy_scoreCOPULA = energy_score,  
                 mstopCOPULA = mstop.bivCopula,
                 TrueKendallRange = TrueKendallTauRange
  )
  
  
  return(output)
  
}

#### Different covariate configuration as well as different functions for some parameters:
# I changed the function in sigma2
sim_JOE_2 <- function(seed, n.train, n.mstop, p, n.test, corr, boost.nu.steplength = 0.1){
  
  
  loss <- function(mu1, mu2, sigma1, sigma2, rho, y, nu2 = NULL){
    
    
    # Margin 1 is ZALG:
    F1 <- pdffz(pZALG(y[,1], mu = mu1, sigma = sigma1))
    pdf1 <- pdffz(dZALG(y[,1], mu = mu1, sigma = sigma1))
    
    # Margin 2 is NBI:
    F2 <- pdffz(pZINBI(y[,2], mu = mu2, sigma = sigma2, nu = nu2))
    pdf2 <- pdffz(dZINBI(y[,2], mu = mu2, sigma = sigma2, nu = nu2))
    
    # Copula parameter 
    thet <- rho
    
    # Minus 1
    CDF1m1 <- pdffz(F1 - pdf1)
    CDF2m1 <- pdffz(F2 - pdf2)
    
    ### Copula Terms
    bit1 <- (1 - F1)^thet
    bit1m1 <- (1 - CDF1m1)^thet
    
    bit2 <- (1 - F2)^thet
    bit2m1 <- (1 - CDF2m1)^thet
    
    
    T1 <- pdffz(  1 - (bit1 + bit2 - bit1*bit2)^(1/thet)  )
    
    T2 <- pdffz(  1 - (bit1m1 + bit2 - bit1m1*bit2)^(1/thet)  ) 
    
    T3 <- pdffz(  1 - (bit1 + bit2m1 - bit1*bit2m1)^(1/thet)   )
    
    T4 <- pdffz(  1 - (bit1m1 + bit2m1 - bit1m1*bit2m1)^(1/thet)  ) 
    
    L_Total <- pdffz( T1 - T2 - T3 + T4 )
    
    # Negative log-likelihood
    return(- log( L_Total ) )
    
  }
  
  data.gen.bivdisc <- function(mu1, mu2, sigma1, sigma2, theta, nu2 = NULL){
    
    y1y2 <- matrix(0, ncol = 2, nrow = length(mu1))
    
    for(i in 1:length(mu1)){
      
      paramlist1 <- list(mu = mu1[i], 
                         sigma = sigma1[i])
      
      
      paramlist2 <- list(mu = mu2[i], 
                         sigma = sigma2[i],
                         nu = nu2[i])
      
      copObj <- copula::joeCopula(param = theta[i])
      
      copThing <- mvdc(copObj, margins = c("ZALG", "ZINBI"), paramMargins = list(paramlist1, paramlist2))
      
      y1y2[i,] <- rMvdc(copThing, n = 1)
    }
    
    
    dat <- data.frame(y1y2[,1], y1y2[,2])
    
    return(dat)
  }
  
  data.gen.bivdisc_energyscore <- function(n, mu1, mu2, sigma1, sigma2, theta, nu = NULL){
    
    
    paramlist1 <- list(mu = mu1, sigma = sigma1)
    
    paramlist2 <- list(mu = mu2, sigma = sigma2, nu = nu) 
    
    copObj <- copula::joeCopula(param = theta)
    
    copThing <- mvdc(copObj, margins = c("ZALG", "ZINBI"), paramMargins = list(paramlist1, paramlist2))
    
    y1y2 <- rMvdc(copThing, n = n)
    
    dat <- data.frame(y1y2[,1], y1y2[,2])
    
    
    return(dat)
  }
  
  
  
  
  ## Coefficients:
  # x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
  mu1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  mu2 = c(0, 0, 0, 0, 0,  0, 0, 0, 0, 0) 
  
  # x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
  sigma1 <- c(0, 0, 0, 0, 0,  0, 0, 0, 0, 0) 
  sigma2 <- c(0, 0, 0, 0, 0,  0, 0, 0, 0, 0) 
  
  nu2 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0) 
  
  # x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
  or  = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  
  ### function for mu1
  nl_function_mu1 <- function(x){
    
    
    nl <- sqrt(x)*x^2 - 2*cos(3*x)
    
    return(nl*0.5)
    
  }
  
  ### function for sigma1
  nl_function_sigma1 <- function(x){
    
    func <- -40*(x^(3/2) - x^(4/3))
    
    func <- 2 * func
    
    return(func)
  }
  
  nl_func_sigma1 <- function(x){
    
    func <- -(log(x)) + cos(x)
    
    func <- 2 * func
    
    return(func)
  }
  
  ## function for mu2
  nl_function_mu2 <- function(x){
    
    nl <- -0.7*exp(x^2) + exp(x^0.4) - 0*log(x+0.2) 
    
    return(nl)
  }
  
  nl_function_sigma2 <- function(x){
    
    return( -1.5*(1.5*cos(2 * x) - 0*log(x + 0.15) + 3*tanh(1*x) ) )
  }
  
  nl_function_nu2 <- function(x){
    
    func <- -(sin(x) - exp(x)^2) 
    
    #func <- func + 2*log(x + 0.01) - 3*log(x^3 + 0.01)
    
    func <- 0.7*func
    
    return(func)
  }
  
  nl_function_rho <- function(x){
    
    #return(2*sin(x*4))
    return(2*sin(x*4))
  }
  
  TrueBeta <-  vector('list')
  TrueBeta$mu1 <- mu1
  TrueBeta$mu2 <- mu2
  
  TrueBeta$sigma1 <- sigma1
  TrueBeta$sigma2 <- sigma2
  
  TrueBeta$nu2 <- nu2
  
  TrueBeta$or <- or
  
  set.seed(seed)
  
  
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
  
  ############################################################################################################################################
  ############################################################################################################################################
  # - test data
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
  
  ############################################################################################################################################
  # GENERATION OF DISTRIBUTION PARAMETERS AND STUFF 
  # predictor
  
  ##### TRAIN DATA
  train.eta.mu1_center <- (x.train %*% mu1)
  train.eta.mu1_center <- train.eta.mu1_center + nl_function_mu1(x1)
  
  ##### TEST DATA
  test.eta.mu1_center <- (x.test %*% mu1)
  test.eta.mu1_center <- test.eta.mu1_center + nl_function_mu1(x1.test)
  
  # apply response function:
  train.eta.mu1<-  plogis(train.eta.mu1_center) 
  range(train.eta.mu1)
  
  test.eta.mu1<-  plogis(test.eta.mu1_center) 
  ########################################################################################################################
  #### SIGMA FOR ZALG RESPONSE
  train.eta.sigma1 <-  x.train %*% sigma1
  
  train.eta.sigma1 <- train.eta.sigma1 - 2 + nl_function_sigma1(x3) 
  train.eta.sigma1 <- plogis(train.eta.sigma1)
  range(train.eta.sigma1)
  
  ##### TEST DATA
  test.eta.sigma1 <-  x.test %*% sigma1
  test.eta.sigma1 <- test.eta.sigma1 - 2 + nl_function_sigma1(x3.test)
  test.eta.sigma1 <- plogis(test.eta.sigma1)
  ########################################################################################################################
  # predictor and apply response function: 
  train.eta.mu2 <-  x.train %*% mu2
  train.eta.mu2 <-  train.eta.mu2 + nl_function_mu2(x2)
  train.eta.mu2 <-  exp(train.eta.mu2) 
  range(train.eta.mu2)
  
  
  test.eta.mu2 <-  x.test %*% mu2
  test.eta.mu2 <-  test.eta.mu2 + nl_function_mu2(x2.test)
  test.eta.mu2 <-  exp(test.eta.mu2) 
  ########################################################################################################################
  # predictor and apply response function: 
  train.eta.sigma2 <-  x.train %*% sigma2
  train.eta.sigma2 <-  train.eta.sigma2 + nl_function_sigma2(x5) + 3
  train.eta.sigma2 <-  exp(train.eta.sigma2) 
  range(train.eta.sigma2)
  
  test.eta.sigma2 <-  x.test %*% sigma2
  test.eta.sigma2 <-  test.eta.sigma2 + nl_function_sigma2(x5.test) + 3
  test.eta.sigma2 <-  exp(test.eta.sigma2) 
  ########################################################################################################################
  # predictor and apply response function: 
  train.eta.nu2 <-  x.train %*% nu2
  train.eta.nu2 <-  train.eta.nu2 + nl_function_nu2(x1) - 3
  train.eta.nu2 <-  plogis(train.eta.nu2) 
  range(train.eta.nu2)
  
  test.eta.nu2 <-  x.test %*% nu2
  test.eta.nu2 <-  test.eta.nu2 + nl_function_nu2(x1.test) - 3
  test.eta.nu2 <-  plogis(test.eta.nu2) 
  
  ########################################################################################################################
  ### Copula parameter
  train.eta.or_center <- (x.train %*% or)
  train.eta.or_center <- train.eta.or_center + nl_function_rho(x4) 
  
  train.copula.parameter <- (exp(train.eta.or_center) + 1)
  range(train.copula.parameter)
  
  TrueKendallTauRange <- VineCopula::BiCopPar2Tau(6, par = range(train.copula.parameter))  
  
  test.eta.or_center <- (x.test %*% or)
  test.eta.or_center <- test.eta.or_center + nl_function_rho(x4.test)
  
  test.copula.parameter <- (exp(test.eta.or_center) + 1)
  ############################################################################################################################################
  # SAMPLING OF BIVARIATE RESPONSE FROM COPULA
  y.train <- data.gen.bivdisc(mu1 = train.eta.mu1, mu2 = train.eta.mu2, sigma1 = train.eta.sigma1, sigma2 = train.eta.sigma2, theta = train.copula.parameter, nu2 = train.eta.nu2)
  
  colnames(y.train)<- c('y1','y2')
  dat.train <- data.frame(y.train,x.train)
  
  
  y.test <- data.gen.bivdisc(mu1 = test.eta.mu1, mu2 = test.eta.mu2, sigma1 = test.eta.sigma1, sigma2 = test.eta.sigma2, theta = test.copula.parameter, nu2 = test.eta.nu2)
  
  colnames(y.test)<- c('y1','y2')
  dat.test <- data.frame(y.test,x.test)
  
  ############################################################################################################################################
  # FITTIIIINGG
  mstop.bivBern <-  vector('list')
  
  AUCbivBern <-  vector('list')
  BrierbivBern <-  vector('list')
  lik <- vector('list')
  MSEbivDisc <- vector("list")
  
  
  
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
                           sigma1 = model_equation,
                           sigma2 = model_equation,
                           rho = model_equation
  )
  
  
  bivModel_Formula <- list(mu1 = model_equation,
                           mu2 = model_equation,
                           sigma1 = model_equation,
                           sigma2 = model_equation,
                           nu2 = model_equation,
                           rho = model_equation
  )
  
  #############################################
  ################# - COPULA - ################
  #############################################
  # Bivariate copula model
  bivDiscCopula <- gamboostLSS(bivModel_Formula, data = dat.train, 
                               families = Joe_Cop_BivDiscrete(marg1 = "ZALG", 
                                                              marg2 = "ZINBI",
                                                              stabilization = "L2"), 
                               control = boost_control(mstop = 2000, risk = 'oobag',
                                                       nu  = boost.nu.steplength, trace = TRUE), 
                               method = 'noncyclic', weights = weight.mstop)
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 1900){
    mstop(bivDiscCopula) <- 3000
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 2900){
    mstop(bivDiscCopula) <- 4000
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 3990){
    bivDiscCopula[8000]
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 7990){
    bivDiscCopula[10000]
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 9990){
    bivDiscCopula[15000]
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 14990){
    bivDiscCopula[20000]
  }
  
  
  
  
  MSTOP_COP <- which.min(risk(bivDiscCopula, merge = T))
  oobag.risk.cop <- risk(bivDiscCopula,merge = T)
  
  # rm(bivDiscCopula)
  # dat.train_biv <- dat.train[weight.mstop == 1, ]
  # 
  # # RE-FIT until OPTIMAL MSTOP
  # bivDiscCopula <- gamboostLSS(bivModel_Formula, data = dat.train_biv, 
  #             families = FGM_Cop_BivDiscrete(marg1 = "ZALG", 
  #                                            marg2 = "ZINBI"),
  #             control = boost_control(mstop = MSTOP_COP, nu  = boost.nu.steplength, trace = TRUE), 
  #             method = 'noncyclic')
  bivDiscCopula <- bivDiscCopula[MSTOP_COP]
  
  
  mstop.bivCopula <-  vector('list')
  mstop.bivCopula$mstop <- MSTOP_COP
  mstop.bivCopula$mu1 <- bivDiscCopula$mu1$mstop()
  mstop.bivCopula$sigma1 <- bivDiscCopula$sigma1$mstop()
  
  mstop.bivCopula$mu2 <- bivDiscCopula$mu2$mstop()
  mstop.bivCopula$sigma2 <- bivDiscCopula$sigma2$mstop()
  mstop.bivCopula$nu2 <- bivDiscCopula$nu2$mstop()
  
  mstop.bivCopula$rho  <- bivDiscCopula$rho$mstop()
  
  # EXTRACT ALL COEFFICIENTS!
  coef.bivDiscCopula <- coef(bivDiscCopula, which = "")
  
  DiscreteMargins_Metrics <- vector("list")
  likBinCopula <- vector('list')
  
  
  ################ Get predictions from each covariate?
  newX <- matrix(rep(seq(from = 0, to = 1, length.out = 500),times = 10), nrow = 500)
  colnames(newX) <-  c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
  
  
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ### Copula predictions / selection
  Copula_MU1_predictions <- list()
  Copula_SIGMA1_predictions <- list()
  
  Copula_MU2_predictions <- list()
  Copula_SIGMA2_predictions <- list()
  Copula_NU2_predictions <- list()
  
  Copula_RHO_predictions <- list()
  
  ### MU 1
  Copula_MU1_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = 1, type = "link" )
  Copula_MU1_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x2", type = "link" )
  Copula_MU1_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x3", type = "link" )
  Copula_MU1_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x4", type = "link" )
  Copula_MU1_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x5", type = "link" )
  Copula_MU1_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x6", type = "link" )
  Copula_MU1_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x7", type = "link" )
  Copula_MU1_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x8", type = "link" )
  Copula_MU1_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x9", type = "link" )
  Copula_MU1_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x10", type = "link" )
  
  #### SIGMA 1
  Copula_SIGMA1_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = 1, type = "link" )
  Copula_SIGMA1_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x2", type = "link" )
  Copula_SIGMA1_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x3", type = "link" )
  Copula_SIGMA1_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x4", type = "link" )
  Copula_SIGMA1_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x5", type = "link" )
  Copula_SIGMA1_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x6", type = "link" )
  Copula_SIGMA1_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x7", type = "link" )
  Copula_SIGMA1_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x8", type = "link" )
  Copula_SIGMA1_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x9", type = "link" )
  Copula_SIGMA1_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x10", type = "link" )
  
  
  #### MU 2
  Copula_MU2_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = 1, type = "link" )
  Copula_MU2_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x2", type = "link" )
  Copula_MU2_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x3", type = "link" )
  Copula_MU2_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x4", type = "link" )
  Copula_MU2_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x5", type = "link" )
  Copula_MU2_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x6", type = "link" )
  Copula_MU2_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x7", type = "link" )
  Copula_MU2_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x8", type = "link" )
  Copula_MU2_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x9", type = "link" )
  Copula_MU2_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x10", type = "link" )
  
  #### SIGMA 2
  Copula_SIGMA2_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = 1, type = "link" )
  Copula_SIGMA2_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x2", type = "link" )
  Copula_SIGMA2_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x3", type = "link" )
  Copula_SIGMA2_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x4", type = "link" )
  Copula_SIGMA2_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x5", type = "link" )
  Copula_SIGMA2_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x6", type = "link" )
  Copula_SIGMA2_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x7", type = "link" )
  Copula_SIGMA2_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x8", type = "link" )
  Copula_SIGMA2_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x9", type = "link" )
  Copula_SIGMA2_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x10", type = "link" )
  
  
  #### NU 2
  Copula_NU2_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = 1, type = "link" )
  Copula_NU2_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x2", type = "link" )
  Copula_NU2_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x3", type = "link" )
  Copula_NU2_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x4", type = "link" )
  Copula_NU2_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x5", type = "link" )
  Copula_NU2_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x6", type = "link" )
  Copula_NU2_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x7", type = "link" )
  Copula_NU2_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x8", type = "link" )
  Copula_NU2_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x9", type = "link" )
  Copula_NU2_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x10", type = "link" )
  
  
  #### RHO
  Copula_RHO_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = 1, type = "link" )
  Copula_RHO_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x2", type = "link" )
  Copula_RHO_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x3", type = "link" )
  Copula_RHO_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x4", type = "link" )
  Copula_RHO_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x5", type = "link" )
  Copula_RHO_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x6", type = "link" )
  Copula_RHO_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x7", type = "link" )
  Copula_RHO_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x8", type = "link" )
  Copula_RHO_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x9", type = "link" )
  Copula_RHO_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x10", type = "link" )
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  
  
  
  # Predict distributional quantities
  predCopula.mu1 <- predict(bivDiscCopula$mu1, newdata = dat.test, type = 'response')
  predCopula.sigma1 <- predict(bivDiscCopula$sigma1, newdata = dat.test, type = 'response')
  
  predCopula.mu2 <- predict(bivDiscCopula$mu2, newdata = dat.test, type = 'response')
  predCopula.sigma2 <- predict(bivDiscCopula$sigma2, newdata = dat.test, type = 'response')
  predCopula.nu2 <- predict(bivDiscCopula$nu2, newdata = dat.test, type = 'response')
  
  
  predCopula.rho <- predict(bivDiscCopula$rho, newdata = dat.test, type = 'response')
  
  ######### COMPUTE PERFORMANCE: 
  
  # MSE (margin 1, then margin 2)
  DiscreteMargins_Metrics$mu1 <- mean((as.numeric(predCopula.mu1)-(as.numeric(dat.test$y1)))^2)
  DiscreteMargins_Metrics$mu2 <- mean((as.numeric(predCopula.mu2)-(as.numeric(dat.test$y2)))^2)
  
  # LOG-SCORE (LOG-LIKELIHOOD) (change to the copula likelihood)
  lik$biv <- sum(loss(mu1 = predCopula.mu1, 
                      mu2 = predCopula.mu2, 
                      sigma1 = predCopula.sigma1, 
                      sigma2 = predCopula.sigma2, 
                      nu2 = predCopula.nu2,
                      rho = predCopula.rho, 
                      y = y.test))
  
  
  ###########################################################################################################################
  ################ - univariat - ############################################################################################
  ###########################################################################################################################
  
  mstop.uni <- vector('list')
  coef.uni  <- vector('list')
  coef.uni_1  <- vector('list')
  
  univariateDiscreteEquation_1 <- list(mu = formula(y1 ~  bbs(x1) + 
                                                      bbs(x2) + 
                                                      bbs(x3) + 
                                                      bbs(x4) + 
                                                      bbs(x5) + 
                                                      bbs(x6) +
                                                      bbs(x7) + 
                                                      bbs(x8) + 
                                                      bbs(x9) + 
                                                      bbs(x10)),
                                       sigma = formula(y1 ~  bbs(x1) + 
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
  
  
  univariateDiscreteEquation_2 <- list(mu = formula(y2 ~  bbs(x1) + 
                                                      bbs(x2) + 
                                                      bbs(x3) + 
                                                      bbs(x4) + 
                                                      bbs(x5) + 
                                                      bbs(x6) +
                                                      bbs(x7) + 
                                                      bbs(x8) + 
                                                      bbs(x9) + 
                                                      bbs(x10)),
                                       sigma = formula(y2 ~  bbs(x1) + 
                                                         bbs(x2) + 
                                                         bbs(x3) + 
                                                         bbs(x4) + 
                                                         bbs(x5) + 
                                                         bbs(x6) +
                                                         bbs(x7) + 
                                                         bbs(x8) + 
                                                         bbs(x9) + 
                                                         bbs(x10),
                                       ),
                                       nu = formula(y2 ~  bbs(x1) + 
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
  
  
  
  
  # univariateDiscreteEquation_1 <- list(mu = formula(y1 ~  bols(x1) + 
  #                                                     bols(x2) + 
  #                                                     bols(x3) + 
  #                                                     bols(x4) + 
  #                                                     bols(x5) + 
  #                                                     bols(x6) +
  #                                                     bols(x7) + 
  #                                                     bols(x8) + 
  #                                                     bols(x9) + 
  #                                                     bols(x10)),
  #                                      sigma = formula(y1 ~  bols(x1) + 
  #                                                        bols(x2) + 
  #                                                        bols(x3) + 
  #                                                        bols(x4) + 
  #                                                        bols(x5) + 
  #                                                        bols(x6) +
  #                                                        bols(x7) + 
  #                                                        bols(x8) + 
  #                                                        bols(x9) + 
  #                                                        bols(x10))
  # )
  # 
  # 
  # univariateDiscreteEquation_2 <- list(mu = formula(y2 ~  bols(x1) + 
  #                                                     bols(x2) + 
  #                                                     bols(x3) + 
  #                                                     bols(x4) + 
  #                                                     bols(x5) + 
  #                                                     bols(x6) +
  #                                                     bols(x7) + 
  #                                                     bols(x8) + 
  #                                                     bols(x9) + 
  #                                                     bols(x10)),
  #                                      sigma = formula(y2 ~  bols(x1) + 
  #                                                        bols(x2) + 
  #                                                        bols(x3) + 
  #                                                        bols(x4) + 
  #                                                        bols(x5) + 
  #                                                        bols(x6) +
  #                                                        bols(x7) + 
  #                                                        bols(x8) + 
  #                                                        bols(x9) + 
  #                                                        bols(x10)),
  #                                      nu = formula(y2 ~  bols(x1) + 
  #                                                     bols(x2) + 
  #                                                     bols(x3) + 
  #                                                     bols(x4) + 
  #                                                     bols(x5) + 
  #                                                     bols(x6) +
  #                                                     bols(x7) + 
  #                                                     bols(x8) + 
  #                                                     bols(x9) + 
  #                                                     bols(x10))
  # )
  # 
  
  # - margin 1
  dat.train.mu1 = dat.train[, !(names(dat.train) %in% "y2")]
  dat.test.mu1  = dat.test[ , !(names(dat.test) %in% "y2")]
  
  glm.uni.mu1 <- gamboostLSS(univariateDiscreteEquation_1, 
                             data = dat.train.mu1, 
                             families = as.families(fname = "ZALG", stabilization = "L2"), 
                             control = boost_control(risk = 'oobag', nu = boost.nu.steplength, mstop = 2000, trace = T), 
                             method = "noncyclic",
                             weights = weight.mstop)
  
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  if(mstop.uni$mu1 >= 1990){
    glm.uni.mu1[5000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  if(mstop.uni$mu1 >= 4990){
    glm.uni.mu1[10000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  if(mstop.uni$mu1 >= 9990){
    glm.uni.mu1[20000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  # dat.train_mu1 <- dat.train.mu1[weight.mstop == 1, ]
  # 
  # glm.uni.mu1 <- gamboostLSS(univariateDiscreteEquation_1, 
  #                            data = dat.train_mu1, 
  #                            families = as.families(fname = "ZALG"), 
  #                         control = boost_control(mstop = mstop.uni$mu1, nu  = boost.nu.steplength, trace = T), 
  #                         method = "noncyclic")
  glm.uni.mu1 <- glm.uni.mu1[mstop.uni$mu1]
  
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 1  
  coef.uni$mu1 <- coef(glm.uni.mu1, which = '')
  
  DiscreteMargins_Metrics$Univariate_mu1 <- mean((as.numeric(predict(glm.uni.mu1$mu, newdata = dat.test.mu1, type = "response")) - (as.numeric(dat.test.mu1$y1)))^2)
  
  # - margin 2
  dat.train.mu2 =  dat.train[, !(names(dat.train) %in% "y1")]
  dat.test.mu2 = dat.test[, !(names(dat.test) %in% "y1")]
  
  glm.uni.mu2 <- gamboostLSS(univariateDiscreteEquation_2, 
                             data = dat.train.mu2, 
                             families = as.families(fname = "ZINBI", stabilization = "L2"), 
                             control = boost_control(risk = 'oobag', nu = boost.nu.steplength, trace = TRUE, mstop = 2000), 
                             method = "noncyclic",
                             weights = weight.mstop)
  
  
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  if(mstop.uni$mu2 >= 1990){
    glm.uni.mu2[5000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  if(mstop.uni$mu2 >= 4990){
    glm.uni.mu2[10000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  if(mstop.uni$mu2 >= 9990){
    glm.uni.mu2[20000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  # dat.train_mu2 <- dat.train.mu2[weight.mstop == 1, ]
  # 
  # glm.uni.mu2 <- gamboostLSS(univariateDiscreteEquation_2, 
  #                            data = dat.train_mu2, 
  #                            families = as.families(fname = "ZINBI"), 
  #                         control = boost_control(mstop = mstop.uni$mu2, nu  = boost.nu.steplength, trace=TRUE), method = "noncyclic")
  glm.uni.mu2 <- glm.uni.mu2[mstop.uni$mu2]
  
  mstop.uni$Margin1 <- mstop(glm.uni.mu1)
  
  mstop.uni$Margin2 <- mstop(glm.uni.mu2)
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 2  
  coef.uni$mu2 <- coef(glm.uni.mu2,which = '')
  
  DiscreteMargins_Metrics$Univariate_mu2 <- mean((as.numeric(predict(glm.uni.mu2$mu, newdata = dat.test.mu2, type = "response"))-(as.numeric(dat.test.mu2$y2)))^2)
  
  # Likelihood
  pred.mu1.uni <- as.numeric(predict(glm.uni.mu1$mu, newdata = dat.test.mu1, type = "response"))
  pred.sigma1.uni <- as.numeric(predict(glm.uni.mu1$sigma, newdata = dat.test.mu1, type = "response"))
  
  pred.mu2.uni <- as.numeric(predict(glm.uni.mu2$mu, newdata = dat.test.mu2, type = "response"))
  pred.sigma2.uni <- as.numeric(predict(glm.uni.mu2$sigma, newdata = dat.test.mu2, type = "response"))
  pred.nu2.uni <- as.numeric(predict(glm.uni.mu2$nu, newdata = dat.test.mu2, type = "response"))
  
  
  mu1.uni.loglik <- sum(-dZALG(x = dat.test.mu1$y1, mu = pred.mu1.uni, sigma = pred.sigma1.uni, log = T))
  mu2.uni.loglik <- sum(-dZINBI(x = dat.test.mu2$y2, mu = pred.mu2.uni, sigma = pred.sigma2.uni, nu = pred.nu2.uni, log = T))
  
  lik$uni<- mu1.uni.loglik + mu2.uni.loglik
  
  
  n.ges = c(n.train, n.mstop, n.test)
  pred.ges <- NULL#list(pred.mu1, pred.mu2, pred.or)
  pred.uni <- list(mu1 = pred.mu1.uni, 
                   sigma1 = pred.sigma1.uni,
                   mu2 = pred.mu2.uni,
                   sigma2 = pred.sigma2.uni,
                   nu2 = pred.nu2.uni
  )
  
  pred.cop <- list(mu1 = predCopula.mu1, 
                   sigma1 = predCopula.sigma1,
                   mu2 = predCopula.mu2,
                   sigma2 = predCopula.sigma2,
                   nu2 = predCopula.nu2,
                   rho = predCopula.rho)
  
  
  ##### Predictions:
  UNIVARIATE_MU1_predictions <- vector("list")
  UNIVARIATE_SIGMA1_predictions <- vector("list")
  
  UNIVARIATE_MU2_predictions <- vector("list")
  UNIVARIATE_SIGMA2_predictions <- vector("list")
  UNIVARIATE_NU2_predictions <- vector("list")
  
  
  #### MU 1
  UNIVARIATE_MU1_predictions$x1 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = 1, type = "link" )
  UNIVARIATE_MU1_predictions$x2 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x2", type = "link" )
  UNIVARIATE_MU1_predictions$x3 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x3", type = "link" )
  UNIVARIATE_MU1_predictions$x4 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x4", type = "link" )
  UNIVARIATE_MU1_predictions$x5 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x5", type = "link" )
  UNIVARIATE_MU1_predictions$x6 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x6", type = "link" )
  UNIVARIATE_MU1_predictions$x7 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x7", type = "link" )
  UNIVARIATE_MU1_predictions$x8 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x8", type = "link" )
  UNIVARIATE_MU1_predictions$x9 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x9", type = "link" )
  UNIVARIATE_MU1_predictions$x10 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x10", type = "link" )
  
  #### SIGMA 1
  UNIVARIATE_SIGMA1_predictions$x1 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = 1, type = "link" )
  UNIVARIATE_SIGMA1_predictions$x2 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x2", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x3 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x3", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x4 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x4", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x5 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x5", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x6 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x6", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x7 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x7", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x8 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x8", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x9 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x9", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x10 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x10", type = "link" )
  
  
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
  
  
  #### NU 2
  UNIVARIATE_NU2_predictions$x1 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = 1, type = "link" )
  UNIVARIATE_NU2_predictions$x2 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x2", type = "link" )
  UNIVARIATE_NU2_predictions$x3 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x3", type = "link" )
  UNIVARIATE_NU2_predictions$x4 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x4", type = "link" )
  UNIVARIATE_NU2_predictions$x5 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x5", type = "link" )
  UNIVARIATE_NU2_predictions$x6 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x6", type = "link" )
  UNIVARIATE_NU2_predictions$x7 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x7", type = "link" )
  UNIVARIATE_NU2_predictions$x8 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x8", type = "link" )
  UNIVARIATE_NU2_predictions$x9 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x9", type = "link" )
  UNIVARIATE_NU2_predictions$x10 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x10", type = "link" )
  
  
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
    sample_uni <- data.gen.bivdisc_energyscore(n = 1000, 
                                               mu1 = pred.mu1.uni[i], 
                                               mu2 = pred.mu2.uni[i], 
                                               sigma1 = pred.sigma1.uni[i], 
                                               sigma2 = pred.sigma2.uni[i],
                                               nu = check_extremes_nuZINBI(pred.nu2.uni[i]),
                                               theta = 1)
    
    pred_sample_uni[1,] <- sample_uni[,1]
    pred_sample_uni[2,] <- sample_uni[,2]
    
    es_uni[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_uni) 
    
    
    # copula
    sample_cop <- data.gen.bivdisc_energyscore(n = 1000, 
                                               mu1 = predCopula.mu1[i],  
                                               mu2 = predCopula.mu2[i], 
                                               sigma1 = predCopula.sigma1[i], 
                                               sigma2 = predCopula.sigma2[i],
                                               nu = check_extremes_nuZINBI(predCopula.nu2[i]),
                                               theta = predCopula.rho[i])
    
    pred_sample_cop[1, ] <- sample_cop[,1]
    pred_sample_cop[2, ] <- sample_cop[,2]
    
    es_cop[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_cop) 
    
    
  }
  
  energy_score <- list()
  #energy_score$biv <- mean(es_biv)
  energy_score$uni <- mean(es_uni)
  energy_score$cop <- mean(es_cop)
  
  
  
  ########################################### We gather: mstop, predictions, MSE margin 1, margin 2, loglikelihood, energy score.
  output <- list(TrueBeta = TrueBeta, 
                 n = n.ges, 
                 p  = p, 
                 Likelihood = lik,
                 predict.uni = pred.uni, 
                 predict.cop = pred.cop,
                 energy_score = energy_score,  
                 ####
                 Coefficients.uni = coef.uni, 
                 mstop.uni = mstop.uni,
                 ###
                 DiscreteMarginMetrics = DiscreteMargins_Metrics,
                 ###
                 Copula_Predictions = list(MU1 = Copula_MU1_predictions, 
                                           SIGMA1 = Copula_SIGMA1_predictions,
                                           MU2 = Copula_MU2_predictions,
                                           SIGMA2 = Copula_SIGMA2_predictions,
                                           NU2 = Copula_NU2_predictions,
                                           RHO = Copula_RHO_predictions),
                 ###
                 Univariate_Predictions = list(MU1 = UNIVARIATE_MU1_predictions, 
                                               SIGMA1 = UNIVARIATE_SIGMA1_predictions,
                                               MU2 = UNIVARIATE_MU2_predictions, 
                                               SIGMA2 = UNIVARIATE_SIGMA2_predictions,
                                               NU2 = UNIVARIATE_NU2_predictions),
                 ###
                 CoefficientsCOPULA = coef.bivDiscCopula, 
                 oobag.riskCOPULA = oobag.risk.cop, 
                 energy_scoreCOPULA = energy_score,  
                 mstopCOPULA = mstop.bivCopula,
                 TrueKendallRange = TrueKendallTauRange
  )
  
  
  return(output)
  
}


#### FLEXIBLE P 
sim_JOE_FLEXP <- function(seed, n.train, n.mstop, p, n.test, corr, boost.nu.steplength = 0.1){
  
  
  loss <- function(mu1, mu2, sigma1, sigma2, rho, y, nu2 = NULL){
    
    
    # Margin 1 is ZALG:
    F1 <- pdffz(pZALG(y[,1], mu = mu1, sigma = sigma1))
    pdf1 <- pdffz(dZALG(y[,1], mu = mu1, sigma = sigma1))
    
    # Margin 2 is NBI:
    F2 <- pdffz(pZINBI(y[,2], mu = mu2, sigma = sigma2, nu = nu2))
    pdf2 <- pdffz(dZINBI(y[,2], mu = mu2, sigma = sigma2, nu = nu2))
    
    # Copula parameter 
    thet <- rho
    
    # Minus 1
    CDF1m1 <- pdffz(F1 - pdf1)
    CDF2m1 <- pdffz(F2 - pdf2)
    
    ### Copula Terms
    bit1 <- (1 - F1)^thet
    bit1m1 <- (1 - CDF1m1)^thet
    
    bit2 <- (1 - F2)^thet
    bit2m1 <- (1 - CDF2m1)^thet
    
    
    T1 <- pdffz(  1 - (bit1 + bit2 - bit1*bit2)^(1/thet)  )
    
    T2 <- pdffz(  1 - (bit1m1 + bit2 - bit1m1*bit2)^(1/thet)  ) 
    
    T3 <- pdffz(  1 - (bit1 + bit2m1 - bit1*bit2m1)^(1/thet)   )
    
    T4 <- pdffz(  1 - (bit1m1 + bit2m1 - bit1m1*bit2m1)^(1/thet)  ) 
    
    L_Total <- pdffz( T1 - T2 - T3 + T4 )
    
    # Negative log-likelihood
    return(- log( L_Total ) )
    
  }
  
  data.gen.bivdisc <- function(mu1, mu2, sigma1, sigma2, theta, nu2 = NULL){
    
    y1y2 <- matrix(0, ncol = 2, nrow = length(mu1))
    
    for(i in 1:length(mu1)){
      
      paramlist1 <- list(mu = mu1[i], 
                         sigma = sigma1[i])
      
      
      paramlist2 <- list(mu = mu2[i], 
                         sigma = sigma2[i],
                         nu = nu2[i])
      
      copObj <- copula::joeCopula(param = theta[i])
      
      copThing <- mvdc(copObj, margins = c("ZALG", "ZINBI"), paramMargins = list(paramlist1, paramlist2))
      
      y1y2[i,] <- rMvdc(copThing, n = 1)
    }
    
    
    dat <- data.frame(y1y2[,1], y1y2[,2])
    
    return(dat)
  }
  
  data.gen.bivdisc_energyscore <- function(n, mu1, mu2, sigma1, sigma2, theta, nu = NULL){
    
    
    paramlist1 <- list(mu = mu1, sigma = sigma1)
    
    paramlist2 <- list(mu = mu2, sigma = sigma2, nu = nu) 
    
    copObj <- copula::joeCopula(param = theta)
    
    copThing <- mvdc(copObj, margins = c("ZALG", "ZINBI"), paramMargins = list(paramlist1, paramlist2))
    
    y1y2 <- rMvdc(copThing, n = n)
    
    dat <- data.frame(y1y2[,1], y1y2[,2])
    
    
    return(dat)
  }
  
  
  
  
  ## Coefficients:
  # x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
  mu1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  mu2 = c(0, 0, 0, 0, 0,  0, 0, 0, 0, 0) 
  
  # x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
  sigma1 <- c(0, 0, 0, 0, 0,  0, 0, 0, 0, 0) 
  sigma2 <- c(0, 0, 0, 0, 0,  0, 0, 0, 0, 0) 
  
  nu2 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0) 
  
  # x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
  or  = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  
  mu1 <- c( 0, 0, 0, 0, 0, 0, rep(0, p-6)) 
  sigma1 <- c( 0, 0, 0, 0, 0, 0, rep(0, p-6)) 
  mu2 <- c( 0, 0, 0, 0, 0, 0, rep(0, p-6)) 
  sigma2 <- c( 0, 0, 0, 0, 0, 0, rep(0, p-6)) 
  nu2 <- c( 0, 0, 0, 0, 0, 0, rep(0, p-6)) 
  or <- c( 0, 0, 0, 0, 0, 0, rep(0, p-6)) 
  
  
  ### function for mu1
  nl_function_mu1 <- function(x){
    
    
    nl <- sqrt(x)*x^2 - 2*cos(3*x)
    
    return(nl*0.5)
    
  }
  
  ### function for sigma1
  nl_function_sigma1 <- function(x){
    
    func <- -40*(x^(3/2) - x^(4/3))
    
    func <- 2 * func
    
    return(func)
  }
  
  nl_func_sigma1 <- function(x){
    
    func <- -(log(x)) + cos(x)
    
    func <- 2 * func
    
    return(func)
  }
  
  ## function for mu2
  nl_function_mu2 <- function(x){
    
    nl <- -0.7*exp(x^2) + exp(x^0.4) - 0*log(x+0.2) 
    
    return(nl)
  }
  
  nl_function_sigma2 <- function(x){
    
    return( -1.5*(1.5*cos(2 * x) - 0*log(x + 0.15) + 3*tanh(1*x) ) )
  }
  
  nl_function_nu2 <- function(x){
    
    func <- -(sin(x) - exp(x)^2) 
    
    #func <- func + 2*log(x + 0.01) - 3*log(x^3 + 0.01)
    
    func <- 0.7*func
    
    return(func)
  }
  
  nl_function_rho <- function(x){
    
    #return(2*sin(x*4))
    return(2*sin(x*4))
  }
  
  TrueBeta <-  vector('list')
  TrueBeta$mu1 <- mu1
  TrueBeta$mu2 <- mu2
  
  TrueBeta$sigma1 <- sigma1
  TrueBeta$sigma2 <- sigma2
  
  TrueBeta$nu2 <- nu2
  
  TrueBeta$or <- or
  
  set.seed(seed)
  
  
  n <- n.train + n.mstop
  weight.mstop <- c(rep(1, times = n.train),rep(0, times = n.mstop)) 
  
  # - training data
  #x.train <- rmvnorm(n = n, mean = rep(0,p), sigma =  toeplitz(sapply(seq(0,p-1), function(x) corr^x)))
  # x1  <- runif(n, 0, 1)
  # x2 <- runif(n, 0, 1)
  # x3 <- runif(n, 0, 1)
  # x4 <- runif(n, 0, 1)
  # x5 <- runif(n, 0, 1)
  # x6 <- runif(n, 0, 1)
  # x7 <- runif(n, 0, 1)
  # x8 <- runif(n, 0, 1)
  # x9 <- runif(n, 0, 1) 
  # x10 <- runif(n, 0, 1)
  # 
  # x.train <- cbind(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
  
  BigXMatrix_Train <- matrix(runif(n * p, 0, 1), nrow = n, ncol = p, byrow = T)
  x.train <- BigXMatrix_Train
  
  ############################################################################################################################################
  ############################################################################################################################################
  # - test data
  #x.test <- rmvnorm(n = n.test, mean = rep(0,p), sigma = toeplitz(sapply(seq(0,p-1), function(x) corr^x)))
  # x1.test <- runif(n.test, 0, 1)
  # x2.test <- runif(n.test, 0, 1)
  # x3.test <- runif(n.test, 0, 1)
  # x4.test <- runif(n.test, 0, 1)
  # x5.test <- runif(n.test, 0, 1)
  # x6.test <- runif(n.test, 0, 1)
  # x7.test <- runif(n.test, 0, 1)
  # x8.test <- runif(n.test, 0, 1)
  # x9.test <- runif(n.test, 0, 1) 
  # x10.test <- runif(n.test, 0, 1)
  # 
  # x.test <- cbind(x1.test, x2.test, x3.test, x4.test, x5.test, x6.test, x7.test, x8.test, x9.test, x10.test)
  # 
  # colnames(x.test) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
  BigXMatrix_Test <- matrix(runif(n.test * p, 0, 1), nrow = n.test, ncol = p, byrow = T)
  x.test <- BigXMatrix_Test
  
  ############################################################################################################################################
  # GENERATION OF DISTRIBUTION PARAMETERS AND STUFF 
  # predictor
  
  ##### TRAIN DATA
  train.eta.mu1_center <- (x.train %*% mu1)
  train.eta.mu1_center <- train.eta.mu1_center + nl_function_mu1(x.train[,1])
  
  ##### TEST DATA
  test.eta.mu1_center <- (x.test %*% mu1)
  test.eta.mu1_center <- test.eta.mu1_center + nl_function_mu1(x.test[,1])
  
  # apply response function:
  train.eta.mu1<-  plogis(train.eta.mu1_center) 
  range(train.eta.mu1)
  
  test.eta.mu1<-  plogis(test.eta.mu1_center) 
  ########################################################################################################################
  #### SIGMA FOR ZALG RESPONSE
  train.eta.sigma1 <-  x.train %*% sigma1
  
  train.eta.sigma1 <- train.eta.sigma1 - 2 + nl_function_sigma1(x.train[,3]) 
  train.eta.sigma1 <- plogis(train.eta.sigma1)
  range(train.eta.sigma1)
  
  ##### TEST DATA
  test.eta.sigma1 <-  x.test %*% sigma1
  test.eta.sigma1 <- test.eta.sigma1 - 2 + nl_function_sigma1(x.test[,3])
  test.eta.sigma1 <- plogis(test.eta.sigma1)
  ########################################################################################################################
  # predictor and apply response function: 
  train.eta.mu2 <-  x.train %*% mu2
  train.eta.mu2 <-  train.eta.mu2 + nl_function_mu2(x.train[,2])
  train.eta.mu2 <-  exp(train.eta.mu2) 
  range(train.eta.mu2)
  
  
  test.eta.mu2 <-  x.test %*% mu2
  test.eta.mu2 <-  test.eta.mu2 + nl_function_mu2(x.test[,2])
  test.eta.mu2 <-  exp(test.eta.mu2) 
  ########################################################################################################################
  # predictor and apply response function: 
  train.eta.sigma2 <-  x.train %*% sigma2
  train.eta.sigma2 <-  train.eta.sigma2 + nl_function_sigma2(x.train[,5]) + 3
  train.eta.sigma2 <-  exp(train.eta.sigma2) 
  range(train.eta.sigma2)
  
  test.eta.sigma2 <-  x.test %*% sigma2
  test.eta.sigma2 <-  test.eta.sigma2 + nl_function_sigma2(x.test[,5]) + 3
  test.eta.sigma2 <-  exp(test.eta.sigma2) 
  ########################################################################################################################
  # predictor and apply response function: 
  train.eta.nu2 <-  x.train %*% nu2
  train.eta.nu2 <-  train.eta.nu2 + nl_function_nu2(x.train[,1]) - 3
  train.eta.nu2 <-  plogis(train.eta.nu2) 
  range(train.eta.nu2)
  
  test.eta.nu2 <-  x.test %*% nu2
  test.eta.nu2 <-  test.eta.nu2 + nl_function_nu2(x.test[,1]) - 3
  test.eta.nu2 <-  plogis(test.eta.nu2) 
  
  ########################################################################################################################
  ### Copula parameter
  train.eta.or_center <- (x.train %*% or)
  train.eta.or_center <- train.eta.or_center + nl_function_rho(x.train[,4]) 
  
  train.copula.parameter <- (exp(train.eta.or_center) + 1)
  range(train.copula.parameter)
  
  TrueKendallTauRange <- VineCopula::BiCopPar2Tau(6, par = range(train.copula.parameter))  
  
  test.eta.or_center <- (x.test %*% or)
  test.eta.or_center <- test.eta.or_center + nl_function_rho(x.test[,4])
  
  test.copula.parameter <- (exp(test.eta.or_center) + 1)
  ############################################################################################################################################
  # SAMPLING OF BIVARIATE RESPONSE FROM COPULA
  y.train <- data.gen.bivdisc(mu1 = train.eta.mu1, mu2 = train.eta.mu2, sigma1 = train.eta.sigma1, sigma2 = train.eta.sigma2, theta = train.copula.parameter, nu2 = train.eta.nu2)
  
  colnames(y.train)<- c('y1','y2')
  dat.train <- data.frame(y.train,x.train)
  
  
  y.test <- data.gen.bivdisc(mu1 = test.eta.mu1, mu2 = test.eta.mu2, sigma1 = test.eta.sigma1, sigma2 = test.eta.sigma2, theta = test.copula.parameter, nu2 = test.eta.nu2)
  
  colnames(y.test)<- c('y1','y2')
  dat.test <- data.frame(y.test,x.test)
  
  ############################################################################################################################################
  # FITTIIIINGG
  mstop.bivBern <-  vector('list')
  
  AUCbivBern <-  vector('list')
  BrierbivBern <-  vector('list')
  lik <- vector('list')
  MSEbivDisc <- vector("list")
  
  
  
  # # Predict distributional quantities
  
  # model_equation <- formula(cbind(y1,y2) ~ bbs(x1, knots = 20, degree = 3, difference = 2) + 
  #                             bbs(x2, knots = 20, degree = 3, difference = 2) + 
  #                             bbs(x3, knots = 20, degree = 3, difference = 2) + 
  #                             bbs(x4, knots = 20, degree = 3, difference = 2) + 
  #                             bbs(x5, knots = 20, degree = 3, difference = 2) + 
  #                             bbs(x6, knots = 20, degree = 3, difference = 2) +
  #                             bbs(x7, knots = 20, degree = 3, difference = 2) + 
  #                             bbs(x8, knots = 20, degree = 3, difference = 2) + 
  #                             bbs(x9, knots = 20, degree = 3, difference = 2) + 
  #                             bbs(x10, knots = 20, degree = 3, difference = 2))
  # 
  # model_equation <- formula(cbind(y1,y2) ~ bbs(x1) + 
  #                             bbs(x2) + 
  #                             bbs(x3) + 
  #                             bbs(x4) + 
  #                             bbs(x5) + 
  #                             bbs(x6) +
  #                             bbs(x7) + 
  #                             bbs(x8) + 
  #                             bbs(x9) + 
  #                             bbs(x10))
  
  model_equation <- formula(cbind(y1,y2) ~ .)
  
  bivModel_Formula <- list(mu1 = model_equation,
                           mu2 = model_equation,
                           sigma1 = model_equation,
                           sigma2 = model_equation,
                           rho = model_equation
  )
  
  
  bivModel_Formula <- list(mu1 = model_equation,
                           mu2 = model_equation,
                           sigma1 = model_equation,
                           sigma2 = model_equation,
                           nu2 = model_equation,
                           rho = model_equation
  )
  
  #############################################
  ################# - COPULA - ################
  #############################################
  # Bivariate copula model
  bivDiscCopula <- gamboostLSS(bivModel_Formula, data = dat.train, 
                               families = Joe_Cop_BivDiscrete(marg1 = "ZALG", 
                                                              marg2 = "ZINBI",
                                                              stabilization = "L2"), 
                               control = boost_control(mstop = 2000, risk = 'oobag',
                                                       nu  = boost.nu.steplength, trace = TRUE), 
                               method = 'noncyclic', weights = weight.mstop)
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 1900){
    mstop(bivDiscCopula) <- 3000
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 2900){
    mstop(bivDiscCopula) <- 4000
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 3990){
    bivDiscCopula[8000]
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 7990){
    bivDiscCopula[10000]
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 9990){
    bivDiscCopula[15000]
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 14990){
    bivDiscCopula[20000]
  }
  
  
  
  
  MSTOP_COP <- which.min(risk(bivDiscCopula, merge = T))
  oobag.risk.cop <- risk(bivDiscCopula,merge = T)
  
  # rm(bivDiscCopula)
  # dat.train_biv <- dat.train[weight.mstop == 1, ]
  # 
  # # RE-FIT until OPTIMAL MSTOP
  # bivDiscCopula <- gamboostLSS(bivModel_Formula, data = dat.train_biv, 
  #             families = FGM_Cop_BivDiscrete(marg1 = "ZALG", 
  #                                            marg2 = "ZINBI"),
  #             control = boost_control(mstop = MSTOP_COP, nu  = boost.nu.steplength, trace = TRUE), 
  #             method = 'noncyclic')
  bivDiscCopula <- bivDiscCopula[MSTOP_COP]
  
  
  mstop.bivCopula <-  vector('list')
  mstop.bivCopula$mstop <- MSTOP_COP
  mstop.bivCopula$mu1 <- bivDiscCopula$mu1$mstop()
  mstop.bivCopula$sigma1 <- bivDiscCopula$sigma1$mstop()
  
  mstop.bivCopula$mu2 <- bivDiscCopula$mu2$mstop()
  mstop.bivCopula$sigma2 <- bivDiscCopula$sigma2$mstop()
  mstop.bivCopula$nu2 <- bivDiscCopula$nu2$mstop()
  
  mstop.bivCopula$rho  <- bivDiscCopula$rho$mstop()
  
  # EXTRACT ALL COEFFICIENTS!
  coef.bivDiscCopula <- coef(bivDiscCopula, which = "")
  
  DiscreteMargins_Metrics <- vector("list")
  likBinCopula <- vector('list')
  
  
  ################ Get predictions from each covariate?
  #newX <- matrix(rep(seq(from = 0, to = 1, length.out = 500),times = 10), nrow = 500)
  #colnames(newX) <-  c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
  newX <- matrix(rep(seq(from = 0, to = 1, length.out = 150), times = p), nrow = 150)
  
  
  colnames(newX) <- colnames(x.train)
  
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ### Copula predictions / selection
  Copula_MU1_predictions <- list()
  Copula_SIGMA1_predictions <- list()
  
  Copula_MU2_predictions <- list()
  Copula_SIGMA2_predictions <- list()
  Copula_NU2_predictions <- list()
  
  Copula_RHO_predictions <- list()
  
  
  
  Copula_MU1_predictions <- sapply(1:ncol(x.train), function(i) predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = i, type = "link"))
  Copula_SIGMA1_predictions <- sapply(1:ncol(x.train), function(i) predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = i, type = "link"))
  
  
  Copula_MU2_predictions <- sapply(1:ncol(x.train), function(i) predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = i, type = "link"))
  Copula_SIGMA2_predictions <- sapply(1:ncol(x.train), function(i) predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = i, type = "link"))
  Copula_NU2_predictions <- sapply(1:ncol(x.train), function(i) predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = i, type = "link"))
  
  Copula_RHO_predictions <- sapply(1:ncol(x.train), function(i) predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = i, type = "link"))
  
  
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  
  
  
  # Predict distributional quantities
  predCopula.mu1 <- predict(bivDiscCopula$mu1, newdata = dat.test, type = 'response')
  predCopula.sigma1 <- predict(bivDiscCopula$sigma1, newdata = dat.test, type = 'response')
  
  predCopula.mu2 <- predict(bivDiscCopula$mu2, newdata = dat.test, type = 'response')
  predCopula.sigma2 <- predict(bivDiscCopula$sigma2, newdata = dat.test, type = 'response')
  predCopula.nu2 <- predict(bivDiscCopula$nu2, newdata = dat.test, type = 'response')
  
  
  predCopula.rho <- predict(bivDiscCopula$rho, newdata = dat.test, type = 'response')
  
  ######### COMPUTE PERFORMANCE: 
  
  # MSE (margin 1, then margin 2)
  DiscreteMargins_Metrics$mu1 <- mean((as.numeric(predCopula.mu1)-(as.numeric(dat.test$y1)))^2)
  DiscreteMargins_Metrics$mu2 <- mean((as.numeric(predCopula.mu2)-(as.numeric(dat.test$y2)))^2)
  
  # LOG-SCORE (LOG-LIKELIHOOD) (change to the copula likelihood)
  lik$biv <- sum(loss(mu1 = predCopula.mu1, 
                      mu2 = predCopula.mu2, 
                      sigma1 = predCopula.sigma1, 
                      sigma2 = predCopula.sigma2, 
                      nu2 = predCopula.nu2,
                      rho = predCopula.rho, 
                      y = y.test))
  
  
  ###########################################################################################################################
  ################ - univariat - ############################################################################################
  ###########################################################################################################################
  
  mstop.uni <- vector('list')
  coef.uni  <- vector('list')
  coef.uni_1  <- vector('list')
  
  # univariateDiscreteEquation_1 <- list(mu = formula(y1 ~  bbs(x1) + 
  #                                                     bbs(x2) + 
  #                                                     bbs(x3) + 
  #                                                     bbs(x4) + 
  #                                                     bbs(x5) + 
  #                                                     bbs(x6) +
  #                                                     bbs(x7) + 
  #                                                     bbs(x8) + 
  #                                                     bbs(x9) + 
  #                                                     bbs(x10)),
  #                                      sigma = formula(y1 ~  bbs(x1) + 
  #                                                        bbs(x2) + 
  #                                                        bbs(x3) + 
  #                                                        bbs(x4) + 
  #                                                        bbs(x5) + 
  #                                                        bbs(x6) +
  #                                                        bbs(x7) + 
  #                                                        bbs(x8) + 
  #                                                        bbs(x9) + 
  #                                                        bbs(x10))
  # )
  # 
  # 
  # univariateDiscreteEquation_2 <- list(mu = formula(y2 ~  bbs(x1) + 
  #                                                     bbs(x2) + 
  #                                                     bbs(x3) + 
  #                                                     bbs(x4) + 
  #                                                     bbs(x5) + 
  #                                                     bbs(x6) +
  #                                                     bbs(x7) + 
  #                                                     bbs(x8) + 
  #                                                     bbs(x9) + 
  #                                                     bbs(x10)),
  #                                      sigma = formula(y2 ~  bbs(x1) + 
  #                                                        bbs(x2) + 
  #                                                        bbs(x3) + 
  #                                                        bbs(x4) + 
  #                                                        bbs(x5) + 
  #                                                        bbs(x6) +
  #                                                        bbs(x7) + 
  #                                                        bbs(x8) + 
  #                                                        bbs(x9) + 
  #                                                        bbs(x10),
  #                                      ),
  #                                      nu = formula(y2 ~  bbs(x1) + 
  #                                                     bbs(x2) + 
  #                                                     bbs(x3) + 
  #                                                     bbs(x4) + 
  #                                                     bbs(x5) + 
  #                                                     bbs(x6) +
  #                                                     bbs(x7) + 
  #                                                     bbs(x8) + 
  #                                                     bbs(x9) + 
  #                                                     bbs(x10))
  # )
  
  
  univariateDiscreteEquation_1 <- formula(y1 ~ .)
  univariateDiscreteEquation_2 <- formula(y2 ~ .)
  
  
  # - margin 1
  dat.train.mu1 = dat.train[, !(names(dat.train) %in% "y2")]
  dat.test.mu1  = dat.test[ , !(names(dat.test) %in% "y2")]
  
  glm.uni.mu1 <- gamboostLSS(univariateDiscreteEquation_1, 
                             data = dat.train.mu1, 
                             families = as.families(fname = "ZALG", stabilization = "L2"), 
                             control = boost_control(risk = 'oobag', nu = boost.nu.steplength, mstop = 2000, trace = T), 
                             method = "noncyclic",
                             weights = weight.mstop)
  
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  if(mstop.uni$mu1 >= 1990){
    glm.uni.mu1[5000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  if(mstop.uni$mu1 >= 4990){
    glm.uni.mu1[10000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  if(mstop.uni$mu1 >= 9990){
    glm.uni.mu1[20000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  # dat.train_mu1 <- dat.train.mu1[weight.mstop == 1, ]
  # 
  # glm.uni.mu1 <- gamboostLSS(univariateDiscreteEquation_1, 
  #                            data = dat.train_mu1, 
  #                            families = as.families(fname = "ZALG"), 
  #                         control = boost_control(mstop = mstop.uni$mu1, nu  = boost.nu.steplength, trace = T), 
  #                         method = "noncyclic")
  glm.uni.mu1 <- glm.uni.mu1[mstop.uni$mu1]
  
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 1  
  coef.uni$mu1 <- coef(glm.uni.mu1, which = '')
  
  DiscreteMargins_Metrics$Univariate_mu1 <- mean((as.numeric(predict(glm.uni.mu1$mu, newdata = dat.test.mu1, type = "response")) - (as.numeric(dat.test.mu1$y1)))^2)
  
  # - margin 2
  dat.train.mu2 =  dat.train[, !(names(dat.train) %in% "y1")]
  dat.test.mu2 = dat.test[, !(names(dat.test) %in% "y1")]
  
  glm.uni.mu2 <- gamboostLSS(univariateDiscreteEquation_2, 
                             data = dat.train.mu2, 
                             families = as.families(fname = "ZINBI", stabilization = "L2"), 
                             control = boost_control(risk = 'oobag', nu = boost.nu.steplength, trace = TRUE, mstop = 2000), 
                             method = "noncyclic",
                             weights = weight.mstop)
  
  
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  if(mstop.uni$mu2 >= 1990){
    glm.uni.mu2[5000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  if(mstop.uni$mu2 >= 4990){
    glm.uni.mu2[10000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  if(mstop.uni$mu2 >= 9990){
    glm.uni.mu2[20000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  # dat.train_mu2 <- dat.train.mu2[weight.mstop == 1, ]
  # 
  # glm.uni.mu2 <- gamboostLSS(univariateDiscreteEquation_2, 
  #                            data = dat.train_mu2, 
  #                            families = as.families(fname = "ZINBI"), 
  #                         control = boost_control(mstop = mstop.uni$mu2, nu  = boost.nu.steplength, trace=TRUE), method = "noncyclic")
  glm.uni.mu2 <- glm.uni.mu2[mstop.uni$mu2]
  
  mstop.uni$Margin1 <- mstop(glm.uni.mu1)
  
  mstop.uni$Margin2 <- mstop(glm.uni.mu2)
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 2  
  coef.uni$mu2 <- coef(glm.uni.mu2,which = '')
  
  DiscreteMargins_Metrics$Univariate_mu2 <- mean((as.numeric(predict(glm.uni.mu2$mu, newdata = dat.test.mu2, type = "response"))-(as.numeric(dat.test.mu2$y2)))^2)
  
  # Likelihood
  pred.mu1.uni <- as.numeric(predict(glm.uni.mu1$mu, newdata = dat.test.mu1, type = "response"))
  pred.sigma1.uni <- as.numeric(predict(glm.uni.mu1$sigma, newdata = dat.test.mu1, type = "response"))
  
  pred.mu2.uni <- as.numeric(predict(glm.uni.mu2$mu, newdata = dat.test.mu2, type = "response"))
  pred.sigma2.uni <- as.numeric(predict(glm.uni.mu2$sigma, newdata = dat.test.mu2, type = "response"))
  pred.nu2.uni <- as.numeric(predict(glm.uni.mu2$nu, newdata = dat.test.mu2, type = "response"))
  
  
  mu1.uni.loglik <- sum(-dZALG(x = dat.test.mu1$y1, mu = pred.mu1.uni, sigma = pred.sigma1.uni, log = T))
  mu2.uni.loglik <- sum(-dZINBI(x = dat.test.mu2$y2, mu = pred.mu2.uni, sigma = pred.sigma2.uni, nu = pred.nu2.uni, log = T))
  
  lik$uni<- mu1.uni.loglik + mu2.uni.loglik
  
  
  n.ges = c(n.train, n.mstop, n.test)
  pred.ges <- NULL#list(pred.mu1, pred.mu2, pred.or)
  pred.uni <- list(mu1 = pred.mu1.uni, 
                   sigma1 = pred.sigma1.uni,
                   mu2 = pred.mu2.uni,
                   sigma2 = pred.sigma2.uni,
                   nu2 = pred.nu2.uni
  )
  
  pred.cop <- list(mu1 = predCopula.mu1, 
                   sigma1 = predCopula.sigma1,
                   mu2 = predCopula.mu2,
                   sigma2 = predCopula.sigma2,
                   nu2 = predCopula.nu2,
                   rho = predCopula.rho)
  
  
  ##### Predictions:
  UNIVARIATE_MU1_predictions <- vector("list")
  UNIVARIATE_SIGMA1_predictions <- vector("list")
  
  UNIVARIATE_MU2_predictions <- vector("list")
  UNIVARIATE_SIGMA2_predictions <- vector("list")
  UNIVARIATE_NU2_predictions <- vector("list")
  
  
  UNIVARIATE_MU1_predictions <- sapply(1:p, function(i) predict(glm.uni.mu1$mu, data.frame(newX), which  = i, type = "link"))
  UNIVARIATE_SIGMA1_predictions <- sapply(1:p, function(i) predict(glm.uni.mu1$sigma, data.frame(newX), which  = i, type = "link"))
  
  UNIVARIATE_MU2_predictions <- sapply(1:p, function(i) predict(glm.uni.mu2$mu, data.frame(newX), which  = i, type = "link"))
  UNIVARIATE_SIGMA2_predictions <- sapply(1:p, function(i) predict(glm.uni.mu2$sigma, data.frame(newX), which  = i, type = "link"))
  UNIVARIATE_NU2_predictions <- sapply(1:p, function(i) predict(glm.uni.mu2$nu, data.frame(newX), which  = i, type = "link"))
  
  
  
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
    sample_uni <- data.gen.bivdisc_energyscore(n = 1000, 
                                               mu1 = pred.mu1.uni[i], 
                                               mu2 = pred.mu2.uni[i], 
                                               sigma1 = pred.sigma1.uni[i], 
                                               sigma2 = pred.sigma2.uni[i],
                                               nu = check_extremes_nuZINBI(pred.nu2.uni[i]),
                                               theta = 1)
    
    pred_sample_uni[1,] <- sample_uni[,1]
    pred_sample_uni[2,] <- sample_uni[,2]
    
    es_uni[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_uni) 
    
    
    # copula
    sample_cop <- data.gen.bivdisc_energyscore(n = 1000, 
                                               mu1 = predCopula.mu1[i],  
                                               mu2 = predCopula.mu2[i], 
                                               sigma1 = predCopula.sigma1[i], 
                                               sigma2 = predCopula.sigma2[i],
                                               nu = check_extremes_nuZINBI(predCopula.nu2[i]),
                                               theta = predCopula.rho[i])
    
    pred_sample_cop[1, ] <- sample_cop[,1]
    pred_sample_cop[2, ] <- sample_cop[,2]
    
    es_cop[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_cop) 
    
    
  }
  
  energy_score <- list()
  energy_score$uni <- mean(es_uni)
  energy_score$cop <- mean(es_cop)
  
  
  
  ########################################### We gather: mstop, predictions, MSE margin 1, margin 2, loglikelihood, energy score.
  output <- list(TrueBeta = TrueBeta, 
                 n = n.ges, 
                 p  = p, 
                 Likelihood = lik,
                 predict.uni = pred.uni, 
                 predict.cop = pred.cop,
                 energy_score = energy_score,  
                 ####
                 Coefficients.uni = coef.uni, 
                 mstop.uni = mstop.uni,
                 ###
                 DiscreteMarginMetrics = DiscreteMargins_Metrics,
                 ###
                 Copula_Predictions = list(MU1 = Copula_MU1_predictions, 
                                           SIGMA1 = Copula_SIGMA1_predictions,
                                           MU2 = Copula_MU2_predictions,
                                           SIGMA2 = Copula_SIGMA2_predictions,
                                           NU2 = Copula_NU2_predictions,
                                           RHO = Copula_RHO_predictions),
                 ###
                 Univariate_Predictions = list(MU1 = UNIVARIATE_MU1_predictions, 
                                               SIGMA1 = UNIVARIATE_SIGMA1_predictions,
                                               MU2 = UNIVARIATE_MU2_predictions, 
                                               SIGMA2 = UNIVARIATE_SIGMA2_predictions,
                                               NU2 = UNIVARIATE_NU2_predictions),
                 ###
                 CoefficientsCOPULA = coef.bivDiscCopula, 
                 oobag.riskCOPULA = oobag.risk.cop, 
                 energy_scoreCOPULA = energy_score,  
                 mstopCOPULA = mstop.bivCopula,
                 TrueKendallRange = TrueKendallTauRange
  )
  
  
  return(output)
  
}


#### USING ROTATED JOE COPULA
sim_JOE180 <- function(seed, n.train, n.mstop, p, n.test, corr, boost.nu.steplength = 0.1){
  
  
  loss <- function(mu1, mu2, sigma1, sigma2, rho, y, nu2 = NULL){
    
    
    # Margin 1 is ZALG:
    F1 <- pdffz(pZALG(y[,1], mu = mu1, sigma = sigma1))
    pdf1 <- pdffz(dZALG(y[,1], mu = mu1, sigma = sigma1))
    
    # Margin 2 is NBI:
    F2 <- pdffz(pZINBI(y[,2], mu = mu2, sigma = sigma2, nu = nu2))
    pdf2 <- pdffz(dZINBI(y[,2], mu = mu2, sigma = sigma2, nu = nu2))
    
    # Copula parameter 
    thet <- rho
    
    # Minus 1
    CDF1m1 <- pdffz(F1 - pdf1)
    CDF2m1 <- pdffz(F2 - pdf2)
    
    ############### Rotation terms:
    F1_Tilde <- pdffz( 1 - F1 )
    F2_Tilde <- pdffz( 1 - F2 )
    
    CDF1m1_Tilde <- pdffz( 1 - CDF1m1 )
    CDF2m1_Tilde <- pdffz( 1 - CDF2m1 )
    
        
    ### Copula Terms
    bit1 <- (1 - F1)^thet
    bit1m1 <- (1 - CDF1m1)^thet
    
    bit2 <- (1 - F2)^thet
    bit2m1 <- (1 - CDF2m1)^thet
    
    
    T1 <- pdffz( F1_Tilde + F2_Tilde - 1 + (  1 - (bit1 + bit2 - bit1*bit2)^(1/thet) ) )
    
    T2 <- pdffz( CDF1m1_Tilde + F2_Tilde - 1 + (  1 - (bit1m1 + bit2 - bit1m1*bit2)^(1/thet) ) ) 
    
    T3 <- pdffz( F1_Tilde + CDF2m1_Tilde - 1 + ( 1 - (bit1 + bit2m1 - bit1*bit2m1)^(1/thet) )  )
    
    T4 <- pdffz( CDF1m1_Tilde + CDF2m1_Tilde - 1 + ( 1 - (bit1m1 + bit2m1 - bit1m1*bit2m1)^(1/thet) ) ) 
    
    
    L_Total <- pdffz( T1 - T2 - T3 + T4 )
    
    # Negative log-likelihood
    return(- log( L_Total ) )
    
  }
  
  data.gen.bivdisc <- function(mu1, mu2, sigma1, sigma2, theta, nu2 = NULL){
    
    y1y2 <- matrix(0, ncol = 2, nrow = length(mu1))
    
    for(i in 1:length(mu1)){
      
      paramlist1 <- list(mu = mu1[i], 
                         sigma = sigma1[i])
      
      
      paramlist2 <- list(mu = mu2[i], 
                         sigma = sigma2[i],
                         nu = nu2[i])
      
      CopulaSamples <- VineCopula::BiCopSim(1, family = 16, par = theta[i])
      
      
      
      y1y2[i, 1] <- qZALG(p = CopulaSamples[,1], mu = mu1[i], sigma = sigma1[i])
      
      y1y2[i, 2] <- qZINBI(p = CopulaSamples[,2], mu = mu2[i], sigma = sigma2[i], nu = nu2[i])
      
      #copObj <- copula::joeCopula(param = theta[i])
      #copThing <- mvdc(copObj, margins = c("ZALG", "ZINBI"), paramMargins = list(paramlist1, paramlist2))
      #y1y2[i,] <- rMvdc(copThing, n = 1)
    }
    
    
    dat <- data.frame(y1y2[,1], y1y2[,2])
    
    return(dat)
  }
  
  data.gen.bivdisc_energyscore <- function(n, mu1, mu2, sigma1, sigma2, theta, nu = NULL){
    
    y1y2 <- matrix(0, nrow = n, ncol = 2)
    
    paramlist1 <- list(mu = mu1, sigma = sigma1)
    
    paramlist2 <- list(mu = mu2, sigma = sigma2, nu = nu) 
    
    # copObj <- copula::joeCopula(param = theta)
    # 
    # copThing <- mvdc(copObj, margins = c("ZALG", "ZINBI"), paramMargins = list(paramlist1, paramlist2))
    # 
    # y1y2 <- rMvdc(copThing, n = n)
    CopulaSamples <- VineCopula::BiCopSim(n, family = 16, par = theta[i])
    
    
    
    y1y2[, 1] <- qZALG(p = CopulaSamples[,1], mu = mu1, sigma = sigma1)
    
    y1y2[, 2] <- qZINBI(p = CopulaSamples[,2], mu = mu2, sigma = sigma2, nu = nu)
    
    dat <- data.frame(y1y2[,1], y1y2[,2])
    
    
    return(dat)
  }
  
  
  
  
  ## Coefficients:
  # x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
  mu1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  mu2 = c(0, 0, 0, 0, 0,  0, 0, 0, 0, 0) 
  
  # x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
  sigma1 <- c(0, 0, 0, 0, 0,  0, 0, 0, 0, 0) 
  sigma2 <- c(0, 0, 0, 0, 0,  0, 0, 0, 0, 0) 
  
  nu2 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0) 
  
  # x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
  or  = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  
  ### function for mu1
  nl_function_mu1 <- function(x){
    
    
    nl <- sqrt(x)*x^2 - 2*cos(3*x)
    
    return(nl*0.5)
    
  }
  
  ### function for sigma1
  nl_function_sigma1 <- function(x){
    
    func <- -40*(x^(3/2) - x^(4/3))
    
    func <- 2 * func
    
    return(func)
  }
  
  nl_func_sigma1 <- function(x){
    
    func <- -(log(x)) + cos(x)
    
    func <- 2 * func
    
    return(func)
  }
  
  ## function for mu2
  nl_function_mu2 <- function(x){
    
    nl <- -0.7*exp(x^2) + exp(x^0.4) - 0*log(x+0.2) 
    
    return(nl)
  }
  
  nl_function_sigma2 <- function(x){
    
    return( -1.5*(1.5*cos(2 * x) + 0.5*log(x)) )
  }
  
  nl_function_nu2 <- function(x){
    
    func <- -(sin(x) - exp(x)^2)
    
    func <- 0.7*func
    
    return(func)
  }
  
  nl_function_rho <- function(x){
    
    #return(2*sin(x*4))
    return(2*sin(x*4))
  }
  
  TrueBeta <-  vector('list')
  TrueBeta$mu1 <- mu1
  TrueBeta$mu2 <- mu2
  
  TrueBeta$sigma1 <- sigma1
  TrueBeta$sigma2 <- sigma2
  
  TrueBeta$nu2 <- nu2
  
  TrueBeta$or <- or
  
  set.seed(seed)
  
  
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
  
  ############################################################################################################################################
  ############################################################################################################################################
  # - test data
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
  
  ############################################################################################################################################
  # GENERATION OF DISTRIBUTION PARAMETERS AND STUFF 
  # predictor
  
  ##### TRAIN DATA
  train.eta.mu1_center <- (x.train %*% mu1)
  train.eta.mu1_center <- train.eta.mu1_center + nl_function_mu1(x1)
  
  ##### TEST DATA
  test.eta.mu1_center <- (x.test %*% mu1)
  test.eta.mu1_center <- test.eta.mu1_center + nl_function_mu1(x1.test)
  
  # apply response function:
  train.eta.mu1<-  plogis(train.eta.mu1_center) 
  range(train.eta.mu1)
  
  test.eta.mu1<-  plogis(test.eta.mu1_center) 
  ########################################################################################################################
  #### SIGMA FOR ZALG RESPONSE
  train.eta.sigma1 <-  x.train %*% sigma1
  
  train.eta.sigma1 <- train.eta.sigma1 - 2 + nl_function_sigma1(x3) 
  train.eta.sigma1 <- plogis(train.eta.sigma1)
  range(train.eta.sigma1)
  
  ##### TEST DATA
  test.eta.sigma1 <-  x.test %*% sigma1
  test.eta.sigma1 <- test.eta.sigma1 - 2 + nl_function_sigma1(x3.test)
  test.eta.sigma1 <- plogis(test.eta.sigma1)
  ########################################################################################################################
  # predictor and apply response function: 
  train.eta.mu2 <-  x.train %*% mu2
  train.eta.mu2 <-  train.eta.mu2 + nl_function_mu2(x2)
  train.eta.mu2 <-  exp(train.eta.mu2) 
  range(train.eta.mu2)
  
  
  test.eta.mu2 <-  x.test %*% mu2
  test.eta.mu2 <-  test.eta.mu2 + nl_function_mu2(x2.test)
  test.eta.mu2 <-  exp(test.eta.mu2) 
  ########################################################################################################################
  # predictor and apply response function: 
  train.eta.sigma2 <-  x.train %*% sigma2
  train.eta.sigma2 <-  train.eta.sigma2 + nl_function_sigma2(x4)
  train.eta.sigma2 <-  exp(train.eta.sigma2) 
  range(train.eta.sigma2)
  
  test.eta.sigma2 <-  x.test %*% sigma2
  test.eta.sigma2 <-  test.eta.sigma2 +  nl_function_sigma2(x4.test)
  test.eta.sigma2 <-  exp(test.eta.sigma2) 
  ########################################################################################################################
  # predictor and apply response function: 
  train.eta.nu2 <-  x.train %*% nu2
  train.eta.nu2 <-  train.eta.nu2 + nl_function_nu2(x1) - 3
  train.eta.nu2 <-  plogis(train.eta.nu2) 
  range(train.eta.nu2)
  
  test.eta.nu2 <-  x.test %*% nu2
  test.eta.nu2 <-  test.eta.nu2 + nl_function_nu2(x1.test) - 3
  test.eta.nu2 <-  plogis(test.eta.nu2) 
  
  ########################################################################################################################
  ### Copula parameter
  train.eta.or_center <- (x.train %*% or)
  train.eta.or_center <- train.eta.or_center + nl_function_rho(x3) 
  
  train.copula.parameter <- (exp(train.eta.or_center) + 1)
  range(train.copula.parameter)
  
  TrueKendallTauRange <- VineCopula::BiCopPar2Tau(6, par = range(train.copula.parameter))  
  
  test.eta.or_center <- (x.test %*% or)
  test.eta.or_center <- test.eta.or_center + nl_function_rho(x3.test)
  
  test.copula.parameter <- (exp(test.eta.or_center) + 1)
  ############################################################################################################################################
  # SAMPLING OF BIVARIATE RESPONSE FROM COPULA
  y.train <- data.gen.bivdisc(mu1 = train.eta.mu1, mu2 = train.eta.mu2, sigma1 = train.eta.sigma1, sigma2 = train.eta.sigma2, theta = train.copula.parameter, nu2 = train.eta.nu2)
  
  colnames(y.train)<- c('y1','y2')
  dat.train <- data.frame(y.train,x.train)
  
  
  y.test <- data.gen.bivdisc(mu1 = test.eta.mu1, mu2 = test.eta.mu2, sigma1 = test.eta.sigma1, sigma2 = test.eta.sigma2, theta = test.copula.parameter, nu2 = test.eta.nu2)
  
  colnames(y.test)<- c('y1','y2')
  dat.test <- data.frame(y.test,x.test)
  
  ############################################################################################################################################
  # FITTIIIINGG
  mstop.bivBern <-  vector('list')
  
  AUCbivBern <-  vector('list')
  BrierbivBern <-  vector('list')
  lik <- vector('list')
  MSEbivDisc <- vector("list")
  
  
  
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
                           sigma1 = model_equation,
                           sigma2 = model_equation,
                           rho = model_equation
  )
  
  
  bivModel_Formula <- list(mu1 = model_equation,
                           mu2 = model_equation,
                           sigma1 = model_equation,
                           sigma2 = model_equation,
                           nu2 = model_equation,
                           rho = model_equation
  )
  
  #############################################
  ################# - COPULA - ################
  #############################################
  # Bivariate copula model
  bivDiscCopula <- gamboostLSS(bivModel_Formula, data = dat.train, 
                               families = Joe180_Cop_BivDiscrete(marg1 = "ZALG", 
                                                                 marg2 = "ZINBI",
                                                                 stabilization = "L2"), 
                               control = boost_control(mstop = 2000, risk = 'oobag',
                                                       nu  = boost.nu.steplength, trace = TRUE), 
                               method = 'noncyclic', weights = weight.mstop)
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 1900){
    mstop(bivDiscCopula) <- 3000
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 2900){
    mstop(bivDiscCopula) <- 4000
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 3990){
    bivDiscCopula[8000]
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 7990){
    bivDiscCopula[10000]
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 9990){
    bivDiscCopula[15000]
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 14990){
    bivDiscCopula[20000]
  }
  
  
  
  
  MSTOP_COP <- which.min(risk(bivDiscCopula, merge = T))
  oobag.risk.cop <- risk(bivDiscCopula,merge = T)
  
  # rm(bivDiscCopula)
  # dat.train_biv <- dat.train[weight.mstop == 1, ]
  # 
  # # RE-FIT until OPTIMAL MSTOP
  # bivDiscCopula <- gamboostLSS(bivModel_Formula, data = dat.train_biv, 
  #             families = FGM_Cop_BivDiscrete(marg1 = "ZALG", 
  #                                            marg2 = "ZINBI"),
  #             control = boost_control(mstop = MSTOP_COP, nu  = boost.nu.steplength, trace = TRUE), 
  #             method = 'noncyclic')
  bivDiscCopula <- bivDiscCopula[MSTOP_COP]
  
  
  mstop.bivCopula <-  vector('list')
  mstop.bivCopula$mstop <- MSTOP_COP
  mstop.bivCopula$mu1 <- bivDiscCopula$mu1$mstop()
  mstop.bivCopula$sigma1 <- bivDiscCopula$sigma1$mstop()
  
  mstop.bivCopula$mu2 <- bivDiscCopula$mu2$mstop()
  mstop.bivCopula$sigma2 <- bivDiscCopula$sigma2$mstop()
  mstop.bivCopula$nu2 <- bivDiscCopula$nu2$mstop()
  
  mstop.bivCopula$rho  <- bivDiscCopula$rho$mstop()
  
  # EXTRACT ALL COEFFICIENTS!
  coef.bivDiscCopula <- coef(bivDiscCopula, which = "")
  
  DiscreteMargins_Metrics <- vector("list")
  likBinCopula <- vector('list')
  
  
  ################ Get predictions from each covariate?
  newX <- matrix(rep(seq(from = 0, to = 1, length.out = 500),times = 10), nrow = 500)
  colnames(newX) <-  c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
  
  
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ### Copula predictions / selection
  Copula_MU1_predictions <- list()
  Copula_SIGMA1_predictions <- list()
  
  Copula_MU2_predictions <- list()
  Copula_SIGMA2_predictions <- list()
  Copula_NU2_predictions <- list()
  
  Copula_RHO_predictions <- list()
  
  ### MU 1
  Copula_MU1_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = 1, type = "link" )
  Copula_MU1_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x2", type = "link" )
  Copula_MU1_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x3", type = "link" )
  Copula_MU1_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x4", type = "link" )
  Copula_MU1_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x5", type = "link" )
  Copula_MU1_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x6", type = "link" )
  Copula_MU1_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x7", type = "link" )
  Copula_MU1_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x8", type = "link" )
  Copula_MU1_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x9", type = "link" )
  Copula_MU1_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x10", type = "link" )
  
  #### SIGMA 1
  Copula_SIGMA1_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = 1, type = "link" )
  Copula_SIGMA1_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x2", type = "link" )
  Copula_SIGMA1_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x3", type = "link" )
  Copula_SIGMA1_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x4", type = "link" )
  Copula_SIGMA1_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x5", type = "link" )
  Copula_SIGMA1_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x6", type = "link" )
  Copula_SIGMA1_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x7", type = "link" )
  Copula_SIGMA1_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x8", type = "link" )
  Copula_SIGMA1_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x9", type = "link" )
  Copula_SIGMA1_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x10", type = "link" )
  
  
  #### MU 2
  Copula_MU2_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = 1, type = "link" )
  Copula_MU2_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x2", type = "link" )
  Copula_MU2_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x3", type = "link" )
  Copula_MU2_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x4", type = "link" )
  Copula_MU2_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x5", type = "link" )
  Copula_MU2_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x6", type = "link" )
  Copula_MU2_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x7", type = "link" )
  Copula_MU2_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x8", type = "link" )
  Copula_MU2_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x9", type = "link" )
  Copula_MU2_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x10", type = "link" )
  
  #### SIGMA 2
  Copula_SIGMA2_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = 1, type = "link" )
  Copula_SIGMA2_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x2", type = "link" )
  Copula_SIGMA2_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x3", type = "link" )
  Copula_SIGMA2_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x4", type = "link" )
  Copula_SIGMA2_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x5", type = "link" )
  Copula_SIGMA2_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x6", type = "link" )
  Copula_SIGMA2_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x7", type = "link" )
  Copula_SIGMA2_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x8", type = "link" )
  Copula_SIGMA2_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x9", type = "link" )
  Copula_SIGMA2_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x10", type = "link" )
  
  
  #### NU 2
  Copula_NU2_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = 1, type = "link" )
  Copula_NU2_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x2", type = "link" )
  Copula_NU2_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x3", type = "link" )
  Copula_NU2_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x4", type = "link" )
  Copula_NU2_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x5", type = "link" )
  Copula_NU2_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x6", type = "link" )
  Copula_NU2_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x7", type = "link" )
  Copula_NU2_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x8", type = "link" )
  Copula_NU2_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x9", type = "link" )
  Copula_NU2_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x10", type = "link" )
  
  
  #### RHO
  Copula_RHO_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = 1, type = "link" )
  Copula_RHO_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x2", type = "link" )
  Copula_RHO_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x3", type = "link" )
  Copula_RHO_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x4", type = "link" )
  Copula_RHO_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x5", type = "link" )
  Copula_RHO_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x6", type = "link" )
  Copula_RHO_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x7", type = "link" )
  Copula_RHO_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x8", type = "link" )
  Copula_RHO_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x9", type = "link" )
  Copula_RHO_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x10", type = "link" )
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  
  
  
  # Predict distributional quantities
  predCopula.mu1 <- predict(bivDiscCopula$mu1, newdata = dat.test, type = 'response')
  predCopula.sigma1 <- predict(bivDiscCopula$sigma1, newdata = dat.test, type = 'response')
  
  predCopula.mu2 <- predict(bivDiscCopula$mu2, newdata = dat.test, type = 'response')
  predCopula.sigma2 <- predict(bivDiscCopula$sigma2, newdata = dat.test, type = 'response')
  predCopula.nu2 <- predict(bivDiscCopula$nu2, newdata = dat.test, type = 'response')
  
  
  predCopula.rho <- predict(bivDiscCopula$rho, newdata = dat.test, type = 'response')
  
  ######### COMPUTE PERFORMANCE: 
  
  # MSE (margin 1, then margin 2)
  DiscreteMargins_Metrics$mu1 <- mean((as.numeric(predCopula.mu1)-(as.numeric(dat.test$y1)))^2)
  DiscreteMargins_Metrics$mu2 <- mean((as.numeric(predCopula.mu2)-(as.numeric(dat.test$y2)))^2)
  
  # LOG-SCORE (LOG-LIKELIHOOD) (change to the copula likelihood)
  lik$biv <- sum(loss(mu1 = predCopula.mu1, 
                      mu2 = predCopula.mu2, 
                      sigma1 = predCopula.sigma1, 
                      sigma2 = predCopula.sigma2, 
                      nu2 = predCopula.nu2,
                      rho = predCopula.rho, 
                      y = y.test))
  
  
  ###########################################################################################################################
  ################ - univariat - ############################################################################################
  ###########################################################################################################################
  
  mstop.uni <- vector('list')
  coef.uni  <- vector('list')
  coef.uni_1  <- vector('list')
  
  univariateDiscreteEquation_1 <- list(mu = formula(y1 ~  bbs(x1) + 
                                                      bbs(x2) + 
                                                      bbs(x3) + 
                                                      bbs(x4) + 
                                                      bbs(x5) + 
                                                      bbs(x6) +
                                                      bbs(x7) + 
                                                      bbs(x8) + 
                                                      bbs(x9) + 
                                                      bbs(x10)),
                                       sigma = formula(y1 ~  bbs(x1) + 
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
  
  
  univariateDiscreteEquation_2 <- list(mu = formula(y2 ~  bbs(x1) + 
                                                      bbs(x2) + 
                                                      bbs(x3) + 
                                                      bbs(x4) + 
                                                      bbs(x5) + 
                                                      bbs(x6) +
                                                      bbs(x7) + 
                                                      bbs(x8) + 
                                                      bbs(x9) + 
                                                      bbs(x10)),
                                       sigma = formula(y2 ~  bbs(x1) + 
                                                         bbs(x2) + 
                                                         bbs(x3) + 
                                                         bbs(x4) + 
                                                         bbs(x5) + 
                                                         bbs(x6) +
                                                         bbs(x7) + 
                                                         bbs(x8) + 
                                                         bbs(x9) + 
                                                         bbs(x10),
                                       ),
                                       nu = formula(y2 ~  bbs(x1) + 
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
  
  
  
  
  # univariateDiscreteEquation_1 <- list(mu = formula(y1 ~  bols(x1) + 
  #                                                     bols(x2) + 
  #                                                     bols(x3) + 
  #                                                     bols(x4) + 
  #                                                     bols(x5) + 
  #                                                     bols(x6) +
  #                                                     bols(x7) + 
  #                                                     bols(x8) + 
  #                                                     bols(x9) + 
  #                                                     bols(x10)),
  #                                      sigma = formula(y1 ~  bols(x1) + 
  #                                                        bols(x2) + 
  #                                                        bols(x3) + 
  #                                                        bols(x4) + 
  #                                                        bols(x5) + 
  #                                                        bols(x6) +
  #                                                        bols(x7) + 
  #                                                        bols(x8) + 
  #                                                        bols(x9) + 
  #                                                        bols(x10))
  # )
  # 
  # 
  # univariateDiscreteEquation_2 <- list(mu = formula(y2 ~  bols(x1) + 
  #                                                     bols(x2) + 
  #                                                     bols(x3) + 
  #                                                     bols(x4) + 
  #                                                     bols(x5) + 
  #                                                     bols(x6) +
  #                                                     bols(x7) + 
  #                                                     bols(x8) + 
  #                                                     bols(x9) + 
  #                                                     bols(x10)),
  #                                      sigma = formula(y2 ~  bols(x1) + 
  #                                                        bols(x2) + 
  #                                                        bols(x3) + 
  #                                                        bols(x4) + 
  #                                                        bols(x5) + 
  #                                                        bols(x6) +
  #                                                        bols(x7) + 
  #                                                        bols(x8) + 
  #                                                        bols(x9) + 
  #                                                        bols(x10)),
  #                                      nu = formula(y2 ~  bols(x1) + 
  #                                                     bols(x2) + 
  #                                                     bols(x3) + 
  #                                                     bols(x4) + 
  #                                                     bols(x5) + 
  #                                                     bols(x6) +
  #                                                     bols(x7) + 
  #                                                     bols(x8) + 
  #                                                     bols(x9) + 
  #                                                     bols(x10))
  # )
  # 
  
  # - margin 1
  dat.train.mu1 = dat.train[, !(names(dat.train) %in% "y2")]
  dat.test.mu1  = dat.test[ , !(names(dat.test) %in% "y2")]
  
  glm.uni.mu1 <- gamboostLSS(univariateDiscreteEquation_1, 
                             data = dat.train.mu1, 
                             families = as.families(fname = "ZALG", stabilization = "L2"), 
                             control = boost_control(risk = 'oobag', nu = boost.nu.steplength, mstop = 2000, trace = T), 
                             method = "noncyclic",
                             weights = weight.mstop)
  
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  if(mstop.uni$mu1 >= 1990){
    glm.uni.mu1[5000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  if(mstop.uni$mu1 >= 4990){
    glm.uni.mu1[10000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  if(mstop.uni$mu1 >= 9990){
    glm.uni.mu1[20000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  # dat.train_mu1 <- dat.train.mu1[weight.mstop == 1, ]
  # 
  # glm.uni.mu1 <- gamboostLSS(univariateDiscreteEquation_1, 
  #                            data = dat.train_mu1, 
  #                            families = as.families(fname = "ZALG"), 
  #                         control = boost_control(mstop = mstop.uni$mu1, nu  = boost.nu.steplength, trace = T), 
  #                         method = "noncyclic")
  glm.uni.mu1 <- glm.uni.mu1[mstop.uni$mu1]
  
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 1  
  coef.uni$mu1 <- coef(glm.uni.mu1, which = '')
  
  DiscreteMargins_Metrics$Univariate_mu1 <- mean((as.numeric(predict(glm.uni.mu1$mu, newdata = dat.test.mu1, type = "response")) - (as.numeric(dat.test.mu1$y1)))^2)
  
  # - margin 2
  dat.train.mu2 =  dat.train[, !(names(dat.train) %in% "y1")]
  dat.test.mu2 = dat.test[, !(names(dat.test) %in% "y1")]
  
  glm.uni.mu2 <- gamboostLSS(univariateDiscreteEquation_2, 
                             data = dat.train.mu2, 
                             families = as.families(fname = "ZINBI", stabilization = "L2"), 
                             control = boost_control(risk = 'oobag', nu = boost.nu.steplength, trace = TRUE, mstop = 2000), 
                             method = "noncyclic",
                             weights = weight.mstop)
  
  
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  if(mstop.uni$mu2 >= 1990){
    glm.uni.mu2[5000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  if(mstop.uni$mu2 >= 4990){
    glm.uni.mu2[10000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  if(mstop.uni$mu2 >= 9990){
    glm.uni.mu2[20000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  # dat.train_mu2 <- dat.train.mu2[weight.mstop == 1, ]
  # 
  # glm.uni.mu2 <- gamboostLSS(univariateDiscreteEquation_2, 
  #                            data = dat.train_mu2, 
  #                            families = as.families(fname = "ZINBI"), 
  #                         control = boost_control(mstop = mstop.uni$mu2, nu  = boost.nu.steplength, trace=TRUE), method = "noncyclic")
  glm.uni.mu2 <- glm.uni.mu2[mstop.uni$mu2]
  
  mstop.uni$Margin1 <- mstop(glm.uni.mu1)
  
  mstop.uni$Margin2 <- mstop(glm.uni.mu2)
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 2  
  coef.uni$mu2 <- coef(glm.uni.mu2,which = '')
  
  DiscreteMargins_Metrics$Univariate_mu2 <- mean((as.numeric(predict(glm.uni.mu2$mu, newdata = dat.test.mu2, type = "response"))-(as.numeric(dat.test.mu2$y2)))^2)
  
  # Likelihood
  pred.mu1.uni <- as.numeric(predict(glm.uni.mu1$mu, newdata = dat.test.mu1, type = "response"))
  pred.sigma1.uni <- as.numeric(predict(glm.uni.mu1$sigma, newdata = dat.test.mu1, type = "response"))
  
  pred.mu2.uni <- as.numeric(predict(glm.uni.mu2$mu, newdata = dat.test.mu2, type = "response"))
  pred.sigma2.uni <- as.numeric(predict(glm.uni.mu2$sigma, newdata = dat.test.mu2, type = "response"))
  pred.nu2.uni <- as.numeric(predict(glm.uni.mu2$nu, newdata = dat.test.mu2, type = "response"))
  
  
  mu1.uni.loglik <- sum(-dZALG(x = dat.test.mu1$y1, mu = pred.mu1.uni, sigma = pred.sigma1.uni, log = T))
  mu2.uni.loglik <- sum(-dZINBI(x = dat.test.mu2$y2, mu = pred.mu2.uni, sigma = pred.sigma2.uni, nu = pred.nu2.uni, log = T))
  
  lik$uni<- mu1.uni.loglik + mu2.uni.loglik
  
  
  n.ges = c(n.train, n.mstop, n.test)
  pred.ges <- NULL#list(pred.mu1, pred.mu2, pred.or)
  pred.uni <- list(mu1 = pred.mu1.uni, 
                   sigma1 = pred.sigma1.uni,
                   mu2 = pred.mu2.uni,
                   sigma2 = pred.sigma2.uni,
                   nu2 = pred.nu2.uni
  )
  
  pred.cop <- list(mu1 = predCopula.mu1, 
                   sigma1 = predCopula.sigma1,
                   mu2 = predCopula.mu2,
                   sigma2 = predCopula.sigma2,
                   nu2 = predCopula.nu2,
                   rho = predCopula.rho)
  
  
  ##### Predictions:
  UNIVARIATE_MU1_predictions <- vector("list")
  UNIVARIATE_SIGMA1_predictions <- vector("list")
  
  UNIVARIATE_MU2_predictions <- vector("list")
  UNIVARIATE_SIGMA2_predictions <- vector("list")
  UNIVARIATE_NU2_predictions <- vector("list")
  
  
  #### MU 1
  UNIVARIATE_MU1_predictions$x1 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = 1, type = "link" )
  UNIVARIATE_MU1_predictions$x2 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x2", type = "link" )
  UNIVARIATE_MU1_predictions$x3 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x3", type = "link" )
  UNIVARIATE_MU1_predictions$x4 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x4", type = "link" )
  UNIVARIATE_MU1_predictions$x5 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x5", type = "link" )
  UNIVARIATE_MU1_predictions$x6 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x6", type = "link" )
  UNIVARIATE_MU1_predictions$x7 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x7", type = "link" )
  UNIVARIATE_MU1_predictions$x8 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x8", type = "link" )
  UNIVARIATE_MU1_predictions$x9 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x9", type = "link" )
  UNIVARIATE_MU1_predictions$x10 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x10", type = "link" )
  
  #### SIGMA 1
  UNIVARIATE_SIGMA1_predictions$x1 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = 1, type = "link" )
  UNIVARIATE_SIGMA1_predictions$x2 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x2", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x3 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x3", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x4 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x4", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x5 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x5", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x6 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x6", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x7 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x7", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x8 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x8", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x9 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x9", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x10 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x10", type = "link" )
  
  
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
  
  
  #### NU 2
  UNIVARIATE_NU2_predictions$x1 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = 1, type = "link" )
  UNIVARIATE_NU2_predictions$x2 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x2", type = "link" )
  UNIVARIATE_NU2_predictions$x3 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x3", type = "link" )
  UNIVARIATE_NU2_predictions$x4 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x4", type = "link" )
  UNIVARIATE_NU2_predictions$x5 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x5", type = "link" )
  UNIVARIATE_NU2_predictions$x6 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x6", type = "link" )
  UNIVARIATE_NU2_predictions$x7 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x7", type = "link" )
  UNIVARIATE_NU2_predictions$x8 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x8", type = "link" )
  UNIVARIATE_NU2_predictions$x9 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x9", type = "link" )
  UNIVARIATE_NU2_predictions$x10 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x10", type = "link" )
  
  
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
    sample_uni <- data.gen.bivdisc_energyscore(n = 1000, 
                                               mu1 = pred.mu1.uni[i], 
                                               mu2 = pred.mu2.uni[i], 
                                               sigma1 = pred.sigma1.uni[i], 
                                               sigma2 = pred.sigma2.uni[i],
                                               nu = check_extremes_nuZINBI(pred.nu2.uni[i]),
                                               theta = 1)
    
    pred_sample_uni[1,] <- sample_uni[,1]
    pred_sample_uni[2,] <- sample_uni[,2]
    
    es_uni[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_uni) 
    
    
    # copula
    sample_cop <- data.gen.bivdisc_energyscore(n = 1000, 
                                               mu1 = predCopula.mu1[i],  
                                               mu2 = predCopula.mu2[i], 
                                               sigma1 = predCopula.sigma1[i], 
                                               sigma2 = predCopula.sigma2[i],
                                               nu = check_extremes_nuZINBI(predCopula.nu2[i]),
                                               theta = predCopula.rho[i])
    
    pred_sample_cop[1, ] <- sample_cop[,1]
    pred_sample_cop[2, ] <- sample_cop[,2]
    
    es_cop[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_cop) 
    
    
  }
  
  energy_score <- list()
  #energy_score$biv <- mean(es_biv)
  energy_score$uni <- mean(es_uni)
  energy_score$cop <- mean(es_cop)
  
  
  
  ########################################### We gather: mstop, predictions, MSE margin 1, margin 2, loglikelihood, energy score.
  output <- list(TrueBeta = TrueBeta, 
                 n = n.ges, 
                 p  = p, 
                 Likelihood = lik,
                 predict.uni = pred.uni, 
                 predict.cop = pred.cop,
                 energy_score = energy_score,  
                 ####
                 Coefficients.uni = coef.uni, 
                 mstop.uni = mstop.uni,
                 ###
                 DiscreteMarginMetrics = DiscreteMargins_Metrics,
                 ###
                 Copula_Predictions = list(MU1 = Copula_MU1_predictions, 
                                           SIGMA1 = Copula_SIGMA1_predictions,
                                           MU2 = Copula_MU2_predictions,
                                           SIGMA2 = Copula_SIGMA2_predictions,
                                           NU2 = Copula_NU2_predictions,
                                           RHO = Copula_RHO_predictions),
                 ###
                 Univariate_Predictions = list(MU1 = UNIVARIATE_MU1_predictions, 
                                               SIGMA1 = UNIVARIATE_SIGMA1_predictions,
                                               MU2 = UNIVARIATE_MU2_predictions, 
                                               SIGMA2 = UNIVARIATE_SIGMA2_predictions,
                                               NU2 = UNIVARIATE_NU2_predictions),
                 ###
                 CoefficientsCOPULA = coef.bivDiscCopula, 
                 oobag.riskCOPULA = oobag.risk.cop, 
                 energy_scoreCOPULA = energy_score,  
                 mstopCOPULA = mstop.bivCopula,
                 TrueKendallRange = TrueKendallTauRange
  )
  
  
  return(output)
  
}

#### USING FRANK COPULA
sim_FRANK <- function(seed, n.train, n.mstop, p, n.test, corr, boost.nu.steplength = 0.1){
  
  
  loss <- function(mu1, mu2, sigma1, sigma2, rho, y, nu2 = NULL){
    
    
    # Margin 1 is ZALG:
    F1 <- pdffz(pZALG(y[,1], mu = mu1, sigma = sigma1))
    pdf1 <- pdffz(dZALG(y[,1], mu = mu1, sigma = sigma1))
    
    # Margin 2 is NBI:
    F2 <- pdffz(pZINBI(y[,2], mu = mu2, sigma = sigma2, nu = nu2))
    pdf2 <- pdffz(dZINBI(y[,2], mu = mu2, sigma = sigma2, nu = nu2))
    
    # Copula parameter 
    thet <- rho
    
    # Minus 1
    CDF1m1 <- pdffz(F1 - pdf1)
    CDF2m1 <- pdffz(F2 - pdf2)
    
    ### Copula Terms
    bit <- -expm1(-thet) 
    
    T1 <- pdffz( -(1/thet)*log( (bit - (1 - exp(-thet*F1))*(1 - exp(-thet*F2)))/bit ) )
    
    T2 <- pdffz( -(1/thet)*log( (bit - (1 - exp(-thet*CDF1m1))*(1 - exp(-thet*F2)))/bit )  ) 
    
    T3 <- pdffz( -(1/thet)*log( (bit - (1 - exp(-thet*F1))*(1 - exp(-thet*CDF2m1)))/bit ) )
    
    T4 <- pdffz( -(1/thet)*log( (bit - (1 - exp(-thet*CDF1m1))*(1 - exp(-thet*CDF2m1)))/bit ) ) 
    
    L_Total <- pdffz( T1 - T2 - T3 + T4 )
    
    # Negative log-likelihood
    return(- log( L_Total ) )
    
  }
  
  data.gen.bivdisc <- function(mu1, mu2, sigma1, sigma2, theta, nu2 = NULL){
    
    y1y2 <- matrix(0, ncol = 2, nrow = length(mu1))
    
    for(i in 1:length(mu1)){
      
      paramlist1 <- list(mu = mu1[i], 
                         sigma = sigma1[i])
      
      
      paramlist2 <- list(mu = mu2[i], 
                         sigma = sigma2[i],
                         nu = nu2[i])
      
      copObj <- copula::frankCopula(param = theta[i])
      
      copThing <- mvdc(copObj, margins = c("ZALG", "ZINBI"), paramMargins = list(paramlist1, paramlist2))
      
      y1y2[i,] <- rMvdc(copThing, n = 1)
    }
    
    
    dat <- data.frame(y1y2[,1], y1y2[,2])
    
    return(dat)
  }
  
  data.gen.bivdisc_energyscore <- function(n, mu1, mu2, sigma1, sigma2, theta, nu = NULL){
    
    
    paramlist1 <- list(mu = mu1, sigma = sigma1)
    
    paramlist2 <- list(mu = mu2, sigma = sigma2, nu = nu) 
    
    copObj <- copula::frankCopula(param = theta)
    
    copThing <- mvdc(copObj, margins = c("ZALG", "ZINBI"), paramMargins = list(paramlist1, paramlist2))
    
    y1y2 <- rMvdc(copThing, n = n)
    
    dat <- data.frame(y1y2[,1], y1y2[,2])
    
    
    return(dat)
  }
  
  
  
  
  ## Coefficients:
  # x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
  mu1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  mu2 = c(0, 0, 0, 0, 0,  0, 0, 0, 0, 0) 
  
  # x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
  sigma1 <- c(0, 0, 0, 0, 0,  0, 0, 0, 0, 0) 
  sigma2 <- c(0, 0, 0, 0, 0,  0, 0, 0, 0, 0) 
  
  nu2 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0) 
  
  # x1, x2, x3, x4, x5, x6, x7, x8, x9, x10
  or  = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  
  ### function for mu1
  nl_function_mu1 <- function(x){
    
    
    nl <- sqrt(x)*x^2 - 2*cos(3*x)
    
    return(nl*0.5)
    
  }
  
  ### function for sigma1
  nl_function_sigma1 <- function(x){
    
    func <- -40*(x^(3/2) - x^(4/3))
    
    func <- 2 * func
    
    return(func)
  }
  
  nl_func_sigma1 <- function(x){
    
    func <- -(log(x)) + cos(x)
    
    func <- 2 * func
    
    return(func)
  }
  
  ## function for mu2
  nl_function_mu2 <- function(x){
    
    nl <- -0.7*exp(x^2) + exp(x^0.4) - 0*log(x+0.2) 
    
    return(nl)
  }
  
  nl_function_sigma2 <- function(x){
    
    #return( -1.5*(1.5*cos(2 * x) + 0.5*log(x)) )
    return( -1.5*(1.5*cos(2 * x) - 0*log(x + 0.15) + 3*tanh(1*x) ) )
  }
  
 
  
  nl_function_nu2 <- function(x){
    
    func <- -(sin(x) - exp(x)^2)
    
    func <- 0.7*func
    
    return(func)
  }
  
  nl_function_rho <- function(x){
    
    return(30*sin(x*4))
    #return( 2*(11*sin(x*4) + log(x^1.5))  )
  }
  
  TrueBeta <-  vector('list')
  TrueBeta$mu1 <- mu1
  TrueBeta$mu2 <- mu2
  
  TrueBeta$sigma1 <- sigma1
  TrueBeta$sigma2 <- sigma2
  
  TrueBeta$nu2 <- nu2
  
  TrueBeta$or <- or
  
  set.seed(seed)
  
  
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
  
  ############################################################################################################################################
  ############################################################################################################################################
  # - test data
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
  
  ############################################################################################################################################
  # GENERATION OF DISTRIBUTION PARAMETERS AND STUFF 
  # predictor
  
  ##### TRAIN DATA
  train.eta.mu1_center <- (x.train %*% mu1)
  train.eta.mu1_center <- train.eta.mu1_center + nl_function_mu1(x1)
  
  ##### TEST DATA
  test.eta.mu1_center <- (x.test %*% mu1)
  test.eta.mu1_center <- test.eta.mu1_center + nl_function_mu1(x1.test)
  
  # apply response function:
  train.eta.mu1<-  plogis(train.eta.mu1_center) 
  range(train.eta.mu1)
  
  test.eta.mu1<-  plogis(test.eta.mu1_center) 
  ########################################################################################################################
  #### SIGMA FOR ZALG RESPONSE
  train.eta.sigma1 <-  x.train %*% sigma1
  
  train.eta.sigma1 <- train.eta.sigma1 - 2 + nl_function_sigma1(x3) 
  train.eta.sigma1 <- plogis(train.eta.sigma1)
  range(train.eta.sigma1)
  
  ##### TEST DATA
  test.eta.sigma1 <-  x.test %*% sigma1
  test.eta.sigma1 <- test.eta.sigma1 - 2 + nl_function_sigma1(x3.test)
  test.eta.sigma1 <- plogis(test.eta.sigma1)
  ########################################################################################################################
  # predictor and apply response function: 
  train.eta.mu2 <-  x.train %*% mu2
  train.eta.mu2 <-  train.eta.mu2 + nl_function_mu2(x2)
  train.eta.mu2 <-  exp(train.eta.mu2) 
  range(train.eta.mu2)
  
  
  test.eta.mu2 <-  x.test %*% mu2
  test.eta.mu2 <-  test.eta.mu2 + nl_function_mu2(x2.test)
  test.eta.mu2 <-  exp(test.eta.mu2) 
  ########################################################################################################################
  # predictor and apply response function: 
  train.eta.sigma2 <-  x.train %*% sigma2
  train.eta.sigma2 <-  train.eta.sigma2 + nl_function_sigma2(x4)
  train.eta.sigma2 <-  exp(train.eta.sigma2) 
  range(train.eta.sigma2)
  
  test.eta.sigma2 <-  x.test %*% sigma2
  test.eta.sigma2 <-  test.eta.sigma2 +  nl_function_sigma2(x4.test)
  test.eta.sigma2 <-  exp(test.eta.sigma2) 
  ########################################################################################################################
  # predictor and apply response function: 
  train.eta.nu2 <-  x.train %*% nu2
  train.eta.nu2 <-  train.eta.nu2 + nl_function_nu2(x1) - 3
  train.eta.nu2 <-  plogis(train.eta.nu2) 
  range(train.eta.nu2)
  
  test.eta.nu2 <-  x.test %*% nu2
  test.eta.nu2 <-  test.eta.nu2 + nl_function_nu2(x1.test) - 3
  test.eta.nu2 <-  plogis(test.eta.nu2) 
  
  ########################################################################################################################
  ### Copula parameter
  train.eta.or_center <- (x.train %*% or)
  train.eta.or_center <- train.eta.or_center + nl_function_rho(x3) 
  
  train.copula.parameter <- train.eta.or_center
  range(train.copula.parameter)
  
  TrueKendallTauRange <- VineCopula::BiCopPar2Tau(5, par = range(train.copula.parameter))  
  
  test.eta.or_center <- (x.test %*% or)
  test.eta.or_center <- test.eta.or_center + nl_function_rho(x3.test)
  
  test.copula.parameter <- test.eta.or_center
  ############################################################################################################################################
  # SAMPLING OF BIVARIATE RESPONSE FROM COPULA
  y.train <- data.gen.bivdisc(mu1 = train.eta.mu1, mu2 = train.eta.mu2, sigma1 = train.eta.sigma1, sigma2 = train.eta.sigma2, theta = train.copula.parameter, nu2 = train.eta.nu2)
  
  colnames(y.train)<- c('y1','y2')
  dat.train <- data.frame(y.train,x.train)
  
  
  y.test <- data.gen.bivdisc(mu1 = test.eta.mu1, mu2 = test.eta.mu2, sigma1 = test.eta.sigma1, sigma2 = test.eta.sigma2, theta = test.copula.parameter, nu2 = test.eta.nu2)
  
  colnames(y.test)<- c('y1','y2')
  dat.test <- data.frame(y.test,x.test)
  
  ############################################################################################################################################
  # FITTIIIINGG
  mstop.bivBern <-  vector('list')
  
  AUCbivBern <-  vector('list')
  BrierbivBern <-  vector('list')
  lik <- vector('list')
  MSEbivDisc <- vector("list")
  
  
  
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
                           sigma1 = model_equation,
                           sigma2 = model_equation,
                           rho = model_equation
  )
  
  
  bivModel_Formula <- list(mu1 = model_equation,
                           mu2 = model_equation,
                           sigma1 = model_equation,
                           sigma2 = model_equation,
                           nu2 = model_equation,
                           rho = model_equation
  )
  
  #############################################
  ################# - COPULA - ################
  #############################################
  # Bivariate copula model
  bivDiscCopula <- gamboostLSS(bivModel_Formula, data = dat.train, 
                               families = Frank_Cop_BivDiscrete(marg1 = "ZALG", 
                                                              marg2 = "ZINBI",
                                                              stabilization = "L2"), 
                               control = boost_control(mstop = 2000, risk = 'oobag',
                                                       nu  = boost.nu.steplength, trace = TRUE), 
                               method = 'noncyclic', weights = weight.mstop)
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 1900){
    mstop(bivDiscCopula) <- 3000
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 2900){
    mstop(bivDiscCopula) <- 4000
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 3990){
    bivDiscCopula[8000]
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 7990){
    bivDiscCopula[10000]
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 9990){
    bivDiscCopula[15000]
  }
  
  MSTOP_COP <- which.min(risk(bivDiscCopula,merge = T))
  
  if(MSTOP_COP >= 14990){
    bivDiscCopula[20000]
  }
  
  
  
  
  MSTOP_COP <- which.min(risk(bivDiscCopula, merge = T))
  oobag.risk.cop <- risk(bivDiscCopula,merge = T)
  
  # rm(bivDiscCopula)
  # dat.train_biv <- dat.train[weight.mstop == 1, ]
  # 
  # # RE-FIT until OPTIMAL MSTOP
  # bivDiscCopula <- gamboostLSS(bivModel_Formula, data = dat.train_biv, 
  #             families = FGM_Cop_BivDiscrete(marg1 = "ZALG", 
  #                                            marg2 = "ZINBI"),
  #             control = boost_control(mstop = MSTOP_COP, nu  = boost.nu.steplength, trace = TRUE), 
  #             method = 'noncyclic')
  bivDiscCopula <- bivDiscCopula[MSTOP_COP]
  
  
  mstop.bivCopula <-  vector('list')
  mstop.bivCopula$mstop <- MSTOP_COP
  mstop.bivCopula$mu1 <- bivDiscCopula$mu1$mstop()
  mstop.bivCopula$sigma1 <- bivDiscCopula$sigma1$mstop()
  
  mstop.bivCopula$mu2 <- bivDiscCopula$mu2$mstop()
  mstop.bivCopula$sigma2 <- bivDiscCopula$sigma2$mstop()
  mstop.bivCopula$nu2 <- bivDiscCopula$nu2$mstop()
  
  mstop.bivCopula$rho  <- bivDiscCopula$rho$mstop()
  
  # EXTRACT ALL COEFFICIENTS!
  coef.bivDiscCopula <- coef(bivDiscCopula, which = "")
  
  DiscreteMargins_Metrics <- vector("list")
  likBinCopula <- vector('list')
  
  
  ################ Get predictions from each covariate?
  newX <- matrix(rep(seq(from = 0, to = 1, length.out = 500),times = 10), nrow = 500)
  colnames(newX) <-  c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
  
  
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ### Copula predictions / selection
  Copula_MU1_predictions <- list()
  Copula_SIGMA1_predictions <- list()
  
  Copula_MU2_predictions <- list()
  Copula_SIGMA2_predictions <- list()
  Copula_NU2_predictions <- list()
  
  Copula_RHO_predictions <- list()
  
  ### MU 1
  Copula_MU1_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = 1, type = "link" )
  Copula_MU1_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x2", type = "link" )
  Copula_MU1_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x3", type = "link" )
  Copula_MU1_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x4", type = "link" )
  Copula_MU1_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x5", type = "link" )
  Copula_MU1_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x6", type = "link" )
  Copula_MU1_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x7", type = "link" )
  Copula_MU1_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x8", type = "link" )
  Copula_MU1_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x9", type = "link" )
  Copula_MU1_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu1", which = "x10", type = "link" )
  
  #### SIGMA 1
  Copula_SIGMA1_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = 1, type = "link" )
  Copula_SIGMA1_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x2", type = "link" )
  Copula_SIGMA1_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x3", type = "link" )
  Copula_SIGMA1_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x4", type = "link" )
  Copula_SIGMA1_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x5", type = "link" )
  Copula_SIGMA1_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x6", type = "link" )
  Copula_SIGMA1_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x7", type = "link" )
  Copula_SIGMA1_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x8", type = "link" )
  Copula_SIGMA1_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x9", type = "link" )
  Copula_SIGMA1_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma1", which = "x10", type = "link" )
  
  
  #### MU 2
  Copula_MU2_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = 1, type = "link" )
  Copula_MU2_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x2", type = "link" )
  Copula_MU2_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x3", type = "link" )
  Copula_MU2_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x4", type = "link" )
  Copula_MU2_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x5", type = "link" )
  Copula_MU2_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x6", type = "link" )
  Copula_MU2_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x7", type = "link" )
  Copula_MU2_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x8", type = "link" )
  Copula_MU2_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x9", type = "link" )
  Copula_MU2_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "mu2", which = "x10", type = "link" )
  
  #### SIGMA 2
  Copula_SIGMA2_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = 1, type = "link" )
  Copula_SIGMA2_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x2", type = "link" )
  Copula_SIGMA2_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x3", type = "link" )
  Copula_SIGMA2_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x4", type = "link" )
  Copula_SIGMA2_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x5", type = "link" )
  Copula_SIGMA2_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x6", type = "link" )
  Copula_SIGMA2_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x7", type = "link" )
  Copula_SIGMA2_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x8", type = "link" )
  Copula_SIGMA2_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x9", type = "link" )
  Copula_SIGMA2_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "sigma2", which = "x10", type = "link" )
  
  
  #### NU 2
  Copula_NU2_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = 1, type = "link" )
  Copula_NU2_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x2", type = "link" )
  Copula_NU2_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x3", type = "link" )
  Copula_NU2_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x4", type = "link" )
  Copula_NU2_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x5", type = "link" )
  Copula_NU2_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x6", type = "link" )
  Copula_NU2_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x7", type = "link" )
  Copula_NU2_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x8", type = "link" )
  Copula_NU2_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x9", type = "link" )
  Copula_NU2_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "nu2", which = "x10", type = "link" )
  
  
  #### RHO
  Copula_RHO_predictions$x1 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = 1, type = "link" )
  Copula_RHO_predictions$x2 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x2", type = "link" )
  Copula_RHO_predictions$x3 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x3", type = "link" )
  Copula_RHO_predictions$x4 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x4", type = "link" )
  Copula_RHO_predictions$x5 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x5", type = "link" )
  Copula_RHO_predictions$x6 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x6", type = "link" )
  Copula_RHO_predictions$x7 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x7", type = "link" )
  Copula_RHO_predictions$x8 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x8", type = "link" )
  Copula_RHO_predictions$x9 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x9", type = "link" )
  Copula_RHO_predictions$x10 <- predict(bivDiscCopula, data.frame(newX), parameter = "rho", which = "x10", type = "link" )
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  ######################################################################################################################################################################################################################################################
  
  
  
  # Predict distributional quantities
  predCopula.mu1 <- predict(bivDiscCopula$mu1, newdata = dat.test, type = 'response')
  predCopula.sigma1 <- predict(bivDiscCopula$sigma1, newdata = dat.test, type = 'response')
  
  predCopula.mu2 <- predict(bivDiscCopula$mu2, newdata = dat.test, type = 'response')
  predCopula.sigma2 <- predict(bivDiscCopula$sigma2, newdata = dat.test, type = 'response')
  predCopula.nu2 <- predict(bivDiscCopula$nu2, newdata = dat.test, type = 'response')
  
  
  predCopula.rho <- predict(bivDiscCopula$rho, newdata = dat.test, type = 'response')
  
  ######### COMPUTE PERFORMANCE: 
  
  # MSE (margin 1, then margin 2)
  DiscreteMargins_Metrics$mu1 <- mean((as.numeric(predCopula.mu1)-(as.numeric(dat.test$y1)))^2)
  DiscreteMargins_Metrics$mu2 <- mean((as.numeric(predCopula.mu2)-(as.numeric(dat.test$y2)))^2)
  
  # LOG-SCORE (LOG-LIKELIHOOD) (change to the copula likelihood)
  lik$biv <- sum(loss(mu1 = predCopula.mu1, 
                      mu2 = predCopula.mu2, 
                      sigma1 = predCopula.sigma1, 
                      sigma2 = predCopula.sigma2, 
                      nu2 = predCopula.nu2,
                      rho = predCopula.rho, 
                      y = y.test))
  
  
  ###########################################################################################################################
  ################ - univariat - ############################################################################################
  ###########################################################################################################################
  
  mstop.uni <- vector('list')
  coef.uni  <- vector('list')
  coef.uni_1  <- vector('list')
  
  univariateDiscreteEquation_1 <- list(mu = formula(y1 ~  bbs(x1) + 
                                                      bbs(x2) + 
                                                      bbs(x3) + 
                                                      bbs(x4) + 
                                                      bbs(x5) + 
                                                      bbs(x6) +
                                                      bbs(x7) + 
                                                      bbs(x8) + 
                                                      bbs(x9) + 
                                                      bbs(x10)),
                                       sigma = formula(y1 ~  bbs(x1) + 
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
  
  
  univariateDiscreteEquation_2 <- list(mu = formula(y2 ~  bbs(x1) + 
                                                      bbs(x2) + 
                                                      bbs(x3) + 
                                                      bbs(x4) + 
                                                      bbs(x5) + 
                                                      bbs(x6) +
                                                      bbs(x7) + 
                                                      bbs(x8) + 
                                                      bbs(x9) + 
                                                      bbs(x10)),
                                       sigma = formula(y2 ~  bbs(x1) + 
                                                         bbs(x2) + 
                                                         bbs(x3) + 
                                                         bbs(x4) + 
                                                         bbs(x5) + 
                                                         bbs(x6) +
                                                         bbs(x7) + 
                                                         bbs(x8) + 
                                                         bbs(x9) + 
                                                         bbs(x10),
                                       ),
                                       nu = formula(y2 ~  bbs(x1) + 
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
  
  
  
  
  # univariateDiscreteEquation_1 <- list(mu = formula(y1 ~  bols(x1) + 
  #                                                     bols(x2) + 
  #                                                     bols(x3) + 
  #                                                     bols(x4) + 
  #                                                     bols(x5) + 
  #                                                     bols(x6) +
  #                                                     bols(x7) + 
  #                                                     bols(x8) + 
  #                                                     bols(x9) + 
  #                                                     bols(x10)),
  #                                      sigma = formula(y1 ~  bols(x1) + 
  #                                                        bols(x2) + 
  #                                                        bols(x3) + 
  #                                                        bols(x4) + 
  #                                                        bols(x5) + 
  #                                                        bols(x6) +
  #                                                        bols(x7) + 
  #                                                        bols(x8) + 
  #                                                        bols(x9) + 
  #                                                        bols(x10))
  # )
  # 
  # 
  # univariateDiscreteEquation_2 <- list(mu = formula(y2 ~  bols(x1) + 
  #                                                     bols(x2) + 
  #                                                     bols(x3) + 
  #                                                     bols(x4) + 
  #                                                     bols(x5) + 
  #                                                     bols(x6) +
  #                                                     bols(x7) + 
  #                                                     bols(x8) + 
  #                                                     bols(x9) + 
  #                                                     bols(x10)),
  #                                      sigma = formula(y2 ~  bols(x1) + 
  #                                                        bols(x2) + 
  #                                                        bols(x3) + 
  #                                                        bols(x4) + 
  #                                                        bols(x5) + 
  #                                                        bols(x6) +
  #                                                        bols(x7) + 
  #                                                        bols(x8) + 
  #                                                        bols(x9) + 
  #                                                        bols(x10)),
  #                                      nu = formula(y2 ~  bols(x1) + 
  #                                                     bols(x2) + 
  #                                                     bols(x3) + 
  #                                                     bols(x4) + 
  #                                                     bols(x5) + 
  #                                                     bols(x6) +
  #                                                     bols(x7) + 
  #                                                     bols(x8) + 
  #                                                     bols(x9) + 
  #                                                     bols(x10))
  # )
  # 
  
  # - margin 1
  dat.train.mu1 = dat.train[, !(names(dat.train) %in% "y2")]
  dat.test.mu1  = dat.test[ , !(names(dat.test) %in% "y2")]
  
  glm.uni.mu1 <- gamboostLSS(univariateDiscreteEquation_1, 
                             data = dat.train.mu1, 
                             families = as.families(fname = "ZALG", stabilization = "L2"), 
                             control = boost_control(risk = 'oobag', nu = boost.nu.steplength, mstop = 2000, trace = T), 
                             method = "noncyclic",
                             weights = weight.mstop)
  
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  if(mstop.uni$mu1 >= 1990){
    glm.uni.mu1[5000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  if(mstop.uni$mu1 >= 4990){
    glm.uni.mu1[10000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  if(mstop.uni$mu1 >= 9990){
    glm.uni.mu1[20000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1, merge = TRUE))
  
  # dat.train_mu1 <- dat.train.mu1[weight.mstop == 1, ]
  # 
  # glm.uni.mu1 <- gamboostLSS(univariateDiscreteEquation_1, 
  #                            data = dat.train_mu1, 
  #                            families = as.families(fname = "ZALG"), 
  #                         control = boost_control(mstop = mstop.uni$mu1, nu  = boost.nu.steplength, trace = T), 
  #                         method = "noncyclic")
  glm.uni.mu1 <- glm.uni.mu1[mstop.uni$mu1]
  
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 1  
  coef.uni$mu1 <- coef(glm.uni.mu1, which = '')
  
  DiscreteMargins_Metrics$Univariate_mu1 <- mean((as.numeric(predict(glm.uni.mu1$mu, newdata = dat.test.mu1, type = "response")) - (as.numeric(dat.test.mu1$y1)))^2)
  
  # - margin 2
  dat.train.mu2 =  dat.train[, !(names(dat.train) %in% "y1")]
  dat.test.mu2 = dat.test[, !(names(dat.test) %in% "y1")]
  
  glm.uni.mu2 <- gamboostLSS(univariateDiscreteEquation_2, 
                             data = dat.train.mu2, 
                             families = as.families(fname = "ZINBI", stabilization = "L2"), 
                             control = boost_control(risk = 'oobag', nu = boost.nu.steplength, trace = TRUE, mstop = 2000), 
                             method = "noncyclic",
                             weights = weight.mstop)
  
  
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  if(mstop.uni$mu2 >= 1990){
    glm.uni.mu2[5000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  if(mstop.uni$mu2 >= 4990){
    glm.uni.mu2[10000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  if(mstop.uni$mu2 >= 9990){
    glm.uni.mu2[20000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2, merge = TRUE))
  
  # dat.train_mu2 <- dat.train.mu2[weight.mstop == 1, ]
  # 
  # glm.uni.mu2 <- gamboostLSS(univariateDiscreteEquation_2, 
  #                            data = dat.train_mu2, 
  #                            families = as.families(fname = "ZINBI"), 
  #                         control = boost_control(mstop = mstop.uni$mu2, nu  = boost.nu.steplength, trace=TRUE), method = "noncyclic")
  glm.uni.mu2 <- glm.uni.mu2[mstop.uni$mu2]
  
  mstop.uni$Margin1 <- mstop(glm.uni.mu1)
  
  mstop.uni$Margin2 <- mstop(glm.uni.mu2)
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 2  
  coef.uni$mu2 <- coef(glm.uni.mu2,which = '')
  
  DiscreteMargins_Metrics$Univariate_mu2 <- mean((as.numeric(predict(glm.uni.mu2$mu, newdata = dat.test.mu2, type = "response"))-(as.numeric(dat.test.mu2$y2)))^2)
  
  # Likelihood
  pred.mu1.uni <- as.numeric(predict(glm.uni.mu1$mu, newdata = dat.test.mu1, type = "response"))
  pred.sigma1.uni <- as.numeric(predict(glm.uni.mu1$sigma, newdata = dat.test.mu1, type = "response"))
  
  pred.mu2.uni <- as.numeric(predict(glm.uni.mu2$mu, newdata = dat.test.mu2, type = "response"))
  pred.sigma2.uni <- as.numeric(predict(glm.uni.mu2$sigma, newdata = dat.test.mu2, type = "response"))
  pred.nu2.uni <- as.numeric(predict(glm.uni.mu2$nu, newdata = dat.test.mu2, type = "response"))
  
  
  mu1.uni.loglik <- sum(-dZALG(x = dat.test.mu1$y1, mu = pred.mu1.uni, sigma = pred.sigma1.uni, log = T))
  mu2.uni.loglik <- sum(-dZINBI(x = dat.test.mu2$y2, mu = pred.mu2.uni, sigma = pred.sigma2.uni, nu = pred.nu2.uni, log = T))
  
  lik$uni<- mu1.uni.loglik + mu2.uni.loglik
  
  
  n.ges = c(n.train, n.mstop, n.test)
  pred.ges <- NULL#list(pred.mu1, pred.mu2, pred.or)
  pred.uni <- list(mu1 = pred.mu1.uni, 
                   sigma1 = pred.sigma1.uni,
                   mu2 = pred.mu2.uni,
                   sigma2 = pred.sigma2.uni,
                   nu2 = pred.nu2.uni
  )
  
  pred.cop <- list(mu1 = predCopula.mu1, 
                   sigma1 = predCopula.sigma1,
                   mu2 = predCopula.mu2,
                   sigma2 = predCopula.sigma2,
                   nu2 = predCopula.nu2,
                   rho = predCopula.rho)
  
  
  ##### Predictions:
  UNIVARIATE_MU1_predictions <- vector("list")
  UNIVARIATE_SIGMA1_predictions <- vector("list")
  
  UNIVARIATE_MU2_predictions <- vector("list")
  UNIVARIATE_SIGMA2_predictions <- vector("list")
  UNIVARIATE_NU2_predictions <- vector("list")
  
  
  #### MU 1
  UNIVARIATE_MU1_predictions$x1 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = 1, type = "link" )
  UNIVARIATE_MU1_predictions$x2 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x2", type = "link" )
  UNIVARIATE_MU1_predictions$x3 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x3", type = "link" )
  UNIVARIATE_MU1_predictions$x4 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x4", type = "link" )
  UNIVARIATE_MU1_predictions$x5 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x5", type = "link" )
  UNIVARIATE_MU1_predictions$x6 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x6", type = "link" )
  UNIVARIATE_MU1_predictions$x7 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x7", type = "link" )
  UNIVARIATE_MU1_predictions$x8 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x8", type = "link" )
  UNIVARIATE_MU1_predictions$x9 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x9", type = "link" )
  UNIVARIATE_MU1_predictions$x10 <- predict(glm.uni.mu1, data.frame(newX), parameter = "mu", which = "x10", type = "link" )
  
  #### SIGMA 1
  UNIVARIATE_SIGMA1_predictions$x1 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = 1, type = "link" )
  UNIVARIATE_SIGMA1_predictions$x2 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x2", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x3 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x3", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x4 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x4", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x5 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x5", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x6 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x6", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x7 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x7", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x8 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x8", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x9 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x9", type = "link" )
  UNIVARIATE_SIGMA1_predictions$x10 <- predict(glm.uni.mu1, data.frame(newX), parameter = "sigma", which = "x10", type = "link" )
  
  
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
  
  
  #### NU 2
  UNIVARIATE_NU2_predictions$x1 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = 1, type = "link" )
  UNIVARIATE_NU2_predictions$x2 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x2", type = "link" )
  UNIVARIATE_NU2_predictions$x3 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x3", type = "link" )
  UNIVARIATE_NU2_predictions$x4 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x4", type = "link" )
  UNIVARIATE_NU2_predictions$x5 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x5", type = "link" )
  UNIVARIATE_NU2_predictions$x6 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x6", type = "link" )
  UNIVARIATE_NU2_predictions$x7 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x7", type = "link" )
  UNIVARIATE_NU2_predictions$x8 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x8", type = "link" )
  UNIVARIATE_NU2_predictions$x9 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x9", type = "link" )
  UNIVARIATE_NU2_predictions$x10 <- predict(glm.uni.mu2, data.frame(newX), parameter = "nu", which = "x10", type = "link" )
  
  
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
    sample_uni <- data.gen.bivdisc_energyscore(n = 1000, 
                                               mu1 = pred.mu1.uni[i], 
                                               mu2 = pred.mu2.uni[i], 
                                               sigma1 = pred.sigma1.uni[i], 
                                               sigma2 = pred.sigma2.uni[i],
                                               nu = check_extremes_nuZINBI(pred.nu2.uni[i]),
                                               theta = sqrt(.Machine$double.eps))
    
    pred_sample_uni[1,] <- sample_uni[,1]
    pred_sample_uni[2,] <- sample_uni[,2]
    
    es_uni[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_uni) 
    
    
    # copula
    sample_cop <- data.gen.bivdisc_energyscore(n = 1000, 
                                               mu1 = predCopula.mu1[i],  
                                               mu2 = predCopula.mu2[i], 
                                               sigma1 = predCopula.sigma1[i], 
                                               sigma2 = predCopula.sigma2[i],
                                               nu = check_extremes_nuZINBI(predCopula.nu2[i]),
                                               theta = predCopula.rho[i])
    
    pred_sample_cop[1, ] <- sample_cop[,1]
    pred_sample_cop[2, ] <- sample_cop[,2]
    
    es_cop[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_cop) 
    
    
  }
  
  energy_score <- list()
  #energy_score$biv <- mean(es_biv)
  energy_score$uni <- mean(es_uni)
  energy_score$cop <- mean(es_cop)
  
  
  
  ########################################### We gather: mstop, predictions, MSE margin 1, margin 2, loglikelihood, energy score.
  output <- list(TrueBeta = TrueBeta, 
                 n = n.ges, 
                 p  = p, 
                 Likelihood = lik,
                 predict.uni = pred.uni, 
                 predict.cop = pred.cop,
                 energy_score = energy_score,  
                 ####
                 Coefficients.uni = coef.uni, 
                 mstop.uni = mstop.uni,
                 ###
                 DiscreteMarginMetrics = DiscreteMargins_Metrics,
                 ###
                 Copula_Predictions = list(MU1 = Copula_MU1_predictions, 
                                           SIGMA1 = Copula_SIGMA1_predictions,
                                           MU2 = Copula_MU2_predictions,
                                           SIGMA2 = Copula_SIGMA2_predictions,
                                           NU2 = Copula_NU2_predictions,
                                           RHO = Copula_RHO_predictions),
                 ###
                 Univariate_Predictions = list(MU1 = UNIVARIATE_MU1_predictions, 
                                               SIGMA1 = UNIVARIATE_SIGMA1_predictions,
                                               MU2 = UNIVARIATE_MU2_predictions, 
                                               SIGMA2 = UNIVARIATE_SIGMA2_predictions,
                                               NU2 = UNIVARIATE_NU2_predictions),
                 ###
                 CoefficientsCOPULA = coef.bivDiscCopula, 
                 oobag.riskCOPULA = oobag.risk.cop, 
                 energy_scoreCOPULA = energy_score,  
                 mstopCOPULA = mstop.bivCopula,
                 TrueKendallRange = TrueKendallTauRange
  )
  
  
  return(output)
  
}
