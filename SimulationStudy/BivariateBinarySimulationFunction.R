#########################################
library("VGAM")
library("gamboostLSS")
library("mvtnorm")
library("pROC")
library("scoringRules")
library("GJRM")
#setwd("BOOSTCOPFILES/")
#########################################
### load Copulas
source("Copulas/BivariateBinary/Copula_Gaussian_BivBin.R")
#source("Copulas/BivariateBinary/Copula_Gaussian2_BivBin.R")
source("Copulas/BivariateBinary/Copula_Clayton_BivBin.R")
source("Copulas/BivariateBinary/Copula_Gumbel_BivBin.R")
#source("Copulas/BivariateBinary/Copula_Frank_BivBin.R")
source("Copulas/BivariateBinary/Copula_Joe_BivBin.R")
source("Copulas/BivariateBinary/Copula_Joe90_BivBin.R")

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



sim <- function(seed, n.train, n.mstop, p, n.test, corr, boost.nu = 0.1){
  
  
  loss <- function(mu1, mu2, rho, y){
    
    cdf1 <- mu1
    
    cdf2 <- mu2
    
    rho <- check_extremes(rho)
    
    #p11 <-  pdffz(VGAM:::pbinorm(qnorm(cdf1), qnorm(cdf2), cov12 = rho))
    p11 <-  pdffz(VineCopula::BiCopCDF(family = 1, u1 = cdf1, u2 = cdf2, par = rho))
    
    return(- ( y[,1]*y[,2] * log( p11 ) +
                 y[,1]*(1-y[,2]) * log( pdffz((cdf1 - p11)) ) +
                 (1-y[,1])*y[,2] * log( pdffz((cdf2 - p11)) ) +
                 (1-y[,1])*(1-y[,2]) * log( pdffz((1 - cdf1 - cdf2 + p11)) ) 
    ) )
    
  }
  
  data.gen.bivbin <- function(FAM, mu1, mu2, theta){
    
    
    u1u2 <- VineCopula::BiCopSim(length(mu1), family = FAM, par = theta)
    
    y1 <- qbinom(p = u1u2[,1], size = 1, prob = mu1)
    
    y2 <- qbinom(p = u1u2[,2], size = 1, prob = mu2)
    
    
    dat <- data.frame(y1, y2)
    
    return(dat)
  }
  
  data.gen.bivbin_energyscore <- function(n, FAM, mu1, mu2, theta){
    
    
    u1u2 <- VineCopula::BiCopSim(n, family = FAM, par = theta)
    
    y1 <- qbinom(p = u1u2[,1], size = 1, prob = mu1)
    
    y2 <- qbinom(p = u1u2[,2], size = 1, prob = mu2)
    
    
    dat <- data.frame(y1, y2)
    
    return(dat)
  }
  
  
  
  ## mu 1
  mu1 <- c(0, -1, +0.5, 1, 0, -0.5, rep(0,p-6))
  
  ### mu 2
  mu2 <- c(0.5, -1, 0.75, 0, 0, rep(0,p-5)) 
  
  # copula parameter
  or  <- c(-0.5, +0.5, -1.5, 1.5, 0, 0, 0, 0,rep(0,p-8))
  
  
  TrueBeta <-  vector('list')
  TrueBeta$mu1 <- mu1
  TrueBeta$mu2 <- mu2
  TrueBeta$or <- or
  
  set.seed(seed)
  
  n <- n.train + n.mstop
  weight.mstop <- c(rep(1, times = n.train),rep(0, times = n.mstop)) 
  
  
  ############################################# ############################################# #############################################
  # - training covariates
  x.train <- rmvnorm(n = n, mean = rep(0,p), sigma =  toeplitz(sapply(seq(0,p-1), function(x) corr^x)))
  
  # - test covariates
  x.test <- rmvnorm(n = n.test, mean = rep(0,p), sigma = toeplitz(sapply(seq(0,p-1), function(x) corr^x)))
  
  # x.train <- apply(x.train, 2, pnorm)
  # x.test <- apply(x.test, 2, pnorm)
   
  # GENERATION OF DISTRIBUTION PARAMETERS AND STUFF 
  # predictor
  train.eta.mu1_center <-  (x.train %*% mu1)
  
  
  ### MU1
  train.eta.mu1<-  pnorm(train.eta.mu1_center) 
  
  test.eta.mu1 <-  pnorm(x.test %*% mu1)
  
  ### MU2
  train.eta.mu2 <-  x.train %*% mu2
  
  train.eta.mu2 <-  1-exp(-exp(train.eta.mu2)) 
  
  test.eta.mu2 <-  1-exp(-exp(x.test %*% mu2)) 
  
  ### Copula parameter
  train.eta.or_center <-   (x.train %*% or)
  
  train.copula.parameter <- tanh(train.eta.or_center)
  
  test.copula.parameter <- tanh( x.test %*% or)
  
  
  TrueKendallTauRange <- VineCopula::BiCopPar2Tau(family = 1, par = check_extremes(range(train.copula.parameter)))
  range(train.eta.mu1)
  range(train.eta.mu2)
  range(train.copula.parameter)

  plot(train.eta.mu1)
  plot(train.eta.mu2)
  plot(train.copula.parameter)

  
  ############################################# ############################################# #############################################
  ############################################# ############################################# #############################################
  
  # SAMPLING OF BIVARIATE RESPONSE FROM COPULA
  y.train <- data.gen.bivbin(FAM = 1, mu1 = train.eta.mu1, mu2 = train.eta.mu2, theta = train.copula.parameter)
  
  colnames(y.train)<- c('y1','y2')
  dat.train <- data.frame(y.train,x.train)
  
  y.test <- data.gen.bivbin(FAM = 1, mu1 = test.eta.mu1, mu2 = test.eta.mu2, theta = test.copula.parameter)
  
  colnames(y.test)<- c('y1','y2')
  dat.test <- data.frame(y.test,x.test)
  
  ############################################# ############################################# #############################################
  ############################################# ############################################# #############################################
  
  
  # FITTIIIINGG
  mstop.bivBern <-  vector('list')
  
  AUCbivBern <-  vector('list')
  BrierbivBern <-  vector('list')
  lik <- vector('list')
  
  
  #############################################
  ################# - COPULA - ################
  #############################################
  # Bivariate copula model
  bivBinCopula <- glmboostLSS(cbind(y1,y2)~., data = dat.train, 
                              families = Gauss_Cop_BivBinary(marg1 = "PROBIT", 
                                                             marg2 = "CLOGLOG", 
                                                             stabilization = "L2"), 
                              control = boost_control(mstop = 1500, 
                                                      risk = 'oobag', 
                                                      nu = boost.nu, 
                                                      trace = TRUE),
                              method = 'noncyclic', 
                              weights = weight.mstop)
  
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 1490){
    bivBinCopula[3000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 2990){
    bivBinCopula[5000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 4990){
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
  oobag.risk.cop <- risk(bivBinCopula,merge = T)
  
  # rm(bivBinCopula)
  # dat.train_biv <- dat.train[weight.mstop == 1, ]
  
  # # RE-FIT until OPTIMAL MSTOP
  # bivBinCopula <- glmboostLSS(cbind(y1,y2)~., data = dat.train_biv, families = Gauss_Cop_BivBinary(marg1 = "PROBIT", marg2 = "CLOGLOG"), 
  #                             control = boost_control(mstop = MSTOP_COP, nu  = 0.05), method = 'noncyclic')
  bivBinCopula <- bivBinCopula[MSTOP_COP]
  
  
  mstop.bivBinCopula <-  vector('list')
  mstop.bivBinCopula$mstop <- MSTOP_COP
  
  mstop.bivBinCopula$mu1 <- bivBinCopula$mu1$mstop()
  mstop.bivBinCopula$mu2 <- bivBinCopula$mu2$mstop()
  
  mstop.bivBinCopula$rho  <- bivBinCopula$rho$mstop()
  
  # EXTRACT ALL COEFFICIENTS!
  coef.bivBinCopula <- coef(bivBinCopula, which = "")
  
  AUCbivBinCopula <-  vector('list')
  BrierbivBinCopula <-  vector('list')
  likBinCopula <- vector('list')
  
  
  # Predict distributional quantities
  predCopula.mu1 <- predict(bivBinCopula$mu1, newdata = dat.test, type = 'response')
  predCopula.mu2 <- predict(bivBinCopula$mu2, newdata = dat.test, type = 'response')
  predCopula.rho <- predict(bivBinCopula$rho, newdata = dat.test, type = 'response')
  
  ######### COMPUTE PERFORMANCE: 
  # AUC (margin 1, then margin 2)
  AUCbivBinCopula$mu1 <- roc(dat.test$y1, as.numeric(predCopula.mu1))$auc
  AUCbivBinCopula$mu2 <- roc(dat.test$y2, as.numeric(predCopula.mu2))$auc
  
  # Brier score (margin 1, then margin 2)
  BrierbivBinCopula$mu1 <- mean((as.numeric(predCopula.mu1)-(as.numeric(dat.test$y1)))^2)
  BrierbivBinCopula$mu2 <- mean((as.numeric(predCopula.mu2)-(as.numeric(dat.test$y2)))^2)
  
  # LOG-SCORE (LOG-LIKELIHOOD) (change to the copula likelihood)
  lik$biv <- sum(loss(mu1 = predCopula.mu1, mu2 = predCopula.mu2, rho = predCopula.rho, y = y.test))
  
  
  ###########################################################################################################################
  ################ - univariat - ############################################################################################
  ###########################################################################################################################
  
  mstop.uni <- vector('list')
  coef.uni  <- vector('list')
  coef.uni_1  <- vector('list')
  AUC.uni   <- vector('list')
  Brier.uni <- vector('list')
  
  # - mu1
  dat.train.mu1 = dat.train[, !(names(dat.train) %in% "y2")]
  dat.test.mu1  = dat.test[ , !(names(dat.test) %in% "y2")]
  
  glm.uni.mu1 <- glmboost(as.factor(y1) ~. , 
                          data = dat.train.mu1, 
                          family = Binomial(link = "probit", type = "glm"), 
                          control = boost_control(risk = 'oobag', nu = boost.nu, mstop = 1000, trace = T), 
                          weights = weight.mstop)
  
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  if(mstop.uni$mu1 >= 990){
    glm.uni.mu1[3000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 2990){
    glm.uni.mu1[5000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 4990){
    glm.uni.mu1[8000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 7990){
    glm.uni.mu1[10000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 9990){
    glm.uni.mu1[15000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 14990){
    glm.uni.mu1[20000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  
  # dat.train_mu1 <- dat.train.mu1[weight.mstop == 1, ]
  # glm.uni.mu1 <- glmboost(as.factor(y1) ~. , 
  #                         data = dat.train_mu1, 
  #                         family = Binomial(link = "cloglog", type = "glm"), 
  #                         control = boost_control(mstop = mstop.uni$mu1, nu  = boost.nu))
  glm.uni.mu1 <- glm.uni.mu1[mstop.uni$mu1]
  
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 1  
  coef.uni$mu1 <- coef(glm.uni.mu1, which = '')
  
  AUC.uni$mu1 <- roc(dat.test.mu1$y1, as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response")))$auc
  Brier.uni$mu1 <- mean((as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response"))-(as.numeric(dat.test.mu1$y1)))^2)
  
  # - mu2
  dat.train.mu2 =  dat.train[, !(names(dat.train) %in% "y1")]
  dat.test.mu2 = dat.test[, !(names(dat.test) %in% "y1")]
  
  glm.uni.mu2 <- glmboost(as.factor(y2) ~. , 
                          data = dat.train.mu2, 
                          family = Binomial(link = "cloglog", type = "glm"), 
                          control = boost_control(risk = 'oobag', nu = boost.nu, mstop = 1000, trace = T), 
                          weights = weight.mstop)
  
  
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  if(mstop.uni$mu2 >= 990){
    glm.uni.mu2[3000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 2990){
    glm.uni.mu2[5000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 4990){
    glm.uni.mu2[8000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 7990){
    glm.uni.mu2[10000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 9990){
    glm.uni.mu2[15000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  
  if(mstop.uni$mu2 >= 14990){
    glm.uni.mu2[20000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  
  
  #   dat.train_mu2 <- dat.train.mu2[weight.mstop == 1, ]
  # 
  # glm.uni.mu2 <- glmboost(as.factor(y2) ~., data = dat.train_mu2, family = Binomial(type = "glm"), 
  #                         control = boost_control(mstop = mstop.uni$mu1, nu  = 0.05))
  glm.uni.mu2 <- glm.uni.mu2[mstop.uni$mu2]
  
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 2  
  coef.uni$mu2 <- coef(glm.uni.mu2,which = '')
  
  AUC.uni$mu2 <- roc(dat.test.mu2$y2, as.numeric(predict(glm.uni.mu2, newdata = dat.test.mu2,type = "response")))$auc
  Brier.uni$mu2 <- mean((as.numeric(predict(glm.uni.mu2, newdata = dat.test.mu2,type = "response"))-(as.numeric(dat.test.mu2$y2)))^2)
  
  # Likelihood
  pred.mu1.uni <- as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response"))
  pred.mu2.uni <- as.numeric(predict(glm.uni.mu2, newdata = dat.test.mu2,type = "response"))
  
  mu1.uni.loglik <- sum(-dbinom(x = dat.test.mu1$y1, size = 1, prob = pred.mu1.uni, log = T))
  mu2.uni.loglik <- sum(-dbinom(x = dat.test.mu2$y2, size = 1, prob = pred.mu2.uni, log = T))
  
  lik$uni<- mu1.uni.loglik + mu2.uni.loglik
  
  
  
  n.ges = c(n.train, n.mstop, n.test)
  pred.ges <- NULL#list(pred.mu1, pred.mu2, pred.or)
  pred.uni <- list(pred.mu1.uni, pred.mu2.uni)
  pred.cop <- list(predCopula.mu1, predCopula.mu2, predCopula.rho)
  
  
  # ###########################################################################################################################
  # ################ - GJRM -      ############################################################################################
  # ###########################################################################################################################
  # 
  # 
  # gjrm_forms <- list(formula(c(" y1 ~", paste(names(dat.test[,-c(1,2)]), collapse = "+"))),
  #                    formula(c(" y2 ~", paste(names(dat.test[,-c(1,2)]), collapse = "+"))),
  #                    formula(c("  ~", paste(names(dat.test[,-c(1,2)]), collapse = "+"))))
  # 
  # # Bivariate copula model in GJRM
  # bivBinGJRM_catch <- myTryCatch( 
  #   gjrm(gjrm_forms, 
  #        data = dat.train, 
  #        BivD = "N", 
  #        Model = "B",
  #        margins = c("probit", "cloglog"))
  # )
  # 
  # if(is.null(bivBinGJRM_catch$value)){
  #   
  #   # EXTRACT ALL COEFFICIENTS!
  #   coef.GJRM <- NULL
  #   AUC.GJRM <-  vector('list')
  #   Brier.GJRM <-  vector('list')
  #   
  #   
  #   # Predict distributional quantities
  #   predGJRM.mu1 <- NULL
  #   predGJRM.mu2 <- NULL
  #   predGJRM.rho <- NULL
  #   
  #   lik$GJRM <- NULL
  #   
  #   
  #   pred.GJRM <- list(predGJRM.mu1, predGJRM.mu2, predGJRM.rho)
  #   
  #   
  #   
  # }else{
  #   
  #   bivBinGJRM <- bivBinGJRM_catch$value
  #   
  #   # EXTRACT ALL COEFFICIENTS!
  #   coef.GJRM <- coef(bivBinGJRM)
  #   
  #   AUC.GJRM <-  vector('list')
  #   Brier.GJRM <-  vector('list')
  #   
  #   
  #   # Predict distributional quantities
  #   predGJRM.mu1 <- predict(bivBinGJRM, eq = 1, newdata = dat.test, type = 'response')
  #   predGJRM.mu2 <- predict(bivBinGJRM, eq = 2, newdata = dat.test, type = 'response')
  #   predGJRM.rho <- predict(bivBinGJRM, eq = 3, newdata = dat.test, type = 'response')
  #   predGJRM.rho <- tanh(predGJRM.rho)
  #   
  #   ######### COMPUTE PERFORMANCE: 
  #   # AUC (margin 1, then margin 2)
  #   AUC.GJRM$mu1 <- roc(dat.test$y1, as.numeric(predGJRM.mu1))$auc
  #   AUC.GJRM$mu2 <- roc(dat.test$y2, as.numeric(predGJRM.mu2))$auc
  #   
  #   # Brier score (margin 1, then margin 2)
  #   Brier.GJRM$mu1 <- mean((as.numeric(predGJRM.mu1)-(as.numeric(dat.test$y1)))^2)
  #   Brier.GJRM$mu2 <- mean((as.numeric(predGJRM.mu2)-(as.numeric(dat.test$y2)))^2)
  #   
  #   # LOG-SCORE (LOG-LIKELIHOOD) (change to the copula likelihood)
  #   lik$GJRM <- sum(loss(mu1 = predGJRM.mu1, mu2 = predGJRM.mu2, rho = predGJRM.rho, y = y.test))
  #   
  #   pred.GJRM <- list(predGJRM.mu1, predGJRM.mu2, predGJRM.rho)
  #   
  #   
  # }
  
  
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
    pred_sample_gjrm <- matrix(NA, nrow = 2, ncol = 1000)
    
    
    # univariate
    sample_uni <- data.gen.bivbin_energyscore(n = 1000, mu1 = pred.mu1.uni[i], mu2 = pred.mu2.uni[i], FAM = 1, theta = 0)
    
    pred_sample_uni[1,] <- sample_uni[,1]
    pred_sample_uni[2,] <- sample_uni[,2]
    
    es_uni[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_uni) 
    
    
    # copula
    sample_cop <- data.gen.bivbin_energyscore(n = 1000, mu1 = predCopula.mu1[i],  mu2 = predCopula.mu2[i], FAM = 1, theta = check_extremes(predCopula.rho[i]))
    
    pred_sample_cop[1, ] <- sample_cop[,1]
    pred_sample_cop[2, ] <- sample_cop[,2]
    
    es_cop[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_cop) 
    
    # if(is.null(bivBinGJRM_catch$value)){
    #   
    #   es_gjrm <- NULL
    #   
    # }else{
    # ### GJRM: 
    # sample_gjrm <- data.gen.bivbin_energyscore(n = 1000, mu1 = predGJRM.mu1[i],  mu2 = predGJRM.mu2[i], FAM = 1, theta = check_extremes(predGJRM.rho[i]))
    # 
    # pred_sample_gjrm[1, ] <- sample_gjrm[,1]
    # pred_sample_gjrm[2, ] <- sample_gjrm[,2]
    # 
    # es_gjrm[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_gjrm) 
    # 
    # }
  }
  
  energy_score <- list()
  energy_score$uni <- mean(es_uni)
  energy_score$cop <- mean(es_cop)
  
  # if(is.null(bivBinGJRM_catch$value)){
  #   
  #   energy_score$gjrm <- NA
  # 
  # }else{
  #   
  #   energy_score$gjrm <- mean(es_gjrm)
  #   
  # }
  
  
  return(list(TrueBeta = TrueBeta, 
              n = n.ges, 
              p  = p, 
              #Coefficients = coef.bivBern, 
              #oobag.risk = oobag.risk, 
              AUCbivBern = AUCbivBern, 
              BrierbivBern = BrierbivBern, 
              Likelihood = lik,
              predict = pred.ges, 
              predict.uni = pred.uni, 
              predict.cop = pred.cop,
              energy_score = energy_score,  
              mstop = mstop.bivBern, 
              ####
              Coefficients.uni = coef.uni, 
              AUC.uni = AUC.uni, 
              Brier.uni = Brier.uni, 
              mstop.uni = mstop.uni,
              ###
              CoefficientsCOPULA = coef.bivBinCopula, 
              oobag.riskCOPULA = oobag.risk.cop, 
              AUCCOPULA = AUCbivBinCopula, 
              BrierCOPULA = BrierbivBinCopula, 
              LikelihoodCOPULA = likBinCopula,
              energy_scoreCOPULA = energy_score,  
              mstopCOPULA = mstop.bivBinCopula, 
              KendallRange = TrueKendallTauRange
              ###
              # Coefficients = coef.GJRM, 
              # AUCbivBern = AUC.GJRM, 
              # BrierbivBern = Brier.GJRM, 
              # Likelihood = lik$GJRM,
              # energy_score = energy_score$gjrm 
  ))
  
}


sim_NONLINEAR <- function(seed, n.train, n.mstop, p, n.test, corr, boost.nu = 0.1){
  
  
  loss <- function(mu1, mu2, rho, y){
    
    cdf1 <- mu1
    
    cdf2 <- mu2
    
    rho <- check_extremes(rho)
    
    #p11 <-  pdffz(VGAM:::pbinorm(qnorm(cdf1), qnorm(cdf2), cov12 = rho))
    p11 <-  pdffz(VineCopula::BiCopCDF(family = 1, u1 = cdf1, u2 = cdf2, par = rho))
    
    return(- ( y[,1]*y[,2] * log( p11 ) +
                 y[,1]*(1-y[,2]) * log( pdffz((cdf1 - p11)) ) +
                 (1-y[,1])*y[,2] * log( pdffz((cdf2 - p11)) ) +
                 (1-y[,1])*(1-y[,2]) * log( pdffz((1 - cdf1 - cdf2 + p11)) ) 
    ) )
    
  }
  
  data.gen.bivbin <- function(FAM, mu1, mu2, theta){
    
    
    u1u2 <- VineCopula::BiCopSim(length(mu1), family = FAM, par = theta)
    
    y1 <- qbinom(p = u1u2[,1], size = 1, prob = mu1)
    
    y2 <- qbinom(p = u1u2[,2], size = 1, prob = mu2)
    
    
    dat <- data.frame(y1, y2)
    
    return(dat)
  }
  
  data.gen.bivbin_energyscore <- function(n, FAM, mu1, mu2, theta){
    
    
    u1u2 <- VineCopula::BiCopSim(n, family = FAM, par = theta)
    
    y1 <- qbinom(p = u1u2[,1], size = 1, prob = mu1)
    
    y2 <- qbinom(p = u1u2[,2], size = 1, prob = mu2)
    
    
    dat <- data.frame(y1, y2)
    
    return(dat)
  }
  
  
  
  ## mu 1
  mu1 <- c(0, 0, 0, 0, 0, 0, rep(0,p-6))
  
  ### mu 2
  mu2 <- c(0, 0, 0, 0, 0, rep(0,p-5)) 
  
  # copula parameter
  or  <- c(0, 0, 0, 0, 0, 0, 0, 0,rep(0,p-8))
  
  
  
  nl_function_mu1 <- function(x){
    
    
    nl <- sqrt(x)*x^2 - 2*cos(3*x)
    
    return(nl*0.5)
    
  }
  
  nl_function_mu2 <- function(x){
    
    func <- -40*(x^(3/2) - x^(4/3))
    
    func <- 0.4 * func
    
    return(func)
  }
  
  nl_function_rho <- function(x){
    
    #return(2*sin(x*4))
    return(2*sin(x*4))
  }
  
  
  TrueBeta <-  vector('list')
  TrueBeta$mu1 <- mu1
  TrueBeta$mu2 <- mu2
  TrueBeta$or <- or
  
  set.seed(seed)
  
  n <- n.train + n.mstop
  weight.mstop <- c(rep(1, times = n.train),rep(0, times = n.mstop)) 
  
  
  ############################################# ############################################# #############################################
  # - training covariates
  x.train <- rmvnorm(n = n, mean = rep(0,p), sigma =  toeplitz(sapply(seq(0,p-1), function(x) corr^x)))
  
  # - test covariates
  x.test <- rmvnorm(n = n.test, mean = rep(0,p), sigma = toeplitz(sapply(seq(0,p-1), function(x) corr^x)))
  
  x.train <- apply(x.train, 2, pnorm)
  x.test <- apply(x.test, 2, pnorm)
  
  # GENERATION OF DISTRIBUTION PARAMETERS AND STUFF 
  # predictor
  train.eta.mu1_center <-  (x.train %*% mu1) + nl_function_mu1(x.train[,2])
  
  ### MU1
  train.eta.mu1<-  pnorm(train.eta.mu1_center) 
  
  test.eta.mu1 <-  pnorm(x.test %*% mu1 + + nl_function_mu1(x.test[,2]))
  
  ### MU2
  train.eta.mu2 <-  x.train %*% mu2 + nl_function_mu2(x.train[,4])
  
  train.eta.mu2 <-  1-exp(-exp(train.eta.mu2)) 
  
  test.eta.mu2 <-  1-exp(-exp(x.test %*% mu2 + nl_function_mu2(x.test[,4]) )) 
  
  ### Copula parameter
  train.eta.or_center <-   (x.train %*% or) + nl_function_rho(x.train[,1])
  
  train.copula.parameter <- tanh(train.eta.or_center)
  
  test.copula.parameter <- tanh( x.test %*% or + nl_function_rho(x.test[,1]))
  
  
  TrueKendallTauRange <- VineCopula::BiCopPar2Tau(family = 1, par = check_extremes(range(train.copula.parameter)))
  range(train.eta.mu1)
  range(train.eta.mu2)
  range(train.copula.parameter)
  
  plot(train.eta.mu1)
  plot(train.eta.mu2)
  plot(train.copula.parameter)
  
  
  ############################################# ############################################# #############################################
  ############################################# ############################################# #############################################
  
  # SAMPLING OF BIVARIATE RESPONSE FROM COPULA
  y.train <- data.gen.bivbin(FAM = 1, mu1 = train.eta.mu1, mu2 = train.eta.mu2, theta = train.copula.parameter)
  
  colnames(y.train)<- c('y1','y2')
  dat.train <- data.frame(y.train,x.train)
  
  y.test <- data.gen.bivbin(FAM = 1, mu1 = test.eta.mu1, mu2 = test.eta.mu2, theta = test.copula.parameter)
  
  colnames(y.test)<- c('y1','y2')
  dat.test <- data.frame(y.test,x.test)
  
  ############################################# ############################################# #############################################
  ############################################# ############################################# #############################################
  
  
  # FITTIIIINGG
  mstop.bivBern <-  vector('list')
  
  AUCbivBern <-  vector('list')
  BrierbivBern <-  vector('list')
  lik <- vector('list')
  
  
  #############################################
  ################# - COPULA - ################
  #############################################
  # Bivariate copula model
  bivBinCopula <- gamboostLSS(cbind(y1,y2)~., data = dat.train, 
                              families = Gauss_Cop_BivBinary(marg1 = "PROBIT", 
                                                             marg2 = "CLOGLOG", 
                                                             stabilization = "L2"), 
                              control = boost_control(mstop = 1500, 
                                                      risk = 'oobag', 
                                                      nu = boost.nu, 
                                                      trace = TRUE),
                              method = 'noncyclic', 
                              weights = weight.mstop)
  
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 1490){
    bivBinCopula[3000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 2990){
    bivBinCopula[5000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 4990){
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
  oobag.risk.cop <- risk(bivBinCopula,merge = T)
  
  # rm(bivBinCopula)
  # dat.train_biv <- dat.train[weight.mstop == 1, ]
  
  # # RE-FIT until OPTIMAL MSTOP
  # bivBinCopula <- glmboostLSS(cbind(y1,y2)~., data = dat.train_biv, families = Gauss_Cop_BivBinary(marg1 = "PROBIT", marg2 = "CLOGLOG"), 
  #                             control = boost_control(mstop = MSTOP_COP, nu  = 0.05), method = 'noncyclic')
  bivBinCopula <- bivBinCopula[MSTOP_COP]
  
  
  mstop.bivBinCopula <-  vector('list')
  mstop.bivBinCopula$mstop <- MSTOP_COP
  
  mstop.bivBinCopula$mu1 <- bivBinCopula$mu1$mstop()
  mstop.bivBinCopula$mu2 <- bivBinCopula$mu2$mstop()
  
  mstop.bivBinCopula$rho  <- bivBinCopula$rho$mstop()
  
  # EXTRACT ALL COEFFICIENTS!
  coef.bivBinCopula <- coef(bivBinCopula, which = "")
  
  AUCbivBinCopula <-  vector('list')
  BrierbivBinCopula <-  vector('list')
  likBinCopula <- vector('list')
  
  
  # Predict distributional quantities
  predCopula.mu1 <- predict(bivBinCopula$mu1, newdata = dat.test, type = 'response')
  predCopula.mu2 <- predict(bivBinCopula$mu2, newdata = dat.test, type = 'response')
  predCopula.rho <- predict(bivBinCopula$rho, newdata = dat.test, type = 'response')
  
  ######### COMPUTE PERFORMANCE: 
  # AUC (margin 1, then margin 2)
  AUCbivBinCopula$mu1 <- roc(dat.test$y1, as.numeric(predCopula.mu1))$auc
  AUCbivBinCopula$mu2 <- roc(dat.test$y2, as.numeric(predCopula.mu2))$auc
  
  # Brier score (margin 1, then margin 2)
  BrierbivBinCopula$mu1 <- mean((as.numeric(predCopula.mu1)-(as.numeric(dat.test$y1)))^2)
  BrierbivBinCopula$mu2 <- mean((as.numeric(predCopula.mu2)-(as.numeric(dat.test$y2)))^2)
  
  # LOG-SCORE (LOG-LIKELIHOOD) (change to the copula likelihood)
  lik$biv <- sum(loss(mu1 = predCopula.mu1, mu2 = predCopula.mu2, rho = predCopula.rho, y = y.test))
  
  
  ### Copula predictions:
  ################ Get predictions from each covariate?
  newX <- matrix(rep(seq(from = 0, to = 1, length.out = 500), times = p), nrow = 500)
  colnames(newX) <-  colnames(x.train)
  
  Copula_MU1_predictions <- vector("list", length = p)
  Copula_MU2_predictions <- vector("list", length = p)
  Copula_RHO_predictions <- vector("list", length = p)
  
  for(i in 1:p){
  
  Copula_MU1_predictions[[i]] <- predict(bivBinCopula, data.frame(newX), parameter = "mu1", which = i, type = "link" )
  
  Copula_MU2_predictions[[i]] <- predict(bivBinCopula, data.frame(newX), parameter = "mu2", which = i, type = "link" )
  
  Copula_RHO_predictions[[i]] <- predict(bivBinCopula, data.frame(newX), parameter = "rho", which = i, type = "link" )
  
  }
  
  ###########################################################################################################################
  ################ - univariat - ############################################################################################
  ###########################################################################################################################
  
  mstop.uni <- vector('list')
  coef.uni  <- vector('list')
  coef.uni_1  <- vector('list')
  AUC.uni   <- vector('list')
  Brier.uni <- vector('list')
  
  # - mu1
  dat.train.mu1 = dat.train[, !(names(dat.train) %in% "y2")]
  dat.test.mu1  = dat.test[ , !(names(dat.test) %in% "y2")]
  
  glm.uni.mu1 <- gamboost(as.factor(y1) ~. , 
                          data = dat.train.mu1, 
                          family = Binomial(link = "probit", type = "glm"), 
                          control = boost_control(risk = 'oobag', nu = boost.nu, mstop = 1000, trace = T), 
                          weights = weight.mstop)
  
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  if(mstop.uni$mu1 >= 990){
    glm.uni.mu1[3000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 2990){
    glm.uni.mu1[5000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 4990){
    glm.uni.mu1[8000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 7990){
    glm.uni.mu1[10000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 9990){
    glm.uni.mu1[15000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 14990){
    glm.uni.mu1[20000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  
  # dat.train_mu1 <- dat.train.mu1[weight.mstop == 1, ]
  # glm.uni.mu1 <- glmboost(as.factor(y1) ~. , 
  #                         data = dat.train_mu1, 
  #                         family = Binomial(link = "cloglog", type = "glm"), 
  #                         control = boost_control(mstop = mstop.uni$mu1, nu  = boost.nu))
  glm.uni.mu1 <- glm.uni.mu1[mstop.uni$mu1]
  
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 1  
  coef.uni$mu1 <- coef(glm.uni.mu1, which = '')
  
  AUC.uni$mu1 <- roc(dat.test.mu1$y1, as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response")))$auc
  Brier.uni$mu1 <- mean((as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response"))-(as.numeric(dat.test.mu1$y1)))^2)
  
  # - mu2
  dat.train.mu2 =  dat.train[, !(names(dat.train) %in% "y1")]
  dat.test.mu2 = dat.test[, !(names(dat.test) %in% "y1")]
  
  glm.uni.mu2 <- glmboost(as.factor(y2) ~. , 
                          data = dat.train.mu2, 
                          family = Binomial(link = "cloglog", type = "glm"), 
                          control = boost_control(risk = 'oobag', nu = boost.nu, mstop = 1000, trace = T), 
                          weights = weight.mstop)
  
  
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  if(mstop.uni$mu2 >= 990){
    glm.uni.mu2[3000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 2990){
    glm.uni.mu2[5000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 4990){
    glm.uni.mu2[8000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 7990){
    glm.uni.mu2[10000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 9990){
    glm.uni.mu2[15000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  
  if(mstop.uni$mu2 >= 14990){
    glm.uni.mu2[20000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  
  
  #   dat.train_mu2 <- dat.train.mu2[weight.mstop == 1, ]
  # 
  # glm.uni.mu2 <- glmboost(as.factor(y2) ~., data = dat.train_mu2, family = Binomial(type = "glm"), 
  #                         control = boost_control(mstop = mstop.uni$mu1, nu  = 0.05))
  glm.uni.mu2 <- glm.uni.mu2[mstop.uni$mu2]
  
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 2  
  coef.uni$mu2 <- coef(glm.uni.mu2,which = '')
  
  AUC.uni$mu2 <- roc(dat.test.mu2$y2, as.numeric(predict(glm.uni.mu2, newdata = dat.test.mu2,type = "response")))$auc
  Brier.uni$mu2 <- mean((as.numeric(predict(glm.uni.mu2, newdata = dat.test.mu2,type = "response"))-(as.numeric(dat.test.mu2$y2)))^2)
  
  # Likelihood
  pred.mu1.uni <- as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response"))
  pred.mu2.uni <- as.numeric(predict(glm.uni.mu2, newdata = dat.test.mu2,type = "response"))
  
  mu1.uni.loglik <- sum(-dbinom(x = dat.test.mu1$y1, size = 1, prob = pred.mu1.uni, log = T))
  mu2.uni.loglik <- sum(-dbinom(x = dat.test.mu2$y2, size = 1, prob = pred.mu2.uni, log = T))
  
  lik$uni<- mu1.uni.loglik + mu2.uni.loglik
  
  
  
  n.ges = c(n.train, n.mstop, n.test)
  pred.ges <- NULL#list(pred.mu1, pred.mu2, pred.or)
  pred.uni <- list(pred.mu1.uni, pred.mu2.uni)
  pred.cop <- list(predCopula.mu1, predCopula.mu2, predCopula.rho)
  
  #### Univariate predictions:
  Univariate_MU1_predictions <- vector("list", length = p)
  Univariate_MU2_predictions <- vector("list", length = p)
  
  for(i in 1:p){
    
    Univariate_MU1_predictions[[i]] <- predict(glm.uni.mu1, data.frame(newX), which = i, type = "link" )
    
    Univariate_MU2_predictions[[i]] <- predict(glm.uni.mu2, data.frame(newX), which = i, type = "link" )
    
    
  }
  
  
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
    pred_sample_gjrm <- matrix(NA, nrow = 2, ncol = 1000)
    
    
    # univariate
    sample_uni <- data.gen.bivbin_energyscore(n = 1000, mu1 = pred.mu1.uni[i], mu2 = pred.mu2.uni[i], FAM = 1, theta = 0)
    
    pred_sample_uni[1,] <- sample_uni[,1]
    pred_sample_uni[2,] <- sample_uni[,2]
    
    es_uni[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_uni) 
    
    
    # copula
    sample_cop <- data.gen.bivbin_energyscore(n = 1000, mu1 = predCopula.mu1[i],  mu2 = predCopula.mu2[i], FAM = 1, theta = check_extremes(predCopula.rho[i]))
    
    pred_sample_cop[1, ] <- sample_cop[,1]
    pred_sample_cop[2, ] <- sample_cop[,2]
    
    es_cop[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_cop) 
    
   
  }
  
  energy_score <- list()
  energy_score$uni <- mean(es_uni)
  energy_score$cop <- mean(es_cop)
  
  
  
  return(list(TrueBeta = TrueBeta, 
              n = n.ges, 
              p  = p, 
              AUCbivBern = AUCbivBern, 
              BrierbivBern = BrierbivBern, 
              Likelihood = lik,
              predict = pred.ges, 
              predict.uni = pred.uni, 
              predict.cop = pred.cop,
              energy_score = energy_score,  
              mstop = mstop.bivBern, 
              ####
              Coefficients.uni = coef.uni, 
              AUC.uni = AUC.uni, 
              Brier.uni = Brier.uni, 
              mstop.uni = mstop.uni,
              ###
              CoefficientsCOPULA = coef.bivBinCopula, 
              oobag.riskCOPULA = oobag.risk.cop, 
              AUCCOPULA = AUCbivBinCopula, 
              BrierCOPULA = BrierbivBinCopula, 
              LikelihoodCOPULA = likBinCopula,
              energy_scoreCOPULA = energy_score,  
              mstopCOPULA = mstop.bivBinCopula, 
              KendallRange = TrueKendallTauRange,
              ###
              CopulaPredictions = list(MU1 = Copula_MU1_predictions, 
                                       MU2 = Copula_MU2_predictions,
                                       RHO = Copula_RHO_predictions), 
              UnivariatePredictions = list(MU1 = Univariate_MU1_predictions,
                                           MU2 = Univariate_MU2_predictions)
  ))
  
}

sim_GUMBEL <- function(seed, n.train, n.mstop, p, n.test, corr, boost.nu = 0.1){
  
  
  loss <- function(mu1, mu2, rho, y){
    
    cdf1 <- mu1
    
    cdf2 <- mu2
    
    rho <- rho
    
    p11 <-  pdffz(VineCopula::BiCopCDF(u1 = cdf1, u2 = cdf2, family = 4, par = rho))
    
    
    return(- ( y[,1]*y[,2] * log( p11 ) +
                 y[,1]*(1-y[,2]) * log( pdffz((cdf1 - p11)) ) +
                 (1-y[,1])*y[,2] * log( pdffz((cdf2 - p11)) ) +
                 (1-y[,1])*(1-y[,2]) * log( pdffz((1 - cdf1 - cdf2 + p11)) ) 
    ) )
    
  }
  
  data.gen.bivbin <- function(FAM, mu1, mu2, theta){
    
    
    u1u2 <- VineCopula::BiCopSim(length(mu1), family = FAM, par = theta)
    
    y1 <- qbinom(p = u1u2[,1], size = 1, prob = mu1)
    
    y2 <- qbinom(p = u1u2[,2], size = 1, prob = mu2)
    
    
    dat <- data.frame(y1, y2)
    
    return(dat)
  }
  
  data.gen.bivbin_energyscore <- function(n, FAM, mu1, mu2, theta){
    
    
    u1u2 <- VineCopula::BiCopSim(n, family = FAM, par = theta)
    
    y1 <- qbinom(p = u1u2[,1], size = 1, prob = mu1)
    
    y2 <- qbinom(p = u1u2[,2], size = 1, prob = mu2)
    
    
    dat <- data.frame(y1, y2)
    
    return(dat)
  }
  
  
  
  # ## mu 1
  # mu1 <- c(0, -1, +0.5, 1, 0, -0.5, rep(0,p-6))
  mu1 <- c(0, -1, +0.5, 1, 0, -0.5, rep(0,p-6))
  # 
  # ### mu 2
  # mu2 <- c(-0.5, -1, 1.5, 0, 0, rep(0,p-5)) 
  mu2 <- c(0.5, -1, 0.75, 0, 0, rep(0,p-5)) 
  # 
  # # copula parameter
  # or  <- c(0, -2, +0.5, +2.5, 0, 0, 0, 0, rep(0,p-8))
  or  <- c(-0.5, +0.5, -1.5, 1.5, 0, 0, 0, 0,rep(0,p-8))
  or  <- c(-0.5, +1.25, -1, 2, 0, 0, 0, 0,rep(0,p-8))
  
  
  TrueBeta <-  vector('list')
  TrueBeta$mu1 <- mu1
  TrueBeta$mu2 <- mu2
  TrueBeta$or <- or
  
  set.seed(seed)
  
  n <- n.train + n.mstop
  weight.mstop <- c(rep(1, times = n.train),rep(0, times = n.mstop)) 
  
  
  ############################################# ############################################# #############################################
  # - training covariates
  x.train <- rmvnorm(n = n, mean = rep(0,p), sigma =  toeplitz(sapply(seq(0,p-1), function(x) corr^x)))
  
  x.train <- apply(x.train, 2, pnorm)
  
  # - test covariates
  x.test <- rmvnorm(n = n.test, mean = rep(0,p), sigma = toeplitz(sapply(seq(0,p-1), function(x) corr^x)))

  x.test <- apply(x.test, 2, pnorm)

    
  # GENERATION OF DISTRIBUTION PARAMETERS AND STUFF 
  # predictor
  train.eta.mu1_center <-  (x.train %*% mu1)
  
  
  ### MU1
  train.eta.mu1<-  pnorm(train.eta.mu1_center) 
  
  test.eta.mu1 <-  pnorm(x.test %*% mu1)
  
  ### MU2
  train.eta.mu2 <-  x.train %*% mu2
  
  train.eta.mu2 <-  1-exp(-exp(train.eta.mu2)) 
  
  test.eta.mu2 <-  1-exp(-exp(x.test %*% mu2)) 
  
  ### Copula parameter
  train.eta.or_center <-   (x.train %*% or)
  
  train.copula.parameter <- exp(train.eta.or_center) + 1
  
  test.copula.parameter <- exp( x.test %*% or) + 1
  
  
  TrueKendallTauRange <- VineCopula::BiCopPar2Tau(family = 4, par = range(train.copula.parameter))
  #TrueKendallTauRange
  
  #summary(VineCopula::BiCopPar2Tau(family = 4, par = train.copula.parameter))
  # range(train.eta.mu1)
  # range(train.eta.mu2)
  # range(train.copula.parameter)
  # 
  # plot(train.eta.mu1)
  # plot(train.eta.mu2)
  # plot(train.copula.parameter)
  # 
  
  ############################################# ############################################# #############################################
  ############################################# ############################################# #############################################
  
  
  
  y.train <- data.gen.bivbin(FAM = 4, mu1 = train.eta.mu1, mu2 = train.eta.mu2, theta = train.copula.parameter)
  
  colnames(y.train)<- c('y1','y2')
  dat.train <- data.frame(y.train,x.train)
  
  y.test <- data.gen.bivbin(FAM = 4, mu1 = test.eta.mu1, mu2 = test.eta.mu2, theta = test.copula.parameter)
  
  colnames(y.test)<- c('y1','y2')
  dat.test <- data.frame(y.test,x.test)
  
  
  
  # FITTIIIINGG
  mstop.bivBern <-  vector('list')
  
  AUCbivBern <-  vector('list')
  BrierbivBern <-  vector('list')
  lik <- vector('list')
  
  
  #############################################
  ################# - COPULA - ################
  #############################################
  # Bivariate copula model
  bivBinCopula <- glmboostLSS(cbind(y1,y2)~., data = dat.train, 
                              families = Gumbel_Cop_BivBinary(marg1 = "PROBIT", 
                                                             marg2 = "CLOGLOG", 
                                                             stabilization = "L2"), 
                              control = boost_control(mstop = 1000, 
                                                      risk = 'oobag', 
                                                      nu = boost.nu, 
                                                      trace = TRUE),
                              method = 'noncyclic', 
                              weights = weight.mstop)
  
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 990){
    bivBinCopula[3000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 2990){
    bivBinCopula[5000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 4990){
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
  oobag.risk.cop <- risk(bivBinCopula,merge = T)
  
  # rm(bivBinCopula)
  # dat.train_biv <- dat.train[weight.mstop == 1, ]
  
  # # RE-FIT until OPTIMAL MSTOP
  # bivBinCopula <- glmboostLSS(cbind(y1,y2)~., data = dat.train_biv, families = Gauss_Cop_BivBinary(marg1 = "PROBIT", marg2 = "CLOGLOG"), 
  #                             control = boost_control(mstop = MSTOP_COP, nu  = 0.05), method = 'noncyclic')
  bivBinCopula <- bivBinCopula[MSTOP_COP]
  
  
  mstop.bivBinCopula <-  vector('list')
  mstop.bivBinCopula$mstop <- MSTOP_COP
  
  mstop.bivBinCopula$mu1 <- bivBinCopula$mu1$mstop()
  mstop.bivBinCopula$mu2 <- bivBinCopula$mu2$mstop()
  
  mstop.bivBinCopula$rho  <- bivBinCopula$rho$mstop()
  
  # EXTRACT ALL COEFFICIENTS!
  coef.bivBinCopula <- coef(bivBinCopula, which = "")
  
  AUCbivBinCopula <-  vector('list')
  BrierbivBinCopula <-  vector('list')
  likBinCopula <- vector('list')
  
  
  # Predict distributional quantities
  predCopula.mu1 <- predict(bivBinCopula$mu1, newdata = dat.test, type = 'response')
  predCopula.mu2 <- predict(bivBinCopula$mu2, newdata = dat.test, type = 'response')
  predCopula.rho <- predict(bivBinCopula$rho, newdata = dat.test, type = 'response')
  
  ######### COMPUTE PERFORMANCE: 
  # AUC (margin 1, then margin 2)
  AUCbivBinCopula$mu1 <- roc(dat.test$y1, as.numeric(predCopula.mu1))$auc
  AUCbivBinCopula$mu2 <- roc(dat.test$y2, as.numeric(predCopula.mu2))$auc
  
  # Brier score (margin 1, then margin 2)
  BrierbivBinCopula$mu1 <- mean((as.numeric(predCopula.mu1)-(as.numeric(dat.test$y1)))^2)
  BrierbivBinCopula$mu2 <- mean((as.numeric(predCopula.mu2)-(as.numeric(dat.test$y2)))^2)
  
  # LOG-SCORE (LOG-LIKELIHOOD) (change to the copula likelihood)
  lik$biv <- sum(loss(mu1 = predCopula.mu1, mu2 = predCopula.mu2, rho = predCopula.rho, y = y.test))
  
  
  ###########################################################################################################################
  ################ - univariat - ############################################################################################
  ###########################################################################################################################
  
  mstop.uni <- vector('list')
  coef.uni  <- vector('list')
  coef.uni_1  <- vector('list')
  AUC.uni   <- vector('list')
  Brier.uni <- vector('list')
  
  # - mu1
  dat.train.mu1 = dat.train[, !(names(dat.train) %in% "y2")]
  dat.test.mu1  = dat.test[ , !(names(dat.test) %in% "y2")]
  
  glm.uni.mu1 <- glmboost(as.factor(y1) ~. , 
                          data = dat.train.mu1, 
                          family = Binomial(link = "probit", type = "glm"), 
                          control = boost_control(risk = 'oobag', nu = boost.nu, mstop = 1000, trace = T), 
                          weights = weight.mstop)
  
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  if(mstop.uni$mu1 >= 990){
    glm.uni.mu1[3000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 2990){
    glm.uni.mu1[5000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 4990){
    glm.uni.mu1[8000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 7990){
    glm.uni.mu1[10000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 9990){
    glm.uni.mu1[15000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 14990){
    glm.uni.mu1[20000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  
  # dat.train_mu1 <- dat.train.mu1[weight.mstop == 1, ]
  # glm.uni.mu1 <- glmboost(as.factor(y1) ~. , 
  #                         data = dat.train_mu1, 
  #                         family = Binomial(link = "cloglog", type = "glm"), 
  #                         control = boost_control(mstop = mstop.uni$mu1, nu  = boost.nu))
  glm.uni.mu1 <- glm.uni.mu1[mstop.uni$mu1]
  
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 1  
  coef.uni$mu1 <- coef(glm.uni.mu1, which = '')
  
  AUC.uni$mu1 <- roc(dat.test.mu1$y1, as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response")))$auc
  Brier.uni$mu1 <- mean((as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response"))-(as.numeric(dat.test.mu1$y1)))^2)
  
  # - mu2
  dat.train.mu2 =  dat.train[, !(names(dat.train) %in% "y1")]
  dat.test.mu2 = dat.test[, !(names(dat.test) %in% "y1")]
  
  glm.uni.mu2 <- glmboost(as.factor(y2) ~. , 
                          data = dat.train.mu2, 
                          family = Binomial(link = "cloglog", type = "glm"), 
                          control = boost_control(risk = 'oobag', nu = boost.nu, mstop = 1500, trace = T), 
                          weights = weight.mstop)
  
  
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  if(mstop.uni$mu2 >= 1490){
    glm.uni.mu2[3000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 2990){
    glm.uni.mu2[5000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 4990){
    glm.uni.mu2[8000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 7990){
    glm.uni.mu2[10000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 9990){
    glm.uni.mu2[15000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  
  if(mstop.uni$mu2 >= 14990){
    glm.uni.mu2[20000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  
  
  #   dat.train_mu2 <- dat.train.mu2[weight.mstop == 1, ]
  # 
  # glm.uni.mu2 <- glmboost(as.factor(y2) ~., data = dat.train_mu2, family = Binomial(type = "glm"), 
  #                         control = boost_control(mstop = mstop.uni$mu1, nu  = 0.05))
  glm.uni.mu2 <- glm.uni.mu2[mstop.uni$mu2]
  
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 2  
  coef.uni$mu2 <- coef(glm.uni.mu2,which = '')
  
  AUC.uni$mu2 <- roc(dat.test.mu2$y2, as.numeric(predict(glm.uni.mu2, newdata = dat.test.mu2,type = "response")))$auc
  Brier.uni$mu2 <- mean((as.numeric(predict(glm.uni.mu2, newdata = dat.test.mu2,type = "response"))-(as.numeric(dat.test.mu2$y2)))^2)
  
  # Likelihood
  pred.mu1.uni <- as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response"))
  pred.mu2.uni <- as.numeric(predict(glm.uni.mu2, newdata = dat.test.mu2,type = "response"))
  
  mu1.uni.loglik <- sum(-dbinom(x = dat.test.mu1$y1, size = 1, prob = pred.mu1.uni, log = T))
  mu2.uni.loglik <- sum(-dbinom(x = dat.test.mu2$y2, size = 1, prob = pred.mu2.uni, log = T))
  
  lik$uni<- mu1.uni.loglik + mu2.uni.loglik
  
  
  
  n.ges = c(n.train, n.mstop, n.test)
  pred.ges <- NULL#list(pred.mu1, pred.mu2, pred.or)
  pred.uni <- list(pred.mu1.uni, pred.mu2.uni)
  pred.cop <- list(predCopula.mu1, predCopula.mu2, predCopula.rho)
  
 
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
    sample_uni <- data.gen.bivbin_energyscore(n = 1000, mu1 = pred.mu1.uni[i], mu2 = pred.mu2.uni[i], FAM = 4, theta = 1)
    
    pred_sample_uni[1,] <- sample_uni[,1]
    pred_sample_uni[2,] <- sample_uni[,2]
    
    es_uni[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_uni) 
    
    
    # copula
    sample_cop <- data.gen.bivbin_energyscore(n = 1000, mu1 = predCopula.mu1[i],  mu2 = predCopula.mu2[i], FAM = 4, theta = predCopula.rho[i])
    
    pred_sample_cop[1, ] <- sample_cop[,1]
    pred_sample_cop[2, ] <- sample_cop[,2]
    
    es_cop[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_cop) 
    
   
  }
  
  energy_score <- list()
  energy_score$uni <- mean(es_uni)
  energy_score$cop <- mean(es_cop)
  
  # if(is.null(bivBinGJRM_catch$value)){
  #   
  #   energy_score$gjrm <- NA
  # 
  # }else{
  #   
  #   energy_score$gjrm <- mean(es_gjrm)
  #   
  # }
  
  
  return(list(TrueBeta = TrueBeta, 
              n = n.ges, 
              p  = p, 
              #Coefficients = coef.bivBern, 
              #oobag.risk = oobag.risk, 
              AUCbivBern = AUCbivBern, 
              BrierbivBern = BrierbivBern, 
              Likelihood = lik,
              predict = pred.ges, 
              predict.uni = pred.uni, 
              predict.cop = pred.cop,
              energy_score = energy_score,  
              mstop = mstop.bivBern, 
              ####
              Coefficients.uni = coef.uni, 
              AUC.uni = AUC.uni, 
              Brier.uni = Brier.uni, 
              mstop.uni = mstop.uni,
              ###
              CoefficientsCOPULA = coef.bivBinCopula, 
              oobag.riskCOPULA = oobag.risk.cop, 
              AUCCOPULA = AUCbivBinCopula, 
              BrierCOPULA = BrierbivBinCopula, 
              LikelihoodCOPULA = likBinCopula,
              energy_scoreCOPULA = energy_score,  
              mstopCOPULA = mstop.bivBinCopula 
              ###
              # Coefficients = coef.GJRM, 
              # AUCbivBern = AUC.GJRM, 
              # BrierbivBern = Brier.GJRM, 
              # Likelihood = lik$GJRM,
              # energy_score = energy_score$gjrm 
  ))
  
}

sim_CLAYTON <- function(seed, n.train, n.mstop, p, n.test, corr, boost.nu = 0.1){
  
  
  loss <- function(mu1, mu2, rho, y){
    
    cdf1 <- mu1
    
    cdf2 <- mu2
    
    rho <- rho
    
    p11 <-  pdffz(VineCopula::BiCopCDF(u1 = cdf1, u2 = cdf2, family = 3, par = rho))
    
    
    return(- ( y[,1]*y[,2] * log( p11 ) +
                 y[,1]*(1-y[,2]) * log( pdffz((cdf1 - p11)) ) +
                 (1-y[,1])*y[,2] * log( pdffz((cdf2 - p11)) ) +
                 (1-y[,1])*(1-y[,2]) * log( pdffz((1 - cdf1 - cdf2 + p11)) ) 
    ) )
    
  }
  
  data.gen.bivbin <- function(FAM, mu1, mu2, theta){
    
    
    u1u2 <- VineCopula::BiCopSim(length(mu1), family = FAM, par = theta)
    
    y1 <- qbinom(p = u1u2[,1], size = 1, prob = mu1)
    
    y2 <- qbinom(p = u1u2[,2], size = 1, prob = mu2)
    
    
    dat <- data.frame(y1, y2)
    
    return(dat)
  }
  
  data.gen.bivbin_energyscore <- function(n, FAM, mu1, mu2, theta){
    
    
    u1u2 <- VineCopula::BiCopSim(n, family = FAM, par = theta)
    
    y1 <- qbinom(p = u1u2[,1], size = 1, prob = mu1)
    
    y2 <- qbinom(p = u1u2[,2], size = 1, prob = mu2)
    
    
    dat <- data.frame(y1, y2)
    
    return(dat)
  }
  
  
  
  # ## mu 1
  # mu1 <- c(0, -1, +0.5, 1, 0, -0.5, rep(0,p-6))
  mu1 <- c(0, -1, +0.5, 1, 0, -0.5, rep(0,p-6))
  # 
  # ### mu 2
  # mu2 <- c(-0.5, -1, 1.5, 0, 0, rep(0,p-5)) 
  mu2 <- c(0.5, -1, 0.75, 0, 0, rep(0,p-5)) 
  # 
  # # copula parameter
  # or  <- c(0, -2, +0.5, +2.5, 0, 0, 0, 0, rep(0,p-8))
  or  <- c(-0.5, +0.5, -1.5, 1.5, 0, 0, 0, 0,rep(0,p-8))
  or  <- c(-0.5, +1.25, -1, 2, 0, 0, 0, 0,rep(0,p-8))
  
  
  TrueBeta <-  vector('list')
  TrueBeta$mu1 <- mu1
  TrueBeta$mu2 <- mu2
  TrueBeta$or <- or
  
  set.seed(seed)
  
  n <- n.train + n.mstop
  weight.mstop <- c(rep(1, times = n.train),rep(0, times = n.mstop)) 
  
  
  ############################################# ############################################# #############################################
  # - training covariates
  x.train <- rmvnorm(n = n, mean = rep(0,p), sigma =  toeplitz(sapply(seq(0,p-1), function(x) corr^x)))
  
  x.train <- apply(x.train, 2, pnorm)
  
  # - test covariates
  x.test <- rmvnorm(n = n.test, mean = rep(0,p), sigma = toeplitz(sapply(seq(0,p-1), function(x) corr^x)))
  
  x.test <- apply(x.test, 2, pnorm)
  
  
  # GENERATION OF DISTRIBUTION PARAMETERS AND STUFF 
  # predictor
  train.eta.mu1_center <-  (x.train %*% mu1)
  
  
  ### MU1
  train.eta.mu1<-  pnorm(train.eta.mu1_center) 
  
  test.eta.mu1 <-  pnorm(x.test %*% mu1)
  
  ### MU2
  train.eta.mu2 <-  x.train %*% mu2
  
  train.eta.mu2 <-  1-exp(-exp(train.eta.mu2)) 
  
  test.eta.mu2 <-  1-exp(-exp(x.test %*% mu2)) 
  
  ### Copula parameter
  train.eta.or_center <-   (x.train %*% or)
  
  train.copula.parameter <- exp(train.eta.or_center) 
  
  test.copula.parameter <- exp( x.test %*% or) 
  
  
  TrueKendallTauRange <- VineCopula::BiCopPar2Tau(family = 3, par = range(train.copula.parameter))
  #TrueKendallTauRange
  
  #summary(VineCopula::BiCopPar2Tau(family = 4, par = train.copula.parameter))
  # range(train.eta.mu1)
  # range(train.eta.mu2)
  # range(train.copula.parameter)
  # 
  # plot(train.eta.mu1)
  # plot(train.eta.mu2)
  # plot(train.copula.parameter)
  # 
  
  ############################################# ############################################# #############################################
  ############################################# ############################################# #############################################
  
  
  
  y.train <- data.gen.bivbin(FAM = 3, mu1 = train.eta.mu1, mu2 = train.eta.mu2, theta = train.copula.parameter)
  
  colnames(y.train)<- c('y1','y2')
  dat.train <- data.frame(y.train,x.train)
  
  y.test <- data.gen.bivbin(FAM = 3, mu1 = test.eta.mu1, mu2 = test.eta.mu2, theta = test.copula.parameter)
  
  colnames(y.test)<- c('y1','y2')
  dat.test <- data.frame(y.test,x.test)
  
  
  
  # FITTIIIINGG
  mstop.bivBern <-  vector('list')
  
  AUCbivBern <-  vector('list')
  BrierbivBern <-  vector('list')
  lik <- vector('list')
  
  
  #############################################
  ################# - COPULA - ################
  #############################################
  # Bivariate copula model
  bivBinCopula <- glmboostLSS(cbind(y1,y2)~., data = dat.train, 
                              families = Clayton_Cop_BivBinary(marg1 = "PROBIT", 
                                                              marg2 = "CLOGLOG", 
                                                              stabilization = "L2"), 
                              control = boost_control(mstop = 1000, 
                                                      risk = 'oobag', 
                                                      nu = boost.nu, 
                                                      trace = TRUE),
                              method = 'noncyclic', 
                              weights = weight.mstop)
  
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 990){
    bivBinCopula[3000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 2990){
    bivBinCopula[5000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 4990){
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
  oobag.risk.cop <- risk(bivBinCopula,merge = T)
  
  # rm(bivBinCopula)
  # dat.train_biv <- dat.train[weight.mstop == 1, ]
  
  # # RE-FIT until OPTIMAL MSTOP
  # bivBinCopula <- glmboostLSS(cbind(y1,y2)~., data = dat.train_biv, families = Gauss_Cop_BivBinary(marg1 = "PROBIT", marg2 = "CLOGLOG"), 
  #                             control = boost_control(mstop = MSTOP_COP, nu  = 0.05), method = 'noncyclic')
  bivBinCopula <- bivBinCopula[MSTOP_COP]
  
  
  mstop.bivBinCopula <-  vector('list')
  mstop.bivBinCopula$mstop <- MSTOP_COP
  
  mstop.bivBinCopula$mu1 <- bivBinCopula$mu1$mstop()
  mstop.bivBinCopula$mu2 <- bivBinCopula$mu2$mstop()
  
  mstop.bivBinCopula$rho  <- bivBinCopula$rho$mstop()
  
  # EXTRACT ALL COEFFICIENTS!
  coef.bivBinCopula <- coef(bivBinCopula, which = "")
  
  AUCbivBinCopula <-  vector('list')
  BrierbivBinCopula <-  vector('list')
  likBinCopula <- vector('list')
  
  
  # Predict distributional quantities
  predCopula.mu1 <- predict(bivBinCopula$mu1, newdata = dat.test, type = 'response')
  predCopula.mu2 <- predict(bivBinCopula$mu2, newdata = dat.test, type = 'response')
  predCopula.rho <- predict(bivBinCopula$rho, newdata = dat.test, type = 'response')
  
  ######### COMPUTE PERFORMANCE: 
  # AUC (margin 1, then margin 2)
  AUCbivBinCopula$mu1 <- roc(dat.test$y1, as.numeric(predCopula.mu1))$auc
  AUCbivBinCopula$mu2 <- roc(dat.test$y2, as.numeric(predCopula.mu2))$auc
  
  # Brier score (margin 1, then margin 2)
  BrierbivBinCopula$mu1 <- mean((as.numeric(predCopula.mu1)-(as.numeric(dat.test$y1)))^2)
  BrierbivBinCopula$mu2 <- mean((as.numeric(predCopula.mu2)-(as.numeric(dat.test$y2)))^2)
  
  # LOG-SCORE (LOG-LIKELIHOOD) (change to the copula likelihood)
  lik$biv <- sum(loss(mu1 = predCopula.mu1, mu2 = predCopula.mu2, rho = predCopula.rho, y = y.test))
  
  
  ###########################################################################################################################
  ################ - univariat - ############################################################################################
  ###########################################################################################################################
  
  mstop.uni <- vector('list')
  coef.uni  <- vector('list')
  coef.uni_1  <- vector('list')
  AUC.uni   <- vector('list')
  Brier.uni <- vector('list')
  
  # - mu1
  dat.train.mu1 = dat.train[, !(names(dat.train) %in% "y2")]
  dat.test.mu1  = dat.test[ , !(names(dat.test) %in% "y2")]
  
  glm.uni.mu1 <- glmboost(as.factor(y1) ~. , 
                          data = dat.train.mu1, 
                          family = Binomial(link = "probit", type = "glm"), 
                          control = boost_control(risk = 'oobag', nu = boost.nu, mstop = 1000, trace = T), 
                          weights = weight.mstop)
  
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  if(mstop.uni$mu1 >= 990){
    glm.uni.mu1[3000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 2990){
    glm.uni.mu1[5000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 4990){
    glm.uni.mu1[8000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 7990){
    glm.uni.mu1[10000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 9990){
    glm.uni.mu1[15000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 14990){
    glm.uni.mu1[20000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  
  # dat.train_mu1 <- dat.train.mu1[weight.mstop == 1, ]
  # glm.uni.mu1 <- glmboost(as.factor(y1) ~. , 
  #                         data = dat.train_mu1, 
  #                         family = Binomial(link = "cloglog", type = "glm"), 
  #                         control = boost_control(mstop = mstop.uni$mu1, nu  = boost.nu))
  glm.uni.mu1 <- glm.uni.mu1[mstop.uni$mu1]
  
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 1  
  coef.uni$mu1 <- coef(glm.uni.mu1, which = '')
  
  AUC.uni$mu1 <- roc(dat.test.mu1$y1, as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response")))$auc
  Brier.uni$mu1 <- mean((as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response"))-(as.numeric(dat.test.mu1$y1)))^2)
  
  # - mu2
  dat.train.mu2 =  dat.train[, !(names(dat.train) %in% "y1")]
  dat.test.mu2 = dat.test[, !(names(dat.test) %in% "y1")]
  
  glm.uni.mu2 <- glmboost(as.factor(y2) ~. , 
                          data = dat.train.mu2, 
                          family = Binomial(link = "cloglog", type = "glm"), 
                          control = boost_control(risk = 'oobag', nu = boost.nu, mstop = 1500, trace = T), 
                          weights = weight.mstop)
  
  
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  if(mstop.uni$mu2 >= 1490){
    glm.uni.mu2[3000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 2990){
    glm.uni.mu2[5000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 4990){
    glm.uni.mu2[8000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 7990){
    glm.uni.mu2[10000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 9990){
    glm.uni.mu2[15000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  
  if(mstop.uni$mu2 >= 14990){
    glm.uni.mu2[20000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  
  
  #   dat.train_mu2 <- dat.train.mu2[weight.mstop == 1, ]
  # 
  # glm.uni.mu2 <- glmboost(as.factor(y2) ~., data = dat.train_mu2, family = Binomial(type = "glm"), 
  #                         control = boost_control(mstop = mstop.uni$mu1, nu  = 0.05))
  glm.uni.mu2 <- glm.uni.mu2[mstop.uni$mu2]
  
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 2  
  coef.uni$mu2 <- coef(glm.uni.mu2,which = '')
  
  AUC.uni$mu2 <- roc(dat.test.mu2$y2, as.numeric(predict(glm.uni.mu2, newdata = dat.test.mu2,type = "response")))$auc
  Brier.uni$mu2 <- mean((as.numeric(predict(glm.uni.mu2, newdata = dat.test.mu2,type = "response"))-(as.numeric(dat.test.mu2$y2)))^2)
  
  # Likelihood
  pred.mu1.uni <- as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response"))
  pred.mu2.uni <- as.numeric(predict(glm.uni.mu2, newdata = dat.test.mu2,type = "response"))
  
  mu1.uni.loglik <- sum(-dbinom(x = dat.test.mu1$y1, size = 1, prob = pred.mu1.uni, log = T))
  mu2.uni.loglik <- sum(-dbinom(x = dat.test.mu2$y2, size = 1, prob = pred.mu2.uni, log = T))
  
  lik$uni<- mu1.uni.loglik + mu2.uni.loglik
  
  
  
  n.ges = c(n.train, n.mstop, n.test)
  pred.ges <- NULL#list(pred.mu1, pred.mu2, pred.or)
  pred.uni <- list(pred.mu1.uni, pred.mu2.uni)
  pred.cop <- list(predCopula.mu1, predCopula.mu2, predCopula.rho)
  
  
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
    sample_uni <- data.gen.bivbin_energyscore(n = 1000, mu1 = pred.mu1.uni[i], mu2 = pred.mu2.uni[i], FAM = 3, theta = sqrt(.Machine$double.eps) )
    
    pred_sample_uni[1,] <- sample_uni[,1]
    pred_sample_uni[2,] <- sample_uni[,2]
    
    es_uni[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_uni) 
    
    
    # copula
    sample_cop <- data.gen.bivbin_energyscore(n = 1000, mu1 = predCopula.mu1[i],  mu2 = predCopula.mu2[i], FAM = 3, theta = predCopula.rho[i])
    
    pred_sample_cop[1, ] <- sample_cop[,1]
    pred_sample_cop[2, ] <- sample_cop[,2]
    
    es_cop[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_cop) 
    
    
  }
  
  energy_score <- list()
  energy_score$uni <- mean(es_uni)
  energy_score$cop <- mean(es_cop)
  
  # if(is.null(bivBinGJRM_catch$value)){
  #   
  #   energy_score$gjrm <- NA
  # 
  # }else{
  #   
  #   energy_score$gjrm <- mean(es_gjrm)
  #   
  # }
  
  
  return(list(TrueBeta = TrueBeta, 
              n = n.ges, 
              p  = p, 
              #Coefficients = coef.bivBern, 
              #oobag.risk = oobag.risk, 
              AUCbivBern = AUCbivBern, 
              BrierbivBern = BrierbivBern, 
              Likelihood = lik,
              predict = pred.ges, 
              predict.uni = pred.uni, 
              predict.cop = pred.cop,
              energy_score = energy_score,  
              mstop = mstop.bivBern, 
              ####
              Coefficients.uni = coef.uni, 
              AUC.uni = AUC.uni, 
              Brier.uni = Brier.uni, 
              mstop.uni = mstop.uni,
              ###
              CoefficientsCOPULA = coef.bivBinCopula, 
              oobag.riskCOPULA = oobag.risk.cop, 
              AUCCOPULA = AUCbivBinCopula, 
              BrierCOPULA = BrierbivBinCopula, 
              LikelihoodCOPULA = likBinCopula,
              energy_scoreCOPULA = energy_score,  
              mstopCOPULA = mstop.bivBinCopula, 
              KendallRange = TrueKendallTauRange
              ###
              # Coefficients = coef.GJRM, 
              # AUCbivBern = AUC.GJRM, 
              # BrierbivBern = Brier.GJRM, 
              # Likelihood = lik$GJRM,
              # energy_score = energy_score$gjrm 
  ))
  
}

sim_JOE <- function(seed, n.train, n.mstop, p, n.test, corr, boost.nu = 0.1){
  
  
  loss <- function(mu1, mu2, rho, y){
    
    cdf1 <- mu1
    
    cdf2 <- mu2
    
    rho <- rho
    
    p11 <-  pdffz(VineCopula::BiCopCDF(u1 = cdf1, u2 = cdf2, family = 6, par = rho))
    
    
    return(- ( y[,1]*y[,2] * log( p11 ) +
                 y[,1]*(1-y[,2]) * log( pdffz((cdf1 - p11)) ) +
                 (1-y[,1])*y[,2] * log( pdffz((cdf2 - p11)) ) +
                 (1-y[,1])*(1-y[,2]) * log( pdffz((1 - cdf1 - cdf2 + p11)) ) 
    ) )
    
  }
  
  data.gen.bivbin <- function(FAM, mu1, mu2, theta){
    
    
    u1u2 <- VineCopula::BiCopSim(length(mu1), family = FAM, par = theta)
    
    y1 <- qbinom(p = u1u2[,1], size = 1, prob = mu1)
    
    y2 <- qbinom(p = u1u2[,2], size = 1, prob = mu2)
    
    
    dat <- data.frame(y1, y2)
    
    return(dat)
  }
  
  data.gen.bivbin_energyscore <- function(n, FAM, mu1, mu2, theta){
    
    
    u1u2 <- VineCopula::BiCopSim(n, family = FAM, par = theta)
    
    y1 <- qbinom(p = u1u2[,1], size = 1, prob = mu1)
    
    y2 <- qbinom(p = u1u2[,2], size = 1, prob = mu2)
    
    
    dat <- data.frame(y1, y2)
    
    return(dat)
  }
  
  
  
  ## mu 1
  mu1 <- c(0, -1, +0.5, 1, 0, -0.5, rep(0,p-6))
  
  ### mu 2
  mu2 <- c(-0.5, -1, 1.5, 0, 0, rep(0,p-5)) 
  
  # copula parameter
  or  <- c(0, -2, +0.5, +2.5, 0, 0, 0, 0, rep(0,p-8))
  
  
  TrueBeta <-  vector('list')
  TrueBeta$mu1 <- mu1
  TrueBeta$mu2 <- mu2
  TrueBeta$or <- or
  
  set.seed(seed)
  
  n <- n.train + n.mstop
  weight.mstop <- c(rep(1, times = n.train),rep(0, times = n.mstop)) 
  
  
  ############################################# ############################################# #############################################
  # - training covariates
  x.train <- rmvnorm(n = n, mean = rep(0,p), sigma =  toeplitz(sapply(seq(0,p-1), function(x) corr^x)))
  
  x.train <- apply(x.train, 2, pnorm)
  
  # - test covariates
  x.test <- rmvnorm(n = n.test, mean = rep(0,p), sigma = toeplitz(sapply(seq(0,p-1), function(x) corr^x)))
  
  x.test <- apply(x.test, 2, pnorm)
  
  
  # GENERATION OF DISTRIBUTION PARAMETERS AND STUFF 
  # predictor
  train.eta.mu1_center <-  (x.train %*% mu1)
  
  
  ### MU1
  train.eta.mu1<-  pnorm(train.eta.mu1_center) 
  
  test.eta.mu1 <-  pnorm(x.test %*% mu1)
  
  ### MU2
  train.eta.mu2 <-  x.train %*% mu2
  
  train.eta.mu2 <-  1-exp(-exp(train.eta.mu2)) 
  
  test.eta.mu2 <-  1-exp(-exp(x.test %*% mu2)) 
  
  ### Copula parameter
  train.eta.or_center <-   (x.train %*% or)
  
  train.copula.parameter <- exp(train.eta.or_center) + 1
  
  test.copula.parameter <- exp( x.test %*% or) + 1
  
  
  TrueKendallTauRange <- VineCopula::BiCopPar2Tau(family = 6, par = range(train.copula.parameter))
  TrueKendallTauRange
  
  #summary(VineCopula::BiCopPar2Tau(family = 4, par = train.copula.parameter))
  range(train.eta.mu1)
  range(train.eta.mu2)
  range(train.copula.parameter)

  plot(train.eta.mu1)
  plot(train.eta.mu2)
  plot(train.copula.parameter)

  
  ############################################# ############################################# #############################################
  ############################################# ############################################# #############################################
  
  
  
  y.train <- data.gen.bivbin(FAM = 6, mu1 = train.eta.mu1, mu2 = train.eta.mu2, theta = train.copula.parameter)
  
  colnames(y.train)<- c('y1','y2')
  dat.train <- data.frame(y.train,x.train)
  
  y.test <- data.gen.bivbin(FAM = 6, mu1 = test.eta.mu1, mu2 = test.eta.mu2, theta = test.copula.parameter)
  
  colnames(y.test)<- c('y1','y2')
  dat.test <- data.frame(y.test,x.test)
  
  
  
  # FITTIIIINGG
  mstop.bivBern <-  vector('list')
  
  AUCbivBern <-  vector('list')
  BrierbivBern <-  vector('list')
  lik <- vector('list')
  
  
  #############################################
  ################# - COPULA - ################
  #############################################
  # Bivariate copula model
  bivBinCopula <- glmboostLSS(cbind(y1,y2) ~., 
                              data = dat.train, 
                              families = Joe_Cop_BivBinary(marg1 = "PROBIT", 
                                                              marg2 = "CLOGLOG", 
                                                              stabilization = "L2"), 
                              control = boost_control(mstop = 1000, 
                                                      risk = 'oobag', 
                                                      nu = boost.nu, 
                                                      trace = TRUE),
                              method = 'noncyclic', 
                              weights = weight.mstop)
  
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 990){
    bivBinCopula[3000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 2990){
    bivBinCopula[5000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 4990){
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
  oobag.risk.cop <- risk(bivBinCopula,merge = T)
  
  # rm(bivBinCopula)
  # dat.train_biv <- dat.train[weight.mstop == 1, ]
  
  # # RE-FIT until OPTIMAL MSTOP
  # bivBinCopula <- glmboostLSS(cbind(y1,y2)~., data = dat.train_biv, families = Gauss_Cop_BivBinary(marg1 = "PROBIT", marg2 = "CLOGLOG"), 
  #                             control = boost_control(mstop = MSTOP_COP, nu  = 0.05), method = 'noncyclic')
  bivBinCopula <- bivBinCopula[MSTOP_COP]
  
  
  mstop.bivBinCopula <-  vector('list')
  mstop.bivBinCopula$mstop <- MSTOP_COP
  
  mstop.bivBinCopula$mu1 <- bivBinCopula$mu1$mstop()
  mstop.bivBinCopula$mu2 <- bivBinCopula$mu2$mstop()
  
  mstop.bivBinCopula$rho  <- bivBinCopula$rho$mstop()
  
  # EXTRACT ALL COEFFICIENTS!
  coef.bivBinCopula <- coef(bivBinCopula, which = "")
  
  AUCbivBinCopula <-  vector('list')
  BrierbivBinCopula <-  vector('list')
  likBinCopula <- vector('list')
  
  
  # Predict distributional quantities
  predCopula.mu1 <- predict(bivBinCopula$mu1, newdata = dat.test, type = 'response')
  predCopula.mu2 <- predict(bivBinCopula$mu2, newdata = dat.test, type = 'response')
  predCopula.rho <- predict(bivBinCopula$rho, newdata = dat.test, type = 'response')
  
  ######### COMPUTE PERFORMANCE: 
  # AUC (margin 1, then margin 2)
  AUCbivBinCopula$mu1 <- roc(dat.test$y1, as.numeric(predCopula.mu1))$auc
  AUCbivBinCopula$mu2 <- roc(dat.test$y2, as.numeric(predCopula.mu2))$auc
  
  # Brier score (margin 1, then margin 2)
  BrierbivBinCopula$mu1 <- mean((as.numeric(predCopula.mu1)-(as.numeric(dat.test$y1)))^2)
  BrierbivBinCopula$mu2 <- mean((as.numeric(predCopula.mu2)-(as.numeric(dat.test$y2)))^2)
  
  # LOG-SCORE (LOG-LIKELIHOOD) (change to the copula likelihood)
  lik$biv <- sum(loss(mu1 = predCopula.mu1, mu2 = predCopula.mu2, rho = predCopula.rho, y = y.test))
  
  
  ###########################################################################################################################
  ################ - univariat - ############################################################################################
  ###########################################################################################################################
  
  mstop.uni <- vector('list')
  coef.uni  <- vector('list')
  coef.uni_1  <- vector('list')
  AUC.uni   <- vector('list')
  Brier.uni <- vector('list')
  
  # - mu1
  dat.train.mu1 = dat.train[, !(names(dat.train) %in% "y2")]
  dat.test.mu1  = dat.test[ , !(names(dat.test) %in% "y2")]
  
  glm.uni.mu1 <- glmboost(as.factor(y1) ~. , 
                          data = dat.train.mu1, 
                          family = Binomial(link = "probit", type = "glm"), 
                          control = boost_control(risk = 'oobag', nu = boost.nu, mstop = 1000, trace = T), 
                          weights = weight.mstop)
  
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  if(mstop.uni$mu1 >= 990){
    glm.uni.mu1[3000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 2990){
    glm.uni.mu1[5000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 4990){
    glm.uni.mu1[8000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 7990){
    glm.uni.mu1[10000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 9990){
    glm.uni.mu1[15000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 14990){
    glm.uni.mu1[20000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  
  # dat.train_mu1 <- dat.train.mu1[weight.mstop == 1, ]
  # glm.uni.mu1 <- glmboost(as.factor(y1) ~. , 
  #                         data = dat.train_mu1, 
  #                         family = Binomial(link = "cloglog", type = "glm"), 
  #                         control = boost_control(mstop = mstop.uni$mu1, nu  = boost.nu))
  glm.uni.mu1 <- glm.uni.mu1[mstop.uni$mu1]
  
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 1  
  coef.uni$mu1 <- coef(glm.uni.mu1, which = '')
  
  AUC.uni$mu1 <- roc(dat.test.mu1$y1, as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response")))$auc
  Brier.uni$mu1 <- mean((as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response"))-(as.numeric(dat.test.mu1$y1)))^2)
  
  # - mu2
  dat.train.mu2 =  dat.train[, !(names(dat.train) %in% "y1")]
  dat.test.mu2 = dat.test[, !(names(dat.test) %in% "y1")]
  
  glm.uni.mu2 <- glmboost(as.factor(y2) ~. , 
                          data = dat.train.mu2, 
                          family = Binomial(link = "cloglog", type = "glm"), 
                          control = boost_control(risk = 'oobag', nu = boost.nu, mstop = 1500, trace = T), 
                          weights = weight.mstop)
  
  
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  if(mstop.uni$mu2 >= 1490){
    glm.uni.mu2[3000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 2990){
    glm.uni.mu2[5000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 4990){
    glm.uni.mu2[8000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 7990){
    glm.uni.mu2[10000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 9990){
    glm.uni.mu2[15000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  
  if(mstop.uni$mu2 >= 14990){
    glm.uni.mu2[20000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  
  
  #   dat.train_mu2 <- dat.train.mu2[weight.mstop == 1, ]
  # 
  # glm.uni.mu2 <- glmboost(as.factor(y2) ~., data = dat.train_mu2, family = Binomial(type = "glm"), 
  #                         control = boost_control(mstop = mstop.uni$mu1, nu  = 0.05))
  glm.uni.mu2 <- glm.uni.mu2[mstop.uni$mu2]
  
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 2  
  coef.uni$mu2 <- coef(glm.uni.mu2,which = '')
  
  AUC.uni$mu2 <- roc(dat.test.mu2$y2, as.numeric(predict(glm.uni.mu2, newdata = dat.test.mu2,type = "response")))$auc
  Brier.uni$mu2 <- mean((as.numeric(predict(glm.uni.mu2, newdata = dat.test.mu2,type = "response"))-(as.numeric(dat.test.mu2$y2)))^2)
  
  # Likelihood
  pred.mu1.uni <- as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response"))
  pred.mu2.uni <- as.numeric(predict(glm.uni.mu2, newdata = dat.test.mu2,type = "response"))
  
  mu1.uni.loglik <- sum(-dbinom(x = dat.test.mu1$y1, size = 1, prob = pred.mu1.uni, log = T))
  mu2.uni.loglik <- sum(-dbinom(x = dat.test.mu2$y2, size = 1, prob = pred.mu2.uni, log = T))
  
  lik$uni<- mu1.uni.loglik + mu2.uni.loglik
  
  
  
  n.ges = c(n.train, n.mstop, n.test)
  pred.ges <- NULL#list(pred.mu1, pred.mu2, pred.or)
  pred.uni <- list(pred.mu1.uni, pred.mu2.uni)
  pred.cop <- list(predCopula.mu1, predCopula.mu2, predCopula.rho)
  
  
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
    sample_uni <- data.gen.bivbin_energyscore(n = 1000, mu1 = pred.mu1.uni[i], mu2 = pred.mu2.uni[i], FAM = 4, theta = 1)
    
    pred_sample_uni[1,] <- sample_uni[,1]
    pred_sample_uni[2,] <- sample_uni[,2]
    
    es_uni[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_uni) 
    
    
    # copula
    sample_cop <- data.gen.bivbin_energyscore(n = 1000, mu1 = predCopula.mu1[i],  mu2 = predCopula.mu2[i], FAM = 4, theta = predCopula.rho[i])
    
    pred_sample_cop[1, ] <- sample_cop[,1]
    pred_sample_cop[2, ] <- sample_cop[,2]
    
    es_cop[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_cop) 
    
    
  }
  
  energy_score <- list()
  energy_score$uni <- mean(es_uni)
  energy_score$cop <- mean(es_cop)
  
  # if(is.null(bivBinGJRM_catch$value)){
  #   
  #   energy_score$gjrm <- NA
  # 
  # }else{
  #   
  #   energy_score$gjrm <- mean(es_gjrm)
  #   
  # }
  
  
  return(list(TrueBeta = TrueBeta, 
              n = n.ges, 
              p  = p, 
              #Coefficients = coef.bivBern, 
              #oobag.risk = oobag.risk, 
              AUCbivBern = AUCbivBern, 
              BrierbivBern = BrierbivBern, 
              Likelihood = lik,
              predict = pred.ges, 
              predict.uni = pred.uni, 
              predict.cop = pred.cop,
              energy_score = energy_score,  
              mstop = mstop.bivBern, 
              ####
              Coefficients.uni = coef.uni, 
              AUC.uni = AUC.uni, 
              Brier.uni = Brier.uni, 
              mstop.uni = mstop.uni,
              ###
              CoefficientsCOPULA = coef.bivBinCopula, 
              oobag.riskCOPULA = oobag.risk.cop, 
              AUCCOPULA = AUCbivBinCopula, 
              BrierCOPULA = BrierbivBinCopula, 
              LikelihoodCOPULA = likBinCopula,
              energy_scoreCOPULA = energy_score,  
              mstopCOPULA = mstop.bivBinCopula 
              ###
              # Coefficients = coef.GJRM, 
              # AUCbivBern = AUC.GJRM, 
              # BrierbivBern = Brier.GJRM, 
              # Likelihood = lik$GJRM,
              # energy_score = energy_score$gjrm 
  ))
  
}

sim_JOE90 <- function(seed, n.train, n.mstop, p, n.test, corr, boost.nu = 0.1){
  
  
  loss <- function(mu1, mu2, rho, y){
    
    cdf1 <- mu1
    
    cdf2 <- mu2
    
    rho <- rho
    
    p11 <-  pdffz(VineCopula::BiCopCDF(u1 = cdf1, u2 = cdf2, family = 6, par = rho))
    
    
    return(- ( y[,1]*y[,2] * log( p11 ) +
                 y[,1]*(1-y[,2]) * log( pdffz((cdf1 - p11)) ) +
                 (1-y[,1])*y[,2] * log( pdffz((cdf2 - p11)) ) +
                 (1-y[,1])*(1-y[,2]) * log( pdffz((1 - cdf1 - cdf2 + p11)) ) 
    ) )
    
  }
  
  data.gen.bivbin <- function(FAM, mu1, mu2, theta){
    
    
    u1u2 <- VineCopula::BiCopSim(length(mu1), family = FAM, par = theta)
    
    y1 <- qbinom(p = u1u2[,1], size = 1, prob = mu1)
    
    y2 <- qbinom(p = u1u2[,2], size = 1, prob = mu2)
    
    
    dat <- data.frame(y1, y2)
    
    return(dat)
  }
  
  data.gen.bivbin_energyscore <- function(n, FAM, mu1, mu2, theta){
    
    
    u1u2 <- VineCopula::BiCopSim(n, family = FAM, par = theta)
    
    y1 <- qbinom(p = u1u2[,1], size = 1, prob = mu1)
    
    y2 <- qbinom(p = u1u2[,2], size = 1, prob = mu2)
    
    
    dat <- data.frame(y1, y2)
    
    return(dat)
  }
  
  
  
  ## mu 1
  mu1 <- c(0, -1, +0.5, 1, 0, -0.5, rep(0,p-6))
  
  ### mu 2
  mu2 <- c(-0.5, -1, 1.5, 0, 0, rep(0,p-5)) 
  
  # copula parameter
  or  <- c(0, -2, +0.5, +2.5, 0, 0, 0, 0, rep(0,p-8))
  
  
  TrueBeta <-  vector('list')
  TrueBeta$mu1 <- mu1
  TrueBeta$mu2 <- mu2
  TrueBeta$or <- or
  
  set.seed(seed)
  
  n <- n.train + n.mstop
  weight.mstop <- c(rep(1, times = n.train),rep(0, times = n.mstop)) 
  
  
  ############################################# ############################################# #############################################
  # - training covariates
  x.train <- rmvnorm(n = n, mean = rep(0,p), sigma =  toeplitz(sapply(seq(0,p-1), function(x) corr^x)))
  
  x.train <- apply(x.train, 2, pnorm)
  
  # - test covariates
  x.test <- rmvnorm(n = n.test, mean = rep(0,p), sigma = toeplitz(sapply(seq(0,p-1), function(x) corr^x)))
  
  x.test <- apply(x.test, 2, pnorm)
  
  
  # GENERATION OF DISTRIBUTION PARAMETERS AND STUFF 
  # predictor
  train.eta.mu1_center <-  (x.train %*% mu1)
  
  
  ### MU1
  train.eta.mu1<-  pnorm(train.eta.mu1_center) 
  
  test.eta.mu1 <-  pnorm(x.test %*% mu1)
  
  ### MU2
  train.eta.mu2 <-  x.train %*% mu2
  
  train.eta.mu2 <-  1-exp(-exp(train.eta.mu2)) 
  
  test.eta.mu2 <-  1-exp(-exp(x.test %*% mu2)) 
  
  ### Copula parameter
  train.eta.or_center <-   (x.train %*% or)
  
  train.copula.parameter <- exp(train.eta.or_center) + 1
  
  test.copula.parameter <- exp( x.test %*% or) + 1
  
  
  TrueKendallTauRange <- VineCopula::BiCopPar2Tau(family = 6, par = range(train.copula.parameter))
  TrueKendallTauRange
  
  #summary(VineCopula::BiCopPar2Tau(family = 4, par = train.copula.parameter))
  range(train.eta.mu1)
  range(train.eta.mu2)
  range(train.copula.parameter)
  
  plot(train.eta.mu1)
  plot(train.eta.mu2)
  plot(train.copula.parameter)
  
  
  ############################################# ############################################# #############################################
  ############################################# ############################################# #############################################
  
  
  
  y.train <- data.gen.bivbin(FAM = 6, mu1 = train.eta.mu1, mu2 = train.eta.mu2, theta = train.copula.parameter)
  
  colnames(y.train)<- c('y1','y2')
  dat.train <- data.frame(y.train,x.train)
  
  y.test <- data.gen.bivbin(FAM = 6, mu1 = test.eta.mu1, mu2 = test.eta.mu2, theta = test.copula.parameter)
  
  colnames(y.test)<- c('y1','y2')
  dat.test <- data.frame(y.test,x.test)
  
  
  
  # FITTIIIINGG
  mstop.bivBern <-  vector('list')
  
  AUCbivBern <-  vector('list')
  BrierbivBern <-  vector('list')
  lik <- vector('list')
  
  
  #############################################
  ################# - COPULA - ################
  #############################################
  # Bivariate copula model
  bivBinCopula <- glmboostLSS(cbind(y1,y2) ~., 
                              data = dat.train, 
                              families = Joe_Cop_BivBinary(marg1 = "PROBIT", 
                                                           marg2 = "CLOGLOG", 
                                                           stabilization = "L2"), 
                              control = boost_control(mstop = 1000, 
                                                      risk = 'oobag', 
                                                      nu = boost.nu, 
                                                      trace = TRUE),
                              method = 'noncyclic', 
                              weights = weight.mstop)
  
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 990){
    bivBinCopula[3000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 2990){
    bivBinCopula[5000]
  }
  
  MSTOP_COP <- which.min(risk(bivBinCopula,merge = T))
  
  if(MSTOP_COP >= 4990){
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
  oobag.risk.cop <- risk(bivBinCopula,merge = T)
  
  # rm(bivBinCopula)
  # dat.train_biv <- dat.train[weight.mstop == 1, ]
  
  # # RE-FIT until OPTIMAL MSTOP
  # bivBinCopula <- glmboostLSS(cbind(y1,y2)~., data = dat.train_biv, families = Gauss_Cop_BivBinary(marg1 = "PROBIT", marg2 = "CLOGLOG"), 
  #                             control = boost_control(mstop = MSTOP_COP, nu  = 0.05), method = 'noncyclic')
  bivBinCopula <- bivBinCopula[MSTOP_COP]
  
  
  mstop.bivBinCopula <-  vector('list')
  mstop.bivBinCopula$mstop <- MSTOP_COP
  
  mstop.bivBinCopula$mu1 <- bivBinCopula$mu1$mstop()
  mstop.bivBinCopula$mu2 <- bivBinCopula$mu2$mstop()
  
  mstop.bivBinCopula$rho  <- bivBinCopula$rho$mstop()
  
  # EXTRACT ALL COEFFICIENTS!
  coef.bivBinCopula <- coef(bivBinCopula, which = "")
  
  AUCbivBinCopula <-  vector('list')
  BrierbivBinCopula <-  vector('list')
  likBinCopula <- vector('list')
  
  
  # Predict distributional quantities
  predCopula.mu1 <- predict(bivBinCopula$mu1, newdata = dat.test, type = 'response')
  predCopula.mu2 <- predict(bivBinCopula$mu2, newdata = dat.test, type = 'response')
  predCopula.rho <- predict(bivBinCopula$rho, newdata = dat.test, type = 'response')
  
  ######### COMPUTE PERFORMANCE: 
  # AUC (margin 1, then margin 2)
  AUCbivBinCopula$mu1 <- roc(dat.test$y1, as.numeric(predCopula.mu1))$auc
  AUCbivBinCopula$mu2 <- roc(dat.test$y2, as.numeric(predCopula.mu2))$auc
  
  # Brier score (margin 1, then margin 2)
  BrierbivBinCopula$mu1 <- mean((as.numeric(predCopula.mu1)-(as.numeric(dat.test$y1)))^2)
  BrierbivBinCopula$mu2 <- mean((as.numeric(predCopula.mu2)-(as.numeric(dat.test$y2)))^2)
  
  # LOG-SCORE (LOG-LIKELIHOOD) (change to the copula likelihood)
  lik$biv <- sum(loss(mu1 = predCopula.mu1, mu2 = predCopula.mu2, rho = predCopula.rho, y = y.test))
  
  
  ###########################################################################################################################
  ################ - univariat - ############################################################################################
  ###########################################################################################################################
  
  mstop.uni <- vector('list')
  coef.uni  <- vector('list')
  coef.uni_1  <- vector('list')
  AUC.uni   <- vector('list')
  Brier.uni <- vector('list')
  
  # - mu1
  dat.train.mu1 = dat.train[, !(names(dat.train) %in% "y2")]
  dat.test.mu1  = dat.test[ , !(names(dat.test) %in% "y2")]
  
  glm.uni.mu1 <- glmboost(as.factor(y1) ~. , 
                          data = dat.train.mu1, 
                          family = Binomial(link = "probit", type = "glm"), 
                          control = boost_control(risk = 'oobag', nu = boost.nu, mstop = 1000, trace = T), 
                          weights = weight.mstop)
  
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  if(mstop.uni$mu1 >= 990){
    glm.uni.mu1[3000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 2990){
    glm.uni.mu1[5000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 4990){
    glm.uni.mu1[8000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 7990){
    glm.uni.mu1[10000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 9990){
    glm.uni.mu1[15000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  if(mstop.uni$mu1 >= 14990){
    glm.uni.mu1[20000]
  }
  mstop.uni$mu1 <- which.min(risk(glm.uni.mu1))
  
  
  # dat.train_mu1 <- dat.train.mu1[weight.mstop == 1, ]
  # glm.uni.mu1 <- glmboost(as.factor(y1) ~. , 
  #                         data = dat.train_mu1, 
  #                         family = Binomial(link = "cloglog", type = "glm"), 
  #                         control = boost_control(mstop = mstop.uni$mu1, nu  = boost.nu))
  glm.uni.mu1 <- glm.uni.mu1[mstop.uni$mu1]
  
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 1  
  coef.uni$mu1 <- coef(glm.uni.mu1, which = '')
  
  AUC.uni$mu1 <- roc(dat.test.mu1$y1, as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response")))$auc
  Brier.uni$mu1 <- mean((as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response"))-(as.numeric(dat.test.mu1$y1)))^2)
  
  # - mu2
  dat.train.mu2 =  dat.train[, !(names(dat.train) %in% "y1")]
  dat.test.mu2 = dat.test[, !(names(dat.test) %in% "y1")]
  
  glm.uni.mu2 <- glmboost(as.factor(y2) ~. , 
                          data = dat.train.mu2, 
                          family = Binomial(link = "cloglog", type = "glm"), 
                          control = boost_control(risk = 'oobag', nu = boost.nu, mstop = 1500, trace = T), 
                          weights = weight.mstop)
  
  
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  if(mstop.uni$mu2 >= 1490){
    glm.uni.mu2[3000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 2990){
    glm.uni.mu2[5000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 4990){
    glm.uni.mu2[8000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 7990){
    glm.uni.mu2[10000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  if(mstop.uni$mu2 >= 9990){
    glm.uni.mu2[15000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  
  if(mstop.uni$mu2 >= 14990){
    glm.uni.mu2[20000]
  }
  mstop.uni$mu2 <- which.min(risk(glm.uni.mu2))
  
  
  
  #   dat.train_mu2 <- dat.train.mu2[weight.mstop == 1, ]
  # 
  # glm.uni.mu2 <- glmboost(as.factor(y2) ~., data = dat.train_mu2, family = Binomial(type = "glm"), 
  #                         control = boost_control(mstop = mstop.uni$mu1, nu  = 0.05))
  glm.uni.mu2 <- glm.uni.mu2[mstop.uni$mu2]
  
  
  # EXTRACT COEFFICIENTS AND COMPUTE AUC / BRIER SCORE FOR MARGIN 2  
  coef.uni$mu2 <- coef(glm.uni.mu2,which = '')
  
  AUC.uni$mu2 <- roc(dat.test.mu2$y2, as.numeric(predict(glm.uni.mu2, newdata = dat.test.mu2,type = "response")))$auc
  Brier.uni$mu2 <- mean((as.numeric(predict(glm.uni.mu2, newdata = dat.test.mu2,type = "response"))-(as.numeric(dat.test.mu2$y2)))^2)
  
  # Likelihood
  pred.mu1.uni <- as.numeric(predict(glm.uni.mu1, newdata = dat.test.mu1,type = "response"))
  pred.mu2.uni <- as.numeric(predict(glm.uni.mu2, newdata = dat.test.mu2,type = "response"))
  
  mu1.uni.loglik <- sum(-dbinom(x = dat.test.mu1$y1, size = 1, prob = pred.mu1.uni, log = T))
  mu2.uni.loglik <- sum(-dbinom(x = dat.test.mu2$y2, size = 1, prob = pred.mu2.uni, log = T))
  
  lik$uni<- mu1.uni.loglik + mu2.uni.loglik
  
  
  
  n.ges = c(n.train, n.mstop, n.test)
  pred.ges <- NULL#list(pred.mu1, pred.mu2, pred.or)
  pred.uni <- list(pred.mu1.uni, pred.mu2.uni)
  pred.cop <- list(predCopula.mu1, predCopula.mu2, predCopula.rho)
  
  
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
    sample_uni <- data.gen.bivbin_energyscore(n = 1000, mu1 = pred.mu1.uni[i], mu2 = pred.mu2.uni[i], FAM = 4, theta = 1)
    
    pred_sample_uni[1,] <- sample_uni[,1]
    pred_sample_uni[2,] <- sample_uni[,2]
    
    es_uni[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_uni) 
    
    
    # copula
    sample_cop <- data.gen.bivbin_energyscore(n = 1000, mu1 = predCopula.mu1[i],  mu2 = predCopula.mu2[i], FAM = 4, theta = predCopula.rho[i])
    
    pred_sample_cop[1, ] <- sample_cop[,1]
    pred_sample_cop[2, ] <- sample_cop[,2]
    
    es_cop[i] <- es_sample(y = c(y.test[i,1], y.test[i,2]), dat = pred_sample_cop) 
    
    
  }
  
  energy_score <- list()
  energy_score$uni <- mean(es_uni)
  energy_score$cop <- mean(es_cop)
  
  # if(is.null(bivBinGJRM_catch$value)){
  #   
  #   energy_score$gjrm <- NA
  # 
  # }else{
  #   
  #   energy_score$gjrm <- mean(es_gjrm)
  #   
  # }
  
  
  return(list(TrueBeta = TrueBeta, 
              n = n.ges, 
              p  = p, 
              #Coefficients = coef.bivBern, 
              #oobag.risk = oobag.risk, 
              AUCbivBern = AUCbivBern, 
              BrierbivBern = BrierbivBern, 
              Likelihood = lik,
              predict = pred.ges, 
              predict.uni = pred.uni, 
              predict.cop = pred.cop,
              energy_score = energy_score,  
              mstop = mstop.bivBern, 
              ####
              Coefficients.uni = coef.uni, 
              AUC.uni = AUC.uni, 
              Brier.uni = Brier.uni, 
              mstop.uni = mstop.uni,
              ###
              CoefficientsCOPULA = coef.bivBinCopula, 
              oobag.riskCOPULA = oobag.risk.cop, 
              AUCCOPULA = AUCbivBinCopula, 
              BrierCOPULA = BrierbivBinCopula, 
              LikelihoodCOPULA = likBinCopula,
              energy_scoreCOPULA = energy_score,  
              mstopCOPULA = mstop.bivBinCopula 
              ###
              # Coefficients = coef.GJRM, 
              # AUCbivBern = AUC.GJRM, 
              # BrierbivBern = Brier.GJRM, 
              # Likelihood = lik$GJRM,
              # energy_score = energy_score$gjrm 
  ))
  
}
