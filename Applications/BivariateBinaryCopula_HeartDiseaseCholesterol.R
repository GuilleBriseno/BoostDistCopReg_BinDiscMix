#########################################
library("gamboostLSS")

#source('bivariateBernoulli.R')
source("Copulas/BivariateBinary/Copula_Gaussian_BivBin.R")
#########################################

#load("../BinaryBinary_Application/cholesterol_ischemic_heart/genotype_data_chol_heart.RData")
#ukbb_i25_cholesterol <- read.delim("../BinaryBinary_Application/cholesterol_ischemic_heart/phenotypes_subset_chol_heart.phe")

load("../BinaryBinary_Application/cholesterol_ischemic_heart/genotype_data_chol_heart.RData")
ukbb_i25_cholesterol <- read.delim("../BinaryBinary_Application/cholesterol_ischemic_heart/phenotypes_subset_chol_heart.phe")


ukbb_i25_cholesterol$cholesterol_0 <- 0
ukbb_i25_cholesterol$cholesterol_0[which(ukbb_i25_cholesterol$cholesterol > 6.19)]  <- 1
  
ukbb_i25_cholesterol[, c('FID','IID', 'PC1','PC2', 'PC3', 'PC4','cholesterol','sex','age','bmi')]  <- list(NULL)
genotype_data[,c('FID', 'IID')] <- list(NULL)
genotype_data = data.frame(genotype_data)


### define step length:
boost.nu <- 0.1


table(ukbb_i25_cholesterol) / nrow(ukbb_i25_cholesterol)

sum( table(ukbb_i25_cholesterol) / nrow(ukbb_i25_cholesterol) )


#### Using Annika's seed:
set.seed(1)

ALL <- data.frame(cbind(ukbb_i25_cholesterol, genotype_data))

### Sample observations used for validating the number of fitting iterations : n_mstop, remaining are n_train
samp <- sample(1:30000, 10000)

ALL_n_mstop <- ALL[samp,]

ALL_n_train <- ALL[-samp,]

weights_train <- rep(1, 30000)  

weights_train[samp] <- 0 
  
perc <- table(ALL$I25_cases)


table(ukbb_i25_cholesterol[-samp, ]) / nrow(ukbb_i25_cholesterol[-samp, ])


#### GAUSSIAN COPULA MODEL: 
GaussianCopulaModel <- glmboostLSS(cbind(I25_cases, cholesterol_0) ~ ., 
                                   data = ALL, 
                                   families = Gauss_Cop_BivBinary(marg1 = "LOGIT", 
                                                                  marg2 = "LOGIT", 
                                                                  stabilization = "L2"), 
                                   method = 'noncyclic', 
                                   weights = weights_train,
                                   control = boost_control(mstop = 50, 
                                                           nu = boost.nu, 
                                                           trace = T, 
                                                           risk = 'oobag'))

GaussianCopulaModel <- GaussianCopulaModel[15000]


GaussianCopulaModel <- GaussianCopulaModel[18000]


GaussianCopulaModel <- GaussianCopulaModel[20000]




##### 
which.min(risk(GaussianCopulaModel, merge = T)) ###### using stabilization and final version of offsets: 17 933 !!!!

plot(risk(GaussianCopulaModel, merge = T), type = "l")

risk_all_iterations <- risk(GaussianCopulaModel, merge = T)

mstop_with_L2_stabilisation <- which.min(risk(GaussianCopulaModel, merge = T))


save(risk_all_iterations, 
     mstop_with_L2_stabilisation,
     file="../BinaryBinary_Application/HeartDiseaseCholesterol_BernoulliGaussCopulaResults_withL2Stabilisation.RData")




rm(GaussianCopulaModel)


ALL_FINALFIT <- ALL[weights_train == 1, ]


GaussianCopulaModel <- glmboostLSS(cbind(I25_cases, cholesterol_0) ~ ., 
                                   data = ALL_FINALFIT, 
                                   families = Gauss_Cop_BivBinary(marg1 = "LOGIT", 
                                                                  marg2 = "LOGIT", 
                                                                  stabilization = "L2"), 
                                   method = 'noncyclic', 
                                   control = boost_control(mstop = 100, 
                                                           nu = boost.nu, 
                                                           trace = T))


GaussianCopulaModel <- GaussianCopulaModel[mstop_with_L2_stabilisation]












#  OLD VERSION WITHOUT STABILISATION
# GaussianCopulaModel <- GaussianCopulaModel[5000]
# 
# GaussianCopulaModel <- GaussianCopulaModel[7000]
# 
# GaussianCopulaModel <- GaussianCopulaModel[8000] # up to here: 7874
# 
# GaussianCopulaModel <- GaussianCopulaModel[9000]
# 
# GaussianCopulaModel <- GaussianCopulaModel[10500] ### up to here: 10457 
# 
# GaussianCopulaModel <- GaussianCopulaModel[12000] ### up to here: 11964
# 
# GaussianCopulaModel <- GaussianCopulaModel[14000] ### up to here: 13075 

 
MSTOP_OPT <- which.min(risk(GaussianCopulaModel, merge = T)) # FIRST TRY: 13075 
MSTOP_OPT # 13075 

GaussianCopulaModel <- GaussianCopulaModel[MSTOP_OPT]

tail(risk(GaussianCopulaModel, merge = T), 1) #### FINAL RISK AT OOBAG: 8603.057 


rm(GaussianCopulaModel)


ALL_FINALFIT <- ALL[weights_train == 1, ]


GaussianCopulaModel <- glmboostLSS(cbind(I25_cases, cholesterol_0) ~ ., 
                                   data = ALL_FINALFIT, 
                                   families = Gauss_Cop_BivBinary(marg1 = "LOGIT", 
                                                                  marg2 = "LOGIT", 
                                                                  stabilization = "L2"), 
                                   method = 'noncyclic', 
                                   control = boost_control(mstop = 2000, 
                                                           nu = boost.nu, 
                                                           trace = T))


GaussianCopulaModel <- GaussianCopulaModel[13075]




length(coef(GaussianCopulaModel$mu1))
length(coef(GaussianCopulaModel$mu2))
length(coef(GaussianCopulaModel$rho))

GaussianCopulaModel$mu1$offset
GaussianCopulaModel$mu2$offset
GaussianCopulaModel$rho$offset

range(predict(GaussianCopulaModel$rho, type = "response"))
print(range(VineCopula::BiCopPar2Tau(family = 1, par = predict(GaussianCopulaModel$rho, type = "response"))), digits = 3)
plot(predict(GaussianCopulaModel$rho, type = "response"))

coef(GaussianCopulaModel$rho)[1]


############################## Risk
risk_all_sample <- risk(GaussianCopulaModel, merge = T)

#### Coefficients 
ListOfCoefficients <- coef(GaussianCopulaModel)

############ MSTOP
mstop.bivBern <-  vector('list')
mstop.bivBern$mstop <- MSTOP
mstop.bivBern$mu1 <- GaussianCopulaModel$mu1$mstop()
mstop.bivBern$mu2 <- GaussianCopulaModel$mu2$mstop()
mstop.bivBern$rho <- GaussianCopulaModel$rho$mstop()

#### COEFFICIENTS AGAIN
coef.bivBern <- coef(GaussianCopulaModel)
  
#### Predicted distribution parameters: 
PredictedMU1 <- predict(GaussianCopulaModel$mu1, type = "response")
PredictedMU2 <- predict(GaussianCopulaModel$mu2, type = "response")
PredictedRHO <- predict(GaussianCopulaModel$rho, type = "response")

ListOfPredictions <- list(MU1 = PredictedMU1, 
                          MU2 = PredictedMU2,
                          RHO = PredictedRHO)

BivBernoulliGaussianCopulaResults <- vector("list")

BivBernoulliGaussianCopulaResults$MSTOP <- mstop.bivBern
BivBernoulliGaussianCopulaResults$Predictios <- ListOfPredictions
BivBernoulliGaussianCopulaResults$Coefficients <- ListOfCoefficients
BivBernoulliGaussianCopulaResults$CoefficientsAgain <- coef.bivBern

save(BivBernoulliGaussianCopulaResults, file="../BinaryBinary_Application/HeartDiseaseCholesterol_BernoulliGaussCopulaResults_witStabilisationBACKUPJAN19.RData")

length(BivBernoulliGaussianCopulaResults$Coefficients$mu1)
length(BivBernoulliGaussianCopulaResults$Coefficients$mu2)
length(BivBernoulliGaussianCopulaResults$Coefficients$rho)


sum(names(BivBernoulliGaussianCopulaResults$Coefficients$mu1)[-1] %in% names(BivBernoulliGaussianCopulaResults$Coefficients$mu2)[-1]) +
sum(names(BivBernoulliGaussianCopulaResults$Coefficients$mu1)[-1] %in% names(BivBernoulliGaussianCopulaResults$Coefficients$rho)[-1]) +
sum(names(BivBernoulliGaussianCopulaResults$Coefficients$mu2)[-1] %in% names(BivBernoulliGaussianCopulaResults$Coefficients$rho)[-1])

## Between margin 1 and margin 2 
sum(names(BivBernoulliGaussianCopulaResults$Coefficients$mu1)[-1] %in% names(BivBernoulliGaussianCopulaResults$Coefficients$mu2)[-1]) 

## Between margin 1 and dependence
sum(names(BivBernoulliGaussianCopulaResults$Coefficients$mu1)[-1] %in% names(BivBernoulliGaussianCopulaResults$Coefficients$rho)[-1]) 

## Between margin 2 and dependence
sum(names(BivBernoulliGaussianCopulaResults$Coefficients$mu2)[-1] %in% names(BivBernoulliGaussianCopulaResults$Coefficients$rho)[-1])




##### Visualise the coefficients:
CoeffDataset <- data.frame(values = exp(c(abs(BivBernoulliGaussianCopulaResults$Coefficients$mu1),
                                      abs(BivBernoulliGaussianCopulaResults$Coefficients$mu2), 
                                      abs(BivBernoulliGaussianCopulaResults$Coefficients$rho))), 
                           ind = factor(c(names(BivBernoulliGaussianCopulaResults$Coefficients$mu1), 
                                   names(BivBernoulliGaussianCopulaResults$Coefficients$mu2), 
                                   names(BivBernoulliGaussianCopulaResults$Coefficients$rho))))

CoeffDataset$Parameter <- factor(c(rep("mu1", length(BivBernoulliGaussianCopulaResults$Coefficients$mu1)), 
                                   rep("mu2",  length(BivBernoulliGaussianCopulaResults$Coefficients$mu2)), 
                                   rep("rho",  length(BivBernoulliGaussianCopulaResults$Coefficients$rho))), 
                                 levels = c("mu1", "mu2", "rho"),
                                 labels = c("vartheta[1]^(1)", "vartheta[1]^(2)", "vartheta^(c)"))



ggplot(CoeffDataset, aes(x = ind, y = values, col = ind)) + 
  geom_point(size = 1) + 
  #labs(y = expression( exp( "|" beta "|") ) ) +
  labs(x = "Chromosome") +
  facet_grid(Parameter~., labeller = label_parsed, scales = "free_y") +
  theme(legend.position = "none", 
        axis.text.x = element_blank())
  


