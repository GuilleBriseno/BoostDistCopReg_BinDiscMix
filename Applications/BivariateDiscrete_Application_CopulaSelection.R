library(gamlss)
library(gamlss.dist)
#data("ex3.health")
#library(xtable)
library(gamboostLSS)
library(VineCopula)


setwd("BOOSTCOPFILES/")

### Export the dataset to be used in LIDO: 
### margins were ZALG and ZINBI: 

load("Applications/bivdisc_healthdata.RData")


### Now we have to search for a copula:
source("Copulas/BivariateDiscrete/Copula_Gaussian_BivDisc_p1m1.R") # GAUSSIAN 
source("Copulas/BivariateDiscrete/Copula_FGM2_BivDisc_p1m1.R")     # FGM
source("Copulas/BivariateDiscrete/Copula_AMH_BivDisc_p1m1.R")     # AMH
source("Copulas/BivariateDiscrete/Copula_Clayton_BivDisc_p1m1.R")     # Clayton
source("Copulas/BivariateDiscrete/Copula_Gumbel_BivDisc_p1m1.R")     # Gumbel
source("Copulas/BivariateDiscrete/Copula_Frank_BivDisc_p1m1.R")     # Frank
source("Copulas/BivariateDiscrete/Copula_Joe_BivDisc_p1m1.R")     # Joe


health.data$sex <- as.factor(health.data$sex)


biv_respones_eqs <- list(mu1 = formula(cbind(doctorco, prescrib) ~ bols(sex) + bbs(age) + bbs(income)),
                         sigma1 = formula(cbind(doctorco, prescrib) ~ bols(sex) + bbs(age) + bbs(income)),
                         mu2 = formula(cbind(doctorco, prescrib) ~ bols(sex) + bbs(age) + bbs(income)),
                         sigma2 = formula(cbind(doctorco, prescrib) ~ bols(sex) + bbs(age) + bbs(income)),
                         nu2 = formula(cbind(doctorco, prescrib) ~ bols(sex) + bbs(age) + bbs(income)),
                         rho = formula(cbind(doctorco, prescrib) ~ bols(sex) + bbs(age) + bbs(income))
)


boost.nu.short <- 0.1

boost.stabilization <- "L2"


GaussCop_families <- Gauss_Cop_BivDiscrete(marg1 = "ZALG", marg2 = "ZINBI", stabilization = boost.stabilization)
FGMCop_families <- FGM_Cop_BivDiscrete(marg1 = "ZALG", marg2 = "ZINBI", stabilization = boost.stabilization)
AMHCop_families <- AMH_Cop_BivDiscrete(marg1 = "ZALG", marg2 = "ZINBI", stabilization = boost.stabilization)
ClaytonCop_families <- Clayton_Cop_BivDiscrete(marg1 = "ZALG", marg2 = "ZINBI", stabilization = boost.stabilization)
GumbelCop_families <- Gumbel_Cop_BivDiscrete(marg1 = "ZALG", marg2 = "ZINBI", stabilization = boost.stabilization)
FrankCop_families <- Frank_Cop_BivDiscrete(marg1 = "ZALG", marg2 = "ZINBI", stabilization = boost.stabilization)
JoeCop_families <- Joe_Cop_BivDiscrete(marg1 = "ZALG", marg2 = "ZINBI", stabilization = boost.stabilization)


### Get weights: 25% of the sample
nrow(health.data) * 0.25
# 1298

set.seed(1)
weights.boost <- sample(c(rep(1, 3892), rep(0, 1298)))

###

round(mean(health.data$doctorco), digits = 3)
round(sd(health.data$doctorco), digits = 3)

round(mean(health.data$prescrib), digits = 3)
round(sd(health.data$prescrib), digits = 3)

round(mean(health.data$age), digits = 3)
round(sd(health.data$age), digits = 3)

round(mean(health.data$income), digits = 3)
round(sd(health.data$income), digits = 3)

round(mean(as.numeric(health.data$sex)-1), digits = 3)





############################################################################################################# GAUSS COPULA MODEL: 
GaussCopula_Model <- gamboostLSS(biv_respones_eqs, 
                           data = health.data,
                           families = GaussCop_families, 
                           weights = weights.boost,
                           control = boost_control(mstop = 3000, risk = "oobag",
                                                   nu  = boost.nu.short, 
                                                   trace = TRUE), 
                           method = 'noncyclic')

plot(risk(GaussCopula_Model, merge = T))

MSTOP_GaussCop <- which.min(risk(GaussCopula_Model, merge = T))
#### 1651

# GaussCopula_Model <- gamboostLSS(biv_respones_eqs, 
#                                  data = health.data,
#                                  families = GaussCop_families, 
#                                  weights = weights.boost,
#                                  control = boost_control(mstop = MSTOP_GaussCop, 
#                                                          risk = "oobag",
#                                                          nu  = boost.nu.short, 
#                                                          trace = TRUE), 
#                                  method = 'noncyclic')
GaussCopula_Model <- GaussCopula_Model[MSTOP_GaussCop]


#### oobag risk:
FinalRisk_Gauss <- tail(risk(GaussCopula_Model, merge = T), 1)
##### 2251.092 

############################################################################################################# FGM COPULA MODEL: 
FGMCopula_Model <- gamboostLSS(biv_respones_eqs, 
                                 data = health.data,
                                 families = FGMCop_families, 
                                 weights = weights.boost,
                                 control = boost_control(mstop = 3000, risk = "oobag",
                                                         nu  = boost.nu.short, 
                                                         trace = TRUE), 
                                 method = 'noncyclic')

FGMCopula_Model <- FGMCopula_Model[5000]

plot(risk(FGMCopula_Model, merge = T))



MSTOP_FGMCop <- which.min(risk(FGMCopula_Model, merge = T))
#### 3019


# FGMCopula_Model <- gamboostLSS(biv_respones_eqs, 
#                                data = health.data,
#                                families = FGMCop_families, 
#                                weights = weights.boost,
#                                control = boost_control(mstop = MSTOP_FGMCop, risk = "oobag",
#                                                        nu  = boost.nu.short, 
#                                                        trace = TRUE), 
#                                method = 'noncyclic')
FGMCopula_Model <- FGMCopula_Model[MSTOP_FGMCop]

#### oobag risk:
FinalRisk_FGM <- tail(risk(FGMCopula_Model, merge = T), 1)
#####  2270.774 



############################################################################################################# AMH COPULA MODEL: 
AMHCopula_Model <- gamboostLSS(biv_respones_eqs, 
                               data = health.data,
                               families = AMHCop_families, 
                               weights = weights.boost,
                               control = boost_control(mstop = 3000, risk = "oobag",
                                                       nu  = boost.nu.short, 
                                                       trace = TRUE), 
                               method = 'noncyclic')

AMHCopula_Model <- AMHCopula_Model[5000]


plot(risk(AMHCopula_Model, merge = T))



MSTOP_AMHCop <- which.min(risk(AMHCopula_Model, merge = T))
####  2993


# AMHCopula_Model <- gamboostLSS(biv_respones_eqs, 
#                                data = health.data,
#                                families = AMHCop_families, 
#                                weights = weights.boost,
#                                control = boost_control(mstop = MSTOP_AMHCop, risk = "oobag",
#                                                        nu  = boost.nu.short, 
#                                                        trace = TRUE), 
#                                method = 'noncyclic')
AMHCopula_Model <- AMHCopula_Model[MSTOP_AMHCop]

#### oobag risk:
FinalRisk_AMH <- tail(risk(AMHCopula_Model, merge = T), 1)
#####  2268.231 




############################################################################################################# CLAYTON COPULA MODEL: 
ClaytonCopula_Model <- gamboostLSS(biv_respones_eqs, 
                                 data = health.data,
                                 families = ClaytonCop_families, 
                                 weights = weights.boost,
                                 control = boost_control(mstop = 3000, risk = "oobag",
                                                         nu  = boost.nu.short, 
                                                         trace = TRUE), 
                                 method = 'noncyclic')

ClaytonCopula_Model <- ClaytonCopula_Model[5000]

plot(risk(ClaytonCopula_Model, merge = T))

MSTOP_ClaytonCop <- which.min(risk(ClaytonCopula_Model, merge = T))
### 2787 


# ClaytonCopula_Model <- gamboostLSS(biv_respones_eqs, 
#                                    data = health.data,
#                                    families = ClaytonCop_families, 
#                                    weights = weights.boost,
#                                    control = boost_control(mstop = MSTOP_ClaytonCop, 
#                                                            risk = "oobag",
#                                                            nu  = boost.nu.short, 
#                                                            trace = TRUE), 
#                                    method = 'noncyclic')

ClaytonCopula_Model <- ClaytonCopula_Model[MSTOP_ClaytonCop]


#### oobag risk:
FinalRisk_Clayton <- tail(risk(ClaytonCopula_Model, merge = T), 1)
#### 2251.06 


### RE-FIT USING ALL OBSERVATIONS:
ClaytonCopula_Model_FINAL <- gamboostLSS(biv_respones_eqs,
                                   data = health.data,
                                   families = ClaytonCop_families,
                                   control = boost_control(mstop = MSTOP_ClaytonCop,
                                                           nu  = boost.nu.short,
                                                           trace = TRUE),
                                   method = 'noncyclic')


############################################################################################################# GUMBEL COPULA MODEL: 
GumbelCopula_Model <- gamboostLSS(biv_respones_eqs, 
                                 data = health.data,
                                 families = GumbelCop_families, 
                                 weights = weights.boost,
                                 control = boost_control(mstop = 3000, risk = "oobag",
                                                         nu  = boost.nu.short, 
                                                         trace = TRUE), 
                                 method = 'noncyclic')

plot(risk(GumbelCopula_Model, merge = T))

MSTOP_GumbelCop <- which.min(risk(GumbelCopula_Model, merge = T))
### 1651


# GumbelCopula_Model <- gamboostLSS(biv_respones_eqs, 
#                                   data = health.data,
#                                   families = GumbelCop_families, 
#                                   weights = weights.boost,
#                                   control = boost_control(mstop = MSTOP_GumbelCop, 
#                                                           risk = "oobag",
#                                                           nu  = boost.nu.short, 
#                                                           trace = TRUE), 
#                                   method = 'noncyclic')
GumbelCopula_Model <- GumbelCopula_Model[MSTOP_GumbelCop]

#### oobag risk:
FinalRisk_Gumbel <- tail(risk(GumbelCopula_Model, merge = T), 1)
##### 2255.1 

############################################################################################################# FRANK COPULA MODEL: 
FrankCopula_Model <- gamboostLSS(biv_respones_eqs, 
                                 data = health.data,
                                 families = FrankCop_families, 
                                 weights = weights.boost,
                                 control = boost_control(mstop = 3000, risk = "oobag",
                                                         nu  = boost.nu.short, 
                                                         trace = TRUE), 
                                 method = 'noncyclic')

FrankCopula_Model <- FrankCopula_Model[18000]


plot(risk(FrankCopula_Model, merge = T))

MSTOP_FrankCop <- which.min(risk(FrankCopula_Model, merge = T))
#### 2222


FrankCopula_Model <- gamboostLSS(biv_respones_eqs, 
                                 data = health.data,
                                 families = FrankCop_families, 
                                 weights = weights.boost,
                                 control = boost_control(mstop = MSTOP_FrankCop, 
                                                         risk = "oobag",
                                                         nu  = boost.nu.short, 
                                                         trace = TRUE), 
                                 method = 'noncyclic')

#### oobag risk:
FinalRisk_Frank <- tail(risk(FrankCopula_Model, merge = T), 1)
###### 2248.893


############################################################################################################# JOE COPULA MODEL: 
JoeCopula_Model <- gamboostLSS(biv_respones_eqs, 
                                 data = health.data,
                                 families = JoeCop_families, 
                                 weights = weights.boost,
                                 control = boost_control(mstop = 3000, risk = "oobag",
                                                         nu  = boost.nu.short, 
                                                         trace = TRUE), 
                                 method = 'noncyclic')

JoeCopula_Model <- JoeCopula_Model[5000]


plot(risk(JoeCopula_Model, merge = T))

MSTOP_JoeCop <- which.min(risk(JoeCopula_Model, merge = T))
#### 1743


# JoeCopula_Model <- gamboostLSS(biv_respones_eqs, 
#                                  data = health.data,
#                                  families = JoeCop_families, 
#                                  weights = weights.boost,
#                                  control = boost_control(mstop = MSTOP_FrankCop, 
#                                                          risk = "oobag",
#                                                          nu  = boost.nu.short, 
#                                                          trace = TRUE), 
#                                  method = 'noncyclic')

JoeCopula_Model <- JoeCopula_Model[MSTOP_JoeCop]


#### oobag risk:
FinalRisk_Joe <- tail(risk(JoeCopula_Model, merge = T), 1)
###### 2260.711 


########### Models according to risk: ----- > LATEST VERSION (02.12.2023)

# 1 - Clayton ---- MSTOP ---------- 3026  --------------- 2249.322 
# 2 - Gauss ------ MSTOP ---------- 1801  --------------- 2249.372 
# 3 - Frank ------ MSTOP ---------- 17994 --------------- 2252.313  (no optimal mstop found yet at 18000)
# 4 - Gumbel ----- MSTOP ---------- 1711  --------------- 2253.522
# 5 - Joe -------- MSTOP ---------- 1743  --------------- 2260.711 
# 6 - AMH -------- MSTOP----------- 3358  --------------- 2266.521
# 7 - FGM -------- MSTOP----------- 3430  --------------- 2269.052

# 8 - Bivariate Poisson ---- 2448.313


# ---------------- Summary of the clayton model: 

# Number of boosting iterations (mstop):  mu1 = 592, sigma1 = 491, mu2 = 285, sigma2 = 510, nu2 = 438, rho = 710 
# Step size:  mu1 = 0.1, sigma1 = 0.1, mu2 = 0.1, sigma2 = 0.1, nu2 = 0.1, rho = 0.1 
# 
# Zero-altered Logarithmic distribution
# 
# Zero-Inflated Negative Binomial distribution
# 
# Clayton copula for bivariate discrete / count responses: rho(... link) 
# 
# Selection frequencies:
#   Parameter mu1:
#   bbs(age) bbs(income) 
# 0.6334459   0.3665541 
# Parameter sigma1:
#   bbs(income)    bbs(age)   bols(sex) 
# 0.4541752   0.3625255   0.1832994 
# Parameter mu2:
#   bbs(age) bbs(income)   bols(sex) 
# 0.4526316   0.3789474   0.1684211 
# Parameter sigma2:
#   bbs(age)   bols(sex) bbs(income) 
# 0.5470588   0.3470588   0.1058824 
# Parameter nu2:
#   bols(sex)    bbs(age) bbs(income) 
# 0.6095890   0.2031963   0.1872146 
# Parameter rho:
#   bbs(age)   bols(sex) bbs(income) 
# 0.4859155   0.3788732   0.1352113


# > coef(ClaytonCopula_Model$sigma1)
# $`bols(sex)`
# (Intercept)        sex1 
# 0.5948705  -0.2557150


# > coef(ClaytonCopula_Model$mu2)
# $`bols(sex)`
# (Intercept)        sex1 
# -0.1653868   0.2718806 

# > coef(ClaytonCopula_Model$sigma2)
# $`bols(sex)`
# (Intercept)        sex1 
# -0.2003583  -0.6865052 

# > coef(ClaytonCopula_Model$nu2)
# $`bols(sex)`
# (Intercept)        sex1 
# 0.04821361 -1.16557823 

# > coef(ClaytonCopula_Model$rho)
# $`bols(sex)`
# (Intercept)        sex1 
# 0.3262854  -0.3286593 

par(mfrow = c(1, 2))
plot(ClaytonCopula_Model_FINAL$mu1)
par(mfrow = c(1, 3))
plot(ClaytonCopula_Model_FINAL$sigma1)

plot(ClaytonCopula_Model_FINAL$mu2)
plot(ClaytonCopula_Model_FINAL$sigma2)
plot(ClaytonCopula_Model_FINAL$nu2)

plot(ClaytonCopula_Model_FINAL$rho)


#### Re-fit clayton without weighted observations until optimal mstop:
ClaytonCopula_Model_FULL <- gamboostLSS(biv_respones_eqs,
                                   data = health.data,
                                   families = ClaytonCop_families,
                                   control = boost_control(mstop = MSTOP_ClaytonCop,
                                                           nu  = boost.nu.short,
                                                           trace = TRUE),
                                   method = 'noncyclic')

summary(ClaytonCopula_Model_FINAL)

# Selection frequencies:
#   Parameter mu1:
#   bbs(age) bbs(income) 
# 0.60199     0.39801 
# Parameter sigma1:
#   bbs(age) bbs(income)   bols(sex) 
# 0.4449761   0.3181818   0.2368421 
# Parameter mu2:
#   bbs(age) bbs(income)   bols(sex) 
# 0.5680934   0.2529183   0.1789883 
# Parameter sigma2:
#   bbs(age)   bols(sex) bbs(income) 
# 0.5599284   0.2826476   0.1574240 
# Parameter nu2:
#   bols(sex)    bbs(age) bbs(income) 
# 0.70702179  0.28087167  0.01210654 
# Parameter rho:
#   bbs(age) bols(sex) 
# 0.564433  0.435567 


# > coef(ClaytonCopula_Model_FULL$sigma1)
# $`bols(sex)`
# (Intercept)        sex1 
# 0.7867547  -0.2288249

# > coef(ClaytonCopula_Model_FULL$mu2)
# $`bols(sex)`
# (Intercept)        sex1 
# -0.1695536   0.2720868

# > coef(ClaytonCopula_Model_FULL$sigma2)
# $`bols(sex)`
# (Intercept)        sex1 
# -0.1811186  -0.6261954 

# > coef(ClaytonCopula_Model_FULL$nu2)
# $`bols(sex)`
# (Intercept)        sex1 
# 0.06266782 -1.18453106

# > coef(ClaytonCopula_Model_FULL$rho)
# $`bols(sex)`
# (Intercept)        sex1 
# 0.4417316  -0.3901140 



par(mfrow = c(1, 2))
plot(ClaytonCopula_Model_FULL$mu1)
par(mfrow = c(1, 3))
plot(ClaytonCopula_Model_FULL$sigma1)

plot(ClaytonCopula_Model_FULL$mu2)
plot(ClaytonCopula_Model_FULL$sigma2)
plot(ClaytonCopula_Model_FULL$nu2)

plot(ClaytonCopula_Model_FULL$rho)

### Optimal model is clayton copula: Save it, as well as its mstop:
MSTOP_ClaytonCop

FinalRisk_Clayton

ClaytonCopula_Model_FULL
ClaytonCopula_Model

weights.boost # for reproducibility

save(MSTOP_ClaytonCop, FinalRisk_Clayton, 
     ClaytonCopula_Model, ClaytonCopula_Model_FULL, 
     weights.boost, 
     file = "Applications/BivariateDiscrete_ClaytonCopulaModels.RData")



# ---- Stuff for visuals:
# Age is age is years divided by 100.  
# Income: Income in AUS divided by 1000. 
AgeNew <- seq(range(health.data$age)[1], range(health.data$age)[2], length.out = 500)
IncomeNew <- seq(range(health.data$income)[1], range(health.data$income)[2], length.out = 500)

AgeReal <- AgeNew * 100
IncomeReal <- IncomeNew * 1000

NewData <- data.frame(sex = factor(rep(0, 500)), 
                      income = IncomeNew, 
                      age = AgeNew)


PredictionsMU1 <- list(Age = predict(ClaytonCopula_Model$mu1, type = "link", which = "age", newdata = NewData), 
                       Income = predict(ClaytonCopula_Model$mu1, type = "link", which = "income", newdata = NewData))
PredictionsSIGMA1 <- list(Age = predict(ClaytonCopula_Model$sigma1, type = "link", which = "age", newdata = NewData), 
                          Income = predict(ClaytonCopula_Model$sigma1, type = "link", which = "income", newdata = NewData))

PredictionsMU2 <- list(Age = predict(ClaytonCopula_Model$mu2, type = "link", which = "age", newdata = NewData), 
                       Income = predict(ClaytonCopula_Model$mu2, type = "link", which = "income", newdata = NewData))
PredictionsSIGMA2 <- list(Age = predict(ClaytonCopula_Model$sigma2, type = "link", which = "age", newdata = NewData), 
                         Income = predict(ClaytonCopula_Model$sigma2, type = "link", which = "income", newdata = NewData))
PredictionsNU2 <- list(Age = predict(ClaytonCopula_Model$nu2, type = "link", which = "age", newdata = NewData), 
                       Income = predict(ClaytonCopula_Model$nu2, type = "link", which = "income", newdata = NewData))

PredictionsRHO <- list(Age = predict(ClaytonCopula_Model$rho, type = "link", which = "age", newdata = NewData), 
                       Income = predict(ClaytonCopula_Model$rho, type = "link", which = "income", newdata = NewData))

# Grid of parmeters and covariates: 
AGE_PREDS <- data.frame(MU1 = PredictionsMU1$Age, 
                        SIGMA1 = PredictionsSIGMA1$Age,
                        MU2 = PredictionsMU2$Age, 
                        SIGMA2 = PredictionsSIGMA2$Age,
                        NU = PredictionsNU2$Age,
                        RHO = PredictionsRHO$Age)

colnames(AGE_PREDS) <- c("mu1", "sigma1", "mu2", "sigma2", "nu2", "theta")

INCOME_PREDS <- data.frame(MU1 = PredictionsMU1$Income, 
                           SIGMA1 = PredictionsSIGMA1$Income,
                           MU2 = PredictionsMU2$Income, 
                           SIGMA2 = PredictionsSIGMA2$Income,
                           NU = PredictionsNU2$Income,
                           RHO = PredictionsRHO$Income)

colnames(INCOME_PREDS) <- c("mu1", "sigma1", "mu2", "sigma2", "nu2", "theta")

AGE_PREDS <- stack(AGE_PREDS)
INCOME_PREDS <- stack(INCOME_PREDS)

AGE_PREDS$Parameter <- factor(AGE_PREDS$ind, 
                        levels = c("mu1", "sigma1", "mu2", "sigma2", "nu2", "theta"), 
                        labels = c("vartheta[1]^(1)", "vartheta[2]^(1)",
                                   "vartheta[1]^(2)", "vartheta[2]^(2)", "vartheta[3]^(2)", 
                                   "vartheta^(c)"), 
                        order = T
                        )

INCOME_PREDS$Parameter <- factor(INCOME_PREDS$ind, 
                              levels = c("mu1", "sigma1", "mu2", "sigma2", "nu2", "theta"), 
                              labels = c("vartheta[1]^(1)", "vartheta[2]^(1)",
                                         "vartheta[1]^(2)", "vartheta[2]^(2)", "vartheta[3]^(2)", 
                                         "vartheta^(c)"), 
                              order = T
)

AGE_PREDS$Margin <- factor(c(rep("Margin 1", 2 * nrow(PredictionsMU1$Age)), 
                                rep("Margin 2", 3 * nrow(PredictionsMU1$Age)),
                                rep("Copula parameter", 1 * nrow(PredictionsMU1$Age))
), 
labels = c("Margin 1", "Margin 2", "Copula parameter"),
order = T)

INCOME_PREDS$Margin <- factor(c(rep("Margin 1", 2 * nrow(PredictionsMU1$Age)), 
                                rep("Margin 2", 3 * nrow(PredictionsMU1$Age)),
                                rep("Copula parameter", 1 * nrow(PredictionsMU1$Age))
                                ), 
                              labels = c("Margin 1", "Margin 2", "Copula parameter"),
                              order = T)



AGE_PREDS$Age <- rep(AgeReal, 6)
INCOME_PREDS$Income <- rep(IncomeReal, 6)

library(ggplot2)
PLOTAGE <- ggplot(AGE_PREDS, aes(x = Age, y = values)) + 
  geom_line() + 
  geom_hline(yintercept = 0, col = "darkgrey") +
  labs(x = "Age (years)", y = expression(hat(f)(Age)), title = "(a)") + 
  facet_grid( ~ Parameter, labeller = labeller(Parameter = label_parsed)) +
  theme_light() +
  theme( strip.text.x = element_text(size = 20), 
         axis.text = element_text(size = 8),
         plot.title = element_text(hjust = 0.5, size = 23))


library(ggplot2)
PLOTINCOME <- ggplot(INCOME_PREDS, aes(x = Income, y = values)) + 
  geom_line() + 
  geom_hline(yintercept = 0, col = "darkgrey") +
  labs(x = "Income (AUD)", y = expression(hat(f)(Income)), title = "(b)") + 
  facet_grid( ~ Parameter, labeller = labeller(Parameter = label_parsed)) +
  theme_light() +
  theme( strip.text.x = element_text(size = 20), 
         axis.text = element_text(size = 8), 
         plot.title = element_text(hjust = 0.5, size = 23))


library(patchwork)

PLOTAGE + PLOTINCOME + plot_layout(nrow = 2)


BiCopPar2Tau(family = 3, range(predict(ClaytonCopula_Model_FINAL$rho, type = "response")))

round(BiCopPar2Tau(family = 3, range(predict(ClaytonCopula_Model_FINAL$rho, type = "response"))), digits = 3)


########### Models according to risk: ----- > OLD VERSION

# 1 - Frank ---------------- 2248.893
# 2 - Clayton -------------- 2251.060
# 3 - Gauss ---------------- 2251.092
# 4 - Gumbel --------------- 2255.100
# 5 - Joe ------------------ 2260.711 
# 6 - AMH ------------------ 2268.231 (?)
# 7 - FGM ------------------ 2270.774
# 8 - Bivariate Poisson ---- 2448.313

##### Frank copula model:
summary(FrankCopula_Model)

# Selection frequencies:
      #   Parameter mu1:
#   bbs(age) bbs(income) 
# 0.6225166   0.3774834 


      # Parameter sigma1:
#   bbs(income)    bbs(age)   bols(sex) 
# 0.4275362   0.3876812   0.1847826 


        # Parameter mu2:
#   bbs(age) bbs(income)   bols(sex) 
# 0.5054945   0.2747253   0.2197802 


      # Parameter sigma2:
#   bols(sex)  bbs(age) 
# 0.6051873 0.3948127 


#         Parameter nu2:
#   bols(sex)  bbs(age) 
# 0.7427386 0.2572614 


#         Parameter rho:
#   bbs(age) bbs(income)   bols(sex) 
# 0.67911480  0.23928077  0.08160443 


### plot final model
par(mfrow = c(1, 2))
plot(FrankCopula_Model$mu1)


par(mfrow = c(1, 3))
plot(FrankCopula_Model$sigma1)


par(mfrow = c(1, 3))
plot(FrankCopula_Model$mu2)



par(mfrow = c(1, 2))
plot(FrankCopula_Model$sigma2)


par(mfrow = c(1, 2))
plot(FrankCopula_Model$nu2)


par(mfrow = c(1, 3))
plot(FrankCopula_Model$rho)


#### Coefficients:
# Margin 1 (Zero altered logarithmic)
coef(FrankCopula_Model$mu1)
coef(FrankCopula_Model$sigma1)

# Margin 2 (Zero inflated negative binomial)
coef(FrankCopula_Model$mu2)
coef(FrankCopula_Model$sigma2)
coef(FrankCopula_Model$nu2)

# Copula parameter
coef(FrankCopula_Model$rho)


save(MSTOP_ClaytonCop, MSTOP_FrankCop, MSTOP_GumbelCop, MSTOP_GaussCop, MSTOP_FGMCop, 
     FrankCopula_Model, ClaytonCopula_Model,
     weights.boost,
     file = "Applications/BivCount_Application_CopulaModels.RData")




############# Apply this to the bivariate poisson loss function: 

source("Copulas/RigidBivariateDistributions/bivariatePoisson.R")

### define the family (distribution has no stabilization)
BivariatePoissonFamilies <- PoissonBV()

biv_responses_eqs_BIVPOIS <- biv_respones_eqs$mu1

bivPoisson_Model <- gamboostLSS(biv_responses_eqs_BIVPOIS, 
                                data = health.data, 
                                families = BivariatePoissonFamilies, 
                                weights = weights.boost, 
                                control = boost_control(mstop = 5000, 
                                                        risk = "oobag",
                                                        nu  = boost.nu.short, 
                                                        trace = TRUE), 
                                method = 'noncyclic'
                                )


plot(risk(bivPoisson_Model, merge = T))

MSTOP_BivariatePoisson <- which.min(risk(bivPoisson_Model, merge = T))
#### 


bivPoisson_Model <- gamboostLSS(biv_responses_eqs_BIVPOIS, 
                                data = health.data, 
                                families = BivariatePoissonFamilies, 
                                weights = weights.boost, 
                                control = boost_control(mstop = MSTOP_BivariatePoisson, 
                                                        risk = "oobag",
                                                        nu  = boost.nu.short, 
                                                        trace = TRUE), 
                                method = 'noncyclic'
)

#### oobag risk:
FinalRisk_BivariatePoisson <- tail(risk(bivPoisson_Model, merge = T), 1)
######  2448.313

