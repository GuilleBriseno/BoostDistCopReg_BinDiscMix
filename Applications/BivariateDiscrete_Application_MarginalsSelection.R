### Bivariate Count application: 

### Distribution search using boosting: We find mstop using oobag risk and nval = 25% of the sample. 

load("bivdisc_healthdata.RData")
library(gamlss.dist)
library(gamboostLSS)


health.data$sex <- as.factor(health.data$sex)


boost.nu <- 0.1

boost.stabilization <- "L2"

### MARGIN 1: doctorco

### 1 parameter distribution
param1_eqs <- doctorco ~ bols(sex) + bbs(age) + bbs(income)

### 2 parameter distribution
param2_eqs <- list(mu = doctorco ~ bols(sex) + bbs(age) + bbs(income), 
                   sigma = doctorco ~ bols(sex) + bbs(age) + bbs(income))

### 3 parameter distribution
param3_eqs <- list(mu = doctorco ~ bols(sex) + bbs(age) + bbs(income), 
                   sigma = doctorco ~ bols(sex) + bbs(age) + bbs(income),
                   nu = doctorco ~ bols(sex) + bbs(age) + bbs(income))

### Get weights: 25% of the sample
nrow(health.data) * 0.25
# 1298

set.seed(1)
weights.boost <- sample(c(rep(1, 3892), rep(0, 1298)))

############################################################################################################### POISSON:
DOCTORCO_POISSON <- gamboost(formula = param1_eqs, 
                             family = Poisson(), 
                             data = health.data,
                             weights = weights.boost,
                             control = boost_control(mstop = 5000, nu = boost.nu, risk = "oobag", trace = TRUE))

plot(risk(DOCTORCO_POISSON))

### optimal mstop for POISSON DISTRIBUTION
MSTOP_DOCTORCO_POISSON <- which.min(risk(DOCTORCO_POISSON))
# 121

### Re-fit until mstop:
DOCTORCO_POISSON <- gamboost(formula = param1_eqs, 
                             family = Poisson(), 
                             data = health.data,
                             weights = weights.boost,
                             control = boost_control(mstop = MSTOP_DOCTORCO_POISSON, nu = boost.nu, risk = "oobag", trace = TRUE))


### evaluate out of bag risk:
FinalRisk_POISSON <- tail(risk(DOCTORCO_POISSON), 1)
#   947.8327



############################################################################################################### GEOMETRIC
DOCTORCO_GEOM <- gamboost(param1_eqs, 
                          family = as.families(fname = "GEOM"), 
                          data = health.data,
                          weights = weights.boost,
                          control = boost_control(mstop = 5000, nu = boost.nu, risk = "oobag", trace = TRUE))

plot(risk(DOCTORCO_GEOM))

### optimal mstop for  DISTRIBUTION
MSTOP_DOCTORCO_GEOM <- which.min(risk(DOCTORCO_GEOM))
# 151

DOCTORCO_GEOM <- gamboost(param1_eqs, 
                          family = as.families(fname = "GEOM"), 
                          data = health.data,
                          weights = weights.boost,
                          control = boost_control(mstop = MSTOP_DOCTORCO_GEOM, nu = boost.nu, risk = "oobag", trace = TRUE))

### evaluate out of bag risk:
FinalRisk_GEOM <- tail(risk(DOCTORCO_GEOM), 1)
##  890.1364


################################################################################################################# ZERO-INFLATED POISSON
DOCTORCO_ZIP <- gamboostLSS(param2_eqs, 
                            families = as.families(fname = "ZIP", stabilization = boost.stabilization), 
                            data = health.data,
                            weights = weights.boost,
                            control = boost_control(mstop = 10000, nu = boost.nu, risk = "oobag", trace = TRUE), 
                            method = "noncyclic")

plot(risk(DOCTORCO_ZIP, merge = T))

### optimal mstop for  DISTRIBUTION
MSTOP_DOCTORCO_ZIP <- which.min(risk(DOCTORCO_ZIP, merge = T))
# 4982 -----> extend 
# 5532 # not stabilised


### STABILISED MSTOP: 1884 


DOCTORCO_ZIP <- gamboostLSS(param2_eqs, 
                            families = as.families(fname = "ZIP", stabilization = boost.stabilization), 
                            data = health.data,
                            weights = weights.boost,
                            control = boost_control(mstop = MSTOP_DOCTORCO_ZIP, nu = boost.nu, risk = "oobag", trace = TRUE), 
                            method = "noncyclic")

### evaluate out of bag risk:
FinalRisk_ZIP <- tail(risk(DOCTORCO_ZIP, merge = T), 1)
##  911.19

### stabilised: 913.053

################################################################################################################ ZERO ALTERED LOGARITHMIC
DOCTORCO_ZALG <- gamboostLSS(param2_eqs, 
                             families = as.families(fname = "ZALG", stabilization = boost.stabilization), 
                             data = health.data,
                             weights = weights.boost,
                             control = boost_control(mstop = 5000, nu = boost.nu, risk = "oobag", trace = TRUE), 
                             method = "noncyclic")

plot(risk(DOCTORCO_ZALG, merge = T))

### optimal mstop  DISTRIBUTION
MSTOP_DOCTORCO_ZALG <- which.min(risk(DOCTORCO_ZALG, merge = T))
# 193

### stabilised: 84


DOCTORCO_ZALG <- gamboostLSS(param2_eqs, 
                             families = as.families(fname = "ZALG", stabilization = boost.stabilization), 
                             data = health.data,
                             weights = weights.boost,
                             control = boost_control(mstop = MSTOP_DOCTORCO_ZALG, nu = boost.nu, risk = "oobag", trace = TRUE), 
                             method = "noncyclic")

### evaluate out of bag risk:
FinalRisk_ZALG <- tail(risk(DOCTORCO_ZALG, merge = T), 1)
##  877.5308 

### stabilised: 877.528


################################################################################################################ NEGATIVE BINOMIAL
DOCTORCO_NBI <- gamboostLSS(param2_eqs, 
                            families = as.families(fname = "NBI", stabilization = boost.stabilization), 
                            data = health.data,
                            weights = weights.boost,
                            control = boost_control(mstop = 5000, nu = boost.nu, risk = "oobag", trace = TRUE), 
                            method = "noncyclic")

plot(risk(DOCTORCO_NBI, merge = T))


### optimal mstop for  DISTRIBUTION
MSTOP_DOCTORCO_NBI <- which.min(risk(DOCTORCO_NBI, merge = T))
# 1033 

### stabilised: 346


DOCTORCO_NBI <- gamboostLSS(param2_eqs, 
                            families = as.families(fname = "NBI", stabilization = boost.stabilization), 
                            data = health.data,
                            weights = weights.boost,
                            control = boost_control(mstop = MSTOP_DOCTORCO_NBI, nu = boost.nu, risk = "oobag", trace = TRUE), 
                            method = "noncyclic")

### evaluate out of bag risk:
FinalRisk_NBI <- tail(risk(DOCTORCO_NBI, merge = T), 1)
##  884.3616


### stabilised: 884.0531



################################################################################################################ ZERO-INFLATED NEGATIVE BINOMIAL
DOCTORCO_ZINBI <- gamboostLSS(param3_eqs, 
                              families = as.families(fname = "ZINBI", stabilization = boost.stabilization), 
                              data = health.data,
                              weights = weights.boost,
                              control = boost_control(mstop = 5000, nu = boost.nu, risk = "oobag", trace = TRUE), 
                              method = "noncyclic")

plot(risk(DOCTORCO_ZINBI, merge = T))

### optimal mstop for POISSON DISTRIBUTION
MSTOP_DOCTORCO_ZINBI <- which.min(risk(DOCTORCO_ZINBI, merge = T))
# 2861 

#### stabilised: 293


######
DOCTORCO_ZINBI <- gamboostLSS(param3_eqs, 
                              families = as.families(fname = "ZINBI", stabilization = boost.stabilization), 
                              data = health.data,
                              weights = weights.boost,
                              control = boost_control(mstop = MSTOP_DOCTORCO_ZINBI, nu = boost.nu, risk = "oobag", trace = TRUE), 
                              method = "noncyclic")

### evaluate out of bag risk:
FinalRisk_ZINBI <- tail(risk(DOCTORCO_ZINBI, merge = T), 1)
##  887.8697


### ### stabilised: 887.7289



################################################################################################################ ZERO-ALTERED NEGATIVE BINOMIAL
DOCTORCO_ZANBI <- gamboostLSS(param3_eqs, 
                              families = as.families(fname = "ZANBI", stabilization = boost.stabilization), 
                              data = health.data,
                              weights = weights.boost,
                              control = boost_control(mstop = 5000, nu = boost.nu, risk = "oobag", trace = TRUE), 
                              method = "noncyclic")

plot(risk(DOCTORCO_ZANBI, merge = T))

### optimal mstop for POISSON DISTRIBUTION
MSTOP_DOCTORCO_ZANBI <- which.min(risk(DOCTORCO_ZANBI, merge = T))
# 876

#### stabilised: 284


DOCTORCO_ZANBI <- gamboostLSS(param3_eqs, 
                              families = as.families(fname = "ZANBI", stabilization = boost.stabilization), 
                              data = health.data,
                              weights = weights.boost,
                              control = boost_control(mstop = MSTOP_DOCTORCO_ZANBI, nu = boost.nu, risk = "oobag", trace = TRUE), 
                              method = "noncyclic")

### evaluate out of bag risk:
FinalRisk_ZANBI <- tail(risk(DOCTORCO_ZANBI, merge = T), 1)
##  883.8037

### stabilised: 882.5501


c(FinalRisk_POISSON, FinalRisk_GEOM, FinalRisk_ZALG, FinalRisk_ZIP, FinalRisk_NBI, FinalRisk_ZINBI, FinalRisk_ZANBI)


##########################################################################################
####### Out of bag risks in order (smallest to largest): ( stabilization)
##########################################################################################
# 1 ZERO ALTERED LOGARITHMIC (best)
# 2 ZERO ALTERED NEGATIVE BINOMIAL                    
# 3 NEGATIVE BINOMIAL                             
# 4 ZERO INFLATED NEGATIVE BINOMIAL                 
# 5 GEOMETRIC
# 6 ZERO-INFLATED POISSON 
# 7 POISSON (worst)

##########################################################################################










#####################################################################################################
##################################################################################################### SECOND MARGIN
#####################################################################################################  PRESCRIB
#####################################################################################################



### MARGIN 2: prescrib
### 1 parameter distribution
param1_eqs <- prescrib ~ bols(sex) + bbs(age) + bbs(income)

### 2 parameter distribution
param2_eqs <- list(mu = prescrib ~ bols(sex) + bbs(age) + bbs(income), 
                   sigma = prescrib ~ bols(sex) + bbs(age) + bbs(income))

### 3 parameter distribution
param3_eqs <- list(mu = prescrib ~ bols(sex) + bbs(age) + bbs(income), 
                   sigma = prescrib ~ bols(sex) + bbs(age) + bbs(income),
                   nu = prescrib ~ bols(sex) + bbs(age) + bbs(income))

############################################################################################################### POISSON:
PRESCRIB_POISSON <- gamboost(formula = param1_eqs, 
                             family = Poisson(), 
                             data = health.data,
                             weights = weights.boost,
                             control = boost_control(mstop = 5000, nu = boost.nu, risk = "oobag", trace = TRUE))

plot(risk(PRESCRIB_POISSON))

### optimal mstop for POISSON DISTRIBUTION
MSTOP_PRESCRIB_POISSON <- which.min(risk(PRESCRIB_POISSON))
# 257

### Re-fit until mstop:
PRESCRIB_POISSON <- gamboost(formula = param1_eqs, 
                             family = Poisson(), 
                             data = health.data,
                             weights = weights.boost,
                             control = boost_control(mstop = MSTOP_PRESCRIB_POISSON, nu = boost.nu, risk = "oobag", trace = TRUE))


### evaluate out of bag risk:
tail(risk(PRESCRIB_POISSON), 1)
#   1580.917



############################################################################################################### GEOMETRIC
PRESCRIB_GEOM <- gamboost(param1_eqs, 
                          family = as.families(fname = "GEOM"), 
                          data = health.data,
                          weights = weights.boost,
                          control = boost_control(mstop = 5000, nu = boost.nu, risk = "oobag", trace = TRUE))

plot(risk(PRESCRIB_GEOM))

### optimal mstop for  DISTRIBUTION
MSTOP_PRESCRIB_GEOM <- which.min(risk(PRESCRIB_GEOM))
# 209

PRESCRIB_GEOM <- gamboost(param1_eqs, 
                          family = as.families(fname = "GEOM"), 
                          data = health.data,
                          weights = weights.boost,
                          control = boost_control(mstop = MSTOP_PRESCRIB_GEOM, nu = boost.nu, risk = "oobag", trace = TRUE))

### evaluate out of bag risk:
tail(risk(PRESCRIB_GEOM), 1)
##  1461.826


################################################################################################################# ZERO-INFLATED POISSON
PRESCRIB_ZIP <- gamboostLSS(param2_eqs, 
                            families = as.families(fname = "ZIP", stabilization = boost.stabilization), 
                            data = health.data,
                            weights = weights.boost,
                            control = boost_control(mstop = 5000, nu = boost.nu, risk = "oobag", trace = TRUE), 
                            method = "noncyclic")

plot(risk(PRESCRIB_ZIP, merge = T))

### optimal mstop for  DISTRIBUTION
MSTOP_PRESCRIB_ZIP <- which.min(risk(PRESCRIB_ZIP, merge = T))
# 1201 


PRESCRIB_ZIP <- gamboostLSS(param2_eqs, 
                            families = as.families(fname = "ZIP", stabilization = boost.stabilization), 
                            data = health.data,
                            weights = weights.boost,
                            control = boost_control(mstop = MSTOP_PRESCRIB_ZIP, nu = boost.nu, risk = "oobag", trace = TRUE), 
                            method = "noncyclic")

### evaluate out of bag risk:
tail(risk(PRESCRIB_ZIP, merge = T), 1)
## 1466.793



################################################################################################################ ZERO ALTERED LOGARITHMIC
PRESCRIB_ZALG <- gamboostLSS(param2_eqs, 
                             families = as.families(fname = "ZALG", stabilization = boost.stabilization), 
                             data = health.data,
                             weights = weights.boost,
                             control = boost_control(mstop = 5000, nu = boost.nu, risk = "oobag", trace = TRUE), 
                             method = "noncyclic")

plot(risk(PRESCRIB_ZALG, merge = T))

### optimal mstop  DISTRIBUTION
MSTOP_PRESCRIB_ZALG <- which.min(risk(PRESCRIB_ZALG, merge = T))
# 440


PRESCRIB_ZALG <- gamboostLSS(param2_eqs, 
                             families = as.families(fname = "ZALG", stabilization = boost.stabilization), 
                             data = health.data,
                             weights = weights.boost,
                             control = boost_control(mstop = MSTOP_PRESCRIB_ZALG, nu = boost.nu, risk = "oobag", trace = TRUE), 
                             method = "noncyclic")

### evaluate out of bag risk:
tail(risk(PRESCRIB_ZALG, merge = T), 1)
##  1471.163 


################################################################################################################ NEGATIVE BINOMIAL
PRESCRIB_NBI <- gamboostLSS(param2_eqs, 
                            families = as.families(fname = "NBI", stabilization = boost.stabilization), 
                            data = health.data,
                            weights = weights.boost,
                            control = boost_control(mstop = 5000, nu = boost.nu, risk = "oobag", trace = TRUE), 
                            method = "noncyclic")

plot(risk(PRESCRIB_NBI, merge = T))


### optimal mstop for POISSON DISTRIBUTION
MSTOP_PRESCRIB_NBI <- which.min(risk(PRESCRIB_NBI, merge = T))
# 627

PRESCRIB_NBI <- gamboostLSS(param2_eqs, 
                            families = as.families(fname = "NBI", stabilization = boost.stabilization), 
                            data = health.data,
                            weights = weights.boost,
                            control = boost_control(mstop = MSTOP_PRESCRIB_NBI, nu = boost.nu, risk = "oobag", trace = TRUE), 
                            method = "noncyclic")

### evaluate out of bag risk:
tail(risk(PRESCRIB_NBI, merge = T), 1)
##  1452.925 



################################################################################################################ ZERO-INFLATED NEGATIVE BINOMIAL
PRESCRIB_ZINBI <- gamboostLSS(param3_eqs, 
                              families = as.families(fname = "ZINBI", stabilization = boost.stabilization), 
                              data = health.data,
                              weights = weights.boost,
                              control = boost_control(mstop = 5000, nu = boost.nu, risk = "oobag", trace = TRUE), 
                              method = "noncyclic")

PRESCRIB_ZINBI <- PRESCRIB_ZINBI[6000]


plot(risk(PRESCRIB_ZINBI, merge = T))

### optimal mstop for POISSON DISTRIBUTION
MSTOP_PRESCRIB_ZINBI <- which.min(risk(PRESCRIB_ZINBI, merge = T))
# 3111


######
PRESCRIB_ZINBI <- gamboostLSS(param3_eqs, 
                              families = as.families(fname = "ZINBI", stabilization = boost.stabilization), 
                              data = health.data,
                              weights = weights.boost,
                              control = boost_control(mstop = MSTOP_PRESCRIB_ZINBI, nu = boost.nu, risk = "oobag", trace = TRUE), 
                              method = "noncyclic")

### evaluate out of bag risk:
tail(risk(PRESCRIB_ZINBI, merge = T), 1)
## 1442.541 




################################################################################################################ ZERO-ALTERED NEGATIVE BINOMIAL
PRESCRIB_ZANBI <- gamboostLSS(param3_eqs, 
                              families = as.families(fname = "ZANBI", stabilization = boost.stabilization), 
                              data = health.data,
                              weights = weights.boost,
                              control = boost_control(mstop = 5000, nu = boost.nu, risk = "oobag", trace = TRUE), 
                              method = "noncyclic")

plot(risk(PRESCRIB_ZANBI, merge = T))

### optimal mstop for POISSON DISTRIBUTION
MSTOP_PRESCRIB_ZANBI <- which.min(risk(PRESCRIB_ZANBI, merge = T))
# 455

PRESCRIB_ZANBI <- gamboostLSS(param3_eqs, 
                              families = as.families(fname = "ZANBI", stabilization = boost.stabilization), 
                              data = health.data,
                              weights = weights.boost,
                              control = boost_control(mstop = MSTOP_PRESCRIB_ZANBI, nu = boost.nu, risk = "oobag", trace = TRUE), 
                              method = "noncyclic")

### evaluate out of bag risk:
tail(risk(PRESCRIB_ZANBI, merge = T), 1)
##  1443.229 




##########################################################################################
####### Out of bag risks in order (smallest to largest): ( stabilization)
##########################################################################################
# 1 ZERO INFLATED NEGATIVE BINOMIAL   (best)              
# 2 ZERO ALTERED NEGATIVE BINOMIAL                    
# 3 NEGATIVE BINOMIAL                             
# 4 ZERO ALTERED LOGARITHMIC 
# 5 GEOMETRIC
# 6 ZERO-INFLATED POISSON 
# 7 POISSON (worst)

##########################################################################################


#### Final marginal distributions for the copula model are the 
# (DOCTORCO) Zero altered logarithmic / ZALG and 
# (PRESCRIB) Zero inflated negative binomial / ZINBI
# distributions 

