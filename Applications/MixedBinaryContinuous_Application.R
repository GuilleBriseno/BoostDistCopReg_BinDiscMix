#load("../MixedBinCont_Application_undernutrition/data/india.Rdata")

india <- read.table("../MixedBinCont_Application_undernutrition/Files for application/data/india.raw", header = T)


ls()
head(india)
dim(india)

###
#library(expectreg)
#data(india)
#data(india.bnd)

#india.bnd
#unique(india$distH)

## load Copula functions
source("Copulas/Mixed_BinaryContinuous/Copula_Clayton270_NEW_BinCont.R")
source("Copulas/Mixed_BinaryContinuous/Copula_Gaussian_BinCont.R")

### Check the covariates used in the model by nadja:
library(gamboostLSS)
#library(R2BayesX)
#library(BayesX)

IndiaGRA <- BayesX::read.gra("../MixedBinCont_Application_undernutrition/Files for application/data/india_dist_sort.gra")




IndiaBND_FROMNK <- BayesX::read.bnd("../MixedBinCont_Application_undernutrition/Files for application/data/india_dist_sort.bnd")

# testfrombnd <- bnd2gra(india.bnd)
# 
# dim(IndiaGRA)
# dim(testfrombnd)
# 
# identical(sort(attributes(IndiaGRA)$dimnames[[1]]), sort(attributes(testfrombnd)$dimnames[[1]]))
# identical(sort(attributes(IndiaGRA)$dimnames[[2]]), sort(attributes(testfrombnd)$dimnames[[2]]))
# 
# colSums(IndiaGRA)
# testfrombnd
# 
# colnames(IndiaGRA)


### Load bnd file:
#IndiaBnd <- read.bnd(file = "../MixedBinCont_Application_undernutrition/data/india_dist_sort.bnd")
#data("india.bnd", package = "expectreg")

#IndiaBnd <- india.bnd

# #### I'm not super sure if the districts here match the ones in the map: 
# # 440 unique values
# uniqueDistricts_expectregObject <- attr(IndiaBnd, "regions")
# sort(as.numeric(uniqueDistricts_expectregObject))
# length(uniqueDistricts_expectregObject)
# 
# 
# uniqueDistricts_fromData <- sort(unique(india$distH))
# length(uniqueDistricts_fromData)




# ### get the levels of those with less than 40 observations. These will always be in the sample. 
districtsTable <- as.data.frame(table(india$distH)) # get the table


lessThan40 <- districtsTable$Var1[which(districtsTable$Freq < 40)] # districts with less than 40 observations
moreThan40 <- districtsTable$Var1[which(districtsTable$Freq >= 40)] # districts with at least 40 observations

### initialise all weights at 1 (all observations would be in the sample )
boost.weights <- rep(1, nrow(india))

## the idea is that we check whether the district of the current observation is one of those that has less than 40 observations.
# if it is, then the weight of the observation will be left at 1, otherwise it will be sampled from (0,1) with the probabilities below
set.seed(1)
for(i in 1:nrow(india)){
  
  if(india$distH[i] %in% moreThan40){
    
    boost.weights[i] <- sample(c(0,1), size = 1, prob = c(0.305, 0.695)) #prob = c(0.41, 0.59)) #
    
  }else{
    
    boost.weights[i] <- 1
  } 
}

mean(boost.weights) # around 2/3 rds will be used for fitting(?)        around 75 % of the observations will be used for fitting


### some formatting for the data:
india$csex <- factor(india$csex)
india$distH <- factor(india$distH)

# standardise wasting!
india$wastingC <- ( india$wasting - mean(india$wasting) )/ sd(india$wasting)
#india$wastingMedian <- ( india$wasting - median(india$wasting) )/ sd(india$wasting)


round(mean(india$fever), digits = 3)
round(sd(india$fever), digits = 3)

round(mean(india$wasting), digits = 3)
round(sd(india$wasting), digits = 3)


round(mean((as.numeric(india$csex) - 1)), digits = 3)

round(mean(india$mbmi), digits = 3)
round(sd(india$mbmi), digits = 3)


round(mean(india$breastfeeding), digits = 3)
round(sd(india$breastfeeding), digits = 3)

round(mean(india$cage), digits = 3)
round(sd(india$cage), digits = 3)




## model formulas: feature all covariates used in the study
modelequation <- cbind(fever, wastingC) ~ bols(csex) + bbs(cage) + bbs(breastfeeding) + bbs(mbmi) + bmrf(distH, bnd = IndiaGRA)

biv_equations <- list(mu1 = modelequation, 
                      mu2 = modelequation, 
                      sigma2 = modelequation, 
                      rho = modelequation)

boost.nu.short <- 0.1

boost.stabilisation <- "L2"


### Model: 
Copula_families_Probit_Gaussian <- Clayton270_Cop_BinCont(marg1 = "PROBIT", marg2 = "NORM", stabilization = boost.stabilisation)


Clayton270_Model <- gamboostLSS(biv_equations, 
                                data = india,
                                families = Copula_families_Probit_Gaussian, 
                                weights = boost.weights,
                                control = boost_control(mstop = 35000,
                                                        risk = "oobag",
                                                        nu  = boost.nu.short, 
                                                        trace = TRUE), 
                                method = 'noncyclic')

plot(risk(Clayton270_Model, merge = TRUE))

OPTIMAL_MSTOP <- which.min(risk(Clayton270_Model, merge = TRUE))
### 28006 

#### 27994  using the different bnd object. <--- this was fucked up anyway

# 27846 

Clayton270_Model <- Clayton270_Model[OPTIMAL_MSTOP]



Clayton270_Model_ALLDATA <- gamboostLSS(biv_equations, 
                                        data = india,
                                        families = Copula_families_Probit_Gaussian, 
                                        control = boost_control(mstop = OPTIMAL_MSTOP,
                                                                nu  = boost.nu.short, 
                                                                trace = TRUE), 
                                        method = 'noncyclic')





summary(Clayton270_Model)

summary(Clayton270_Model_ALLDATA)

### out-of-bag-risk:
tail(risk(Clayton270_Model, merge = T), 1)
# 11504.52 

## 11504.18 using the different bnd object


######## 11503.79 final risk with the actual GRA file: 


VineCopula::BiCopPar2Tau(family = 33, range(-predict(Clayton270_Model_ALLDATA$rho, type = "response")))

round(VineCopula::BiCopPar2Tau(family = 33, range(-predict(Clayton270_Model_ALLDATA$rho, type = "response"))), digits = 3)



################## Now lets try to plot these things: 
summary(Clayton270_Model_ALLDATA)

par(mfrow = c(1, 5))
plot(Clayton270_Model_ALLDATA$mu1)

par(mfrow = c(1, 4))
plot(Clayton270_Model_ALLDATA$mu2)


par(mfrow = c(1, 5))
plot(Clayton270_Model_ALLDATA$sigma2)


par(mfrow = c(1, 3))
plot(Clayton270_Model_ALLDATA$rho)




####################### Prepare some plots:  bols(csex) + bbs(cage) + bbs(breastfeeding) + bbs(mbmi) + bmrf(distH, bnd = IndiaBnd)
range(india$cage)
range(india$breastfeeding)
range(india$mbmi)
levels(india$distH)


IndiaPredictions <- data.frame(cage = seq(0, 35, length.out = 250), 
                               breastfeeding = seq(0, 35, length.out = 250), 
                               mbmi = seq(12.45, 39.85, length.out = 250), 
                               distH = factor(rep("5", 250)))


##### mu 1
par(mfrow = c(1, 3))
plot(predict(Clayton270_Model_ALLDATA$mu1, which = "cage", newdata = IndiaPredictions), type = "l", ylim = c(-1, 0.5))
plot(predict(Clayton270_Model_ALLDATA$mu1, which = "breastfeeding", newdata = IndiaPredictions), type = "l", ylim = c(-1, 0.5))
plot(predict(Clayton270_Model_ALLDATA$mu1, which = "mbmi", newdata = IndiaPredictions), type = "l", ylim = c(-1, 0.5))


##### mu 2
par(mfrow = c(1, 3))
plot(predict(Clayton270_Model_ALLDATA$mu2, which = "cage", newdata = IndiaPredictions), type = "l", ylim = c(-0.5, 0.5))
plot(predict(Clayton270_Model_ALLDATA$mu2, which = "breastfeeding", newdata = IndiaPredictions), type = "l", ylim = c(-0.5, 0.5))
plot(predict(Clayton270_Model_ALLDATA$mu2, which = "mbmi", newdata = IndiaPredictions), type = "l", ylim = c(-0.5, 0.5))

#### sigma 2
plot(predict(Clayton270_Model_ALLDATA$sigma2, which = "cage", newdata = IndiaPredictions), type = "l")
plot(predict(Clayton270_Model_ALLDATA$sigma2, which = "breastfeeding", newdata = IndiaPredictions), type = "l")
plot(predict(Clayton270_Model_ALLDATA$sigma2, which = "mbmi", newdata = IndiaPredictions), type = "l")

#### rho 
plot(predict(Clayton270_Model_ALLDATA$rho, which = "cage", newdata = IndiaPredictions), type = "l", ylim = c(-1, 1))
abline(a=0, b = 0, col = "red")
plot(predict(Clayton270_Model_ALLDATA$rho, which = "breastfeeding", newdata = IndiaPredictions), type = "l")
plot(predict(Clayton270_Model_ALLDATA$rho, which = "mbmi", newdata = IndiaPredictions), type = "l")



######### cage  covariate
par(mfrow = c(2,2))
plot(predict(Clayton270_Model_ALLDATA$mu1, which = "cage", newdata = IndiaPredictions), type = "l", ylim = c(-1, 0.5))
plot(predict(Clayton270_Model_ALLDATA$mu2, which = "cage", newdata = IndiaPredictions), type = "l", ylim = c(-1, 0.5))
plot(predict(Clayton270_Model_ALLDATA$sigma2, which = "cage", newdata = IndiaPredictions), type = "l", ylim = c(-1, 0.5))
plot(predict(Clayton270_Model_ALLDATA$rho, which = "cage", newdata = IndiaPredictions), type = "l", ylim = c(-1, 0.5))


######### breastfeeding  covariate
par(mfrow = c(2,2))
plot(predict(Clayton270_Model_ALLDATA$mu1, which = "breastfeeding", newdata = IndiaPredictions), type = "l", ylim = c(-0.25, 0.25))
plot(predict(Clayton270_Model_ALLDATA$mu2, which = "breastfeeding", newdata = IndiaPredictions), type = "l", ylim = c(-0.25, 0.25))
plot(predict(Clayton270_Model_ALLDATA$sigma2, which = "breastfeeding", newdata = IndiaPredictions), type = "l", ylim = c(-0.25, 0.25))
plot(predict(Clayton270_Model_ALLDATA$rho, which = "breastfeeding", newdata = IndiaPredictions), type = "l", ylim = c(-0.25, 0.25))

######### mbmi  covariate
par(mfrow = c(2,2))
plot(predict(Clayton270_Model_ALLDATA$mu1, which = "mbmi", newdata = IndiaPredictions), type = "l", ylim = c(-1, 0.5))
plot(predict(Clayton270_Model_ALLDATA$mu2, which = "mbmi", newdata = IndiaPredictions), type = "l", ylim = c(-1, 0.5))
plot(predict(Clayton270_Model_ALLDATA$sigma2, which = "mbmi", newdata = IndiaPredictions), type = "l", ylim = c(-1, 0.5))
plot(predict(Clayton270_Model_ALLDATA$rho, which = "mbmi", newdata = IndiaPredictions), type = "l", ylim = c(-1, 0.5))

##################
Plot_CAGE <- data.frame(mu1 = predict(Clayton270_Model_ALLDATA$mu1, which = "cage", newdata = IndiaPredictions), 
                        mu2 = predict(Clayton270_Model_ALLDATA$mu2, which = "cage", newdata = IndiaPredictions),
                        sigma2 = predict(Clayton270_Model_ALLDATA$sigma2, which = "cage", newdata = IndiaPredictions), 
                        rho = predict(Clayton270_Model_ALLDATA$rho, which = "cage", newdata = IndiaPredictions))

colnames(Plot_CAGE) <- c("mu1", "mu2", "sigma2", "rho")

Plot_CAGE <- stack(Plot_CAGE)
Plot_CAGE$cage <- rep(IndiaPredictions$cage, 4)
Plot_CAGE$ind <- factor(Plot_CAGE$ind,
                        levels = c("mu1", "mu2", "sigma2", "rho"),
                        labels = c("vartheta[1]^(1)", 
                                   "vartheta[1]^(2)",
                                   "vartheta[2]^(2)",
                                   "vartheta^(c)"))

CAGE_Plot <- ggplot(Plot_CAGE, aes(x = cage, y = values)) +
  geom_line() +
  labs(x = "Child's age (months)", y = expression(hat(f)( paste("Child's age") )), title = "(a)" ) + 
  #lims(y = c(-1, 0.5)) +
  facet_grid(~ind, labeller=labeller(ind = label_parsed )) +
  scale_x_continuous(breaks = c(0,5,10, 15, 20,25,30,35)) + 
  scale_y_continuous(breaks =  c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5), labels = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5)) +
  geom_hline(col = "grey", yintercept = 0) + 
  theme_light() +
  theme(strip.text = element_text(size = 15), 
        axis.text.y = element_text(size = 8), 
        axis.text.x = element_text(size = 8), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = -5, b = 0, l = 0)),
        axis.title.x = element_text(size = 17, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size = 23))

CAGE_Plot

Plot_breastfeed <- data.frame(mu1 = predict(Clayton270_Model_ALLDATA$mu1, which = "breastfeeding", newdata = IndiaPredictions), 
                              mu2 = predict(Clayton270_Model_ALLDATA$mu2, which = "breastfeeding", newdata = IndiaPredictions),
                              sigma2 = predict(Clayton270_Model_ALLDATA$sigma2, which = "breastfeeding", newdata = IndiaPredictions), 
                              rho = predict(Clayton270_Model_ALLDATA$rho, which = "breastfeeding", newdata = IndiaPredictions))


colnames(Plot_breastfeed) <- c("mu1", "mu2", "sigma2", "rho")

Plot_breastfeed <- stack(Plot_breastfeed)
Plot_breastfeed$breastfeeding <- rep(IndiaPredictions$breastfeeding, 4)
Plot_breastfeed$ind <- factor(Plot_breastfeed$ind,
                              levels = c("mu1", "mu2", "sigma2", "rho"),
                              labels = c("vartheta[1]^(1)", 
                                         "vartheta[1]^(2)",
                                         "vartheta[2]^(2)",
                                         "vartheta^(c)"))

Breastfeed_PLOT <- ggplot(Plot_breastfeed, aes(x = breastfeeding, y = values)) +
  geom_line() +
  labs(x = "Months of breastfeeding", y = expression(hat(f)( paste("Breastfeeding") )), title = "(b)" ) + 
  #lims(y = c(-0.15, 0.15)) +
  facet_grid(~ind, labeller=labeller(ind = label_parsed )) +
  scale_x_continuous(breaks = c(0,5,10, 15, 20,25,30,35)) + 
  scale_y_continuous(breaks = c(-0.15, -0.05, 0, 0.05, 0.15), labels = c(-0.15, -0.05, 0, 0.05, 0.15)) +
  geom_hline(col = "grey", yintercept = 0) + 
  theme_light() +
  theme(strip.text = element_text(size = 15), 
        axis.text.y = element_text(size = 8), 
        axis.text.x = element_text(size = 8), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = -5, b = 0, l = 0)),
        axis.title.x = element_text(size = 17, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size = 23))


Breastfeed_PLOT

Plot_mbmi <- data.frame(mu1 = predict(Clayton270_Model_ALLDATA$mu1, which = "mbmi", newdata = IndiaPredictions), 
                        mu2 = predict(Clayton270_Model_ALLDATA$mu2, which = "mbmi", newdata = IndiaPredictions),
                        sigma2 = predict(Clayton270_Model_ALLDATA$sigma2, which = "mbmi", newdata = IndiaPredictions), 
                        rho = predict(Clayton270_Model_ALLDATA$rho, which = "mbmi", newdata = IndiaPredictions))


colnames(Plot_mbmi) <- c("mu1", "mu2", "sigma2", "rho")

Plot_mbmi <- stack(Plot_mbmi)
Plot_mbmi$mbmi <- rep(IndiaPredictions$mbmi, 4)
Plot_mbmi$ind <- factor(Plot_mbmi$ind,
                        levels = c("mu1", "mu2", "sigma2", "rho"),
                        labels = c("vartheta[1]^(1)", 
                                   "vartheta[1]^(2)",
                                   "vartheta[2]^(2)",
                                   "vartheta^(c)"))


MBMI_PLOT <- ggplot(Plot_mbmi, aes(x = mbmi, y = values)) +
  geom_line() +
  labs(x = "Mother's Body-Mass-Index", y = expression(hat(f)( paste("Mother's BMI") )), title = "(c)" ) + 
  #lims(y = c(-1, 2)) +
  facet_grid(~ind, labeller=labeller(ind = label_parsed )) +
  scale_x_continuous(breaks = c(12,15,20,25,30,35, 40)) + 
  scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5,  2), labels = c(-1, -0.5, 0, 0.5, 1, 1.5,  2)) +
  geom_hline(col = "grey", yintercept = 0) + 
  theme_light() +
  theme(strip.text = element_text(size = 15), 
        axis.text.y = element_text(size = 8), 
        axis.text.x = element_text(size = 8), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = -5, b = 0, l = 0)),
        axis.title.x = element_text(size = 17, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size = 23))

MBMI_PLOT


library(patchwork)

CAGE_Plot + Breastfeed_PLOT + MBMI_PLOT + plot_layout(nrow = 3, ncol = 1)


#### NOW THE SPATIAL EFFECT:











################################
R2BayesX::plotmap(map = india.bnd, 
                  cbind(india$distH, predict(Clayton270_Model_ALLDATA$mu1, which = "distH", type = "link")), range = c(-0.6, +0.45))
R2BayesX::plotmap(map = IndiaBND_FROMNK,
                  cbind(india$distH, predict(Clayton270_Model_ALLDATA$mu1, which = "distH", type = "link")), range = c(-0.6, +0.45))


R2BayesX::plotmap(map = india.bnd, 
                  cbind(india$distH, predict(Clayton270_Model_ALLDATA$mu2, which = "distH", type = "link")), range = c(-0.5, 0.7))
R2BayesX::plotmap(map = IndiaBND_FROMNK, 
                  cbind(india$distH, predict(Clayton270_Model_ALLDATA$mu2, which = "distH", type = "link")), range = c(-0.5, 0.7))


R2BayesX::plotmap(map = india.bnd, 
                  cbind(india$distH, predict(Clayton270_Model_ALLDATA$sigma2, which = "distH", type = "link")), range = c(-0.5, 0.65))
R2BayesX::plotmap(map = IndiaBND_FROMNK,
                  cbind(india$distH, predict(Clayton270_Model_ALLDATA$sigma2, which = "distH", type = "link")), range = c(-0.5, 0.65))


R2BayesX::plotmap(map = india.bnd, 
                  cbind(india$distH, predict(Clayton270_Model_ALLDATA$rho, which = "distH", type = "link")), range = c(-0.97, -0.8))
R2BayesX::plotmap(map = IndiaBND_FROMNK, 
                  cbind(india$distH, predict(Clayton270_Model_ALLDATA$rho, which = "distH", type = "link")), range = c(-0.97, -0.8))

length(unique(predict(Clayton270_Model_ALLDATA$rho, which = "distH", type = "link")))


range(predict(Clayton270_Model_ALLDATA$mu1, which = "distH", type = "link"))
range(predict(Clayton270_Model_ALLDATA$mu2, which = "distH", type = "link"))
range(predict(Clayton270_Model_ALLDATA$sigma2, which = "distH", type = "link"))
range(predict(Clayton270_Model_ALLDATA$rho, which = "distH", type = "link"))


length(levels(india$distH))
sort(unique(table((india$distH))))

length(unique(attributes(IndiaBND_FROMNK)$regions))

sort(as.numeric(unique(attributes(IndiaBND_FROMNK)$regions)))



length(unique(attributes(india.bnd)$regions))

sort(as.numeric(unique(attributes(india.bnd)$regions)))

sort(as.numeric(colnames(IndiaGRA)))

sum(colnames(IndiaGRA) == rownames(IndiaGRA))




sum(sort(as.numeric(unique(attributes(india.bnd)$regions))) %in% sort(as.numeric(colnames(IndiaGRA))))


sum(sort(as.numeric(unique(attributes(IndiaBND_FROMNK)$regions))) %in% sort(unique((india$distH))))


############################# Predict distributional quantities: 
### probability of fever: 
Pmargin1 <- predict(Clayton270_Model_ALLDATA$mu1, which = "distH", type = "response")

DistrQuants <- india %>% select(distH, wasting2, fever) %>% mutate(EMargin1 = predict(Clayton270_Model_ALLDATA$mu1, type = "response"), 
                                                                   SDMargin1 = sqrt(predict(Clayton270_Model_ALLDATA$mu1, type = "response")*(1 - predict(Clayton270_Model_ALLDATA$mu1, type = "response"))),
                                                                   EMargin2 = predict(Clayton270_Model_ALLDATA$mu2, type = "response"), 
                                                                   SDMargin2 = predict(Clayton270_Model_ALLDATA$sigma2, type = "response"), 
                                                                   CopulaParameter = -predict(Clayton270_Model_ALLDATA$rho, type = "response"), 
                                                                   KendallTau = VineCopula::BiCopPar2Tau(33, -predict(Clayton270_Model_ALLDATA$rho, type = "response")), 
                                                                   JointProb = VineCopula::BiCopCDF(family = 33, par = -predict(Clayton270_Model_ALLDATA$rho, type = "response"),
                                                                                                    u2 = (1-predict(Clayton270_Model_ALLDATA$mu1, type = "response")), 
                                                                                                    u1 = pnorm(wasting2, 
                                                                                                               mean = predict(Clayton270_Model_ALLDATA$mu2, type = "response"), 
                                                                                                               sd = predict(Clayton270_Model_ALLDATA$sigma2, type = "response")))
                                                                   )


DistH_Summaries <- DistrQuants %>% group_by(distH) %>% summarise(EXM1 = mean(EMargin1), 
                                                                 SDM1 = mean(SDMargin1), 
                                                                 EXM2 = mean(EMargin2), 
                                                                 SDM2 = mean(SDMargin2), 
                                                                 CopParam = mean(CopulaParameter), 
                                                                 Kendall = mean(KendallTau),
                                                                 JointProb = mean(JointProb))


DistH_Summaries_MED <- DistrQuants %>% group_by(distH) %>% summarise(EXM1 = median(EMargin1), 
                                                                 SDM1 = median(SDMargin1), 
                                                                 EXM2 = median(EMargin2), 
                                                                 SDM2 = median(SDMargin2), 
                                                                 CopParam = median(CopulaParameter), 
                                                                 Kendall = median(KendallTau),
                                                                 JointProb = mean(JointProb))





colnames(DistH_Summaries)

apply(DistH_Summaries[,-1], 2, range)
apply(DistH_Summaries_MED[,-1], 2, range)

###### MARGIN 1
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,2)], range = c(0.1, 0.6), col = gray.colors(nrow(DistH_Summaries)))
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,3)], range = c(0.3, 0.5), col = gray.colors(nrow(DistH_Summaries)))

###### MARGIN 2
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,4)], range = c(-0.6, 1.1), col = gray.colors(nrow(DistH_Summaries)))
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,5)], range = c(0.6, 1.85), col = gray.colors(nrow(DistH_Summaries)))

#### Copula parameter: 
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,6)], range = c(-0.23, -0.1), col = gray.colors(nrow(DistH_Summaries)))
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,7)], range = c(-0.1, -0.05), col = gray.colors(nrow(DistH_Summaries)))


### Joint probabilities: 
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,8)], range = c(0.1, 0.65), col = gray.colors(nrow(DistH_Summaries)))



################################ Different palette
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,2)], range = c(0.1, 0.6), col = hcl.colors(nrow(DistH_Summaries), "viridis"))
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,3)], range = c(0.3, 0.5), col = hcl.colors(nrow(DistH_Summaries), "viridis"))

###### MARGIN 2
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,4)], range = c(-0.6, 1.1), col = hcl.colors(nrow(DistH_Summaries), "viridis"))
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,5)], range = c(0.6, 1.85), col = hcl.colors(nrow(DistH_Summaries), "viridis"))

#### Copula parameter: 
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,6)], range = c(-0.23, -0.1), col = hcl.colors(nrow(DistH_Summaries), "viridis"))
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,7)], range = c(-0.1, -0.05), col = hcl.colors(nrow(DistH_Summaries), "viridis"))

### Joint probabilities: 
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,8)], range = c(0.1, 0.65), col = hcl.colors(nrow(DistH_Summaries), "viridis"))




R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,2)], range = c(0.1, 0.6), col = hcl.colors(nrow(DistH_Summaries), "blues"))
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,3)], range = c(0.3, 0.5), col = hcl.colors(nrow(DistH_Summaries), "blues"))
###### MARGIN 2
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,4)], range = c(-0.6, 1.1), col = hcl.colors(nrow(DistH_Summaries), "blues"))
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,5)], range = c(0.6, 1.85), col = hcl.colors(nrow(DistH_Summaries), "blues"))
#### Copula parameter: 
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,6)], range = c(-0.23, -0.1), col = hcl.colors(nrow(DistH_Summaries), "blues"))
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,7)], range = c(-0.1, -0.05), col = hcl.colors(nrow(DistH_Summaries), "blues"))
### Joint probabilities: 
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,8)], range = c(0.1, 0.65), col = hcl.colors(nrow(DistH_Summaries), "blues"))



### Terrain colors
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,2)], range = c(0.1, 0.6), col = terrain.colors(nrow(DistH_Summaries), alpha = 1),
                  main = "(a)")
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,3)], range = c(0.3, 0.5), col = terrain.colors(nrow(DistH_Summaries), alpha = 1),
                  main = "(b)")
###### MARGIN 2
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,4)], range = c(-0.6, 1.1), col = terrain.colors(nrow(DistH_Summaries), alpha = 1),
                  main = "(c)")
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,5)], range = c(0.6, 1.85), col = terrain.colors(nrow(DistH_Summaries), alpha = 1),
                  main = "(d)")
#### Copula parameter: 
# R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,6)], range = c(-0.23, -0.1), col = terrain.colors(nrow(DistH_Summaries), alpha = 1),
#                   main = "(e)")
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,7)], range = c(-0.1, -0.05), col = terrain.colors(nrow(DistH_Summaries), alpha = 1),
                  main = "(e)")
### Joint probabilities: 
R2BayesX::plotmap(map = india.bnd, DistH_Summaries[,c(1,8)], range = c(0.1, 0.65), col = terrain.colors(nrow(DistH_Summaries), alpha = 1),
                  main = "(f)")

######################### 











######################
par(mfrow = c(1, 5))
plot(Gauss_Model$mu1)

par(mfrow = c(1, 4))
plot(Gauss_Model$mu2)


par(mfrow = c(1, 5))
plot(Gauss_Model[OPTIMAL_MSTOP_GAUSS]$sigma2)


par(mfrow = c(1, 5))
plot(Gauss_Model[OPTIMAL_MSTOP_GAUSS]$rho)
########################


par(mfrow = c(2, 5))
plot(Clayton270_Model[OPTIMAL_MSTOP]$mu1)
plot(Gauss_Model[OPTIMAL_MSTOP_GAUSS]$mu1)

par(mfrow = c(2, 5))
plot(Clayton270_Model[OPTIMAL_MSTOP]$mu2)
plot(Gauss_Model[OPTIMAL_MSTOP_GAUSS]$mu2)


par(mfrow = c(2, 5))
plot(Clayton270_Model[OPTIMAL_MSTOP]$sigma2)
plot(Gauss_Model[OPTIMAL_MSTOP_GAUSS]$sigma2)


par(mfrow = c(2, 5))
plot(Clayton270_Model[OPTIMAL_MSTOP]$rho)
plot(Gauss_Model[OPTIMAL_MSTOP_GAUSS]$rho)
















# 
# 
# ## save the model:
# save(OPTIMAL_MSTOP, boost.weights, 
#      Clayton270_Model, 
#      file = "Applications/MixedBinCont_Application_Clayton270Model.RData")
# 
# 
# par(mfrow = c(1, 5))
# plot(Clayton270_Model$mu1)
# 
# par(mfrow = c(1, 5))
# plot(Clayton270_Model$mu2)
# 
# 
# par(mfrow = c(1, 5))
# plot(Clayton270_Model$sigma2)
# 
# 
# par(mfrow = c(1, 4))
# plot(Clayton270_Model$rho)
# 





###############
### Model: 
Copula_families_Probit_Gaussian <- Gauss_Cop_BinCont(marg1 = "PROBIT", marg2 = "NORM", stabilization = boost.stabilisation)


Gauss_Model <- gamboostLSS(biv_equations, 
                                data = india,
                                families = Copula_families_Probit_Gaussian, 
                                weights = boost.weights,
                                control = boost_control(mstop = 10000,
                                                        risk = "oobag",
                                                        nu  = boost.nu.short, 
                                                        trace = TRUE), 
                                method = 'noncyclic')

Gauss_Model <- Gauss_Model[20000]
Gauss_Model <- Gauss_Model[30000]


plot(risk(Gauss_Model, merge = TRUE))

OPTIMAL_MSTOP_GAUSS <- which.min(risk(Gauss_Model, merge = TRUE))
### 19652 


tail(risk(Clayton270_Model[OPTIMAL_MSTOP], merge = T), 1) # 11504.52 
tail(risk(Gauss_Model[OPTIMAL_MSTOP_GAUSS], merge = T), 1) # 11503.77 

