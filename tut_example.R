library(rstan)
schools_dat <- list(J = 8,
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
fit <- stan(file = 'tut_example.stan', data = schools_dat)
print(fit)
plot(fit)
pairs(fit, pars = c("mu", "tau", "lp__"))

la <- extract(fit, permuted = TRUE) # return a list of arrays
mu <- la$mu

### return an array of three dimensions: iterations, chains, parameters
a <- extract(fit, permuted = FALSE)
a
### use S3 functions on stanfit objects
a2 <- as.array(fit)
m <- as.matrix(fit)
d <- as.data.frame(fit)
a2
m
d

rm(list=ls())
library(PoolTestR)
head(SimpleExampleData)
PrevWholeDataset <- PoolPrev(SimpleExampleData,Result,NumInPool)
PrevWholeDataset
PrevByRegion <- PoolPrev(SimpleExampleData, Result, NumInPool, Region)
PrevByRegion
PrevByVillage <- PoolPrev(SimpleExampleData, Result, NumInPool, Village)
PrevByVillage
PrevByYear <- PoolPrev(SimpleExampleData, Result, NumInPool, Year)
PrevByYear
PrevByRegionYear <- PoolPrev(SimpleExampleData, Result, NumInPool, Region, Year)
PrevByRegionYear
Model <- PoolReg(Result ~ Year + Region, SimpleExampleData, NumInPool)
summary(Model)
getPrevalence(Model)
PrevByHier <- HierPoolPrev(SimpleExampleData, Result, NumInPool,
                           c("Village","Site"))
PrevByHier
PrevByYearHier <- HierPoolPrev(SimpleExampleData, Result, NumInPool,
                               c("Village","Site"), Year)
PrevByYearHier
BayesMod <- PoolRegBayes(Result ~ Region + Year,
                         data = SimpleExampleData,
                         poolSize = NumInPool)
summary(BayesMod)
getPrevalence(BayesMod) #Bayesian model
