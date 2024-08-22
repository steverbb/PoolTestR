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


# testing the basic things
library(devtools)
load_all(path = "/Users/gifted/Desktop/ANU/RA/brms")
load_all()
document()
dat <- mgcv::gamSim(1, n = 30, scale = 2)
# test it works for multiple gp terms
fit1 <- brm(y ~ gp(x0+x2, cov = "exp") + x1 + gp(x2) + x3 , data = dat, chains = 2)
summary(fit1)
fit2 <- PoolRegBayes(Result ~ Region + Year + gp(NumInPool) + gp(NumInPool, cov="exp"),
                                 data = SimpleExampleData,
                                 poolSize = NumInPool)
summary(fit2)
##################################################################
# simulation
rm(list=ls())
library(MASS)  # For multivariate normal distribution
library(dplyr)  # For data manipulation
library(devtools)
load_all(path = "/Users/gifted/Desktop/ANU/RA/brms")
load_all()
document()

# Step 1: Define Parameters
p <- 0.01 # prevalence
l <- 1.0  # Example length scale
sigma <- 1.75  # Example sigma

# Step 2: Generate Spatial Grid
set.seed(123)  # For reproducibility
N_points <- 100
lat <- runif(N_points, -2, 2)
lon <- runif(N_points, -2, 2)
coords <- cbind(lat, lon)

# Step 3: Simulate Gaussian Process with exponential kernel
dists <- as.matrix(dist(coords))
K <- sigma^2 * exp(-dists /  l )
det(K)
#View(K)
# K <- sigma^2 * exp(-dists^2 / (2 * l^2)) # exponential quadratic kernel
gp_values <- mvrnorm(1, mu = rep(0, N_points), Sigma = K)
mean(gp_values)
hist(gp_values)

gp_values=gp_values - 6 # adjust the mean
#install.packages("logitnorm")
library(logitnorm)
momentsLogitnorm(-6,sigma = 1.75)

hist(plogis(gp_values))
mean(plogis(gp_values))
# Step 4: Generate Prevalence Values
prevalence <- 1 / (1 + exp(-gp_values)) # inverse logit
prevalence <- prevalence * (p / mean(prevalence))  # Adjust to have mean ~1%
hist(prevalence)
mean(prevalence)

# Step 5: Create Pools
set.seed(2)
pool_data <- data.frame()
for (i in 1:N_points) {
  n_pools <- sample(2:3, 1)
  for (j in 1:n_pools) {
    pool_size <- sample(1:30, 1)
    #Assign Prevalence to Pools
    site_prevalence <- prevalence[i]
    positive = rbinom(1, 1, 1-(1-site_prevalence)^(pool_size) )
    pool_data <- rbind(pool_data, data.frame(Latitude = coords[i, 1], Longitude = coords[i, 2]
                                             , poolsize = pool_size, Positives = positive))
  }
}
head(pool_data)

# Prepare the data for PoolTestR
pool_data <- pool_data %>%
  mutate(id = 1:nrow(pool_data))

# Ensure no invalid data points
if (any(pool_data$Positives < 0 | pool_data$Positives > pool_data$poolsize)) {
  stop("Data contains invalid values for Positives or poolsize.")
}
# Check if any n values are negative or zero
if (any(pool_data$poolsize <= 0)) {
  stop("Pool sizes must be positive integers.")
}

# Step 6: Run PoolRegBayes
# Fit the model using PoolRegBayes
# Define initial values for each chain
init_values <- list(
  list(alpha = 0.01, sigma = 1.75, l = 1.0),
  list(alpha = 0.01, sigma = 1.75, l = 1.0),
  list(alpha = 0.01, sigma = 1.75, l = 1.0),
  list(alpha = 0.01, sigma = 1.75, l = 1.0)
)

fit <- PoolRegBayes(Positives ~ gp(Latitude,Longitude,cov="exp",scale = FALSE), data = pool_data, poolSize = poolsize
                    , control = list(adapt_delta = 0.9)
                    ) # errors
get_prior(fit)
fit$model
summary(fit)

fit <- PoolRegBayes(Positives ~ gp(Latitude,Longitude,cov="exp",scale = FALSE), data = pool_data, poolSize = poolsize
                    , control = list(adapt_delta = 0.9)
) # errors
get_prior(fit)
fit$model
summary(fit)


# next, we could increase the simulation sample size, 200, 500,
# run a few and record the time, time difference,
# see whether the approximation can work well
# measure for the estimation, for each dataset, how the estimation error, will it decrease as sample size increase?
# as k increase, the estimation error will decrease. find a reasonable k, that runs quickly, but we still have a nice estimation
#
# keep the number of divergent minimum, better to have no divergent transitions, adjust the adapt_delta
# record the number of divergent transitions


# the k in the gp terms, corresponding to which term in the paper, and how k changes will affect the approximation
# figure out, does the scale, isometric option: is the length scale the same for different distance, northsouth and eastwest,
#
# next next step, doing multivariate gaussian process,
fit <- PoolRegBayes(Positives ~ gp(Latitude,Longitude,cov="exp",scale = FALSE,k=10), data = pool_data, poolSize = poolsize
                    , control = list(adapt_delta = 0.9)
) # errors
get_prior(fit)
fit$model
summary(fit)
