# -------------------------------------------------------
#
#                    Load libraries
#
# -------------------------------------------------------

# Install packages (if needed)
# install.packages("tidyverse")
# install.packages("jagsUI")
# install.packages("coda")
# install.packages("mcmcplots")
# install.packages("loo")

# Load library
library(tidyverse)
library(jagsUI)
library(coda)
library(mcmcplots)
library(loo)

# Set seed, scientific notation options, and working directory
set.seed(123)
options(scipen = 9999)
setwd(".")

# Setting up cores
Ncores <- parallel::detectCores()
print(Ncores) # Number of available cores
workers <- Ncores * 0.5 # For low background use 80%, for medium use 50% of Ncores

# -------------------------------------------------------
#
#             Variable and Object Definitions
#
# -------------------------------------------------------

# beta.0 = abundance intercept 
# beta.1 = abundance trend estimate
# alpha.0 = prob (on logit scale) of detecting at least one vocalization at a site
#           that is not occupied.
# alpha.1 = additional prob (on logit scale) of detecting at least one vocalization at 
#           a site that is not occupied. 
# omega = mean # of false positive acoustic detections
# p = detection probability of an individual in point count data
# tau.day = precision for random day effect on true vocalization detection rate. 
# a.phi = overdispersion parameter for zero-truncated negative binomial. 
# gamma.1 = random day effect on true vocalization detection rate
# n.days = number of recording days.
# N = latent abundance process
# tp = true positive rate
# p.a = prob of detecting at least one vocalization in an acoustic recording
# v = acoustic vocalization data from clustering algorithm
# y = binary summary of v, takes value 0 if v = 0 and value 1 if v > 0. 
# c = point count data
# R = number of total sites
# J = number of repeat visits for acoustic data at each site
# J.A = max number of repeat visits at each acoustic data site. 
# n.count = number of repeat visits for count data
# R.val = number of sites where validation of acoustic data occurred
# days = variable used to index different recording days. 
# A.times = indexing variable used to determine specific indexes with v > 0.
# K = true number of acosutic vocalizations
# k = true number of manually validated acoustic vocalizations
# n = total number of manually validated acoustic vocalizations 
# sites.a = specific indices of sites where acoustic data were obtained
# R.val = number of validated sites for acoustic data
# Other variables not defined are for computation of Bayesian p-values. 

# -------------------------------------------------------
#
#                    Load Data
#
# -------------------------------------------------------

# BirdNet detections
bnet_dat_all <- read.csv("./Data/Acoustic_Data/NOBO_BirdNETall_2024.csv")

# ARU weather covariates
weather_dat <- read.csv("./Data/Acoustic_Data/ARU_weathercovs.csv")

# Site covariates
site_covs <- read.csv("./Data/Acoustic_Data/ARU_siteCovs.csv")


# -------------------------------------------------------
#
#    Data Wrangling
#
# -------------------------------------------------------


# Subset to 14 days starting at May 26 and ending on July 17. Dates are every 4 days.
date_order <- c("2024-05-26", "2024-05-30", "2024-06-03", "2024-06-07", # in ascending order
                "2024-06-11", "2024-06-15", "2024-06-19", "2024-06-23",
                "2024-06-27", "2024-07-01", "2024-07-05", "2024-07-09",
                "2024-07-13", "2024-07-17")

# Dates and their corresponding occasion numbers
bnet_dat <- bnet_dat_all %>% filter(Date %in% date_order)
head(bnet_dat)

# Adding occasion column
bnet_dat <- bnet_dat %>%
  mutate(Occasion = match(Date, date_order))

# Adding a row count
bnet_dat$Count <- 1


# Initialize the matrix
v <- matrix(0, nrow = 27, ncol = 14)        

# Extract count data
for (i in 1:nrow(bnet_dat)) {
  
  # Extracting plot ID
  site <- bnet_dat$Site_Number[i]
  
  # Extracting Occasion
  occasion <- bnet_dat$Occasion[i]
  
  # Fill in the matrix with the number of individuals
  v[site, occasion] =  v[site, occasion] + bnet_dat$Count[i]
  
} # end loop 

# take a look
print(v)
                             
# ARU at site 24 stopped recording on 6/20. Making Columns 8:14 NA
v[24,8:14] <- NA
print(v)  

# Renaming columns to date Month_day
formatted_dates <- format(as.Date(date_order), "%b_%d")
colnames(v) <- formatted_dates
print(v) 

# Adding rownames
rownames(v) <- as.numeric(1:27)

# Get the site numbers with at least one call
v_mat <- as.matrix(v)
sites.a <- which(rowSums(!is.na(v_mat) & v >= 1) > 0)
print(sites.a)



# Creating a binary detection matrix, values >= 1 become 1, and 0 stays 0
v_df <- as.data.frame(v)
y <- v_df %>% mutate(across(where(is.numeric), ~ if_else(. >= 1, 1, 0)))
y <- as.matrix(y)
rownames(y) <- NULL  
colnames(y) <- NULL 
print(y)

# Total number of sites
R <- as.integer(27) 

# Number of repeat visits for each site
J <- rep(14, R)  



# ---------------------------------------
# Manually validated  
# ---------------------------------------


# Validated Calls
# Do not include sites with no calls 
# only include occasions where at least 1 call was validated for a site
n <- read.csv("./Data/Acoustic_Data/Bnet14day_n.csv", row.names = 1)

# True Calls
# Calls found to be true, same dimension as n
k <- read.csv("./Data/Acoustic_Data/Bnet14day_k.csv", row.names = 1)

# Survey days calls were validated, same dimension as n
val.times <- read.csv("./Data/Acoustic_Data/Bnet14day_val.times.csv", row.names = 1)

# Total number of sites with manually validated data
R.val <- nrow(n)

# How many surveys were validate
J.val <- rep(14, R.val)  
 

# Check
dim(n) # dimensions
dim(k)
dim(val.times)
print(n) # format
print(k)
print(val.times)


# ---------------------------------------
# Covariates 
# ---------------------------------------

# For day random effect on cue production rate
days <- matrix(NA, nrow = 27, ncol = 14)  
for (col in 1:14) {
  days[, col] <- col
}
print(days)




# J.r contains the number of surveys at each acoustic site that contains at least 1 detected vocalization. 
J.r <- apply(v_mat, 1, function(a) {sum(a != 0)})
J.r <- ifelse(is.na(J.r), 0, J.r)
J.r <- as.numeric(J.r)
print(J.r)

# A.times is a R x J matrix that contains the indices of the
# specific surveys at each acoustic site that have at least one detected
# vocalization. This determines the specific components of v[i, j] that
# are used in the zero-truncated Poisson.
A.times <- matrix(NA, R, max(J))
tmp <- apply(v, 1, function(a) {which(a != 0)})
for (i in 1:R) {
  if (length(tmp[[i]]) > 0) {
    A.times[i, 1:length(tmp[[i]])] <- tmp[[i]]
  }
}
print(A.times)



# For Bayesian P-value
R.A <- sum(J.r > 0)
sites.a.v <- which(J.r > 0)
J.A <- max(J)


# Bundle data 
Bnet14.data <- list(R = R, 
                    J = J, 
                    v = v, 
                    y = y, 
                    k = k, 
                    n = n, 
                    sites.a = sites.a, 
                    R.val = R.val, 
                    J.val = J.val, 
                    J.r = J.r, 
                    A.times = A.times, 
                    val.times = val.times, 
                    R.A = R.A, 
                    J.A = J.A, 
                    sites.a.v = sites.a.v, 
                    days = days,
                    n.days = max(J))

# ,
#                     herbPdens = herbPdens,
#                     woodyParea = woodyParea

# Check structure
str(Bnet14.data)


# -------------------------------------------------------
# Model Specifications
# -------------------------------------------------------

# Parameters monitored
params <- c('alpha.0', 
              'alpha.1', 
              'beta.0',
              'beta.1',
              'beta.2',
              'tau', 
              'tau.day', 
              'a.phi', 
              'omega', 
              'bp.y', 
              'bp.v',
              'lambda')

# Initial Values 
inits <- function() {
              list(
                N = rep(1, R), 
                beta.0 = rnorm(1),
                beta.1 = rnorm(1), 
                beta.2 = rnorm(1),
                omega = runif(1, 0, 10), 
                tau.day = runif(1, 0.1, 1),
                tau = runif(1, 0.1, 1),
                tau.p = runif(1, 0.1, 1),
                tau.day.p = runif(1, 0.1, 1),
                alpha.1 = runif(1, 0, 1), 
                a.phi = runif(1, 0, 5)
              )
}#end inits
  
  

# MCMC
n.iter = 1000 # 200000
n.burnin = 100# 60000 
n.chains = 3
n.thin = 50



# -------------------------------------------------------
# Acoustic Model
# -------------------------------------------------------
cat(" model {
  
  # -------------------------------------------
  # Priors 
  # -------------------------------------------
  beta.0 ~ dnorm(0, .01)
  beta.1 ~ dnorm(0, .01)
  beta.2 ~ dnorm(0, .01)
  alpha.0 <- logit(mu.alpha)
  mu.alpha ~ dunif(0, 1)
  alpha.1 ~ dunif(0, 1000) # Constrained to be positive
  omega ~ dunif(0, 1000)
  tau.day ~ dgamma(.01, .01)
  a.phi ~ dunif(0, 100)

  for (i in 1:n.days) {
    gamma.1[i] ~ dnorm(0, tau.day)
  }

  for (i in 1:R) {
    for (j in 1:J.A) {
      phi[i, j] ~ dgamma(a.phi, a.phi)
    }
  }
  
  # -------------------------------------------
  # Likelihood and process model 
  # -------------------------------------------
  for (i in 1:R) {
  
    # Abundance Model
    log(lambda[i]) <- beta.0 # + beta.1 * herbPdens[i] + beta.2 * woodyParea[i]
    N[i] ~ dpois(lambda[i])
    
    # Detection Model
    logit(p.a[i]) <- alpha.0 + alpha.1 * N[i]
    
    # Acoustic Data 
    for (j in 1:J[i]) {
    
    # Vocalization Model
      log(delta[i, j]) <- gamma.1[days[i, j]]
      y[i, j] ~ dbin(p.a[i], 1)
      tp[i, j] <- delta[i, j] * N[i] / (delta[i, j] * N[i] + omega)
      
      # Posterior predictive checks for Bayesian P-value
      y.pred[i, j] ~ dbin(p.a[i], 1)
      resid.y[i, j] <- pow(pow(y[i, j], 0.5) - pow(p.a[i], 0.5), 2)
      resid.y.pred[i, j] <- pow(pow(y.pred[i, j], 0.5) - pow(p.a[i], 0.5), 2)
    } # j
    
    for (j in 1:J.r[i]) {
      v[i, A.times[i, j]] ~ dpois((delta[i, A.times[i, j]] * N[i] + omega) * phi[i, A.times[i, j]] * y[i, A.times[i, j]]) T(1, )
      v.pred[i, j] ~ dpois((delta[i, A.times[i, j]] * N[i] + omega) * phi[i, A.times[i, j]] * y[i, A.times[i, j]]) T(1, )
      mu.v[i, j] <- ((delta[i, A.times[i, j]] * N[i] + omega) * phi[i, A.times[i, j]]) / (1 - exp(-1 * ((delta[i, A.times[i, j]] * N[i] + omega) * phi[i, A.times[i, j]])))
      resid.v[i, j] <- pow(pow(v[i, A.times[i, j]], 0.5) - pow(mu.v[i, j], 0.5), 2)
      resid.v.pred[i, j] <- pow(pow(v.pred[i, j], 0.5) - pow(mu.v[i, j], 0.5), 2)
    } # j
  } # i
  
  # ------------------------------------------- 
  # Manual validation 
  # -------------------------------------------
  for (i in 1:R.val) {
    for (j in 1:J.val[i]) {
      K[i, j] ~ dbin(tp[sites.a[i], j], v[sites.a[i], val.times[i, j]])
      #k[i, val.times[i, j]] ~ dhyper(K[i, j], v[sites.a[i], val.times[i, j]] - K[i, j], n[i, val.times[i, j]], 1)
      k[i, val.times[i, j]] ~ dhyper(K[i, j], max(0, v[sites.a[i], val.times[i, j]] - K[i, j]), n[i, val.times[i, j]], 1)

    } # j
  } # i
  
  # -------------------------------------------
  # Bayesian P-value
  # -------------------------------------------
  for (i in 1:R.A) {
    tmp.v[i] <- sum(resid.v[sites.a.v[i], 1:J.r[sites.a.v[i]]])
    tmp.v.pred[i] <- sum(resid.v.pred[sites.a.v[i], 1:J.r[sites.a.v[i]]])
  }
  fit.y <- sum(resid.y[sites.a, 1:J.A])
  fit.y.pred <- sum(resid.y.pred[sites.a, 1:J.A])
  fit.v <- sum(tmp.v[1:R.A])
  fit.v.pred <- sum(tmp.v.pred[1:R.A])
  bp.y <- step(fit.y.pred - fit.y)
  bp.v <- step(fit.v.pred - fit.v)
}
", fill=TRUE, file="./jags_models/ARU_mod.txt")
# ------------End Model-------------



 


# Fit Model
fm <- jags(data = Bnet14.data, 
                    inits = inits, 
                    parameters.to.save = params, 
                    model.file = "./jags_models/ARU_mod.txt", 
                    n.iter = n.iter,
                    n.burnin = n.burnin,
                    n.chains = n.chains, 
                    n.thin = n.thin,
                    parallel = TRUE,
                    n.cores = workers,
                    DIC = TRUE) 
                    
 

# Model summary
print(fm, digits = 3)

# Rhat
fm.0$Rhat

# Trace plots
mcmcplot(fm$samples)

# Bayesian P value
cat("Bayesian p-value =", fm.0$summary["bp.v",1], "\n")

# Average number of false positives detections
cat("False positives =", fm.0$summary["omega",1], "\n")
