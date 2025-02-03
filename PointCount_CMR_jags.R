# -------------------------------------------------------
#
#                    Load libraries
#
# -------------------------------------------------------

# Install packages (if needed)
# install.packages("nimble")
# install.packages("coda")

# Load library
library(tidyverse)
library(jagsUI)
library(coda)
library(mcmcplots)


# Set seed, scientific notation options, and working directory
set.seed(123)
options(scipen = 9999)
setwd(".")

# -------------------------------------------------------
#
#                    Load Data
#
# -------------------------------------------------------

# Load in capture data
pc_dat <- read.csv("./Data/Point_Count_Data/NOBO_PC_Summer2024data.csv")

# Load in site covariates
site_covs <- read.csv("./Data/Point_Count_Data/PointCount_siteCovs.csv")

# -------------------------------------------------------
#
#                   Data Wrangling
#
# -------------------------------------------------------

# Remove rows with NA values
pc_dat <- na.omit(pc_dat)

# Creating a day of year column
pc_dat$Date <- mdy(pc_dat$Date)
pc_dat$DOY <- yday(pc_dat$Date)
pc_dat$DOYscaled <- ((pc_dat$DOY - 1) / 365)


# Add a unique id column based on Point Count number, NOBO number, and Day of Year
pc_dat$UniqueID <- paste(pc_dat$PointNum, pc_dat$Survey, pc_dat$NOBOnum,  sep = "_")

# Subset the capture-mark data
pc_CMR <- pc_dat[, c( "UniqueID" ,"PointNum", "Survey", "NOBOnum",
                      "int1", "int2", "int3", "int4"),
                 drop = FALSE]


# Reorder pc_CMR by PointNum, Survey
pc_CMR <- pc_CMR %>%
  arrange(PointNum, Survey, NOBOnum)

# Take a look
head(pc_CMR)

# Subset detection intervals
y <- as.matrix(pc_CMR[,c("int1","int2","int3", "int4")])

# Number of time intervals
J <- 4

# number of sites
nsites <- 10

# Number of individuals detected
nind <- nrow(y)

# Superpopulation 
M <- 400

# Data Augmentation
y <- rbind(y, matrix(0, nrow=(M-nind), ncol=4))

# Number of "pseudo-groups" added
nz <- M - nind

# Make site (Point_Individual) into an integer
site <- as.numeric(factor(pc_CMR$PointNum))

# Add in NA for M
site <- c(site, rep(NA, M-nind))

# Take a look
print(site)

# Extract site covariates
herbPdens = as.vector(scale(site_covs[,c("herb_Pdens")])) 
woodymnParea = as.vector(site_covs[,c("woody_mnParea")]) 
mnElev = as.vector(scale(site_covs[,c("mnElev")])) 
woody_Npatches  = as.vector(scale(site_covs[,c("woody_Npatches")])) 
woody_mnFocal5mRadi = as.vector(site_covs[,c("woody_mnFocal5mRadi")]) 
optim_mnFocal5mRadi = as.vector(site_covs[,c("optim_mnFocal5mRadi")])  
 

# Extract detection covariates
obsvr <- as.matrix(as.numeric(as.factor(pc_dat[,"Observer"])))  # Convert names to numeric indices
temp <- as.matrix(as.numeric(pc_dat[,"Temp.deg.F"]))
wind <- as.matrix(as.numeric(as.factor(pc_dat[,"Wind.Beau.Code"]))) 
sky <- as.matrix(as.numeric(as.factor(pc_dat[,"Sky.Beau.Code"])))
doy <- as.matrix(pc_dat[,"DOYscaled"])

# Transforming VegDens to stacked
merged_df <- merge(site_covs, pc_CMR, by = "PointNum")
vegDens = as.matrix(merged_df[,c("vegDens50m")]) # --- match into detection matrix currently by site. needs to be by detection (dim)


head(obsvr)
head(temp)
head(wind)
head(sky)
head(doy)


# Bundle data for JAGs
data <- list(J = J,
             M = M,
             nsites = nsites,
             y = y,
             herbPdens = herbPdens,
             woodymnParea = woodymnParea,
             mnElev = mnElev,
             woody_Npatches = woody_Npatches,
             woody_mnFocal5mRadi = woody_Npatches,
             optim_mnFocal5mRadi =optim_mnFocal5mRadi,
             obsvr = obsvr,
             temp = temp,
             wind = wind,
             sky = sky,
             doy = doy,
             vegDens = vegDens,
             group = site)
 
# Take a look
str(data)

 


# ------------------------------------------------------------------------------
#
#                             Model Fitting
#
# ------------------------------------------------------------------------------

## Distributions
# library(ggfortify)
# ggdistribution(dunif, seq(0, 1, 0.001), min = 0, max = 1) # p0
# ggdistribution(dnorm, seq(0, 0.01, 0.0001)) # alpha's and beta's 
# ggdistribution(dgamma, seq(0, 5, 0.01), shape = 0.1, rate = 0.1) # tau


# -------------------------------------------------------
#
#                  Detection Models 
#
# -------------------------------------------------------




# -------------------------------------------------------
# Detection Model 0: Null Model
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  beta0 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0  
    probs[s] <- lambda[s] / sum(lambda[])
  }

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }
  
  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Detection model
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + S.raneff[group[i]] 
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file="./jags_models/CMR_detmod0.txt")
# ------------End Model-------------


# Parameters monitored
detparams.0 <- c("lambda",
              "p0", 
              "alpha0", 
              "beta0", 
              "psi",
              "S.raneff",
              "tau",
              "sigma",
              "p_Bayes")


# Initial values
detinits.0 <- function() {
  list (p0 = runif(1),
        beta0 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
detfm.0 <- jags(data = data, 
             parameters.to.save = detparams.0,
             inits = detinits.0, 
             model.file = "./jags_models/CMR_detmod0.txt",
             n.iter = 250000,
             n.burnin = 15000,
             n.chains = 3, 
             n.thin = 10,
             parallel = TRUE,
             n.cores = 8,
             DIC = TRUE)  



# Rhat
detfm.0$Rhat

# Model summary
print(detfm.0, digits = 3)

# Trace plots
mcmcplot(detfm.0$samples)

# Bayesian P value
cat("Bayesian p-value =", detfm.0$summary["p_Bayes",1], "\n")





# -------------------------------------------------------
# Detection Model 1: Observer
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 
    probs[s] <- lambda[s] / sum(lambda[])
  }

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }
  
  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Detection model
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 * obsvr[group[i],1]  + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_detmod1.txt")
# ------------End Model-------------


# Parameters monitored
detparams.1 <- c("lambda",
              "p0", 
              "alpha0", 
              "alpha1", 
              "beta0", 
              "psi",
              "S.raneff",
              "tau",
              "sigma",
              "p_Bayes")


# Initial values
detinits.1 <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
detfm.1 <- jags(data = data, 
             parameters.to.save = detparams.1,
             inits = detinits.1, 
             model.file = "./jags_models/CMR_detmod1.txt",
             n.iter = 250000,
             n.burnin = 15000,
             n.chains = 3, 
             n.thin = 10,
             parallel = TRUE,
             n.cores = 8,
             DIC = TRUE)  



# Rhat
detfm.1$Rhat

# Model summary
print(detfm.1, digits = 3)

# Trace plots
mcmcplot(detfm.1$samples)

# Bayesian P value
cat("Bayesian p-value =", detfm.1$summary["p_Bayes",1], "\n")



# -------------------------------------------------------
# Detection Model 2: Temperature
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0  
    probs[s] <- lambda[s] / sum(lambda[])
  }

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }
  
  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Detection model
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 * temp[group[i],1]  + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_detmod2.txt")
# ------------End Model-------------


# Parameters monitored
detparams.2 <- c("lambda",
              "p0", 
              "alpha0", 
              "alpha1", 
              "beta0", 
              "psi",
              "S.raneff",
              "tau",
              "sigma",
              "p_Bayes")


# Initial values
detinits.2 <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
detfm.2 <- jags(data = data, 
             parameters.to.save = detparams.2,
             inits = detinits.2, 
             model.file = "./jags_models/CMR_detmod2.txt",
             n.iter = 250000,
             n.burnin = 15000,
             n.chains = 3, 
             n.thin = 10,
             parallel = TRUE,
             n.cores = 8,
             DIC = TRUE)  



# Rhat
detfm.2$Rhat

# Model summary
print(detfm.2, digits = 3)

# Trace plots
mcmcplot(detfm.2$samples)

# Bayesian P value
cat("Bayesian p-value =", detfm.2$summary["p_Bayes",1], "\n")





# -------------------------------------------------------
# Detection Model 3: Wind
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 
    probs[s] <- lambda[s] / sum(lambda[])
  }

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }
  
  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Detection model
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 * wind[group[i],1]  + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_detmod3.txt")
# ------------End Model-------------


# Parameters monitored
detparams.3 <- c("lambda",
              "p0", 
              "alpha0", 
              "alpha1", 
              "beta0", 
              "psi",
              "S.raneff",
              "tau",
              "sigma",
              "p_Bayes")


# Initial values
detinits.3 <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
detfm.3 <- jags(data = data, 
             parameters.to.save = detparams.3,
             inits = detinits.3, 
             model.file = "./jags_models/CMR_detmod3.txt",
             n.iter = 250000,
             n.burnin = 15000,
             n.chains = 3, 
             n.thin = 10,
             parallel = TRUE,
             n.cores = 8,
             DIC = TRUE)  



# Rhat
detfm.3$Rhat

# Model summary
print(detfm.3, digits = 3)

# Trace plots
mcmcplot(detfm.3$samples)

# Bayesian P value
cat("Bayesian p-value =", detfm.3$summary["p_Bayes",1], "\n")



# -------------------------------------------------------
# Detection Model 4: Sky
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 
    probs[s] <- lambda[s] / sum(lambda[])
  }

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }
  
  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Detection model
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 * sky[group[i],1]  + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/detCMR_mod4.txt")
# ------------End Model-------------


# Parameters monitored
detparams.4 <- c("lambda",
              "p0", 
              "alpha0", 
              "alpha1", 
              "beta0", 
              "psi",
              "S.raneff",
              "tau",
              "sigma",
              "p_Bayes")


# Initial values
detinits.4 <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
detfm.4 <- jags(data = data, 
             parameters.to.save = detparams.4,
             inits = detinits.4, 
             model.file = "./jags_models/detCMR_mod4.txt",
             n.iter = 250000,
             n.burnin = 15000,
             n.chains = 3, 
             n.thin = 10,
             parallel = TRUE,
             n.cores = 8,
             DIC = TRUE)  



# Rhat
detfm.4$Rhat

# Model summary
print(detfm.4, digits = 4)

# Trace plots
mcmcplot(detfm.4$samples)

# Bayesian P value
cat("Bayesian p-value =", detfm.4$summary["p_Bayes",1], "\n")


# -------------------------------------------------------
# Detection Model 5: Day of Year
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0  
    probs[s] <- lambda[s] / sum(lambda[])
  }

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }
  
  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Detection model
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 * doy[group[i],1]  + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_detmod5.txt")
# ------------End Model-------------


# Parameters monitored
detparams.5 <- c("lambda",
              "p0", 
              "alpha0", 
              "alpha1", 
              "beta0", 
              "beta1",
              "beta2",
              "psi",
              "S.raneff",
              "tau",
              "sigma",
              "p_Bayes")


# Initial values
detinits.5 <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
detfm.5 <- jags(data = data, 
             parameters.to.save = detparams.5,
             inits = detinits.5, 
             model.file = "./jags_models/CMR_detmod5.txt",
             n.iter = 250000,
             n.burnin = 15000,
             n.chains = 3, 
             n.thin = 10,
             parallel = TRUE,
             n.cores = 8,
             DIC = TRUE)  



# Rhat
detfm.5$Rhat

# Model summary
print(detfm.5, digits = 4)

# Trace plots
mcmcplot(detfm.5$samples)

# Bayesian P value
cat("Bayesian p-value =", detfm.5$summary["p_Bayes",1], "\n")



# -------------------------------------------------------
# Detection Model 6: Vegetation Density
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0  
    probs[s] <- lambda[s] / sum(lambda[])
  }

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }
  
  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Detection model
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 * vegDens[group[i],1]  + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_detmod6.txt")
# ------------End Model-------------


# Parameters monitored
detparams.6 <- c("lambda",
                 "p0", 
                 "alpha0", 
                 "alpha1", 
                 "beta0", 
                 "beta1",
                 "beta2",
                 "psi",
                 "S.raneff",
                 "tau",
                 "sigma",
                 "p_Bayes")


# Initial values
detinits.6 <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
detfm.6 <- jags(data = data, 
                parameters.to.save = detparams.6,
                inits = detinits.6, 
                model.file = "./jags_models/CMR_detmod6.txt",
                n.iter = 250000,
                n.burnin = 15000,
                n.chains = 3, 
                n.thin = 10,
                parallel = TRUE,
                n.cores = 8,
                DIC = TRUE)  



# Rhat
detfm.6$Rhat

# Model summary
print(detfm.6, digits = 4)

# Trace plots
mcmcplot(detfm.6$samples)

# Bayesian P value
cat("Bayesian p-value =", detfm.6$summary["p_Bayes",1], "\n")


# -------------------------------------------------------
# Detection Model 7: Wind + Vegetation Density
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  alpha1 ~ dnorm(0, 0.01)
  alpha2 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0  
    probs[s] <- lambda[s] / sum(lambda[])
  }

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }
  
  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Detection model
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 * vegDens[group[i],1] + alpha2 * wind[group[i],1]  + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_detmod7.txt")
# ------------End Model-------------


# Parameters monitored
detparams.7 <- c("lambda",
                 "p0", 
                 "alpha0", 
                 "alpha1",
                 "alpha2",
                 "beta0", 
                 "beta1",
                 "beta2",
                 "psi",
                 "S.raneff",
                 "tau",
                 "sigma",
                 "p_Bayes")


# Initial values
detinits.7 <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        alpha2 = 0,
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
detfm.7 <- jags(data = data, 
                parameters.to.save = detparams.7,
                inits = detinits.7, 
                model.file = "./jags_models/CMR_detmod7.txt",
                n.iter = 250000,
                n.burnin = 15000,
                n.chains = 3, 
                n.thin = 10,
                parallel = TRUE,
                n.cores = 8,
                DIC = TRUE)  



# Rhat
detfm.7$Rhat

# Model summary
print(detfm.7, digits = 4)

# Trace plots
mcmcplot(detfm.7$samples)

# Bayesian P value
cat("Bayesian p-value =", detfm.7$summary["p_Bayes",1], "\n")




# -------------------------------------------------------
# Detection Model Ranking
# -------------------------------------------------------

# Bayesian P value
cat("Bayesian p-value =", detfm.0$summary["p_Bayes",1], "\n")
cat("Bayesian p-value =", detfm.1$summary["p_Bayes",1], "\n")
cat("Bayesian p-value =", detfm.2$summary["p_Bayes",1], "\n")
cat("Bayesian p-value =", detfm.3$summary["p_Bayes",1], "\n")
cat("Bayesian p-value =", detfm.4$summary["p_Bayes",1], "\n")
cat("Bayesian p-value =", detfm.5$summary["p_Bayes",1], "\n")
cat("Bayesian p-value =", detfm.6$summary["p_Bayes",1], "\n")
cat("Bayesian p-value =", detfm.7$summary["p_Bayes",1], "\n")

# Combine the DIC values
detDIC_values <- c(detfm.0$DIC, 
                   detfm.1$DIC, 
                   detfm.2$DIC, 
                   detfm.3$DIC, 
                   detfm.4$DIC, 
                   detfm.5$DIC,
                   detfm.6$DIC,
                   detfm.7$DIC)

detDIC_df <- data.frame(Model = c("detfm.0", 
                                  "detfm.1", 
                                  "detfm.2", 
                                  "detfm.3", 
                                  "detfm.4", 
                                  "detfm.5",
                                  "detfm.6",
                                  "detfm.7",),
                        DIC = detDIC_values)

# Rank the values from lowest to highest
detDIC_df <- detDIC_df[order(detDIC_df$DIC), ]
print(detDIC_df)



# -------------------------------------------------------
#
#                  Abundance Models 
#
# -------------------------------------------------------




# -------------------------------------------------------
# Abundance Model 0: Null Model
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 
    probs[s] <- lambda[s] / sum(lambda[])
  }

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }
  
  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Detection model
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  wind[group[i],1] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_abundmod0.txt")
# ------------End Model-------------


# Parameters monitored
abundparams.0 <- c("lambda",
              "p0", 
              "alpha0", 
              "alpha1", 
              "beta0", 
              "psi",
              "S.raneff",
              "tau",
              "sigma",
              "p_Bayes")


# Initial values
abundinits.0 <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
abundfm.0 <- jags(data = data, 
             parameters.to.save = abundparams.0,
             inits = abundinits.0, 
             model.file = "./jags_models/CMR_abundmod0.txt",
             n.iter = 250000,
             n.burnin = 15000,
             n.chains = 3, 
             n.thin = 10,
             parallel = TRUE,
             n.cores = 8,
             DIC = TRUE)  



# Rhat
abundfm.0$Rhat

# Model summary
print(abundfm.0, digits = 4)

# Trace plots
mcmcplot(abundfm.0$samples)

# Bayesian P value
cat("Bayesian p-value =", abundfm.0$summary["p_Bayes",1], "\n")


# -------------------------------------------------------
# Abundance Model 1: herbPdens
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1 * herbPdens[s]  
    probs[s] <- lambda[s] / sum(lambda[])
  }

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }
  
  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Detection model
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  wind[group[i],1]  + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_abundmod1.txt")
# ------------End Model-------------


# Parameters monitored
abundparams.1 <- c("lambda",
              "p0", 
              "alpha0", 
              "alpha1", 
              "beta0", 
              "beta1",
              "psi",
              "S.raneff",
              "tau",
              "sigma",
              "p_Bayes")


# Initial values
abundinits.1 <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        beta1 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
abundfm.1 <- jags(data = data, 
             parameters.to.save = abundparams.1,
             inits = abundinits.1, 
             model.file = "./jags_models/CMR_abundmod1.txt",
             n.iter = 250000,
             n.burnin = 15000,
             n.chains = 3, 
             n.thin = 10,
             parallel = TRUE,
             n.cores = 8,
             DIC = TRUE)  



# Rhat
abundfm.1$Rhat

# Model summary
print(abundfm.1, digits = 4)

# Trace plots
mcmcplot(abundfm.1$samples)

# Bayesian P value
cat("Bayesian p-value =", abundfm.1$summary["p_Bayes",1], "\n")



# -------------------------------------------------------
# Abundance Model 2: herbPdens + woodymnParea
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1 * herbPdens[s] + beta2 * woodymnParea[s]
    probs[s] <- lambda[s] / sum(lambda[])
  }

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }
  
  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Detection model
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  wind[group[i],1] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_abundmod2.txt")
# ------------End Model-------------


# Parameters monitored
abundparams.2 <- c("lambda",
                   "p0", 
                   "alpha0", 
                   "alpha1", 
                   "beta0", 
                   "beta1",
                   "beta2",
                   "psi",
                   "S.raneff",
                   "tau",
                   "sigma",
                   "p_Bayes")


# Initial values
abundinits.2 <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
abundfm.2 <- jags(data = data, 
                  parameters.to.save = abundparams.2,
                  inits = abundinits.2, 
                  model.file = "./jags_models/CMR_abundmod2.txt",
                  n.iter = 250000,
                  n.burnin = 15000,
                  n.chains = 3, 
                  n.thin = 10,
                  parallel = TRUE,
                  n.cores = 8,
                  DIC = TRUE)  



# Rhat
abundfm.2$Rhat

# Model summary
print(abundfm.2, digits = 4)

# Trace plots
mcmcplot(abundfm.2$samples)

# Bayesian P value
cat("Bayesian p-value =", abundfm.2$summary["p_Bayes",1], "\n")


# -------------------------------------------------------
# Abundance Model 3: herbPdens + woodymnParea + mnElev
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  beta3 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1 * herbPdens[s] + beta2 * woodymnParea[s] + beta3 * mnElev[s]
    probs[s] <- lambda[s] / sum(lambda[])
  }

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }
  
  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Detection model
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 * wind[group[i],1]  + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_abundmod3.txt")
# ------------End Model-------------


# Parameters monitored
abundparams.3 <- c("lambda",
                   "p0", 
                   "alpha0", 
                   "alpha1", 
                   "beta0", 
                   "beta1",
                   "beta2",
                   "beta3",
                   "psi",
                   "S.raneff",
                   "tau",
                   "sigma",
                   "p_Bayes")


# Initial values
abundinits.3 <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        beta3 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
abundfm.3 <- jags(data = data, 
                  parameters.to.save = abundparams.3,
                  inits = abundinits.3, 
                  model.file = "./jags_models/CMR_abundmod3.txt",
                  n.iter = 250000,
                  n.burnin = 15000,
                  n.chains = 3, 
                  n.thin = 10,
                  parallel = TRUE,
                  n.cores = 8,
                  DIC = TRUE)  



# Rhat
abundfm.3$Rhat

# Model summary
print(abundfm.3, digits = 4)

# Trace plots
mcmcplot(abundfm.3$samples)

# Bayesian P value
cat("Bayesian p-value =", abundfm.3$summary["p_Bayes",1], "\n")



# -------------------------------------------------------
# Abundance Model 4: herbPdens + woody_Npatches
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1 * herbPdens[s] + beta2 * woody_Npatches[s]
    probs[s] <- lambda[s] / sum(lambda[])
  }

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }
  
  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Detection model
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 * wind[group[i],1]  + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_abundmod4.txt")
# ------------End Model-------------


# Parameters monitored
abundparams.4 <- c("lambda",
                   "p0", 
                   "alpha0", 
                   "alpha1", 
                   "beta0", 
                   "beta1",
                   "beta2",
                   "psi",
                   "S.raneff",
                   "tau",
                   "sigma",
                   "p_Bayes")


# Initial values
abundinits.4 <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
abundfm.4 <- jags(data = data, 
                  parameters.to.save = abundparams.4,
                  inits = abundinits.4, 
                  model.file = "./jags_models/CMR_abundmod4.txt",
                  n.iter = 250000,
                  n.burnin = 15000,
                  n.chains = 3, 
                  n.thin = 10,
                  parallel = TRUE,
                  n.cores = 8,
                  DIC = TRUE)  



# Rhat
abundfm.4$Rhat

# Model summary
print(abundfm.4, digits = 4)

# Trace plots
mcmcplot(abundfm.4$samples)

# Bayesian P value
cat("Bayesian p-value =", abundfm.4$summary["p_Bayes",1], "\n")




# -------------------------------------------------------
# Abundance Model 5: herbPdens + woody_mnFocal5mRadi
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1 * herbPdens[s] + beta2 * woody_mnFocal5mRadi[s]
    probs[s] <- lambda[s] / sum(lambda[])
  }

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }
  
  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Detection model
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 * wind[group[i],1]  + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_abundmod5.txt")
# ------------End Model-------------


# Parameters monitored
abundparams.5 <- c("lambda",
                   "p0", 
                   "alpha0", 
                   "alpha1", 
                   "beta0", 
                   "beta1",
                   "beta2",
                   "psi",
                   "S.raneff",
                   "tau",
                   "sigma",
                   "p_Bayes")


# Initial values
abundinits.5 <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
abundfm.5 <- jags(data = data, 
                  parameters.to.save = abundparams.5,
                  inits = abundinits.5, 
                  model.file = "./jags_models/CMR_abundmod5.txt",
                  n.iter = 250000,
                  n.burnin = 15000,
                  n.chains = 3, 
                  n.thin = 10,
                  parallel = TRUE,
                  n.cores = 8,
                  DIC = TRUE)  



# Rhat
abundfm.5$Rhat

# Model summary
print(abundfm.5, digits = 4)

# Trace plots
mcmcplot(abundfm.5$samples)

# Bayesian P value
cat("Bayesian p-value =", abundfm.5$summary["p_Bayes",1], "\n")




# -------------------------------------------------------
# Abundance Model 6: herbPdens + optim_mnFocal5mRadi
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1 * herbPdens[s] + beta2 * optim_mnFocal5mRadi[s]
    probs[s] <- lambda[s] / sum(lambda[])
  }

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }
  
  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Detection model
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 * wind[group[i],1]  + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_abundmod6.txt")
# ------------End Model-------------


# Parameters monitored
abundparams.6 <- c("lambda",
                   "p0", 
                   "alpha0", 
                   "alpha1", 
                   "beta0", 
                   "beta1",
                   "beta2",
                   "psi",
                   "S.raneff",
                   "tau",
                   "sigma",
                   "p_Bayes")


# Initial values
abundinits.6 <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
abundfm.6 <- jags(data = data, 
                  parameters.to.save = abundparams.6,
                  inits = abundinits.6, 
                  model.file = "./jags_models/CMR_abundmod6.txt",
                  n.iter = 250000,
                  n.burnin = 15000,
                  n.chains = 3, 
                  n.thin = 10,
                  parallel = TRUE,
                  n.cores = 8,
                  DIC = TRUE)  



# Rhat
abundfm.6$Rhat

# Model summary
print(abundfm.6, digits = 4)

# Trace plots
mcmcplot(abundfm.6$samples)

# Bayesian P value
cat("Bayesian p-value =", abundfm.6$summary["p_Bayes",1], "\n")





# -------------------------------------------------------
# Abundance Model 7: herbPdens  + mnElev
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1 * herbPdens[s] + beta2 * mnElev[s]
    probs[s] <- lambda[s] / sum(lambda[])
  }

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }
  
  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Detection model
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 * wind[group[i],1]  + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_abundmod7.txt")
# ------------End Model-------------


# Parameters monitored
abundparams.7 <- c("lambda",
                   "p0", 
                   "alpha0", 
                   "alpha1", 
                   "beta0", 
                   "beta1",
                   "beta2",
                   "psi",
                   "S.raneff",
                   "tau",
                   "sigma",
                   "p_Bayes")


# Initial values
abundinits.7 <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
abundfm.7 <- jags(data = data, 
                  parameters.to.save = abundparams.7,
                  inits = abundinits.7, 
                  model.file = "./jags_models/CMR_abundmod7.txt",
                  n.iter = 250000,
                  n.burnin = 15000,
                  n.chains = 3, 
                  n.thin = 10,
                  parallel = TRUE,
                  n.cores = 8,
                  DIC = TRUE)  



# Rhat
abundfm.7$Rhat

# Model summary
print(abundfm.7, digits = 4)

# Trace plots
mcmcplot(abundfm.7$samples)

# Bayesian P value
cat("Bayesian p-value =", abundfm.7$summary["p_Bayes",1], "\n")



# -------------------------------------------------------
# Abundance Model Ranking
# -------------------------------------------------------

# Combine the DIC values
abundDIC_values <- c(abundfm.0$DIC, 
                     abundfm.1$DIC, 
                     abundfm.2$DIC, 
                     abundfm.3$DIC, 
                     abundfm.4$DIC, 
                     abundfm.5$DIC,
                     abundfm.6$DIC)

abundDIC_df <- data.frame(Model = c("abundfm.0", 
                                    "abundfm.1", 
                                    "abundfm.2", 
                                    "abundfm.3", 
                                    "abundfm.4", 
                                    "abundfm.5",
                                    "abundfm.6"),
                          DIC = abundDIC_values)

# Rank the values from lowest to highest
abundDIC_df <- abundDIC_df[order(abundDIC_df$DIC), ]
print(abundDIC_df)

# -------------------------------------------------------
#
#   Estimating Abundance 
#
# -------------------------------------------------------


# Combine chains
combined_chains <- as.mcmc(do.call(rbind, fm.4$samples))

# Extract lambda estimates
lambda_columns <- grep("^lambda\\[", colnames(combined_chains))
lambda_samples <- combined_chains[ ,lambda_columns]

# mean abundance 
lambda_tot <- rowSums(lambda_samples)

# Area in hectares
area <- pi*(200^2)/10000

# Getting density
dens_df <- as.data.frame(lambda_tot/area)

# Summarize by row
colnames(dens_df)[1] <- "Density"
dens_df[,2] <- "PC CMR"
colnames(dens_df)[2] <- "Model"
dens_df <- dens_df[, c("Model", "Density")]# Switch the order of columns
head(dens_df)

# Calculate the 95% Credible Interval
ci_bounds <- quantile(dens_df$Density, probs = c(0.025, 0.975))


# Subset the data frame to 95% CI
dens_df <- subset(dens_df, Density >= ci_bounds[1] & Density <= ci_bounds[2])


# Plot
ggplot(dens_df, aes(x = Model, y = Density, fill = Model)) +
  geom_violin() +
  geom_boxplot(aes(x = Model, y = Density),
               width = 0.2, position = position_dodge(width = 0.8)) +
  labs(
    title = "",
    x = "Model",
    y = "Density (N/hectare)") +
  scale_y_continuous(limits = c(0, 10),
                     breaks = seq(0, 10, by = 5),
                     labels = scales::comma) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )


# total density
mean_dens <- mean(dens_df$Density)
LCI_dens <- min(dens_df$Density)
HCI_dens <- max(dens_df$Density)

print(mean_dens)
print(LCI_dens)
print(HCI_dens)

# total abundance
mean_dens * 1096.698
LCI_dens * 1096.698
HCI_dens * 1096.698


# Export Density data frame
saveRDS(dens_df, "./Data/Fitted_Models/PC_CMR_Dens.rds")

















