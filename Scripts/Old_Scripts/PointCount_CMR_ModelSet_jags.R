# -------------------------------------------------------
#
#                    Load libraries
#
# -------------------------------------------------------

# Install packages (if needed)
install.packages("tidyverse")
install.packages("jagsUI")
install.packages("coda")
install.packages("mcmcplots")
install.packages("loo")

# Load library
library(tidyverse)
library(jagsUI)
library(coda)
library(mcmcplots)
library(loo)

# Check JAGs Version
# Latest version as of (Feb 5, 2025) JAGS 4.3.1.  
# Download here: https://sourceforge.net/projects/mcmc-jags/
rjags::jags.version() 

# Set seed, scientific notation options, and working directory
set.seed(123)
options(scipen = 9999)
setwd(".")

# Setting up cores
Ncores <- parallel::detectCores()
print(Ncores) # number of available cores
workers <- Ncores * 0.5 # may need to change, for low background use 80% num, for medium use 50%

# -------------------------------------------------------
#
#                    Load Data
#
# -------------------------------------------------------

# Load in capture data
pc_dat <- read.csv("./NOBO_PC_Summer2024data.csv")

# Load in site covariates
site_covs <- read.csv("./PointCount_siteCovs.csv")

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

# Extract and scale site covariates for X.abund
X.abund <- site_covs[,-c(1:4,25)] 
X.abund$woody_lrgPInx <- scale(X.abund$woody_lrgPInx)
X.abund$herb_lrgPInx  <- scale(X.abund$herb_lrgPInx)
X.abund$woody_AggInx <- scale(X.abund$woody_AggInx)
X.abund$herb_AggInx  <- scale(X.abund$herb_AggInx)
X.abund$woody_EdgDens <- scale(X.abund$woody_EdgDens)
X.abund$herb_EdgDens  <- scale(X.abund$herb_EdgDens)
X.abund$woody_Pdens <- scale(X.abund$woody_Pdens)
X.abund$herb_Pdens  <- scale(X.abund$herb_Pdens)
X.abund$woody_Npatches <- scale(X.abund$woody_Npatches)
X.abund$herb_Npatches  <- scale(X.abund$herb_Npatches)
X.abund$mnElev  <- scale(X.abund$mnElev)
print(X.abund)

# Extract and scale site covariates for X.det
X.det <- pc_dat[,c(2,7:9,18)]
X.det$Observer <- as.numeric(as.factor(X.det$Observer))
X.det$Temp.deg.F <- as.numeric(X.det$Temp.deg.F)
X.det$Wind.Beau.Code <- as.numeric(as.factor(X.det$Wind.Beau.Code)) 
X.det$Sky.Beau.Code <- as.numeric(as.factor(X.det$Sky.Beau.Code))
merged_df <- merge(site_covs, pc_CMR, by = "PointNum") # VegDens to stacked 
X.det$vegDens = as.matrix(merged_df[,c("vegDens50m")]) 
print(X.det)


# Bundle data for JAGs
data <- list(J = J,
             M = M,
             nsites = nsites,
             y = y,
             X.abund = X.abund,
             X.det = X.det,
             group = site)
 
# Check structure 
str(data)

 
# ------------------------------------------------------------------------------
#
#                             Model Fitting
#
# ------------------------------------------------------------------------------

# -------------------------------------------------------
# Model Specifications
# -------------------------------------------------------

## Distributions
# library(ggfortify)
# ggdistribution(dunif, seq(0, 1, 0.001), min = 0, max = 1) # p0
# ggdistribution(dnorm, seq(0, 0.01, 0.0001)) # alpha's and beta's 
# ggdistribution(dgamma, seq(0, 5, 0.01), shape = 0.1, rate = 0.1) # tau

# Parameters monitored
params <- c("lambda",
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
             "p_Bayes",
             "log_lik")
            


# Initial values
 inits  <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        beta3 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
)}

# MCMC 
n.iter = 100#300000
n.burnin =0 #20000
n.chains = 3 
n.thin = 10


# -------------------------------------------------------
# Model 0: Null Model
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod0.txt")
# ------------End Model-------------

# Fit Model
fm.0 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod0.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")

# -------------------------------------------------------
# Model 1: woody prp
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,1]
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod1.txt")
# ------------End Model-------------

# Fit Model
fm.1 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod1.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")

# -------------------------------------------------------
# Model 2: woody mnParea
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,4] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod2.txt")
# ------------End Model-------------

# Fit Model
fm.2 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod2.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")

# -------------------------------------------------------
# Model 3: woody ClmInx
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,6]
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod3.txt")
# ------------End Model-------------

# Fit Model
fm.3 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod3.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 4: woody shp inx
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,8] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod4.txt")
# ------------End Model-------------

# Fit Model
fm.4 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod4.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 5: woody lrg p 
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,10]
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod5.txt")
# ------------End Model-------------

# Fit Model
fm.5 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod5.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 6: woody ag inx
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,12] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod6.txt")
# ------------End Model-------------

# Fit Model
fm.6 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod6.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 7: wood edge dens
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,14] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod7.txt")
# ------------End Model-------------

# Fit Model
fm.7 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod7.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 8: woody p dens
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,16]
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod8.txt")
# ------------End Model-------------

# Fit Model
fm.8 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod8.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 9: woody n patches
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,18] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod9.txt")
# ------------End Model-------------

# Fit Model
fm.9 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod9.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 10: woody focal
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,20] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod10.txt")
# ------------End Model-------------

# Fit Model
fm.10 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod10.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 11: herp prp
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,2] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod11.txt")
# ------------End Model-------------

# Fit Model
fm.11 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod11.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 12: herb mnParea
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,4] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod12.txt")
# ------------End Model-------------

# Fit Model
fm.12 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod12.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 13: herb ClmInx 
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,7] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod13.txt")
# ------------End Model-------------

# Fit Model
fm.13 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod13.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 14: herp shp inx
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,9]
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod14.txt")
# ------------End Model-------------

# Fit Model
fm.14 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod14.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 15: herb lrg p
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,11] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod15.txt")
# ------------End Model-------------

# Fit Model
fm.15 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod15.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 16: herb ag inx 
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,13]
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod16.txt")
# ------------End Model-------------

# Fit Model
fm.16 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod16.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 17: herb edge dens
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,15] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod17.txt")
# ------------End Model-------------

# Fit Model
fm.17 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod17.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 18: herb p dens
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,17] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod18.txt")
# ------------End Model-------------

# Fit Model
fm.18 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod18.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 19: herb n patches
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,19] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod19.txt")
# ------------End Model-------------

# Fit Model
fm.19 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod19.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 20: open prp
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,3]
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod20.txt")
# ------------End Model-------------

# Fit Model
fm.20 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod20.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 21: elev 
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,21] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod21.txt")
# ------------End Model-------------

# Fit Model
fm.21 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod21.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 22: optim focal
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,22]
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod22.txt")
# ------------End Model-------------

# Fit Model
fm.22 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod22.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 23:   woody_EdgDens + herb_prp
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,14] + beta2*X.abund[s,2]
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod23.txt")
# ------------End Model-------------

# Fit Model
fm.23 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod23.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 24: woody_mnFocal30m + herb_prp
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,20] + beta2*X.abund[s,2]
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod24.txt")
# ------------End Model-------------

# Fit Model
fm.24 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod24.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 25: herb_EdgDens + woody_mnParea
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,15] + beta2*X.abund[s,4]
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod25.txt")
# ------------End Model-------------

# Fit Model
fm.25 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod25.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 26: woody_lrgPInx + herb_mnParea
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,10] + beta2*X.abund[s,4]
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod26.txt")
# ------------End Model-------------

# Fit Model
fm.26 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod26.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 27: woody_AggInx + herb_mnParea
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,12] + beta2*X.abund[s,4] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod27.txt")
# ------------End Model-------------

# Fit Model
fm.27 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod27.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 28:  woody_mnFocal30m + herb_mnParea
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,20] + beta2*X.abund[s,4] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod28.txt")
# ------------End Model-------------

# Fit Model
fm.28 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod28.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 29:  woody_lrgPInx + herb_ClmIdx  
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,10] + beta2*X.abund[s,7] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod29.txt")
# ------------End Model-------------

# Fit Model
fm.29 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod29.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 30: woody_AggInx + herb_ShpInx  
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,12] + beta2*X.abund[s,9]
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod30.txt")
# ------------End Model-------------

# Fit Model
fm.30 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod30.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 31: woody_mnFocal30m + herb_ShpInx  
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,20] + beta2*X.abund[s,9] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod31.txt")
# ------------End Model-------------

# Fit Model
fm.31 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod31.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 32: herb_AggInx + woody_lrgPInx 
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,12] + beta2*X.abund[s,10] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod32.txt")
# ------------End Model-------------

# Fit Model
fm.32 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod32.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 33: herb_Pdens + woody_lrgPInx 
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,17] + beta2*X.abund[s,10]
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod33.txt")
# ------------End Model-------------

# Fit Model
fm.33 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod33.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 34: herb_Npatches + woody_lrgPInx  
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,19] + beta2*X.abund[s,10] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod34.txt")
# ------------End Model-------------

# Fit Model
fm.34 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod34.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 35: woody_AggInx + herb_lrgPInx  
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,12] + beta2*X.abund[s,11]  
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod35.txt")
# ------------End Model-------------

# Fit Model
fm.35 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod35.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 36: herb_AggInx + woody_AggInx  
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,13] + beta2*X.abund[s,12] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod36.txt")
# ------------End Model-------------

# Fit Model
fm.36 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod36.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 37: herb_Pdens + woody_AggInx 
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,17] + beta2*X.abund[s,12] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod37.txt")
# ------------End Model-------------

# Fit Model
fm.37 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod37.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 38: herb_Npatches + woody_AggInx  
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,19] + beta2*X.abund[s,12] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod38.txt")
# ------------End Model-------------

# Fit Model
fm.38 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod38.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 39: woody_mnFocal30m + herb_AggInx  
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,20] + beta2*X.abund[s,13] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod39.txt")
# ------------End Model-------------

# Fit Model
fm.39 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod39.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 40: woody_Pdens + herb_EdgDens   
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,16] + beta2*X.abund[s,15]
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod40.txt")
# ------------End Model-------------

# Fit Model
fm.40 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod40.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 41: woody_Npatches + herb_EdgDens 
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,18] + beta2*X.abund[s,15] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod41.txt")
# ------------End Model-------------

# Fit Model
fm.41 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod41.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 42: woody_mnFocal30m + herb_Pdens   
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,20] + beta2*X.abund[s,17] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod42.txt")
# ------------End Model-------------

# Fit Model
fm.42 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod42.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 43: woody_mnFocal30m + herb_Npatches  
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,20] + beta2*X.abund[s,19] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod43.txt")
# ------------End Model-------------

# Fit Model
fm.43 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod43.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")
# -------------------------------------------------------
# Model 44: optim_mnFocal30m + herb_EdgDens   
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
    log(lambda[s]) <- beta0 + beta1*X.abund[s,22] + beta2*X.abund[s,15] 
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
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod44.txt")
# ------------End Model-------------

# Fit Model
fm.44 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_mod44.txt",
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_JAGs.RData")


# -------------------------------------------------------
# Ranking Models
# -------------------------------------------------------

# Total number of models - 44 + 1 (fm.0)
total_fits <- 3#45  # Adjust as needed
waic_values <- numeric(total_fits)  # For WAIC values
fitnames <- character(total_fits)   # For model names
for (i in 0:(total_fits - 1)) {
  model_name <- paste0("fm.", i)
  
  # Extract log-likelihood samples
  log_lik_array <- get(model_name)$sims.list$log_lik # Extract log-likelihood samples
  log_lik_matrix <- apply(log_lik_array, c(1, 2), sum)  # Summing across J
  
  # Compute WAIC
  waic_result <- loo::waic(log_lik_matrix)
  
  # Store WAIC values
  waic_values[i + 1] <- waic_result$estimates[3,1]
  fitnames[i + 1] <- model_name  # Store model name
}

# Combine into data frame
WAIC_df <- data.frame(Model = fitnames, WAIC = waic_values)

# Order by WAIC (ascending order, better fit first)
WAIC_df <- WAIC_df[order(WAIC_df$WAIC),]

# Print results
print(WAIC_df)

# Best model?
bm <- get(WAIC_df[1,1]) 
  
# Best model fit. P-value = 0.5 means good fit, = 1 or 0 is a poor fit
cat("Bayesian p-value =", bm$summary["p_Bayes",1], "\n")

# Check convergence
bm$Rhat # Rhat: less than 1.1 means good convergence
mcmcplot(bm$samples)# Visually inspect trace plots

# Model summary
summary(bm$samples)

# Save Environment
save.image(file = "./CMR_JAGs.RData")

# -------------------------------------------------------
#
#   Estimating Abundance 
#
# -------------------------------------------------------

# Combine chains
combined_chains <- as.mcmc(do.call(rbind, bm$samples))

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
dens_df <- dens_df[, c("Model", "Density")] # Switch the order of columns
head(dens_df)

# Calculate the 95% Credible Interval
ci_bounds <- quantile(dens_df$Density, probs = c(0.025, 0.975))

# Subset the data frame to 95% CI
dens_df <- subset(dens_df, Density >= ci_bounds[1] & Density <= ci_bounds[2])


# Violin plot
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
print(mean(dens_df$Density))
print(min(dens_df$Density))
print(max(dens_df$Density))


# Total abundance
mean_dens * 1096.698
LCI_dens * 1096.698
HCI_dens * 1096.698


# Save Environment
save.image(file = "CMR_JAGs.RData")


# End Script