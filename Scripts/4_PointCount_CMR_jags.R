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

# Check JAGs Version
# Latest version as of (5 Feb 2025) JAGS 4.3.1.  
# Download here: https://sourceforge.net/projects/mcmc-jags/
rjags::jags.version() 

# Set seed, scientific notation options, and working directory
set.seed(123)
options(scipen = 9999)
setwd(".")

# Setting up cores
Ncores <- parallel::detectCores()
print(Ncores) # Number of available cores
workers <- Ncores  * 0.5 # For low background use 80%, for medium use 50% of Ncores

# Source custom function for checking Rhat values > 1.1
source("./Scripts/Rhat_check_function.R")

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
print(nind)
M <- nind + 100 # ensure that M is larger than number detected
print(M)

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
X.abund <- site_covs[,-c(1:4,25:26)] 
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

# Area surveyed 
area <- pi * (200^2) / 4046.86  # in acres


# Bundle data for JAGs
data <- list(J = J,
             M = M,
             nsites = nsites,
             y = y,
             X.abund = as.matrix(X.abund),
             X.det = as.matrix(X.det),
             group = site,
             area = area)
 
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
            "N",
            "N_tot",
            "p0", 
            'p',
            "alpha0", 
            "alpha1", 
            "beta0", 
            "beta1",
            "beta2",
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
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
)}

# MCMC 
n.iter = 300000
n.burnin = 20000
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
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }

  # Abundance model 
  for(s in 1:nsites){
    log(lambda[s]) <- beta0
    probs[s] <- lambda[s] / sum(lambda[])
    N[s] <- sum(group[] == s)  # Estimated abundance at site s
  }


  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Detection model = intercept + wind + site random effect (for pseudoreplication)
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
  
  # Total estimated population
  N_tot <- sum(z[])
  
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

# -------------------------------------------------------
# Model 1: Woody_EdgDens + herb_prp  
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

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }

  # Abundance model 
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1*X.abund[s, 14] + beta2*X.abund[s, 2]
    probs[s] <- lambda[s] / sum(lambda[])
    N[s] <- sum(group[] == s)  # Estimated abundance at site s
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
  
  # Total estimated population
  N_tot <- sum(z[])
  
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

# -------------------------------------------------------
# Model 2: woody_shrubmnFocal30m + herb_prp  
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

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }

  # Abundance model 
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1*X.abund[s, 21] + beta2*X.abund[s, 2]
    probs[s] <- lambda[s] / sum(lambda[])
    N[s] <- sum(group[] == s)  # Estimated abundance at site s
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
  
  # Total estimated population
  N_tot <- sum(z[])
  
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


# -------------------------------------------------------
# Model 3: woody_lrgPInx + herb_ClmIdx 
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

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }

  # Abundance model 
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1*X.abund[s, 10] + beta2*X.abund[s, 7]
    probs[s] <- lambda[s] / sum(lambda[])
    N[s] <- sum(group[] == s)  # Estimated abundance at site s
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
  
  # Total estimated population
  N_tot <- sum(z[])
  
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

# -------------------------------------------------------
# Model 4: woody_AggInx + herb_ShpInx 
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

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }

  # Abundance model 
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1*X.abund[s, 12] + beta2*X.abund[s, 9]
    probs[s] <- lambda[s] / sum(lambda[])
    N[s] <- sum(group[] == s)  # Estimated abundance at site s
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
  
  # Total estimated population
  N_tot <- sum(z[])
  
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

# -------------------------------------------------------
# Model 5: woody_shrubmnFocal30m - herb_lrgPInx 
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

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }

  # Abundance model 
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1*X.abund[s, 21] + beta2*X.abund[s, 11]
    probs[s] <- lambda[s] / sum(lambda[])
    N[s] <- sum(group[] == s)  # Estimated abundance at site s
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
  
  # Total estimated population
  N_tot <- sum(z[])
  
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


# -------------------------------------------------------
# Check Convergence Models
# -------------------------------------------------------

# Rhat values
check_rhat(fm.0$Rhat, threshold = 1.1)
check_rhat(fm.1$Rhat, threshold = 1.1)
check_rhat(fm.2$Rhat, threshold = 1.1)
check_rhat(fm.3$Rhat, threshold = 1.1)
check_rhat(fm.4$Rhat, threshold = 1.1)
check_rhat(fm.5$Rhat, threshold = 1.1)

# Trace plots
# mcmcplot(fm.0$samples) 
# mcmcplot(fm.1$samples) 
# mcmcplot(fm.2$samples)
# mcmcplot(fm.3$samples) 
# mcmcplot(fm.4$samples) 
# mcmcplot(fm.5$samples)

# Save Environment
save.image(file = "./Data/Model_Environments/CMR_bm_JAGs.RData")

# -------------------------------------------------------
# Ranking Models
# -------------------------------------------------------

# Total number of models = 5 fits + 1 null
total_fits <- 6  # Adjust as needed
waic_values <- numeric(total_fits)  # WAIC values
fitnames <- character(total_fits)   # Model names
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

# Model summary
summary(bm$samples)


# -------------------------------------------------------
#
#   Beta Estimates and Covariate Effects 
#
# -------------------------------------------------------

# Combine chains
combined_chains <- as.mcmc(do.call(rbind, fm.2$samples))

# -------------------------------------------------------
# Beta Estimates
# -------------------------------------------------------

# Extract beta estimates
beta0_samples <- combined_chains[, "beta0"]
beta1_samples <- combined_chains[, "beta1"]
beta2_samples <- combined_chains[, "beta2"]

# Means
beta0 <- mean(beta0_samples)
beta1 <- mean(beta1_samples)
beta2 <- mean(beta2_samples)
 
# Credible Intervals
# beta0_CI_lower <- quantile(beta0_samples, probs = 0.025)
# beta0_CI_upper <- quantile(beta0_samples, probs = 0.975)
# 
# beta1_CI_lower <- quantile(beta1_samples, probs = 0.025)
# beta1_CI_upper <- quantile(beta1_samples, probs = 0.975)
# 
# beta2_CI_lower <- quantile(beta2_samples, probs = 0.025)
# beta2_CI_upper <- quantile(beta2_samples, probs = 0.975)


# Compute 95% CI for each beta
beta_df <- data.frame(
  value = c(beta0_samples, beta1_samples, beta2_samples),
  parameter = rep(c("beta0", "beta1", "beta2"), each = length(beta0_samples))
) %>%
  group_by(parameter) %>%
  filter(value >= quantile(value, 0.025) & value <= quantile(value, 0.975))  # Keep only values within 95% CI

# Add model
beta_df$Model <- "PC CMR"

# Plot density
ggplot(beta_df, aes(x = value, fill = parameter)) +
  geom_density(alpha = 0.5) +
  labs(title = "Posterior Density Plots for Beta Estimates", x = "Estimate", y = "Density") +
  theme_minimal()

# Create violin plot
ggplot(beta_df, aes(x = parameter, y = value, fill = parameter)) +
  geom_violin(alpha = 0.5, trim = TRUE) +  # Violin plot with smoothing
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +  # Line at y = 0
  labs(title = "Violin Plots for Beta Estimates", x
       = "Parameter", 
       y = "Estimate") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")  # Nice color scheme

# Export beta dataframe
saveRDS(beta_df, "./Data/Fitted_Models/PC_CMR_beta_df.rds")

# -------------------------------------------------------
# Covariate Effects
# -------------------------------------------------------
 
herbCov_name <- "herb_prp"
woodyCov_name <- 'Woody_shrub_mnFocal30m'

# Create a prediction of covariate values
herb_cov_pred_vals <- seq(min(site_covs[, herbCov_name]), max(site_covs[, herbCov_name]), length.out = 1000)
woody_cov_pred_vals <- seq(min(site_covs[, woodyCov_name]), max(site_covs[, woodyCov_name]), length.out = 1000)

# Scale to have a mean of 0
herb_cov_pred_vals_scaled <- (herb_cov_pred_vals - mean(site_covs[, herbCov_name])) / sd(site_covs[, herbCov_name])
woody_cov_pred_vals_scaled <- (woody_cov_pred_vals - mean(site_covs[, woodyCov_name])) / sd(site_covs[, woodyCov_name])

# Matrices for storing predictions
herb_preds <- matrix(NA, nrow = length(beta1_samples), ncol = length(herb_cov_pred_vals_scaled))
woody_preds <- matrix(NA, nrow = length(beta2_samples), ncol = length(woody_cov_pred_vals_scaled))

# Generate predictions
for (i in 1:length(beta0_samples)) {
  herb_preds[i, ] <- beta0_samples[i] + beta1_samples[i] * herb_cov_pred_vals_scaled
  woody_preds[i, ] <- beta0_samples[i] + beta2_samples[i] * woody_cov_pred_vals_scaled
}

# Calculate credible intervals
herb_preds_CI_lower <- apply(herb_preds, 2, quantile, probs = 0.025)
herb_preds_CI_upper <- apply(herb_preds, 2, quantile, probs = 0.975)

woody_preds_CI_lower <- apply(woody_preds, 2, quantile, probs = 0.025)
woody_preds_CI_upper <- apply(woody_preds, 2, quantile, probs = 0.975)

# Calculate mean predictions
herb_preds_mean <- apply(herb_preds, 2, mean)
woody_preds_mean <- apply(woody_preds, 2, mean)

# Combine into a single data frame
herbaceous_data <- data.frame(
  herb_cov_pred_vals_scaled = herb_cov_pred_vals_scaled,
  herb_preds_mean = herb_preds_mean,
  herb_preds_CI_lower = herb_preds_CI_lower,
  herb_preds_CI_upper = herb_preds_CI_upper
)

woody_data <- data.frame(
  woody_cov_pred_vals_scaled = woody_cov_pred_vals_scaled,
  woody_preds_mean = woody_preds_mean,
  woody_preds_CI_lower = woody_preds_CI_lower,
  woody_preds_CI_upper = woody_preds_CI_upper
)



# Plot Herbaceous Patch
herbcovEff_plot <- ggplot(herbaceous_data, aes(x = herb_cov_pred_vals_scaled, y = herb_preds_mean)) +
                  geom_line(color = "lightgreen", linewidth = 1.5) +  # Line plot
                  geom_ribbon(aes(ymin = herb_preds_CI_lower, ymax = herb_preds_CI_upper), 
                              fill = rgb(0.2, 0.6, 0.2, 0.2), alpha = 0.5) +  # CI shading
                  labs(x = "Herbaceous Covariate", 
                       y = "Predicted Effect", 
                       title = "Predicted Effect of Herbaceous Proportion") +
                  theme_minimal() +
                  theme(panel.grid = element_blank())
# View
herbcovEff_plot

# Export                
ggsave(plot = herbcovEff_plot, "Figures/CMR_HerbCovEffect_plot.jpeg",  
       width = 8, height = 5, dpi = 300) 


# Plot Woody Largest Patch Index
woodycovEff_plot <- ggplot(woody_data, aes(x = woody_cov_pred_vals_scaled, y = woody_preds_mean)) +
                  geom_line(color = "forestgreen", linewidth = 1.5) +  # Line plot
                  geom_ribbon(aes(ymin = woody_preds_CI_lower, 
                                  ymax = woody_preds_CI_upper), 
                              fill = rgb(0.2, 0.6, 0.2, 0.2), alpha = 0.5) +  # CI shading
                  labs(x = "Woody Covariate", 
                       y = "", 
                       title = "Predicted Effect of Shrub Focal") +
                  theme_minimal() +
                  theme(panel.grid = element_blank())
# View
woodycovEff_plot

# Export                
ggsave(plot = woodycovEff_plot, "Figures/CMR_WoodyCovEffect_plot.jpeg",  
       width = 8, height = 5, dpi = 300) 

# -------------------------------------------------------
#
#   Estimating Abundance 
#
# -------------------------------------------------------


# Extracting Abundance
Ntot_samples <- combined_chains[ ,"N_tot"]

# Ntotal is the abundance based on 10 point counts at a radius of 200m.
# To correct for density, Ntotal needs to be divided by 10 * area surveyed
area <- pi * (250^2) / 4046.86  # Area in acres
dens_samples <- Ntot_samples / (area * 10)

# Create data frame for density
dens_df <- data.frame(Model = rep("PC CMR", length(dens_samples)), Density = dens_samples)
colnames(dens_df)[2] <- "Density"
head(dens_df)

# Calculate the mean and 95% Credible Interval
dens_summary <- dens_df %>%
  group_by(Model) %>%
  summarise(
    Mean = mean(Density),
    Lower_CI = quantile(Density, 0.025),
    Upper_CI = quantile(Density, 0.975)
  )

# Subset the data within the 95% credible interval
dens_df <- dens_df[dens_df$Density >= dens_summary$Lower_CI 
                   & dens_df$Density <= dens_summary$Upper_CI, ]


# Plot density violin
ggplot(dens_df, aes(x = Model, y = Density, fill = Model)) + 
  geom_violin(trim = FALSE, alpha = 0.6, adjust = 5) +  # Adjust bandwidth for smoothing
  labs(x = "Model", y = "Density (N/acre)") +
  scale_fill_manual(values = c("PC CMR" = "orange", 
                               "PC HDS" = "purple", 
                               "AV Bnet" = "blue")) +  # Custom colors
  scale_y_continuous(limits = c(0, 0.5),
                     breaks = seq(0, 0.5, by = 0.25),
                     labels = scales::comma) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),  
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
        panel.grid = element_blank(),
        legend.position = "none")  # Removes legend




# Total density
print(dens_summary)




# Getting total abundance
abund_summary <- dens_summary
abund_summary[,2:4] <- abund_summary[,2:4] * 2710

# Plot Abundance - Violin
abund_df <- dens_df
abund_df$Density <- abund_df$Density * 2710

ggplot(abund_df, aes(x = Model, y = Density, fill = Model)) + 
  geom_violin(trim = FALSE, alpha = 0.6, adjust = 5) +  # Adjust bandwidth for smoothing
  labs(x = "Model", y = "Density (N/acre)") +
  scale_fill_manual(values = c("PC CMR" = "orange", 
                               "PC HDS" = "purple", 
                               "AV Bnet" = "blue")) +  # Custom colors
  scale_y_continuous(limits = c(0, 1000),
                     breaks = seq(0, 1000, by = 100),
                     labels = scales::comma) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),  
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
        panel.grid = element_blank(),
        legend.position = "none")  # Removes legend



# Total abundance
abund_summary
  


# Export density dataframe
saveRDS(dens_df, "./Data/Fitted_Models/PC_CMR_dens_df.rds")
saveRDS(dens_summary, "./Data/Fitted_Models/PC_CMR_dens_summary.rds")
saveRDS(abund_summary, "./Data/Fitted_Models/PC_CMR_abund_summary.rds")

# Save Environment
save.image(file = "./Data/Model_Environments/CMR_bm_JAGs.RData")

# End Script