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
pc_dat$DOY_sin <- sin(2 * pi * pc_dat$DOY / 365)


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
woodyParea = as.vector(site_covs[,c("woody_Parea")]) 

# Extract detection covariates
obsvr <- as.matrix(as.numeric(as.factor(pc_dat[,"Observer"])))  # Convert names to numeric indices
temp <- as.matrix(as.numeric(pc_dat[,"Temp.deg.F"]))
wind <- as.matrix(as.numeric(as.factor(pc_dat[,"Wind.Beau.Code"]))) 
sky <- as.matrix(as.numeric(as.factor(pc_dat[,"Sky.Beau.Code"])))
doy <- as.matrix(pc_dat[,"DOY_sin"])

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
             woodyParea = woodyParea,
             obsvr = obsvr,
             temp = temp,
             wind = wind,
             sky = sky,
             doy = doy,
             group = site)
 
# Take a look
str(data)

 


# -------------------------------------------------------
#
#                   Model Fitting
#
# -------------------------------------------------------

## Distributions
# library(ggfortify)
# ggdistribution(dunif, seq(0, 1, 0.001), min = 0, max = 1) # p0
# ggdistribution(dnorm, seq(0, 0.01, 0.0001)) # alpha's and beta's 
# ggdistribution(dgamma, seq(0, 5, 0.01), shape = 0.1, rate = 0.1) # tau


# -------------------------------------------------------
# Null Model
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1 * herbPdens[s] + beta2 * woodyParea[s]
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
", fill=TRUE, file="CMR_mod0.txt")
# ------------End Model-------------


# Parameters monitored
params.0 <- c("lambda",
              "p0", 
              "alpha0", 
              "beta0", 
              "beta1",
              "beta2",
              "psi",
              "S.raneff",
              "tau",
              "sigma",
              "p_Bayes")


# Initial values
inits.0 <- function() {
  list (p0 = runif(1),
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
fm.0 <- jags(data = data, 
             parameters.to.save = params.0,
             inits = inits.0, 
             model.file = "CMR_mod0.txt",
             n.iter = 250000,
             n.burnin = 15000,
             n.chains = 3, 
             n.thin = 10,
             parallel = TRUE,
             n.cores = 8,
             DIC = TRUE)  



# Rhat
fm.0$Rhat

# Model summary
print(fm.0, digits = 3)

# Trace plots
mcmcplot(fm.0$samples)

# Bayesian P value
cat("Bayesian p-value =", fm.0$summary["p_Bayes",1], "\n")





# -------------------------------------------------------
# Model 1: Observer
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  #alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1 * herbPdens[s] + beta2 * woodyParea[s]
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
", fill=TRUE, file="CMR_mod1.txt")
# ------------End Model-------------


# Parameters monitored
params.1 <- c("lambda",
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
inits.1 <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
fm.1 <- jags(data = data, 
             parameters.to.save = params.1,
             inits = inits.1, 
             model.file = "CMR_mod1.txt",
             n.iter = 250000,
             n.burnin = 15000,
             n.chains = 3, 
             n.thin = 10,
             parallel = TRUE,
             n.cores = 8,
             DIC = TRUE)  



# Rhat
fm.1$Rhat

# Model summary
print(fm.1, digits = 3)

# Trace plots
mcmcplot(fm.1$samples)

# Bayesian P value
cat("Bayesian p-value =", fm.1$summary["p_Bayes",1], "\n")



# -------------------------------------------------------
# Model 2: Temperature
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  #alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1 * herbPdens[s] + beta2 * woodyParea[s]
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
", fill=TRUE, file="CMR_mod2.txt")
# ------------End Model-------------


# Parameters monitored
params.2 <- c("lambda",
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
inits.2 <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
fm.2 <- jags(data = data, 
             parameters.to.save = params.2,
             inits = inits.2, 
             model.file = "CMR_mod2.txt",
             n.iter = 250000,
             n.burnin = 15000,
             n.chains = 3, 
             n.thin = 10,
             parallel = TRUE,
             n.cores = 8,
             DIC = TRUE)  



# Rhat
fm.2$Rhat

# Model summary
print(fm.2, digits = 3)

# Trace plots
mcmcplot(fm.2$samples)

# Bayesian P value
cat("Bayesian p-value =", fm.2$summary["p_Bayes",1], "\n")





# -------------------------------------------------------
# Model 3: Wind
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  #alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1 * herbPdens[s] + beta2 * woodyParea[s]
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
", fill=TRUE, file="CMR_mod3.txt")
# ------------End Model-------------


# Parameters monitored
params.3 <- c("lambda",
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
inits.3 <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
fm.3 <- jags(data = data, 
             parameters.to.save = params.3,
             inits = inits.3, 
             model.file = "CMR_mod3.txt",
             n.iter = 250000,
             n.burnin = 15000,
             n.chains = 3, 
             n.thin = 10,
             parallel = TRUE,
             n.cores = 8,
             DIC = TRUE)  



# Rhat
fm.3$Rhat

# Model summary
print(fm.3, digits = 3)

# Trace plots
mcmcplot(fm.3$samples)

# Bayesian P value
cat("Bayesian p-value =", fm.3$summary["p_Bayes",1], "\n")



# -------------------------------------------------------
# Model 3: Sky
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  #alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1 * herbPdens[s] + beta2 * woodyParea[s]
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
", fill=TRUE, file="CMR_mod4.txt")
# ------------End Model-------------


# Parameters monitored
params.4 <- c("lambda",
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
inits.4 <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
fm.4 <- jags(data = data, 
             parameters.to.save = params.4,
             inits = inits.4, 
             model.file = "CMR_mod4.txt",
             n.iter = 250000,
             n.burnin = 15000,
             n.chains = 3, 
             n.thin = 10,
             parallel = TRUE,
             n.cores = 8,
             DIC = TRUE)  



# Rhat
fm.4$Rhat

# Model summary
print(fm.4, digits = 4)

# Trace plots
mcmcplot(fm.4$samples)

# Bayesian P value
cat("Bayesian p-value =", fm.4$summary["p_Bayes",1], "\n")


# -------------------------------------------------------
# Model 4: sin Day of Year
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  #alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1 * herbPdens[s] + beta2 * woodyParea[s]
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
", fill=TRUE, file="CMR_mod5.txt")
# ------------End Model-------------


# Parameters monitored
params.5 <- c("lambda",
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
inits.5 <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
fm.5 <- jags(data = data, 
             parameters.to.save = params.5,
             inits = inits.5, 
             model.file = "CMR_mod5.txt",
             n.iter = 250000,
             n.burnin = 15000,
             n.chains = 3, 
             n.thin = 10,
             parallel = TRUE,
             n.cores = 8,
             DIC = TRUE)  



# Rhat
fm.5$Rhat

# Model summary
print(fm.5, digits = 4)

# Trace plots
mcmcplot(fm.5$samples)

# Bayesian P value
cat("Bayesian p-value =", fm.5$summary["p_Bayes",1], "\n")





# -------------------------------------------------------
#
#   Estimating Abundance 
#
# -------------------------------------------------------


# Combine chains
combined_chains <- as.mcmc(do.call(rbind, fm.rcov$samples))

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
    y = "Density (N/hectare") +
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

















