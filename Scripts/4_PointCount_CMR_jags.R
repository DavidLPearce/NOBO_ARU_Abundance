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
        beta0 = 1,
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
# Model 1: Woody N patches
# -------------------------------------------------------
cat("
model {

  # ---------------------------------
  # Abundance Priors
  # ---------------------------------
  beta0 ~ dnorm(0, 10)
  beta1 ~ dnorm(0, 10)
  beta2 ~ dnorm(0, 5)
  
  # ---------------------------------
  # Detection Priors
  # ---------------------------------
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0)) 
  alpha1 ~ dnorm(0, 0.01)
  
  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  

  # Site-level random effect for pseudoreplication
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }
  
  # ---------------------------------
  # Individual Encounter/Presence Derived 
  # ---------------------------------
  psi <- sum(lambda[])/M
  
  # -------------------------------------------
  #
  # Likelihood and Process Model 
  #
  # -------------------------------------------

  # ---------------------------------
  # Abundance Model
  # ---------------------------------

  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1 * X.abund[s, 2] 
    
    # Individual site probability
    probs[s] <- lambda[s] / sum(lambda[])
    
    # Estimated abundance at site s
    N[s] <- sum(group[] == s)  
  }

  # ---------------------------------
  # Presence Model 
  # ---------------------------------

  for(i in 1:M){
  
    # Group == site membership
    group[i] ~ dcat(probs[])  
    
    # Presence: Data augmentation variables
    z[i] ~ dbern(psi)         
    
  # ---------------------------------
  # Detection Model
  # ---------------------------------

    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
  # ---------------------------------
  # Posterior Predictive checks
  # ---------------------------------

    # Generate replicated data
    y_rep[i,j] ~ dbern(pz[i,j]) 
    
    # Discrepancy for observed data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j]) 
    
    # Discrepancy for replicated data  
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  
  
  # ---------------------------------
  # Log-Likelihood
  # ---------------------------------
 
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    
    
    } # End J
  } # End M
  
  # ---------------------------------
  # Total Abundance
  # ---------------------------------
 
  N_tot <- sum(z[])
  
  # ---------------------------------
  # Bayesian p-value
  # ---------------------------------
  
  # Sum of discrepancies for observed data
  sum_obs <- sum(discrepancy_obs[,])  
  
  # Sum of discrepancies for replicated data
  sum_rep <- sum(discrepancy_rep[,])  
  
  # Proportion of times replicated > observed
  p_Bayes <- step(sum_rep - sum_obs)  
  
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

# Check Convergence
check_rhat(fm.1$Rhat, threshold = 1.1) # Rhat values
# mcmcplot(fm.1$samples) # Trace plots

# -------------------------------------------------------
#
#   Beta Estimates and Covariate Effects 
#
# -------------------------------------------------------

# Combine chains
combined_chains <- as.mcmc(do.call(rbind, fm.1$samples))

# -------------------------------------------------------
# Beta Estimates
# -------------------------------------------------------

# Extract beta estimates
beta0_samples <- combined_chains[, "beta0"]
beta1_samples <- combined_chains[, "beta1"]
 

# Means
beta0 <- mean(beta0_samples)
beta1 <- mean(beta1_samples)
 

# Compute 95% CI for each beta
beta_df <- data.frame(
  value = c(beta0_samples, beta1_samples),
  parameter = rep(c("beta0", "beta1"), each = length(beta0_samples))
) %>%
  group_by(parameter) %>%
  filter(value >= quantile(value, 0.025) & value <= quantile(value, 0.975))  # Keep only values within 95% CI

# Add model
beta_df$Model <- "PC CMR"

# Export beta dataframe
saveRDS(beta_df, "./Data/Fitted_Models/PC_CMR_beta_df.rds")

# -------------------------------------------------------
# Covariate Effects
# -------------------------------------------------------

# Set covariate name 
Cov_name <- "herb_mnParea"

# Create a prediction of covariate values
cov_pred_vals <- seq(min(site_covs[, Cov_name]), max(site_covs[, Cov_name]), length.out = 100)

# Matrices for storing predictions
cov_preds <- matrix(NA, nrow = length(beta0_samples), ncol = length(cov_pred_vals))

# Generate predictions
for (i in 1:length(beta0_samples)) {
  cov_preds[i, ] <- beta0_samples[i] + beta1_samples[i] * cov_pred_vals
}

# Calculate mean predictions
cov_preds_mean <- apply(cov_preds, 2, mean)# mean
cov_preds_LCI <- apply(cov_preds, 2, quantile, probs = 0.025) # LCI
cov_preds_HCI <- apply(cov_preds, 2, quantile, probs = 0.975) # HCI

# Combine into a single data frame
cov_pred_df <- data.frame(
  cov_pred_vals = cov_pred_vals,
  cov_preds_mean = cov_preds_mean,
  cov_preds_LCI = cov_preds_LCI,
  cov_preds_HCI = cov_preds_HCI
)



# Plot Cov effect
ggplot(cov_pred_df, aes(x = cov_pred_vals, y = cov_preds_mean)) +
  geom_line(color = "black", linewidth = 1.5) +   
  geom_ribbon(aes(ymin = cov_preds_LCI, 
                  ymax = cov_preds_HCI), 
              fill = "forestgreen", alpha = 0.3) +
  labs(x = "Covariate Value", 
       y = "Effect Estimate", 
       title = paste0("Predicted Effect of ", Cov_name)) +
  theme_minimal() +
  theme(panel.grid = element_blank())


# -------------------------------------------------------
#   Estimating Abundance 
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
#save.image(file = "./Data/Model_Environments/CMR_bm_JAGs.RData")

# End Script