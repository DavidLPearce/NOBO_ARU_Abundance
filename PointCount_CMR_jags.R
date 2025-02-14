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
workers <- 2 # Ncores  * 0.6 # For low background use 80%, for medium use 50% of Ncores

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

# Area surveyed 
area <- pi * (200^2) / 4046.86  # in acres


# Bundle data for JAGs
data <- list(J = J,
             M = M,
             nsites = nsites,
             y = y,
             X.abund = X.abund,
             X.det = X.det,
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
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
)}

# MCMC 
n.iter = 300000
n.burnin = 20000
n.chains = 3 
n.thin = 10


# -------------------------------------------------------
# Model 33: herb_Pdens  
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

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }

  # Abundance model 
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1*X.abund[s,17] 
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
save.image(file = "./CMR_bm_JAGs.RData")

# Trace plots
mcmcplot(fm.33$samples)

# Rhat
fm.33$Rhat

# Model summary
print(fm.33, digits = 2)

# # Bayes p-value
# cat("Bayesian p-value =", fm.33$summary["p_Bayes",1], "\n") 


# WAIC
# log_lik_array <- fm.33$sims.list$log_lik  # Extract log-likelihood samples
# log_lik_matrix <- apply(log_lik_array, c(1, 2), sum)  # Summing across J
# waic_result <- loo::waic(log_lik_matrix)
# print(waic_result)

# # Leave One Out
# loo_result <- loo::loo(log_lik_matrix)
# print(loo_result)

# -------------------------------------------------------
#
#   Covariate Effects 
#
# -------------------------------------------------------

# Combine chains
combined_chains <- as.mcmc(do.call(rbind, fm.33$samples))
# 
# # Extract beta estimates
# beta0_samples <- combined_chains[, "beta0"]
# beta1_samples <- combined_chains[, "beta1"]
# beta2_samples <- combined_chains[, "beta2"]
# 
# # Means
# beta0 <- mean(beta0_samples)
# beta1 <- mean(beta1_samples)
# beta2 <- mean(beta2_samples)
# 
# # # Credible Intervals
# # beta0_CI_lower <- quantile(beta0_samples, probs = 0.025)
# # beta0_CI_upper <- quantile(beta0_samples, probs = 0.975)
# # 
# # beta1_CI_lower <- quantile(beta1_samples, probs = 0.025)
# # beta1_CI_upper <- quantile(beta1_samples, probs = 0.975)
# # 
# # beta2_CI_lower <- quantile(beta2_samples, probs = 0.025)
# # beta2_CI_upper <- quantile(beta2_samples, probs = 0.975)
# 
# 
# # Create a prediction of covariate values
# herb_Pdens_pred_vals <- seq(min(site_covs[,'herb_Pdens']), max(site_covs[,'herb_Pdens']), length.out = 100)
# woody_lrgPInx_pred_vals <- seq(min(site_covs[,'woody_lrgPInx']), max(site_covs[,'woody_lrgPInx']), length.out = 100)
# 
# # Scale to have a mean of 0
# herb_Pdens_pred_vals_scaled <- (herb_Pdens_pred_vals - mean(site_covs[,'herb_Pdens'])) / sd(site_covs[,'herb_Pdens'])
# woody_lrgPInx_pred_vals_scaled <- (woody_lrgPInx_pred_vals - mean(site_covs[,'woody_lrgPInx'])) / sd(site_covs[,'woody_lrgPInx'])
# 
# # Matrices for storing predictions
# herbPdens_preds <- matrix(NA, nrow = length(beta1_samples), ncol = length(herb_Pdens_pred_vals_scaled))
# woodylrgPInx_preds <- matrix(NA, nrow = length(beta2_samples), ncol = length(woody_lrgPInx_pred_vals_scaled))
# 
# # Generate predictions  
# for (i in 1:length(beta0_samples)) {
#   herbPdens_preds[i, ] <- beta0_samples[i] + beta1_samples[i] * herb_Pdens_pred_vals_scaled
#   woodylrgPInx_preds[i, ] <- beta0_samples[i] + beta2_samples[i] * woody_lrgPInx_pred_vals_scaled
# }
# 
# # Calculate credible intervals
# herbPdens_CI_lower <- apply(herbPdens_preds, 2, quantile, probs = 0.025)
# herbPdens_CI_upper <- apply(herbPdens_preds, 2, quantile, probs = 0.975)
# 
# woodylrgPInx_CI_lower <- apply(woodylrgPInx_preds, 2, quantile, probs = 0.025)
# woodylrgPInx_CI_upper <- apply(woodylrgPInx_preds, 2, quantile, probs = 0.975)
# 
# # Calculate mean predictions
# herbPdens_mean <- apply(herbPdens_preds, 2, mean)
# woodylrgPInx_mean <- apply(woodylrgPInx_preds, 2, mean)
# 
# # Plot Herbaceous Patch
# plot(herb_Pdens_pred_vals_scaled, herbPdens_mean, type = "l", col = "lightgreen", lwd = 2, 
#      xlab = "Herbaceous Patch Density", ylab = "Predicted Effect",
#      main = "Effect of Herbaceous Patch Density")
# polygon(c(herb_Pdens_pred_vals_scaled, rev(herb_Pdens_pred_vals_scaled)),
#         c(herbPdens_CI_lower, rev(herbPdens_CI_upper)),
#         col = rgb(0.2, 0.6, 0.2, 0.2), border = NA)  # Add CI shading
# 
# # Plot Woody Largest Patch Index
# plot(woody_lrgPInx_pred_vals_scaled, woodylrgPInx_mean, type = "l", col = "forestgreen", lwd = 2, 
#      xlab = "Woody Largest Patch Index", ylab = "Predicted Effect",
#      main = "Effect of Woody Largest Patch Index")
# polygon(c(woody_lrgPInx_pred_vals_scaled, rev(woody_lrgPInx_pred_vals_scaled)),
#         c(woodylrgPInx_CI_lower, rev(woodylrgPInx_CI_upper)),
#         col = rgb(0.2, 0.6, 0.2, 0.2), border = NA)  # Add CI shading
# 


# -------------------------------------------------------
#
#   Estimating Abundance 
#
# -------------------------------------------------------

 
# Extracting Abundance
Ntot_samples <- combined_chains[ ,"N_tot"]

# Ntotal is the abundance based on 10 point counts at a radius of 200m.
# To correct for density, Ntotal needs to be divided by 10 * area surveyed
area <- pi * (200^2) / 4046.86  # Area in acres
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

# Plot density
ggplot(dens_summary, aes(x = Model)) +
  geom_point(aes(y = Mean), size = 3) +  # Add mean points
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.05) +  # Add error bars for CI
  labs(
    title = "",
    x = "Model",
    y = "Density (N/acre)"
  ) +
  scale_y_continuous(limits = c(0, 0.5),
                     breaks = seq(0, 0.5, by = 0.25),
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



# Total density
print(mean(dens_df$Density))
print(min(dens_df$Density))
print(max(dens_df$Density))



# Getting total abundance
abund_summary <- dens_summary
abund_summary[,2:4] <- abund_summary[,2:4] * 2710



# Plot abundance
ggplot(abund_summary, aes(x = Model)) +
  geom_point(aes(y = Mean), size = 3) +  # Add mean points
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.05) +  # Add error bars for CI
  labs(
    title = "",
    x = "Model",
    y = "Abundance"
  ) +
  scale_y_continuous(limits = c(500, 700),
                     breaks = seq(500, 700, by = 25),
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





# Total abundance
mean(dens_df$Density) * 2710
  


# Save Environment
save.image(file = "./CMR_bm_JAGs.RData")

# Export density dataframe
saveRDS(dens_summary, "./Data/Fitted_Models/PC_CMR_dens_summary.rds")
saveRDS(abund_summary, "./Data/Fitted_Models/PC_CMR_abund_summary.rds")


# End Script