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
library(plotly)
library(jagsUI)
library(coda)
library(mcmcplots)
library(loo)

# Check JAGs Version
# Download here: https://sourceforge.net/projects/mcmc-jags/
rjags::jags.version() # Code written under JAGS 4.3.1. version as of (5 Feb 2025) 

# Set seed, scientific notation options, and working directory
set.seed(123)
options(scipen = 9999)
setwd(".")

# Setting up cores
Ncores <- parallel::detectCores()
print(Ncores) # Number of available cores
workers <- Ncores  * 0.5 # For low background use 80%, for medium use 50% of Ncores
print(workers)

# Source custom function for checking Rhat values > 1.1
source("./Scripts/Rhat_check_function.R")

# Model name object
model_name <- "PC CMR"

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

# Rename SiteID to PointNum for matching
colnames(site_covs)[2] <- "PointNum"


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

# Number of sites
S <- 10

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


# Extract and scale site covariates for X.det
X.det <- pc_dat[,c(2,7:9,18)]
X.det$Observer <- as.numeric(as.factor(X.det$Observer))
X.det$Temp.deg.F <- as.numeric(X.det$Temp.deg.F)
X.det$Wind.Beau.Code <- as.numeric(as.factor(X.det$Wind.Beau.Code)) 
X.det$Sky.Beau.Code <- as.numeric(as.factor(X.det$Sky.Beau.Code))
merged_df <- merge(site_covs, pc_CMR, by = "PointNum") # VegDens to stacked 
X.det$vegDens = as.matrix(merged_df[,c("vegDens50m")])
X.det <- as.matrix(X.det)
print(X.det)

# Format X.abund
X.abund <- as.matrix(site_covs[,-c(1:2)]) # Remove X and site id
print(X.abund)


# Area surveyed 
area <- pi * (200^2) / 4046.86  # in acres


# Bundle data for JAGs
data <- list(J = J,
             M = M,
             S = S,
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

# -------------------
# MCMC Specifications
# -------------------

n.iter = 200000
n.burnin = 60000
n.chains = 3 
n.thin = 10

# Posterior samples
est_post_samps = ((n.iter* n.chains) - n.burnin) / n.thin
print(est_post_samps)


# Parameters monitored
params <- c("lambda",
            "N",
            "N_tot",
            "pz_mean", 
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


# -------------------------------------------------------
# Model 1 Statement
# -------------------------------------------------------
cat("
model {

  # ---------------------------------
  # Abundance Priors
  # ---------------------------------
  beta0 ~ dnorm(0, 10)
  beta1 ~ dnorm(0, 10)
  beta2 ~ dnorm(0, 10)
  beta3 ~ dnorm(0, 10)

  # ---------------------------------
  # Detection Priors
  # ---------------------------------
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0)) 
  alpha1 ~ dnorm(0, 10)
  
  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  

  # Site-level random effect for pseudoreplication
  for(s in 1:S){  
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
  # Abundance Submodel
  # ---------------------------------

  for(s in 1:S){
    log(lambda[s]) <- beta0 + beta1 * X.abund[s, 7] +  beta2 * X.abund[s, 12]
    
    # Individual site probability
    probs[s] <- lambda[s] / sum(lambda[])
    
    # Estimated abundance at site s
    N[s] <- sum(group[] == s)  
  }

  # ---------------------------------
  # Presence Submodel 
  # ---------------------------------

  for(i in 1:M){
  
    # Group == site membership
    group[i] ~ dcat(probs[])  
    
    # Presence: Data augmentation variables
    z[i] ~ dbern(psi)         
    
  # ---------------------------------
  # Detection Submodel
  # ---------------------------------

    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],6] + S.raneff[group[i]]
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
  # Derived Metrics
  # ---------------------------------
 
  # Abundance
  N_tot <- sum(z[])
  
  # Detection
  pz_mean <- mean(pz[,])
  
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
mcmcplot(fm.1$samples) # Trace plots

# Check Fit
cat("Bayesian p-value =", fm.1$summary["p_Bayes",1], "\n")# P-value = 0.5 means good fit, = 1 or 0 is a poor fit

# Model summary
#summary(fm.1$samples)

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
beta2_samples <- combined_chains[, "beta2"]
#beta3_samples <- combined_chains[, "beta3"]

# Means
beta0 <- mean(beta0_samples)
beta1 <- mean(beta1_samples)
beta2 <- mean(beta2_samples)
# beta3 <- mean(beta3_samples)

# Compute 95% CI for each beta
beta_df <- data.frame(
  value = c(beta0_samples, beta1_samples, beta2_samples ), # , beta3_samples
  parameter = rep(c("beta0", "beta1", "beta2" ), each = length(beta0_samples)) # , "beta3"
) %>%
  group_by(parameter) %>%
  filter(value >= quantile(value, 0.025) & value <= quantile(value, 0.975))  # Keep only values within 95% CI

# Add model
beta_df$Model <- model_name

# Plot
ggplot(beta_df, aes(x = parameter, y = value, fill = parameter)) +
  geom_violin(alpha = 0.5, trim = TRUE) +   
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) + 
  labs(title = "Violin Plots for Beta Estimates", x
       = "Parameter", 
       y = "Estimate") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") 




# Export beta dataframe
saveRDS(beta_df, "./Data/Model_Data/PC_CMR_beta_df.rds")

# -------------------------------------------------------
# Covariate Effects
# -------------------------------------------------------


# Set covariate name 
Cov1_name <- "herb_ClmIdx"
Cov2_name <- "woody_AggInx"

# Create a prediction of covariate values
cov1_pred_vals <- seq(min(X.abund[, Cov1_name]), max(X.abund[, Cov1_name]), length.out = 1000)
cov2_pred_vals <- seq(min(X.abund[, Cov2_name]), max(X.abund[, Cov2_name]), length.out = 1000)

# Mean scaling covariates
cov1_scaled <- (X.abund[, Cov1_name] - mean(X.abund[, Cov1_name])) / (max(X.abund[, Cov1_name]) - min(X.abund[, Cov1_name]))
cov2_scaled <- (X.abund[, Cov2_name] - mean(X.abund[, Cov2_name])) / (max(X.abund[, Cov2_name]) - min(X.abund[, Cov2_name]))

# Matrices for storing predictions
cov1_preds <- matrix(NA, nrow = length(beta0_samples), ncol = length(cov1_scaled))
cov2_preds <- matrix(NA, nrow = length(beta0_samples), ncol = length(cov2_scaled))

# # Create a meshgrid of covariate values for interaction predictions
# interaction_grid <- expand.grid(cov1_scaled = cov1_scaled, cov2_scaled = cov2_scaled)
# interaction_grid$interaction_term <- interaction_grid$cov1_scaled * interaction_grid$cov2_scaled
# 
# # Initialize matrix to store predictions
# interaction_preds <- matrix(NA, nrow = length(beta0_samples), ncol = nrow(interaction_grid))
# 

# Generate predictions
for (i in 1:length(beta0_samples)) {
  cov1_preds[i, ] <- beta0_samples[i] + beta1_samples[i] * cov1_scaled # Linear
  cov2_preds[i, ] <- beta0_samples[i] + beta2_samples[i] * cov2_scaled
  
  # interaction_preds[i, ] <- beta0_samples[i] +                       # Interaction
  #   beta1_samples[i] * interaction_grid$cov1_scaled +
  #   beta2_samples[i] * interaction_grid$cov2_scaled +
  #   beta3_samples[i] * interaction_grid$interaction_term
  
}


# Calculate mean predictions
cov1_preds_mean <- apply(cov1_preds, 2, mean)# mean
cov1_preds_LCI <- apply(cov1_preds, 2, quantile, probs = 0.025) # LCI
cov1_preds_HCI <- apply(cov1_preds, 2, quantile, probs = 0.975) # HCI

cov2_preds_mean <- apply(cov2_preds, 2, mean) 
cov2_preds_LCI <- apply(cov2_preds, 2, quantile, probs = 0.025) 
cov2_preds_HCI <- apply(cov2_preds, 2, quantile, probs = 0.975) 

# interaction_preds_mean <- apply(interaction_preds, 2, mean)
# interaction_preds_LCI <- apply(interaction_preds, 2, quantile, probs = 0.025)
# interaction_preds_HCI <- apply(interaction_preds, 2, quantile, probs = 0.975)

# Combine into a single data frame
cov1_pred_df <- data.frame(
  cov1_scaled = cov1_scaled,
  cov1_preds_mean = cov1_preds_mean,
  cov1_preds_LCI = cov1_preds_LCI,
  cov1_preds_HCI = cov1_preds_HCI)

cov2_pred_df <- data.frame(
  cov2_scaled = cov2_scaled,
  cov2_preds_mean = cov2_preds_mean,
  cov2_preds_LCI = cov2_preds_LCI,
  cov2_preds_HCI = cov2_preds_HCI)

# interaction_pred_df <- data.frame(
#   interaction_term = interaction_grid$interaction_term,
#   cov1_scaled = interaction_grid$cov1_scaled,
#   cov2_scaled = interaction_grid$cov2_scaled,
#   interaction_preds_mean = interaction_preds_mean,
#   interaction_preds_LCI = interaction_preds_LCI,
#   interaction_preds_HCI = interaction_preds_HCI)


# Plot effect

# Cov 1
ggplot(cov1_pred_df, aes(x = cov1_scaled, y = cov1_preds_mean)) +
  geom_line(color = "black", linewidth = 1.5) +   
  geom_ribbon(aes(ymin = cov1_preds_LCI, 
                  ymax = cov1_preds_HCI), 
              fill = "forestgreen", alpha = 0.3) +
  labs(x = "Covariate Value", 
       y = "Effect Estimate", 
       title = paste0(model_name, " | Predicted Effect of ", Cov1_name)) +
  theme_minimal() +
  theme(panel.grid = element_blank())

# Cov 2
ggplot(cov2_pred_df, aes(x = cov2_scaled, y = cov2_preds_mean)) +
  geom_line(color = "black", linewidth = 1.5) +   
  geom_ribbon(aes(ymin = cov2_preds_LCI, 
                  ymax = cov2_preds_HCI), 
              fill = "forestgreen", alpha = 0.3) +
  labs(x = "Covariate Value", 
       y = "Effect Estimate", 
       title = paste0(model_name, " | Predicted Effect of ", Cov2_name)) +
  theme_minimal() +
  theme(panel.grid = element_blank())


# # Interactive effects
# plot_ly(interaction_pred_df, 
#         x = ~cov1_scaled, 
#         y = ~cov2_scaled, 
#         z = ~interaction_preds_mean, 
#         type = "contour",
#         colorscale = "Viridis",
#         ncontours = 80,
#         contours = list(showlabels = FALSE), 
#         colorbar = list(title = "Interaction Effect", titlefont = list(size = 14)) 
#         )%>%
#         layout(
#           title = list(
#             text = paste0(model_name, " | Interaction Effect "),
#             font = list(size = 18),  
#             x = 0.5,   
#             xanchor = "center"   
#           ),
#           xaxis = list(
#             title = list(text = "Herbaceous Clumpy Index (scaled)", 
#                          font = list(size = 14))),
#           yaxis = list(
#             title = list(text = "Woody Number of Patches (scaled)", 
#                          font = list(size = 14))))



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
dens_df <- data.frame(Model = rep(model_name, length(dens_samples)), Density = dens_samples)
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

# Plot Abundance - Violin model_name
abund_df <- dens_df
abund_df$Density <- abund_df$Density * 2710

ggplot(abund_df, aes(x = Model, y = Density, fill = Model)) + 
  geom_violin(trim = FALSE, alpha = 0.6, adjust = 5) +
  labs(x = "Model", y = "Density (N/acre)") +
  scale_fill_manual(values = c("PC CMR" = "orange")) +  
  scale_y_continuous(limits = c(0, 1000),
                     breaks = seq(0, 1000, by = 100),
                     labels = scales::comma) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),  
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
        panel.grid = element_blank(),
        legend.position = "none")   



# Total abundance
abund_summary
  


# Export density dataframe
saveRDS(dens_df, "./Data/Model_Data/PC_CMR_dens_df.rds")
saveRDS(dens_summary, "./Data/Model_Data/PC_CMR_dens_summary.rds")
saveRDS(abund_summary, "./Data/Model_Data/PC_CMR_abund_summary.rds")

# Save Environment
#save.image(file = "./Data/Model_Environments/CMR_bm_JAGs.RData")

# End Script