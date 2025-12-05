# Author: David L. Pearce
# Description:
#             TBD

# This code extends code from the following sources: 
#     1. Kery, M. and Royle, J. A. (2020). Applied hierarchical modeling 
#           in ecology: Analysis of distribution, abundance, and species 
#           richness in R and BUGS: Volume 2: Dynamic and advanced models. 
#           Academic Press.


# Citation: 
#      TBD

# -------------------------------------------------------
#
#                    Load libraries
#
# -------------------------------------------------------

# Install packages (if needed)
# install.packages("tidyverse")
# install.packages("jagsUI")
# install.packages("coda")
# install.packages("MCMCvis")

# Load library
library(tidyverse)
library(jagsUI)
library(coda)
library(MCMCvis)

# Check JAGs Version
# Download here: https://sourceforge.net/projects/mcmc-jags/
rjags::jags.version() # Code written under JAGS 4.3.1.  

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
source("./Scripts/Functions/Rhat_check_function.R")

# Model name object
model_name <- "PC MCR"

# -------------------------------------------------------
#
# Variable and Object Definitions ******* Finish this
# 
# -------------------------------------------------------

# beta0 = abundance intercept 

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
colnames(site_covs)[1] <- "PointNum"


# Remove rows with NA values (no birds detected)
pc_dat <- na.omit(pc_dat)

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

# AVERAGE CALL RATE?? - use this and Lutima for AV model gamma0

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


# Extract and scale site covariates for detection
Wind <- pc_dat[, 'Wind.Beau.Code']
Wind <- as.matrix(as.numeric(as.factor(Wind)))
Wind_Lvls = nrow(unique(Wind))
merged_df <- merge(site_covs, pc_CMR, by = "PointNum") # VegDens to stacked 
VegDens <- merged_df[,c("vegDens50m")]
VegDens <- scale(VegDens)

# Inspect
head(Wind)
head(VegDens)

# Format abundance covariates
Herb_COH <- as.matrix(scale(site_covs[,'herb_COH'])) #  herb_ClmIdx
Woody_SPLIT <- as.matrix(scale(site_covs[, 'woody_SPLIT']))#  woody_AggInx

# Inspect
head(Herb_COH)
head(Woody_SPLIT)

# Bundle data for JAGs
data <- list(J = J,
             M = M,
             S = S,
             y = y,
             Wind = Wind,
             Wind_Lvls = Wind_Lvls,
             VegDens = VegDens,
             Herb_COH = Herb_COH,
             Woody_SPLIT = Woody_SPLIT,
             group = site
)

# Check structure 
str(data)


# ------------------------------------------------------------------------------
#
#                             Model Fitting
#
# ------------------------------------------------------------------------------

# -------------------
# MCMC Specifications
# -------------------
# n_iter <- 300000
# n_burnin <- 100000
# n_chains <- 3
# n_thin <- 50
# n_adapt <- 5000

# Test Settings
n_iter <- 20000
n_burnin <- 0
n_chains <- 6
n_thin <- 10
n_adapt <- 5000

# Posterior samples
n_post_samps = ((n_iter * n_chains) - n_burnin) / n_thin
print(n_post_samps)


# ----------------------
# Model Specifications
# ----------------------

# Parameters monitored
params <- c("lambda",
            "lambda_pred",
            "N",
            "N_pred",
            "N_tot",
            "alpha0", 
            "alpha1", 
            "alpha2", 
            "p",
            "beta0", 
            "beta1",
            "beta2",
            "psi",
            'tau_sre',
            "S_RE",
            "tau_det",
            "p_Bayes",
            "sum_obs",
            "sum_rep")



# Initial values
inits  <- function() {
  list (p0 = runif(1),
        alpha1 = rnorm(Wind_Lvls, 0, 1),
        alpha2 = rnorm(1, 0, 1),
        beta0 = rnorm(1, 0, 1),
        beta1 = rnorm(1, 0, 1),
        beta2 = rnorm(1, 0, 1),
        z = c( rep(1,nind), rep(0, (M-nind)))
)}


# -------------------------------------------------------
# Model Statement
# -------------------------------------------------------
cat("
model {

  # ---------------------------------
  # Abundance Priors
  # ---------------------------------
  beta0 ~ dnorm(0, 0.1) # Intercept
  beta1 ~ dnorm(0, 0.1) # Herbaceous cohesion index
  beta2 ~ dnorm(0, 0.1) # Woody splitting index
  
  # # Site random effect
  # tau_lam ~ dgamma(0.001, 0.001)
  # for (s in 1:S) {
  #    S_RE_lam[s] ~ dnorm(beta0, tau_lam)
  # }

  # ---------------------------------
  # Detection Priors
  # ---------------------------------
  p0 ~ dunif(0, 1)
  alpha0 <- log(p0/(1-p0)) 
  
  # Wind is a categorical covariate
  for (w in 1:Wind_Lvls){
  alpha1[w] ~ dnorm(0, 0.1)
  }
  
  # Vegetation Density
  alpha2 ~ dnorm(0, 0.1)
  
  # Site-level random effect for pseudoreplication
  tau_sre ~ dgamma(0.01, 0.01) 
  for(s in 1:S){  
    S_RE[s] ~ dnorm(0, tau_det)  
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
  # Abundance
  # ---------------------------------

  for(s in 1:S){
    log(lambda[s]) <- beta0 + beta1 * Herb_COH[s, 1] +  beta2 * Woody_SPLIT[s, 1]
    
    # Individual site probability
    probs[s] <- lambda[s] / sum(lambda[])
    
    # Estimated abundance at site s
    N[s] <- sum(group[] == s) 
    
  } # End S

  # ---------------------------------
  # Presence 
  # ---------------------------------

  for(i in 1:M){
  
    # Group == site membership
    group[i] ~ dcat(probs[])  
    
    # Presence: Data augmentation variables
    z[i] ~ dbern(psi)         
    
  # ---------------------------------
  # Detection
  # ---------------------------------
 
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1[Wind[group[i], 1]] + alpha2 * VegDens[group[i],1] + S_RE[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    } # End J
  } # End M
  
  # -------------------------------------------
  # Predict Abundance at Unsampled Areas 
  # -------------------------------------------
  
  for (u in 1:U) {

    log(lambda_pred[u]) <- log(rho[u]) + beta0 + beta1 * Herb_COH_pred[u] + beta2 * Woody_SPLIT_pred[u]
    N_pred[u] ~ dpois(lambda_pred[u]) 
 
  }   
  
  # -------------------------------------------
  # PPC and Bayesian P-value
  # -------------------------------------------

  for(i in 1:M){
   for(j in 1:J){
    # Generate replicated data
    y_rep[i,j] ~ dbern(pz[i,j]) 
    
    # Discrepancy for observed data
    # discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j]) 
    discrepancy_obs[i,j] <- pow(y[i,j] - p[i,j], 2)
    
    # Discrepancy for replicated data  
    # discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  
    discrepancy_rep[i,j] <- pow(y_rep[i,j] - p[i,j], 2)

    } # End J
  } # End M
    
    # Sum of discrepancies for observed data
    sum_obs <- sum(discrepancy_obs[,])  
  
    # Sum of discrepancies for replicated data
    sum_rep <- sum(discrepancy_rep[,])  
  
    # Proportion of times replicated > observed
    p_Bayes <- step(sum_rep - sum_obs)

  # ---------------------------------
  # Derived Metrics
  # ---------------------------------
 
  # Abundance
  N_tot <- sum(z[])
  
  # Abundance at unsampled sites
  N_pred_tot <- sum(N_pred[])
  
  # Total abundance
  N_tot <- sum(N[]) + sum(N_pred[])
  
  # Overall site mean abundances
  samp_site_Nmean <-  sum(N[]) / S
  pred_site_Nmean <- sum(N_pred[]) / U
  
  # Expected site mean abundance
  samp_site_lambda_mean <-  sum(lambda[]) / S
  pred_site_lambda_mean <- sum(lambda_pred[]) / U
  
  # Variance
  samp_site_lambda_var <- sum( (lambda[] - samp_site_lambda_mean)^2 ) / (S - 1)
  pred_site_lambda_var <- sum( (lambda_pred[] - pred_site_lambda_mean)^2 ) / (U - 1)

  
}

", fill=TRUE, file = "./jags_models/Model_PC_MCR.txt")
# ------------End Model-------------

# Fit Model
fm1 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/Model_PC_MCR.txt",
             n.iter = n_iter,
             n.burnin = n_burnin,
             n.chains = n_chains, 
             n.thin = n_thin,
             n.adapt = n_adapt,
             parallel = TRUE,
             n.cores = workers,
             verbose = TRUE,
             DIC = FALSE
)

# -------------------------------------------------------
# Check Convergence
# -------------------------------------------------------


# Trace plots
MCMCvis::MCMCtrace(fm1, 
                   params = c("N_tot",
                              "N_samp_tot",
                              "N_pred_tot",
                              "beta0", 
                              "beta1",
                              "beta2",
                              "alpha0", 
                              "alpha1", 
                              "alpha2", 
                              'psi',
                              'tau_sre',
                              "S_RE",
                              "samp_site_mu_var",
                              "pred_site_mu_var",
                              "samp_site_mu_mean",
                              "pred_site_mu_mean",
                              "samp_site_Nmean",
                              "pred_site_Nmean",
                              "lambda",
                              "lambda_pred",
                              "N",
                              "N_pred",

                   ),
                   pdf = T,
                   filename = "TracePlots_PC_MCR.pdf",
                   wd = "./Figures"
)

# Rhat
check_rhat(fm1$Rhat, threshold = 1.1)


# Save model
saveRDS(fm1, "./Data/Model_Data/Fit_Model_PC_MCR.rds")

 
# -------------------------------------------------------
# Posterior Predictive Checks
# -------------------------------------------------------


# ----------------------
# Extract Fit
# ----------------------

# Abundance
fit_y_data <- data.frame(
  Observed = fm1$sims.list$sum_obs, # Observed values
  Predicted = fm1$sims.list$sum_rep,  # Predicted values
  type = rep(c("Observed", "Predicted"), each = length(fm1$sims.list$sum_obs))
)


# ----------------------
# Density Plot
# ----------------------

# Bayes P-value
# P-value = 0.5 means good fit, = 1 or 0 is a poor fit
mn_bpy <- round(mean(fm1$summary["p_Bayes",1]), 2) 

# y
y_PPC_Dens <- ggplot(fit_y_data) +
  geom_density(aes(x = Observed, fill = "Observed"), alpha = 0.5) +   
  geom_density(aes(x = Predicted, fill = "Predicted"), alpha = 0.5) +  
  scale_fill_manual(values = c("Observed" = "blue", "Predicted" = "red")) +  
  labs(title = "", 
       x = "Fit Values", 
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") + 
  annotate("text", x = 100, y = 0.03, label = paste0("Bayes p-value = ", mn_bpy), hjust = 0) 

# View
print(y_PPC_Dens)

# Export                
ggsave(plot = y_PPC_Dens, "./Figures/PPC_BP_PC_MCR.jpeg", width = 8, height = 5, dpi = 300)
dev.off()

# -------------------------------------------------------
#
#                 Posterior Estimates  
#
# -------------------------------------------------------

# Combine chains for posterior inference
combined_chains <- as.mcmc(do.call(rbind, fm1$samples))

# -------------------------------------------------------
# Beta Estimates
# -------------------------------------------------------

# Extract beta estimates
beta0_samples <- combined_chains[, "beta0"]
beta1_samples <- combined_chains[, "beta1"]
beta2_samples <- combined_chains[, "beta2"]

# Compute 95% CI for each beta
beta_df <- data.frame(
  value = c(beta0_samples, beta1_samples, beta2_samples ),
  parameter = rep(c("beta0", "beta1", "beta2" ), each = length(beta0_samples))
) %>%
  group_by(parameter) %>%
  filter(value >= quantile(value, 0.025) & value <= quantile(value, 0.975))  # Keep only values within 95% CI

# Add model
beta_df$Model <- model_name

# Create Summary
beta_summary <- beta_df %>%
  group_by(parameter, Model) %>%
  summarise(
    mean = mean(value),
    LCI = min(value),
    UCI = max(value),
    .groups = "drop"  
  )
# view
print(beta_summary)

# Export beta dataframe and summary
saveRDS(beta_df, "./Data/Model_Data/Beta_df_PC_MCR.rds")
write.csv(beta_summary, "./Data/Model_Data/Beta_Summary_PC_MCR.csv")

# -------------------------------------------------------
# Alpha Estimates
# -------------------------------------------------------

# Extract alpha estimates
alpha0_samples <- combined_chains[, "alpha0"]
alpha1_samples <- combined_chains[, grepl("^alpha1\\[", colnames(combined_chains))]
alpha1_samples_1 <- alpha1_samples[,1]
alpha1_samples_2 <- alpha1_samples[,2]
alpha1_samples_3 <- alpha1_samples[,3] 
alpha1_samples_4 <- alpha1_samples[,4]
alpha1_samples_5 <- alpha1_samples[,5]
alpha2_samples <- combined_chains[, "alpha2"]


# Extract site random effect
sRE_cols  <- grep("^S_RE_det\\[", colnames(combined_chains), value = TRUE)
sRE_samples <- combined_chains[, sRE_cols]
sRE_samples <- rowMeans(sRE_samples) # Row means



# Compute 95% CI  
alpha_df <- data.frame(
  value = c(alpha0_samples, 
            alpha1_samples_1, 
            alpha1_samples_2,
            alpha1_samples_3,
            alpha1_samples_4,
            alpha1_samples_5,
            alpha2_samples,
            sRE_samples),  
  parameter = rep(c("alpha0", 
                    "alpha1_1", 
                    'alpha1_2',
                    'alpha1_3',
                    'alpha1_4',
                    'alpha1_5',
                    "alpha2",
                    "sRE_det"), each = length(alpha0_samples))
) %>%
  group_by(parameter) %>%
  filter(value >= quantile(value, 0.025) & value <= quantile(value, 0.975))  

# Add model
alpha_df$Model <- model_name

# Create summary
alpha_summary <- alpha_df %>%
  group_by(parameter, Model) %>%
  summarise(
    mean = mean(value),
    LCI = min(value),
    UCI = max(value),
    .groups = "drop"  
  )

# view
print(alpha_summary)

# Export alpha summary
saveRDS(alpha_summary, "./Data/Model_Data/Alpha_Summary_PC_MCR.rds")


# -------------------------------------------------------
# Detection probability 
# -------------------------------------------------------

# Extract detection posterior samples
p_samples <- combined_chains[, grepl("^p\\[", colnames(combined_chains))]

# Combine all samples into a single long vector
all_p_samples <- as.vector(p_samples)

# Create summary
p_summary <- data.frame(
  mean = mean(all_p_samples),
  LCI = quantile(all_p_samples, 0.025),
  UCI = quantile(all_p_samples, 0.975)
)

# Add model
p_summary$Model <- model_name

# Add parameter name
p_summary$Parameter <- "detection"

# view
print(p_summary)

# -------------------------------------------------------
#   Estimating Abundance 
# -------------------------------------------------------

# Extract abundance at surveyed sites
N_survyed_cols <- grep("^N\\[", colnames(combined_chains), value = TRUE)
N_survyed_samples <- combined_chains[, N_survyed_cols, drop = FALSE] 
colnames(N_survyed_samples) <- paste0("S", 1:S)

# Extract abundance at unsurveyed sites
N_unsurvyed_cols <- grep("^N_pred\\[", colnames(combined_chains), value = TRUE)
N_unsurvyed_samples <- combined_chains[, N_unsurvyed_cols, drop = FALSE] 
colnames(N_unsurvyed_samples) <- paste0("U", 1:U)
min(N_unsurvyed_samples)

# Combine into one dataframe
N_samples_df <- cbind(N_survyed_samples, N_unsurvyed_samples)
str(N_samples_df)

# Site Summary
site_summary <- data.frame(
  Site = colnames(N_samples_df),
  Mean = apply(N_samples_df, 2, mean),
  LCL  = apply(N_samples_df, 2, quantile, 0.025),
  UCL  = apply(N_samples_df, 2, quantile, 0.975)
)

print(site_summary)

# Total Summary 
total_N_samples <- combined_chains[,'N_tot']
total_summary <- data.frame(
  Mean = mean(total_N_samples),
  LCI  = quantile(total_N_samples, 0.025),
  UCI  = quantile(total_N_samples, 0.975),
  Model = model_name,
  Parameter = "Abundance"
)

print(total_summary)

# Density hectare
dens_summary_ha <- total_summary[,1:3] / 1098
dens_summary_ha$Model <- model_name
dens_summary_ha$Parameter <- "Density (N/ha)"
dens_summary_ha

# Density Acre
dens_summary_ac <- total_summary[,1:3] / 2710
dens_summary_ac$Model <- model_name
dens_summary_ac$Parameter <- "Density (N/ac)"
dens_summary_ac

# Combine with detection and vocal rate
param_summary <- rbind(param_summary, total_summary, dens_summary_ha, dens_summary_ac)
print(param_summary)

# Trim total_abundance_samples to 95% CI
total_N_samples_95CI <- total_N_samples[
  total_N_samples >= quantile(total_N_samples, 0.025) &
    total_N_samples <= quantile(total_N_samples, 0.975)
]

  
# Export alpha summary
saveRDS(abund_95df, "./Data/Model_Data/Abund_df_PC_MCR.rds")
saveRDS(param_summary, "./Data/Model_Data/Param_Summary_PC_MCR.rds")

# ------------ End Script -----------------