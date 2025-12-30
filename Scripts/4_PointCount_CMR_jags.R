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
# install.packages("gridExtra")
# install.packages("jagsUI")
# install.packages("coda")
# install.packages("MCMCvis")

# Load library
library(tidyverse)
library(gridExtra)
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

# Unsampled area covariates
predict_covs <- read.csv("./Data/Point_Count_Data/PC_PredictsiteCovs.csv")

# -------------------------------------------------------
#
#                   Data Wrangling
#
# -------------------------------------------------------

# ----------------------
# Observation Matrix
# ----------------------

# Subset to surveys 2, 3, and 4  
pc_CMR <- pc_dat %>%
  filter(Survey %in% c(2:4))

# Check the first few rows
head(pc_CMR, 10)

# Function to get first detection interval
get_first_detection <- function(int1, int2, int3, int4) {
  if(all(is.na(c(int1, int2, int3, int4)))) {
    return(NA)  
  }
  if(int1 == 1) return(1)
  if(int2 == 1) return(2)
  if(int3 == 1) return(3)
  if(int4 == 1) return(4)
  return(5)  
}
# ----

# Add first detection variable
pc_CMR <- pc_CMR %>%
  rowwise() %>%
  mutate(first_detection = get_first_detection(int1, int2, int3, int4)) %>%
  ungroup()

# Remove rows with NA 
pc_CMR <- pc_CMR %>%
  filter(!is.na(first_detection))

# Add a unique identifier  
pc_CMR <- pc_CMR %>%
  mutate(UniqueID = paste(PointNum, Survey, NOBOnum, sep = "_"))

# Recode Survey to Visit: Survey 2 -> Visit 1, Survey 3 -> Visit 2, Survey 4 -> Visit 3
pc_CMR <- pc_CMR %>%
  mutate(Visit = case_when(
    Survey == 2 ~ 1,
    Survey == 3 ~ 2,
    Survey == 4 ~ 3
  ))

# Inspect
pc_CMR %>%
  select(PointNum, 
         Visit, 
         NOBOnum, 
         int1, int2, int3, int4, 
         first_detection) %>%
  head(15) %>%
  print()

# Define constants
S <- 10  # Number of sites
K <- 3   # Number of survey occasions 
J <- 4   # Number of time intervals

# Create ordered site list
site_list <- sort(unique(pc_CMR$PointNum))

# Check visits at each site 
print(table(pc_CMR$PointNum, pc_CMR$Visit))

# All 10 sites were surveyed twice
all_sites <- 1:S
all_visits <- 1:K

complete_surveys <- expand.grid(
  PointNum = all_sites,
  Visit = all_visits
) %>%
  arrange(PointNum, Visit)

print(complete_surveys)


# Number of individuals detected
nind <- nrow(pc_CMR)

# detection variable
y <- pc_CMR$first_detection

# Create site assignment for each individual
site <- pc_CMR$PointNum

# Create visit assignment for each individual
visit <- pc_CMR$Visit

# Data Augmentation; M = augmented population size
M <- nind + 200

# Augment y
# Detected individuals: 1-4 (their first detection interval)
# Augmented individuals: 5 (never detected)
y_aug <- c(y, rep(5, M - nind))

# Augment site assignment
# Detected individuals: known sites (2-10)
# Augmented individuals: NA (unknown, to be estimated)
site_aug <- c(site, rep(NA, M - nind))

# Augment visit assignment
# Detected individuals: known visits (1 or 2)
# Augmented individuals: NA (unknown, to be estimated)
visit_aug <- c(visit, rep(NA, M - nind))

# ----------------------
# Covariates
# ----------------------

# Convert Wind to numeric
pc_CMR <- pc_CMR %>%
  mutate(Wind = as.numeric(Wind.Beau.Code))

# Get wind cov for every visit 
survey_covariates <- pc_dat %>%
  filter(Survey %in% c(2, 3, 4)) %>%
  mutate(Visit = case_when(
    Survey == 2 ~ 1,
    Survey == 3 ~ 2,
    Survey == 4 ~ 3
  )) %>%
  group_by(PointNum, Visit) %>%
  summarise(
    Wind = first(Wind.Beau.Code),
    Survey = first(Survey),
    Date = first(Date),
    .groups = "drop"
  ) %>%
  arrange(PointNum, Visit)

# Create Wind matrix: site x visit 
Wind_matrix <- matrix(NA, nrow = S, ncol = K)

for(i in 1:nrow(survey_covariates)) {
  site_num <- survey_covariates$PointNum[i]
  visit <- survey_covariates$Visit[i]
  Wind_matrix[site_num, visit] <- survey_covariates$Wind[i]
}

# Add 1 to all Wind values, currently starts at 0
Wind_matrix <- Wind_matrix + 1
print(Wind_matrix)

# Number of levels
Wind_Lvls <- length(unique(as.vector(Wind_matrix)))


# Rename SiteID to PointNum for matching
colnames(site_covs)[1] <- "PointNum"

# Detection covariates  
VegDens <- as.vector(scale(site_covs$vegDens50m))

# Abundance covariates 
Herb_COH <- as.vector(scale(site_covs$herb_COH))
Woody_SPLIT <- as.vector(scale(site_covs$woody_SPLIT))

# ----------------------
# Prepare prediction data
# ----------------------

# Scale prediction covariates 
Herb_COH_pred <- as.numeric(scale(predict_covs[,'herb_COH']))
Woody_SPLIT_pred <- as.numeric(scale(predict_covs[,'woody_SPLIT']))

# Inspect
print(Herb_COH_pred)
print(Woody_SPLIT_pred)

# Calculate area ratios
A_sampled <- pi * (227^2) # ~Area of NOBO home range/scale covs were extracted at
A_unsampled <- predict_covs$area_m2
rho <- A_unsampled / A_sampled

# Number of unsampled sites
U <- nrow(predict_covs)

# ------------------------------------------
# Bundle all data for JAGS
# ------------------------------------------
data <- list(
  # Sample size
  M = M,
  S = S,
  K = K,
  nind = nind,
  # Individual-level data
  y = y_aug,
  site = site_aug,
  visit = visit_aug,
  # Survey-level covariates
  Wind = Wind_matrix,    
  Wind_Lvls = Wind_Lvls,
  # Site-level covariates
  VegDens = VegDens,      
  Herb_COH = Herb_COH,      
  Woody_SPLIT = Woody_SPLIT,
  # Prediction
  U = U,
  Herb_COH_pred = Herb_COH_pred,
  Woody_SPLIT_pred = Woody_SPLIT_pred,
  rho = rho
  
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
n_iter <- 300000
n_burnin <- 100000
n_chains <- 3
n_thin <- 50
n_adapt <- 5000

# Test Settings
# n_iter <- 10000
# n_burnin <- 0
# n_chains <- 6
# n_thin <- 10
# n_adapt <- 5000

# Posterior samples
n_post_samps = ((n_iter * n_chains) - n_burnin) / n_thin
print(n_post_samps)

# ----------------------
# Model Specifications
# ----------------------

# Parameters monitored
params <- c(
  "lambda",
  "lambda_pred",
  "N",
  "N_pred",
  "N_samp_tot",
  "N_pred_tot",
  "N_tot",
  "samp_site_lambda_var",
  "pred_site_lambda_var",
  "samp_site_lambda_mean",
  "pred_site_lambda_mean",
  "samp_site_Nmean",
  "pred_site_Nmean",
  "alpha0", 
  "alpha1", 
  "alpha2", 
  "p",
  "beta0", 
  "beta1",
  "beta2",
  "psi",
  "sigma_sre",
  "S_RE",
  "tau_sre",
  "p_Bayes",
  "sum_obs",
  "sum_rep"
  
)

# Initial values
make_inits <- function() {
  list(
    # Detection parameters
    p0 = runif(1),
    alpha1 = rnorm(Wind_Lvls, 0, 1),
    alpha2 = rnorm(1, 0, 1),
    
    # Abundance parameters
    beta0 = rnorm(1, 0, 1),
    beta1 = rnorm(1, 0, 1),
    beta2 = rnorm(1, 0, 1),
    
    # Data augmentation variable
    z = c(rep(1, nind), rbinom(M - nind, 1, 0.5)),
    
    # Site random effects
    site_re = rnorm(S, 0, 1)
  )
}

# Initial Values for each chain
inits <- lapply(1:n_chains, function(x) make_inits())

# -------------------------------------------------------
# Model Statement
# -------------------------------------------------------
cat("
model {

  # ---------------------------------
  # Abundance Priors 
  # ---------------------------------
  
  # Intercept
  beta0 ~ dnorm(0, 0.1)
  
  # Covariates
  beta1 ~ dnorm(0, 0.1) # Herbaceous cohesion index
  beta2 ~ dnorm(0, 0.1) # Woody splitting index
  
  # ---------------------------------
  # Detection Priors 
  # ---------------------------------
  
  # Baseline detection probability (per interval)
  p0 ~ dunif(0, 1)
  alpha0 <- log(p0 / (1 - p0))   # Logit scale
  
  # Wind effect - categorical
  for(w in 1:Wind_Lvls) {
    alpha1[w] ~ dnorm(0, 0.1)
  }
  
  # Vegetation Density
  alpha2 ~ dnorm(0, 0.1)
  
  # Site random effect for pseudoreplication
  tau_sre ~ dgamma(0.01, 0.01)  
  sigma_sre <- 1 / sqrt(tau_sre) 
  for(s in 1:S){  
    S_RE[s] ~ dnorm(0, tau_sre)  
  }
  
  # Uniform visit probabilities
  for(k in 1:K) {
    visit_probs[k] <- 1 / K
  }
  
  # ---------------------------------
  # Abundance 
  # ---------------------------------
  
  for(s in 1:S) {

    log(lambda[s]) <- beta0 + beta1 * Herb_COH[s] + beta2 * Woody_SPLIT[s]
    
    # Probability of individual site membership
    probs[s] <- lambda[s] / sum(lambda[])
    
    # Estimated abundance 
    N[s] <- sum(z[] * equals(site_mem[], s))
  }
  
  # Inclusion probability
  psi <- sum(lambda[]) / M
  
  # ---------------------------------
  # Individual-Level Model
  # ---------------------------------
  
  for(i in 1:M) {
    
    # Presence/absence 
    z[i] ~ dbern(psi)
    
    # Site membership 
    site_mem[i] ~ dcat(probs[])
    
    # Visit assignment, uniform across K visits
    visit_mem[i] ~ dcat(visit_probs[])
    
    # Detection probability for individual
    logit(p[i]) <- alpha0 + alpha1[Wind[site_mem[i], visit_mem[i]]] + alpha2 * VegDens[site_mem[i]] + S_RE[site_mem[i]]
                   
                    
    # ---------------------------------
    # Cell probabilities
    # ---------------------------------

    ## pi[i,1] = detected in interval 1
    ## pi[i,2] = detected in interval 2 (not detected in 1)
    ## pi[i,3] = detected in interval 3 (not detected in 1 or 2)
    ## pi[i,4] = detected in interval 4 (not detected in 1, 2, or 3)
    ## pi[i,5] = never detected
    
    pi[i,1] <- p[i] * z[i]
    pi[i,2] <- (1 - p[i]) * p[i] * z[i]
    pi[i,3] <- pow(1 - p[i], 2) * p[i] * z[i]
    pi[i,4] <- pow(1 - p[i], 3) * p[i] * z[i]
    pi[i,5] <- 1 - sum(pi[i, 1:4])
    
    # Likelihood: categorical distribution
    y[i] ~ dcat(pi[i, ])
  }
  

  
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
    
    # Generate replicated data 
    y_rep[i] ~ dcat(pi[i,])
    
    for(j in 1:5){
      # Observed data indicator
      y_obs_ind[i,j] <- equals(y[i], j)
      
      # Expected probability
      exp_prob[i,j] <- pi[i,j]
      
      # Freeman-Tukey discrepancy for observed
      ft_obs[i,j] <- pow(sqrt(y_obs_ind[i,j]) - sqrt(exp_prob[i,j]), 2)
      
      # Freeman-Tukey discrepancy for replicated
      y_rep_ind[i,j] <- equals(y_rep[i], j)
      ft_rep[i,j] <- pow(sqrt(y_rep_ind[i,j]) - sqrt(exp_prob[i,j]), 2)
    }
    
    # Sum across categories for this individual
    discrepancy_obs[i] <- sum(ft_obs[i,])
    discrepancy_rep[i] <- sum(ft_rep[i,])
  }
  
  # Total discrepancy
  sum_obs <- sum(discrepancy_obs[])
  sum_rep <- sum(discrepancy_rep[])
  
  # Bayesian p-value
  p_Bayes <- step(sum_rep - sum_obs)
  
  
  # ---------------------------------
  # Derived Quantities
  # ---------------------------------
  
  # Abundance
  N_samp_tot <- sum(z[])
  
  # Abundance at unsampled sites
  N_pred_tot <- sum(N_pred[])
  
  # Total abundance
  N_tot <- sum(N[]) + sum(N_pred[])
  
  # Observed site mean abundances
  samp_site_Nmean <-  sum(N[]) / S
  pred_site_Nmean <- sum(N_pred[]) / U
  
  # Expected site mean abundance
  samp_site_lambda_mean <-  sum(lambda[]) / S
  pred_site_lambda_mean <- sum(lambda_pred[]) / U
  
  # Variance
  samp_site_lambda_var <- sum( (lambda[] - samp_site_lambda_mean)^2 ) / (S - 1)
  pred_site_lambda_var <- sum( (lambda_pred[] - pred_site_lambda_mean)^2 ) / (U - 1)

}
", fill = TRUE, file = "./jags_models/Model_TOD_MCR.txt")
# ------------End Model-------------


# -------------------------------------------------------
# Fit Model
# -------------------------------------------------------

fm1 <- jags(data = data, 
            parameters.to.save = params,
            inits = inits, 
            model.file = "./jags_models/Model_TOD_MCR.txt",
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
                   params = c(
                     "N_tot",
                     "N_samp_tot",
                     "N_pred_tot",
                     "beta0", 
                     "beta1",
                     "beta2",
                     "alpha0", 
                     "alpha1", 
                     "alpha2", 
                     'psi',
                     "sigma_sre",
                     "S_RE",
                     "samp_site_lambda_var",
                     "pred_site_lambda_var",
                     "samp_site_lambda_mean",
                     "pred_site_lambda_mean",
                     "samp_site_Nmean",
                     "pred_site_Nmean",
                     "lambda",
                     "lambda_pred",
                     "N",
                     "N_pred"
                   ),
                   pdf = T,
                   filename = "TracePlots_PC_TOD.pdf",
                   wd = "./Figures"
)

# Rhat
check_rhat(fm1$Rhat, threshold = 1.1)


# Save model
# saveRDS(fm1, "./Data/Model_Data/Fit_Model_PC_MCR.rds")

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
  annotate("text", x = 10, y = 0.05, label = paste0("Bayes p-value = ", mn_bpy), hjust = 0) 

# View
print(y_PPC_Dens)


# -------------------------------------------------------
#
#                 Posterior Estimates  
#
# -------------------------------------------------------

# Combine chains
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
  value = c(beta0_samples, 
            beta1_samples,
            beta2_samples),  
  parameter = rep(c("beta0", 
                    "beta1", 
                    "beta2"), each = length(beta0_samples))
) %>%
  group_by(parameter) %>%
  filter(value >= quantile(value, 0.025) & value <= quantile(value, 0.975))  

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
alpha2_samples <- combined_chains[, "alpha2"]


# Extract site random effect
sRE_cols  <- grep("^S_RE\\[", colnames(combined_chains), value = TRUE)
sRE_samples <- combined_chains[, sRE_cols]
sRE_samples <- rowMeans(sRE_samples) # Row means



# Compute 95% CI  
alpha_df <- data.frame(
  value = c(alpha0_samples, 
            alpha1_samples_1, 
            alpha1_samples_2,
            alpha1_samples_3,
            alpha2_samples,
            sRE_samples),  
  parameter = rep(c("alpha0", 
                    "alpha1_1", 
                    'alpha1_2',
                    'alpha1_3',
                    "alpha2",
                    "S_RE"), each = length(alpha0_samples))
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
#   Estimating Abundance --- set up to extrapolate from surveyed sites N estimate
# -------------------------------------------------------

# # Extract abundance posterior
# Ntot_samples <- combined_chains[ ,"N_samp_tot"]
# 
# # Ntotal is the abundance based on 10 point counts at a radius of 200m.
# # To correct for density, Ntotal needs to be divided by 10 * area surveyed
# area <- pi * (200^2) / 10000  # Area in acres
# dens_samples <- Ntot_samples / (area * 10)
# 
# # Create data frame for density
# dens_df <- data.frame(Model = rep(model_name, length(dens_samples)), Density = dens_samples)
# colnames(dens_df)[2] <- "Density"
# head(dens_df)
# 
# # Calculate the mean and 95% Credible Interval
# dens_summary <- dens_df %>%
#   group_by(Model) %>%
#   summarise(
#     Mean = mean(Density),
#     Lower_CI = quantile(Density, 0.025),
#     Upper_CI = quantile(Density, 0.975)
#   )
# 
# # Subset the data within the 95% credible interval
# dens_df <- dens_df[dens_df$Density >= dens_summary$Lower_CI & dens_df$Density <= dens_summary$Upper_CI, ]
# 
# # Getting total abundance
# abund_summary <- dens_summary
# abund_summary[,2:4] <- abund_summary[,2:4] * 1096
# 
# # Plot Abundance - Violin
# abund_df <- dens_df
# abund_df$Density <- abund_df$Density * 1096
# 
# ggplot(abund_df, aes(x = Model, y = Density, fill = Model)) + 
#   geom_violin(trim = FALSE, alpha = 0.6, adjust = 5) +   
#   labs(x = "Model", y = "Total Abundance ") +
#   scale_fill_manual(values = c("PC CMR" = "orange")) +   
#   scale_y_continuous(limits = c(0, 1000),
#                      breaks = seq(0, 1000, by = 100),
#                      labels = scales::comma) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  
#         axis.title.x = element_text(face = "bold", margin = margin(t = 10)),  
#         axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
#         panel.grid = element_blank(),
#         legend.position = "none")  
# 
# 
# 
# # Total abundance
# abund_summary
# 
# # Clear
# dev.off()
# 
# # Export Dataframes
# saveRDS(dens_df, "./Data/Model_Data/Density_df_PC-CMR.rds")
# saveRDS(dens_summary, "./Data/Model_Data/Density_summary_PC-CMR.rds")
# saveRDS(abund_df, "./Data/Model_Data/Abundance_df_PC-CMR.rds")
# saveRDS(abund_summary, "./Data/Model_Data/Abundance_summary_PC-CMR.rds")
