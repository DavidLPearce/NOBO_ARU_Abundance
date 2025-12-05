# Author: David L. Pearce
# Description:
#             TBD

# This code extends code from the following sources: 
#     1. Chambert, T., Waddle, J. H., Miller, D. A., Walls, S. C., 
#           and Nichols, J. D. (2018b). A new framework for analysing 
#           automated acoustic species detection data: Occupancy estimation 
#           and optimization of recordings post-processing. 
#           Methods in Ecology and Evolution, 9(3):560–570.
#     2. Kery, M. and Royle, J. A. (2020). Applied hierarchical modeling 
#           in ecology: Analysis of distribution, abundance, and species 
#           richness in R and BUGS: Volume 2: Dynamic and advanced models. 
#           Academic Press.
#     3. Doser, J. W., A. O. Finley, A. S. Weed, and E. F. Zipkin. 2021. 
#           Integrating automated acoustic vocalization data and point count 
#           surveys for estimation of bird abundance. 
#           Methods in Ecology and Evolution 12:1040–1049.

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

# Set seed, scientific notation options, and working directory
set.seed(123)
options(scipen = 9999)
setwd(".")

# Setting up cores
Ncores <- parallel::detectCores()
print(Ncores) # Number of available cores
workers <- Ncores * 0.5 # For low background use 80%, for medium use 50% of Ncores
print(workers)

# Source custom function for checking Rhat values > 1.1
source("./Scripts/Functions/Rhat_check_function.R")

# Model name object
model_name <- "AV Bsong"

# -------------------------------------------------------
#
# Variable and Object Definitions ******* Finish this
# 
# -------------------------------------------------------

# beta0 = abundance intercept 
# beta.1 = abundance trend estimate
# alpha0 = prob (on logit scale) of detecting at least one vocalization at a site that is not occupied.
# alpha1 = additional prob (on logit scale) of detecting at least one vocalization at a site that is not occupied. 
# omega = mean # of false positive acoustic detections
# p = detection probability of an individual in point count data
# tau.day = precision for random day effect on true vocalization detection rate. 
# a.phi = overdispersion parameter for zero-truncated negative binomial. 
# gamma.1 = random day effect on true vocalization detection rate
# n_days = number of recording days.
# N = latent abundance process
# tp = true positive rate
# p_a = prob of detecting at least one vocalization in an acoustic recording
# v = acoustic vocalization data from clustering algorithm
# y = binary summary of v, takes value 0 if v = 0 and value 1 if v > 0. 
# c = point count data
# S = number of total sites
# J = number of repeat visits for acoustic data at each site
# J_A = max number of repeat visits at each acoustic data site. 
# S_val = number of sites where validation of acoustic data occurred
# days = variable used to index different recording days. 
# A_times = indexing variable used to determine specific indexes with v > 0.
# K = true number of acoustic vocalizations
# k = true number of manually validated acoustic vocalizations
# n = total number of manually validated acoustic vocalizations 
# sites_a = specific indices of sites where acoustic data were obtained
# R.val = number of validated sites for acoustic data
# Other variables not defined are for computation of Bayesian p-values. 

# -------------------------------------------------------
#
#                    Load Data
#
# -------------------------------------------------------

# BirdNet detections
Bsong_dat <- read.csv("./Data/Acoustic_Data/BirdSong_Classifier/Subset14day_BirdSong_NOBO.csv")

# ARU weather covariates
weather_dat <- read.csv("./Data/Acoustic_Data/ARU_weathercovs.csv")

# Surveyed Site covariates
site_covs <- read.csv("./Data/Acoustic_Data/ARU_siteCovs.csv")

# Unsampled area covariates
predict_covs <- read.csv("./Data/Acoustic_Data/ARU_PredictsiteCovs.csv")

# -------------------------------------------------------
#
#    Data Wrangling
#
# -------------------------------------------------------


# Adding a row count
Bsong_dat$Count <- 1

# Initialize a site by survey matrix
v <- matrix(0, nrow = 27, ncol = 14)        

# Extract count data
for (i in 1:nrow(Bsong_dat)) {
  
  # Extracting plot ID
  site <- Bsong_dat$Site_Number[i]
  
  # Extracting Occasion
  occasion <- Bsong_dat$Occasion[i]
  
  # Fill in the matrix with the number of individuals
  v[site, occasion] =  v[site, occasion] + Bsong_dat$Count[i]
  
} # end loop 

# Take a look
print(v)
sum(v) # Total calls


# Renaming columns to date Month_day
date_order <- c("2024-05-26", "2024-05-30", "2024-06-03", "2024-06-07", 
                "2024-06-11", "2024-06-15", "2024-06-19", "2024-06-23",
                "2024-06-27", "2024-07-01", "2024-07-05", "2024-07-09",
                "2024-07-13", "2024-07-17")

formatted_dates <- format(as.Date(date_order), "%b_%d")
colnames(v) <- formatted_dates

# Adding rownames
rownames(v) <- as.numeric(1:27)

# Get the site numbers with at least one call
v_mat <- as.matrix(v)
sites_a <- which(rowSums(!is.na(v_mat) & v >= 1) > 0)
print(sites_a)

# Creating a binary detection matrix, values >= 1 become 1, and 0 stays 0
v_df <- as.data.frame(v)
y <- v_df %>% mutate(across(where(is.numeric), ~ if_else(. >= 1, 1, 0)))
y <- as.matrix(y)
rownames(y) <- NULL  
colnames(y) <- NULL 
print(y)

# Total number of sites
S <- as.integer(27) 

# Number of repeat visits for each site
J <- rep(14, S)  

# J_r contains the number of surveys at each acoustic site that contains at least 1 detected vocalization. 
J_r <- apply(v_mat, 1, function(a) {sum(a != 0)})
J_r <- ifelse(is.na(J_r), 0, J_r)
J_r <- as.numeric(J_r)
print(J_r)

# A_times is a site by survey matrix that contains the indices of the
# specific surveys at each acoustic site that have at least one detected
# vocalization. This determines the specific components of v[s j] that
# are used in the zero-truncated Poisson vocalization model.
A_times <- matrix(NA, S, max(J))
tmp <- apply(v, 1, function(a) {which(a != 0)})
for (i in 1:S) {
  if (length(tmp[[i]]) > 0) {
    A_times[i, 1:length(tmp[[i]])] <- tmp[[i]]
  }
}
print(A_times)

# ----------------------
# Manually validated  
# ----------------------


# Validated Calls
# Do not include sites with no calls 
# only include occasions where at least 1 call was validated for a site
n <- read.csv("./Data/Acoustic_Data/BirdSong_Classifier/Subset14day_BirdSong_n.csv", row.names = 1)
n <- as.matrix(n)
 
# True Calls
# Calls found to be true, same dimension as n
k <- read.csv("./Data/Acoustic_Data/BirdSong_Classifier/Subset14day_BirdSong_k.csv", row.names = 1)
k <- as.matrix(k)

# Survey days calls were validated, same dimension as n
val_times <- read.csv("./Data/Acoustic_Data/BirdSong_Classifier/Subset14day_BirdSong_valtimes.csv", row.names = 1)
val_times <- as.matrix(val_times)

# Total number of sites with manually validated data
S_val <- nrow(n)

# How many surveys were validate
J_val <- rep(14, S_val)  

# Check dimensions
dim(n) 
dim(k)
dim(val_times)
 
# ----------------------
# Covariates 
# ----------------------

# survey random effect index
days <- matrix(rep(1:14, times = 27), nrow = 27, ncol = 14, byrow = TRUE)

# Format abundance covariates
Herb_COH <- as.matrix(scale(site_covs[,'herb_COH']))  
Woody_SPLIT <- as.matrix(scale(site_covs[, 'woody_SPLIT']))  

# Inspect
head(Herb_COH)
head(Woody_SPLIT)

# Format detection covariates
Wind <- as.matrix(scale(weather_dat[, 'Wind_mph']))
VegDens <- scale(site_covs[,'vegDens50m'])

# ----------------------
# Prepare prediction data
# ----------------------

# Scale prediction covariates using the SAME scaling from training data
pred_Herb_COH <- scale(predict_covs[,'herb_COH'])
pred_Woody_SPLIT <- scale(predict_covs[,'woody_SPLIT'])

# Calculate area ratios
A_sampled <- pi * (227^2)
A_unsampled <- predict_covs$area_m2
rho <- A_unsampled / A_sampled

# Number of unsampled sites
U <- nrow(predict_covs)

# ----------------------
# Bayesian P-value
# ----------------------

S_A <- sum(J_r > 0)
sites_a_v <- which(J_r > 0)
J_A <- max(J)

# ----------------------
# Bundle Data 
# ----------------------

data <- list(S = S, 
             J = J, 
             v = v, 
             y = y,
             n = n,
             k = k, 
             val_times = val_times, 
             sites_a = sites_a, 
             S_val = S_val, 
             J_val = J_val, 
             J_r = J_r, 
             A_times = A_times, 
             S_A = S_A, 
             J_A = J_A, 
             sites_a_v = sites_a_v, 
             n_days = max(J),
             Herb_COH = as.numeric(Herb_COH),
             Woody_SPLIT = as.numeric(Woody_SPLIT),
             Wind = as.numeric(Wind),
             VegDens = as.numeric(VegDens),
             U = U,
             Herb_COH_pred = as.numeric(pred_Herb_COH),
             Woody_SPLIT_pred = as.numeric(pred_Woody_SPLIT),
             rho = rho
)

# Check structure
str(data)


# ---------------------------------------------------------- 
# 
#           Acoustic HM with Validated Calls
# 
# ----------------------------------------------------------

# ----------------------
# MCMC Specifications
# ----------------------
n_iter <- 300000
n_burnin <- 100000
n_chains <- 3
n_thin <- 50
n_adapt <- 5000

# Test Settings
# n_iter <- 20000
# n_burnin <- 0
# n_chains <- 6
# n_thin <- 10
# n_adapt <- 5000

# posterior samples
post_samps = (((n_iter - n_burnin) / n_thin) * n_chains)
print(post_samps)


# ----------------------
# Model Specifications
# ----------------------

# Parameters monitored
params <- c(# Abundance
            'beta0',  
            'beta1',
            'beta2',
            'tau_N',
            'sigma_N',
            'sigma2_N',
            'mu', 
            'N',
            'N_samp_tot',
            'mu_pred',
            'N_pred',
            'N_pred_tot',
            'N_tot',
            'samp_site_mu_var',
            'pred_site_mu_var',
            'samp_site_mu_mean',
            'pred_site_mu_mean',
            'samp_site_Nmean',
            'pred_site_Nmean',
            # Detection
            'p_a', 
            'mu_alpha',
            'alpha0', 
            'alpha1', 
            'alpha2',
            'alpha3',
            #  Vocalization
            'gamma0',
            'tau_jre',
            'sigma_jre',
            'J_RE',
            'omega',
            'delta',
            'phi',
            'r_phi',
            # Posterior Predictive Checks
            'fit_y',   
            'fit_y_pred',
            'fit_v',
            'fit_v_pred',
            'bp_y',
            'bp_v')

# Initial Values 
make_inits <- function() {
  list(
    # Abundance
    N      = rep(1, S),
    beta0  = rnorm(1, 0, 1),
    beta1  = rnorm(1, 0, 1),
    beta2  = rnorm(1, 0, 1),
    # Detection
    alpha1 = runif(1, 0, 1), 
    alpha2 = rnorm(1, 0, 1),
    alpha3 = rnorm(1, 0, 1),
    # Vocalization
    gamma0 = runif(1, log(6), log(15)),
    omega  = runif(1, 0, 0.25)
  )
}

# Initial Values for each chain
inits <- lapply(1:n_chains, function(x) make_inits())


# ----------------------------- 
# Model Statement 
# ----------------------------- 
cat(" model {
  
  # ----------------------
  # Abundance Priors
  # ----------------------
  
  # Intercept
  beta0 ~ dnorm(0, 0.1) 
  
  # Covariate effect
  beta1 ~ dnorm(0, 0.1) # Herbaceous
  beta2 ~ dnorm(0, 0.1) # Woody 
  
  # Underdispersion: variance < mean
  sigma_N ~ dunif(0.1, 1)
  tau_N = 1/pow(sigma_N, 2)
  sigma2_N <- 1/tau_N
  
  # variance (sigma2), standard deviation (sigma), precision (tau) relationship
  #     variance = (stddev)^2 = 1/precision
  #     stddev = 1/sqrt(precision) = sqrt(variance)
  #     precision = 1/stddev^2 = 1/variance
  
  # ------------------------
  # Detection Priors
  # ------------------------
  
  # Intercept
  alpha0 <- logit(mu_alpha)
  mu_alpha ~ dunif(0, 1)
  
  # True individuals
  alpha1 ~ dunif(0, 1000) # Constrained to be positive
  
  # Covariate effect
  alpha2 ~ dnorm(0, 0.1) # Wind
  alpha3 ~ dnorm(0, 0.1) # Vegetation Density
  
  # ------------------------
  # Call Rate Priors
  # ------------------------
  
  # False positive rate
  omega  ~ dunif(0, 1000) # From Doser et al. 
  
  # base vocal rate/30 min: 
  # ~2 calls/10 mins * 3 ten min periods = 6 calls/30 mins
  # must be on the log scale to match the log-link function on delta,
  # making the intercept 6 calls on the natural scale.
  gamma0 ~ dnorm(log(6), 1/pow(0.75, 2)) T(0,)
  
  
  # Survey random effect
  tau_jre ~ dgamma(0.01, 0.01)  
  sigma_jre <- 1 / sqrt(tau_jre)
  for (j in 1:n_days) {
     J_RE[j] ~ dnorm(0, tau_jre)
  }

  # Overdispersion
  r_phi ~ dgamma(0.01, 0.01)
  for (s in 1:S) {
    for (j in 1:J_A) {
      phi[s, j] ~ dgamma(r_phi, r_phi)
    }
  }

  # -------------------------------------------
  #
  # Likelihood and Process Model 
  #
  # -------------------------------------------

  # Site
  for (s in 1:S) {
    
    # ---------------------------------
    # Abundance  
    # ---------------------------------
    
    # ztNormal
    log(mu[s]) <- beta0 + beta1 * Herb_COH[s] +  beta2 * Woody_SPLIT[s]
    N[s] ~ dnorm(mu[s], tau_N) T(0,)


    # Survey
    for (j in 1:J[s]) {
    
    # ---------------------------------
    # Detection   
    # ---------------------------------
    logit(p_a[s, j]) <- alpha0 + alpha1 * N[s]  + alpha2 * Wind[j] + alpha3 * VegDens[s] 

    # ---------------------------------
    # Call rate  
    # ---------------------------------
    
    # Survey Random Effect
    log(delta[s, j]) <- gamma0  + J_RE[j] # + gamma1[Sky[j]]
    
    # ---------------------------------
    # Observations
    # ---------------------------------
    y[s, j] ~ dbin(p_a[s, j], 1)

    # ---------------------------------
    # True Positives 
    # ---------------------------------
    tp[s, j] <- delta[s, j] * N[s] / (delta[s, j] * N[s] + omega)

    # ---------------------------------
    # PPC Abundance  
    # ---------------------------------
    y_pred[s, j] ~ dbin(p_a[s, j], 1)
    resid_y[s, j] <- pow(pow(y[s, j], 0.5) - pow(p_a[s, j], 0.5), 2)
    resid_y_pred[s, j] <- pow(pow(y_pred[s, j], 0.5) - pow(p_a[s, j], 0.5), 2)

    } # End J
    
    # Surveys with Vocalizations
    for (j in 1:J_r[s]) {
  
    # ---------------------------------
    # Vocalizations  
    # ---------------------------------
    
    # Zero Truncated Negative Binomial
    v[s, A_times[s, j]] ~ dpois((delta[s, A_times[s, j]] * N[s] + omega) * phi[s, A_times[s, j]] * y[s, A_times[s, j]]) T(1, )

    # ---------------------------------
    # PPC calls  
    # ---------------------------------
    v_pred[s, j] ~ dpois((delta[s, A_times[s, j]] * N[s] + omega) * phi[s, A_times[s, j]] * y[s, A_times[s, j]]) T(1, )
    mu_v[s, j] <- ((delta[s, A_times[s, j]] * N[s] + omega) * phi[s, A_times[s, j]]) / (1 - exp(-1 * ((delta[s, A_times[s, j]] * N[s] + omega) * phi[s, A_times[s, j]])))
    resid_v[s, j] <- pow(pow(v[s, A_times[s, j]], 0.5) - pow(mu_v[s, j], 0.5), 2)
    resid_v_pred[s, j] <- pow(pow(v_pred[s, j], 0.5) - pow(mu_v[s, j], 0.5), 2)
    
    } # End J_r
  } # End S
  
  # ------------------------------------------- 
  # Manual validation 
  # -------------------------------------------
  for (s in 1:S_val) {
    for (j in 1:J_val[s]) {
      K[s, j] ~ dbin(tp[sites_a[s], j], v[sites_a[s], val_times[s, j]])
      k[s, val_times[s, j]] ~ dhyper(K[s, j], v[sites_a[s], val_times[s, j]] - K[s, j], n[s, val_times[s, j]], 1)
    } # End J
  } # End S
  
  
  # -------------------------------------------
  # Predict Abundance at Unsampled Areas 
  # -------------------------------------------
  for (u in 1:U) {
  
    # ztNormal
    log(mu_pred[u]) <- log(rho[u]) + beta0 + beta1 * Herb_COH_pred[u] + beta2 * Woody_SPLIT_pred[u]
    N_pred[u] ~ dnorm(mu_pred[u], tau_N / rho[u]) T(0,)
 
  }
  
  # -------------------------------------------
  # PPC and Bayesian P-value
  # -------------------------------------------
  for (s in 1:S_A) {
    tmp_v[s] <- sum(resid_v[sites_a_v[s], 1:J_r[sites_a_v[s]]])
    tmp_v_pred[s] <- sum(resid_v_pred[sites_a_v[s], 1:J_r[sites_a_v[s]]])
  }
  fit_y <- sum(resid_y[sites_a, 1:J_A])
  fit_y_pred <- sum(resid_y_pred[sites_a, 1:J_A])
  fit_v <- sum(tmp_v[1:S_A])
  fit_v_pred <- sum(tmp_v_pred[1:S_A])
  bp_y <- step(fit_y_pred - fit_y)
  bp_v <- step(fit_v_pred - fit_v)
  
  # -------------------------------------------
  # Derive Parameters
  # -------------------------------------------
  
  # Abundance at sampled sites
  N_samp_tot <- sum(N[])
  
  # Abundance at unsampled sites
  N_pred_tot <- sum(N_pred[])
  
  # Total abundance
  N_tot <- sum(N[]) + sum(N_pred[])
  
  # Overall site mean abundances
  samp_site_Nmean <-  sum(N[]) / S
  pred_site_Nmean <- sum(N_pred[]) / U
  
  # Expected site mean abundance
  samp_site_mu_mean <-  sum(mu[]) / S
  pred_site_mu_mean <- sum(mu_pred[]) / U
  
  # Variance
  samp_site_mu_var <- sum( (mu[] - samp_site_mu_mean)^2 ) / (S - 1)
  pred_site_mu_var <- sum( (mu_pred[] - pred_site_mu_mean)^2 ) / (U - 1)
  
}
", fill = TRUE, file = "./JAGs_Models/Model_AV_Bsong.txt")
# ------------End Model-------------

# -------------------------------------------------------
# Fit Model
# -------------------------------------------------------

fm1 <- jagsUI::jags(data = data,
                     inits = inits,
                     parameters.to.save = params,
                     model.file = "./JAGs_Models/Model_AV_Bsong.txt",
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
                   params = c('N_tot',
                              'N_samp_tot',
                              'N_pred_tot',
                              'beta0',
                              'beta1',
                              'beta2',
                              'mu_alpha',
                              'alpha0',  
                              'alpha1', 
                              'alpha2',
                              'alpha3',
                              'r_phi',
                              'omega',
                              'gamma0',
                              'tau_jre',
                              'sigma_jre',
                              'J_RE',
                              'tau_N',
                              'sigma_N',
                              'sigma2_N',
                              'samp_site_mu_var',
                              'pred_site_mu_var',
                              'samp_site_mu_mean',
                              'pred_site_mu_mean',
                              'samp_site_Nmean',
                              'pred_site_Nmean',
                              'mu',
                              'mu_pred',
                              'N',
                              'N_pred',
                              'delta',
                              'phi'
                   ),
                   pdf = TRUE,
                   filename = "TracePlots_ARU_Bsong.pdf",
                   wd = "./Figures"
)

# Rhat
check_rhat(fm1$Rhat, threshold = 1.1) 

# Save model
saveRDS(fm1, "./Data/Model_Data/Fit_Model_ARU_Bsong.rds")

# -------------------------------------------------------
# Posterior Predictive Checks
# -------------------------------------------------------

# ----------------------
# Extract Fit
# ----------------------

# Abundance
fit_y_data <- data.frame(
  observed = fm1$sims.list$fit_y,  # Observed values
  predicted = fm1$sims.list$fit_y_pred,  # Predicted values
  type = rep(c("Observed", "Predicted"), each = length(fm1$sims.list$fit_y))
)

# Calls
fit_v_data <- data.frame(
  observed = fm1$sims.list$fit_v,  # Observed values
  predicted = fm1$sims.list$fit_v_pred,  # Predicted values
  type = rep(c("Observed", "Predicted"), each = length(fm1$sims.list$fit_v))
)


# ----------------------
# Density Plot
# ----------------------
 
# Bayes P-value
# P-value = 0.5 means good fit, = 1 or 0 is a poor fit
mn_bpy <- round(mean(fm1$summary["bp_y",1]), 2) 
mn_bpv <- round(mean(fm1$summary["bp_v",1]), 2)

# y
y_PPC_Dens <- ggplot(fit_y_data) +
              geom_density(aes(x = observed, fill = "Observed"), alpha = 0.5) +   
              geom_density(aes(x = predicted, fill = "Predicted"), alpha = 0.5) +  
              scale_fill_manual(values = c("Observed" = "blue", "Predicted" = "red")) +  
              labs(title = "A)", 
                   x = "", 
                   y = "Density") +
              theme_minimal() +
              theme(legend.title = element_blank()) +
              theme(legend.position = "none") + 
              annotate("text", x = 40, y = 0.13, label = paste0("Bayes p-value = ", mn_bpy), hjust = 0)

# v
v_PPC_Dens <- ggplot(fit_v_data) +
              geom_density(aes(x = observed, fill = "Observed"), alpha = 0.5) +   
              geom_density(aes(x = predicted, fill = "Predicted"), alpha = 0.5) +  
              scale_fill_manual(values = c("Observed" = "blue", "Predicted" = "red")) +  
              labs(title = "B)", 
                   x = "Fit Values", 
                   y = "Density") +
              theme_minimal() +
              theme(legend.title = element_blank()) +
              theme(legend.position = "none") + 
              annotate("text", x = 40, y = 0.06, label = paste0("Bayes p-value = ", mn_bpv), hjust = 0)

# Multipanel figure
grid.arrange(y_PPC_Dens, v_PPC_Dens, nrow = 2)

# Save to file
jpeg("Figures/PPC_ARU_Bsong.jpg", width = 10, height = 8, units = "in", res = 300)
grid.arrange(y_PPC_Dens, v_PPC_Dens, nrow = 2)
dev.off()

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
saveRDS(beta_df, "./Data/Model_Data/Beta_df_ARU_Bsong.rds")
write.csv(beta_summary, "./Data/Model_Data/Beta_Summary_ARU_Bsong.csv")


# -------------------------------------------------------
# Alpha Estimates
# -------------------------------------------------------

# Extract alpha estimates
alpha0_samples <- combined_chains[, "alpha0"]
alpha1_samples <- combined_chains[, "alpha1"]
alpha2_samples <- combined_chains[, "alpha2"]
alpha3_samples <- combined_chains[, "alpha3"]

# Compute 95% CI for each
alpha_df <- data.frame(
  value = c(alpha0_samples, 
            alpha1_samples, 
            alpha2_samples,
            alpha3_samples),  
  parameter = rep(c("alpha0", 
                    "alpha1", 
                    "alpha2",
                    "alpha3"), each = length(alpha0_samples))
) %>%
  group_by(parameter) %>%
  filter(value >= quantile(value, 0.025) & 
           value <= quantile(value, 0.975))  

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
saveRDS(alpha_summary, "./Data/Model_Data/Alpha_Summary_ARU_Bsong.rds")


# -------------------------------------------------------
# Detection probability 
# -------------------------------------------------------

# Extract samples
p_samples <- combined_chains[, grepl("^p_a\\[", colnames(combined_chains))]

# Combine
all_p_samples <- as.vector(p_samples)

# Create summary
p_summary <- data.frame(
  Mean = mean(all_p_samples),
  LCI = quantile(all_p_samples, 0.025),
  UCI = quantile(all_p_samples, 0.975)
)

# Add model
p_summary$Model <- model_name

# Add parameter name
p_summary$Parameter <- "Detection"

# view
print(p_summary)



# -------------------------------------------------------
# Vocalization Estimates
# -------------------------------------------------------

# Extract samples
delta_samples <- as.matrix(fm1$samples)[, grepl("^delta\\[", colnames(as.matrix(fm1$samples)))]

# Combine
all_deltas <- as.vector(delta_samples)

# Create summary
delta_summary <- data.frame(
  Mean = mean(all_deltas),
  LCI = quantile(all_deltas, 0.025),
  UCI = quantile(all_deltas, 0.975)
)

# Add model
delta_summary$Model <- model_name

# Add parameter name
delta_summary$Parameter <- "Vocal Rate"

# View
print(delta_summary)

# Combine with detection
param_summary <- rbind(p_summary, delta_summary)

# View
print(param_summary)

# Random Effect on vocal rate
jRE_cols  <- grep("^J_RE\\[", colnames(combined_chains), value = TRUE)
jRE_samples <- combined_chains[, jRE_cols]
jRE_samples <- rowMeans(jRE_samples) # Row means

# Extract gamma0 estimates
gamma0_samples <- combined_chains[, "gamma0"]

# Compute 95% CI for each
gamma_df <- data.frame(
  value = c(gamma0_samples,
            jRE_samples),  
  parameter = rep(c("gamma0",
                    "jRE"), each = length(gamma0_samples))
) %>%
  group_by(parameter) %>%
  filter(value >= quantile(value, 0.025) & 
           value <= quantile(value, 0.975))  

# Add model
gamma_df$Model <- model_name

# Create summary
gamma_summary <- gamma_df %>%
  group_by(parameter, Model) %>%
  summarise(
    mean = mean(value),
    LCI = min(value),
    UCI = max(value),
    .groups = "drop"  
  )

# Exponentiate to get on real scale
gamma_summary[3,3:5] <- exp(gamma_summary[1,3:5])
gamma_summary$Model <- model_name
gamma_summary$parameter <- "exp(gamma0)"

# Take a look
print(gamma_summary)

# Export gamma summary
saveRDS(gamma_summary, "./Data/Model_Data/Gamma_Summary_ARU_Bsong.rds")



# -------------------------------------------------------
#  Predict Abundance 
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


# Export abundance df and parameter summary
saveRDS(total_N_samples_95CI, "./Data/Model_Data/Abund_df_ARU_Bsong.rds")
saveRDS(param_summary, "./Data/Model_Data/Param_Summary_ARU_Bsong.rds")



# ------------ End Script -----------------
