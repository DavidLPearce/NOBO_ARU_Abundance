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
workers <- Ncores * 0.3 # For low background use 80%, for medium use 50% of Ncores
print(workers)

# Source custom function for checking Rhat values > 1.1
source("./Scripts/Rhat_check_function.R")

# Model name object
model_name <- "AV Bnet"

# -------------------------------------------------------
#
#             Variable and Object Definitions ******* Finish this
#
# -------------------------------------------------------

# beta0 = abundance intercept 
# beta.1 = abundance trend estimate
# alpha0 = prob (on logit scale) of detecting at least one vocalization at a site
#           that is not occupied.
# alpha1 = additional prob (on logit scale) of detecting at least one vocalization at 
#           a site that is not occupied. 
# omega = mean # of false positive acoustic detections
# p = detection probability of an individual in point count data
# tau.day = precision for random day effect on true vocalization detection rate. 
# a.phi = overdispersion parameter for zero-truncated negative binomial. 
# gamma.1 = random day effect on true vocalization detection rate
# N = latent abundance process
# tp = true positive rate
# p_a = prob of detecting at least one vocalization in an acoustic recording
# v = acoustic vocalization data from clustering algorithm
# y = binary summary of v, takes value 0 if v = 0 and value 1 if v > 0. 
# c = point count data
# R = number of total sites
# J = number of repeat visits for acoustic data at each site
# J.A = max number of repeat visits at each acoustic data site. 
# n.count = number of repeat visits for count data
# R.val = number of sites where validation of acoustic data occurred
# A_times = indexing variable used to determine specific indexes with v > 0.
# K = true number of acosutic vocalizations
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
bnet_dat <- read.csv("./Data/Acoustic_Data/BirdNET_Classifier/FullDeploy_BirdNET.csv")

# ARU weather covariates
weather_dat <- read.csv("./Data/Acoustic_Data/ARU_weathercovs.csv")

# Site covariates
site_covs <- read.csv("./Data/Acoustic_Data/ARU_siteCovs.csv")

# -------------------------------------------------------
#
#    Data Wrangling
#
# -------------------------------------------------------

# Subset to 14 days starting at May 26 and ending on July 17. Dates are every 4 days.
date_order <- c("2024-05-26", "2024-05-30", "2024-06-03", "2024-06-07", # in ascending order
                "2024-06-11", "2024-06-15", "2024-06-19", "2024-06-23",
                "2024-06-27", "2024-07-01", "2024-07-05", "2024-07-09",
                "2024-07-13", "2024-07-17")

# Dates and their corresponding occasion numbers
bnet_dat <- bnet_dat %>% filter(Date %in% date_order)
head(bnet_dat)
nrow(bnet_dat)

# Adding occasion column
bnet_dat <- bnet_dat %>%
  mutate(Occasion = match(Date, date_order))

# ----------------------------
# Observation Matrix
# ----------------------------

# Adding a row count
bnet_dat$Count <- 1


# Initialize a site by survey matrix
v <- matrix(0, nrow = 27, ncol = 14)        

# Extract count data
for (i in 1:nrow(bnet_dat)) {
  
  # Extracting plot ID
  site <- bnet_dat$Site_Number[i]
  
  # Extracting Occasion
  occasion <- bnet_dat$Occasion[i]
  
  # Fill in the matrix with the number of individuals
  v[site, occasion] =  v[site, occasion] + bnet_dat$Count[i]
  
} # end loop 

# Take a look
print(v)
sum(v) # Total calls

# Renaming columns to date Month_day
formatted_dates <- format(as.Date(date_order), "%b_%d")
colnames(v) <- formatted_dates
print(v) 

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

# A_times is a R x J matrix that contains the indices of the
# specific surveys at each acoustic site that have at least one detected
# vocalization. This determines the specific components of v[i, j] that
# are used in the zero-truncated Poisson.
A_times <- matrix(NA, S, max(J))
tmp <- apply(v, 1, function(a) {which(a != 0)})
for (i in 1:S) {
  if (length(tmp[[i]]) > 0) {
    A_times[i, 1:length(tmp[[i]])] <- tmp[[i]]
  }
}
print(A_times)

# ---------------------------------------
# Manually validated  
# ---------------------------------------

# Validated Calls
# Do not include sites with no calls
# only include occasions where at least 1 call was validated for a site
n <- read.csv("./Data/Acoustic_Data/BirdNET_Classifier/Subset14day_BirdNET_n.csv", row.names = 1)
n <- as.matrix(n)

# True Calls
# Calls found to be true, same dimension as n
k <- read.csv("./Data/Acoustic_Data/BirdNET_Classifier/Subset14day_BirdNET_k.csv", row.names = 1)
k <- as.matrix(k)

# Survey days calls were validated, same dimension as n
val_times <- read.csv("./Data/Acoustic_Data/BirdNET_Classifier/Subset14day_BirdNET_valtimes.csv", row.names = 1)
val_times <- as.matrix(val_times)

# Total number of sites with manually validated data
S_val <- nrow(n)

# How many surveys were validate
J_val <- rep(14, S_val)


# ---------------------------------------
# Covariates 
# ---------------------------------------

# survey random effect index
days <- matrix(rep(1:14, times = 27), nrow = 27, ncol = 14, byrow = TRUE)

# Format abundance covariates
Herb_COH <- as.matrix(scale(site_covs[,'herb_ClmIdx'])) # herb_COH
Woody_SPLIT <- as.matrix(scale(site_covs[, 'woody_AggInx']))# woody_SPLIT

# Inspect
head(Herb_COH)
head(Woody_SPLIT)

# Format detection covariates
Wind <- as.matrix(scale(weather_dat[, 'Wind_mph']))
VegDens <- scale(site_covs[,'vegDens50m'])

# Area surveyed 
area <- pi * (200^2) / 10000  # in hectares
Offset <- rep(area, 27)


# ---------------------------------------
# Bayesian P-value
# ---------------------------------------

S_A <- sum(J_r > 0)
sites_a_v <- which(J_r > 0)
J_A <- max(J)


# ---------------------------------------
# Bundle Data 
# ---------------------------------------

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
             days = days,
             n.days = max(J),
             Herb_COH = Herb_COH,
             Woody_SPLIT = Woody_SPLIT,
             Wind = Wind,
             VegDens = VegDens,
             Offset = area)

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
n_iter = 800000
n_burnin = 200000
n_chains = 3 
n_thin = 10
n_adapt = 5000

# posterior samples
post_samps = (((n_iter - n_burnin) / n_thin) * n_chains)
print(post_samps)

# ----------------------
# Model Specifications
# ----------------------

# Parameters monitored
params <- c('mu', # Abundance
            'N_tot',
            'N',
            'beta0',
            'beta1',
            'beta2',
            'tau_s',
            'S_RE',
            'sigma_mu',
            'p_a',# Detection 
            'alpha0', 
            'alpha1', 
            'alpha2',
            'alpha3',
            'mu_j',#  Vocalization
            'tau_j',
            'J_RE',
            'omega',
            'delta',
            'phi',
            'r_phi',
            'lam_phi',
            'fit_y',# Posterior Predictive Checks
            'fit_y_pred',
            'fit_v',
            'fit_v_pred',
            'bp_y', # Bayes p-value
            'bp_v')



# Initial Values 
make_inits <- function() {
  list(
    N      = rep (1, S),  # Abundance
    beta0  = rnorm(1, 0, 1),
    beta1  = rnorm(1, 0, 1),
    beta2  = rnorm(1, 0, 1),
    alpha1 = runif(1, 0, 1),        # Detection
    alpha2 = rnorm(1, 0, 1),
    alpha3 = rnorm(1, 0, 1),
    omega  = runif(1, 0, 1)         # Vocalization
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
  beta0 ~ dnorm(0, 0.001) 
  
  # Covariate effect
  beta1 ~ dnorm(0, 0.001) # Herbaceous
  beta2 ~ dnorm(0, 0.001) # Woody 
 
  # Site random effect
  tau_s ~ dgamma(0.001, 0.001)
  for (s in 1:S) {
     S_RE[s] ~ dnorm(beta0, tau_s)
  }

  # Underdispersion
  for (s in 1:S) {
  sigma_mu[s] ~ dnorm(0,1) T(0,)
  }
  
  # ------------------------
  # Detection Priors
  # ------------------------
  
  # Intercept
  alpha0 <- logit(mu_alpha) # Constrains alpha0 to be between 0 and 1 on the logit scale (propability)
  mu_alpha ~ dunif(0, 1)
  
  # True individuals
  alpha1 ~ dunif(0, 1000) # Constrained to be positive
  
  # Covariate effect
  alpha2 ~ dnorm(0, 0.001) # Wind
  alpha3 ~ dnorm(0, 0.001) # Vegetation Density
  
  # ------------------------
  # Call Rate Priors
  # ------------------------
  
  # False positive rate
  omega  ~ dunif(0, 1000) # From Doser et al.  
  
  # Survey random effect
  tau_j ~ dgamma(0.001, 0.001)
  for (j in 1:n.days) {
     J_RE[j] ~ dnorm(0, tau_j)
  }

  # Overdispersion
  r_phi ~ dgamma(0.001, 0.001)
  lam_phi ~ dgamma(0.001, 0.001)
  for (s in 1:S) {
    for (j in 1:J_A) {
      phi[s, j] ~ dgamma(r_phi, lam_phi)
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
    # Abundance Submodel  
    # ---------------------------------
    
    # Normal 
    log(mu[s]) <- S_RE[s] + beta1 * Herb_COH[s, 1] +  beta2 * Woody_SPLIT[s, 1]
    N[s] ~ dnorm(mu[s], sigma_mu[s])
    
    # Log-Normal
    # mu[s] <-  S_RE[s] + beta1 * Herb_COH[s, 1] +  beta2 * Woody_SPLIT[s, 1]
    # N[s] ~ dlnorm(mu[s], sigma_mu[s]) 

    # Survey
    for (j in 1:J[s]) {

    # ---------------------------------
    # Detection Submodel  
    # ---------------------------------
    logit(p_a[s, j]) <- alpha0 + alpha1 * N[s]  + alpha2 * Wind[j,1] + alpha3 * VegDens[s, 1] 

    # ---------------------------------
    # Call rate Submodel  
    # ---------------------------------
    
    # Survey Random Effect
    log(delta[s, j]) <- J_RE[j]

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
  
  # Abundance
  N_tot <- sum(N[])

}
", fill = TRUE, file = "./JAGs_Models/AV_Bnet_Model.txt")
# ------------End Model-------------

# -------------------------------------------------------
# Fit Model
# -------------------------------------------------------

fm1 <- jagsUI::jags(data = data,
                    inits = inits,
                    parameters.to.save = params,
                    model.file = "./JAGs_Models/AV_Bnet_Model.txt",
                    n.iter = n_iter,
                    n.burnin = n_burnin,
                    n.chains = n_chains,
                    n.thin = n_thin,
                    n.adapt = n_adapt,
                    parallel = TRUE,
                    n.cores = workers,
                    verbose = TRUE,
                    DIC = FALSE)


# -------------------------------------------------------
# Check Convergence
# -------------------------------------------------------

# Rhat
check_rhat(fm1$Rhat, threshold = 1.1) 

# Trace plots
MCMCvis::MCMCtrace(fm1, 
                   params = c('mu',  
                              'N_tot',
                              'N',
                              'beta0',
                              'beta1',
                              'beta2',
                              'tau_s',
                              'S_RE',
                              'sigma_mu',
                              'alpha0',  
                              'alpha1', 
                              'alpha2',
                              'alpha3',
                              'mu_j', 
                              'tau_j',
                              'J_RE',
                              'omega',
                              'delta',
                              'phi',
                              'r_phi',
                              'lam_phi'
                   ),
                   pdf = T,
                   filename = "ARU_Bnet_TracePlots.pdf",
                   wd = "./Figures"
)


# -------------------------------------------------------
# Posterior Predictive Checks
# -------------------------------------------------------

# ----------------------
# Extract Residuals
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
  annotate("text", x = 30, y = 0.13, label = paste0("Bayes p-value = ", mn_bpy), hjust = 0)

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
  annotate("text", x = 10, y = 0.10, label = paste0("Bayes p-value = ", mn_bpv), hjust = 0)


# Multipanel figure
grid.arrange(y_PPC_Dens, v_PPC_Dens, nrow = 2)

# Save to file
jpeg("Figures/PPC_BP_ARU-Bnet.jpg", width = 10, height = 8, units = "in", res = 300)
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
saveRDS(beta_df, "./Data/Model_Data/ARU_Bnet_beta_df.rds")
write.csv(beta_summary, "./Figures/ARU_Bnet_BetaSummary.csv")



# -------------------------------------------------------
# Alpha Estimates
# -------------------------------------------------------

# Extract alpha estimates
alpha0_samples <- combined_chains[, "alpha0"]
alpha1_samples <- combined_chains[, "alpha1"]
alpha2_samples <- combined_chains[, "alpha2"]
alpha3_samples <- combined_chains[, "alpha3"]

# Extract site random effect
jRE_samples <- combined_chains[, c("J_RE[1]", "J_RE[2]","J_RE[3]","J_RE[4]","J_RE[5]",
                                   "J_RE[6]", "J_RE[7]", "J_RE[8]", "J_RE[9]", "J_RE[10]",
                                   "J_RE[11]", "J_RE[12]", "J_RE[13]", "J_RE[14]"
)]

jRE_samples <- rowMeans(jRE_samples) # Row means

# Compute 95% CI for each beta
alpha_df <- data.frame(
  value = c(alpha0_samples, 
            alpha1_samples, 
            alpha2_samples,
            alpha3_samples,
            jRE_samples),  
  parameter = rep(c("alpha0", 
                    "alpha1", 
                    "alpha2",
                    "alpha3",
                    "jRE"), each = length(alpha0_samples))
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
    .groups = "drop"   # optional: ungroups after summarising
  )

# view
print(alpha_summary)

# Export alpha summary
saveRDS(alpha_summary, "./Data/Model_Data/ARU_Bnet_alpha_summary.rds")


# -------------------------------------------------------
# Detection probability 
# -------------------------------------------------------

# Extract samples
p_samples <- combined_chains[, grepl("^p_a\\[", colnames(combined_chains))]

# Combine
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
# Vocalization Estimates
# -------------------------------------------------------

# Extract samples
delta_samples <- as.matrix(fm1$samples)[, grepl("^delta\\[", colnames(as.matrix(fm1$samples)))]

# Combine
all_deltas <- as.vector(delta_samples)

# Create summary
delta_summary <- data.frame(
  mean = mean(all_deltas),
  LCI = quantile(all_deltas, 0.025),
  UCI = quantile(all_deltas, 0.975)
)

# Add model
delta_summary$Model <- model_name

# Add parameter name
delta_summary$Parameter <- "vocal rate"

# View
print(delta_summary)

# Combine with detection
param_summary <- rbind(p_summary, delta_summary)

# View
print(param_summary)

# -------------------------------------------------------
#  Estimating Abundance 
# -------------------------------------------------------

# Extract abundance posterior
Ntot_samples <- combined_chains[ ,"N_tot"]

# Ntotal is the abundance based on 27 acoustic sites at a radius of 200m.

# Area in hectares of a 200m radius circle
area <- pi * (200^2) / 10000  # Area in hectares

# Calculate density (individuals per hectare)
dens_samples <- Ntot_samples / (area * 27)

# Create data frame for density
dens_df <- data.frame(Model = rep(model_name, length(dens_samples)), Density = dens_samples)
colnames(dens_df)[2] <- "Density"
head(dens_df)

# Abundance estimates
abund_df <- dens_df # posterior estimates
abund_df$Density <- abund_df$Density * 1096  

abund_summary <- abund_df %>% # summary
  group_by(Model) %>%
  summarise(
    mean = mean(Density),
    LCI = quantile(Density, 0.025),
    UCI = quantile(Density, 0.975)
  )

# Add model
abund_summary$Model <- model_name

# Add parameter name
abund_summary$Parameter <- "abundance"

# View
print(abund_summary)

# Combine with detection
param_summary <- rbind(param_summary, abund_summary)

# View
print(param_summary)

# Trim abund df to 95% CI
abund_95df <- abund_df %>%
  left_join(abund_summary, by = "Model") %>%
  filter(Density >= LCI & Density <= UCI) %>%
  select(-LCI, -UCI)


# Export abundance df and parameter summary
saveRDS(abund_95df, "./Data/Model_Data/ARU_Bnet_abund_df.rds")
saveRDS(param_summary, "./Data/Model_Data/ARU_Bnet_param_summary.rds")

# ------------ End Script -----------------