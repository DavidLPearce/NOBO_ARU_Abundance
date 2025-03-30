# Author: David L. Pearce
# Description:
#             TBD

# This code extends code from the following sources: 
#     1. Chambert, T., Waddle, J. H., Miller, D. A., Walls, S. C., 
#           and Nichols, J. D. (2018). A new framework for analysing 
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
# install.packages("nimble")
# install.packages("coda")
# install.packages("mcmcplots")
# install.packages("COMPoissonReg")

# Load library
library(tidyverse)
library(plotly)
library(nimble)
library(coda)
library(mcmcplots)
library(COMPoissonReg)

# Set seed, scientific notation options, and working directory
set.seed(123)
options(scipen = 9999)
setwd(".")

# Source custom function for checking Rhat values > 1.1
source("./Scripts/Custom_Functions/Rhat_check_function.R")

# Source, assign, and register hypergeometric distribution for NIMBLE
source("./Scripts/Custom_Functions/dhyper_nimble.R")
source("./Scripts/Custom_Functions/rhyper_nimble.R")
source("./Scripts/Custom_Functions/register_hypergeometric_distribution.R")

# Source, assign, and register Conway-Maxwell Poisson distribution for NIMBLE
source("./Scripts/Custom_Functions/dCMPois_nimble.R")
source("./Scripts/Custom_Functions/rCMPois_nimble.R")
source("./Scripts/Custom_Functions/register_Conway-Maxwell_Poisson_distribution.R")


# Model name object
model_name <- "AV Wolfe"

# -------------------------------------------------------
#
# Variable and Object Definitions *** Not finished
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
# n.days = number of recording days.
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

# Wolfe et al. detections
wolfe_dat <- read.csv("./Data/Acoustic_Data/NOBO_Wolfe_14day2024.csv")

# ARU weather covariates
weather_dat <- read.csv("./Data/Acoustic_Data/ARU_weathercovs.csv")

# Site covariates
site_covs <- read.csv("./Data/Acoustic_Data/ARU_siteCovs.csv")

# -------------------------------------------------------
#
#    Data Wrangling
#
# -------------------------------------------------------


# Adding a row count
wolfe_dat$Count <- 1

# Initialize a site by survey matrix
v <- matrix(0, nrow = 27, ncol = 14)        

# Extract count data
for (i in 1:nrow(wolfe_dat)) {
  
  # Extracting plot ID
  site <- wolfe_dat$Site_Number[i]
  
  # Extracting Occasion
  occasion <- wolfe_dat$Occasion[i]
  
  # Fill in the matrix with the number of individuals
  v[site, occasion] =  v[site, occasion] + wolfe_dat$Count[i]
  
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
sites_a <- as.vector(sites_a)
print(sites_a)

# Creating a binary detection matrix, values >= 1 become 1, and 0 stays 0
v_df <- as.data.frame(v)
y <- v_df %>% mutate(across(where(is.numeric), ~ if_else(. >= 1, 1, 0)))
y <- as.matrix(y)
rownames(y) <- NULL  
colnames(y) <- NULL 
print(y)

# ARU at site 24 stopped recording on 6/20. Making Columns 8:14 NA
# v[24,8:14] <- NA
# print(v)

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
n <- read.csv("./Data/Acoustic_Data/Wolfe14day_n.csv", row.names = 1)
n <- as.matrix(n)

# True Calls
# Calls found to be true, same dimension as n
k <- read.csv("./Data/Acoustic_Data/Wolfe14day_k.csv", row.names = 1)
k <- as.matrix(k)

# Survey days calls were validated, same dimension as n
val_times <- read.csv("./Data/Acoustic_Data/Wolfe14day_val_times.csv", row.names = 1)
val_times <- as.matrix(val_times)

# Total number of sites with manually validated data
S_val <- nrow(n)

# How many surveys were validate
J_val <- rep(14, S_val)  
J_val <- as.vector(unname(J_val))

# Check dimensions
dim(n) 
dim(k)
dim(val_times)



# ----------------------
# Covariates 
# ----------------------

# survey random effect index
days <- matrix(rep(1:14, times = 27), nrow = 27, ncol = 14, byrow = TRUE)

# Format X_abund
# X_abund <- as.matrix(site_covs[,-c(1:2)]) # Remove X and site id
# print(X_abund)

# Remove site, lat, long
X_abund <- site_covs[,-c(1:4)]

# Mean and center scale values 
# scale() also works, but stating to center or scale (scale(object, center = TRUE, scale = TRUE)) causes issues. 
# So, manually doing it to avoid what ever issue scale() was causing.
X_abund$woody_prp <- (X_abund$woody_prp - mean(X_abund$woody_prp, na.rm = TRUE)) / sd(X_abund$woody_prp, na.rm = TRUE) 
X_abund$herb_prp <- (X_abund$herb_prp - mean(X_abund$herb_prp, na.rm = TRUE)) / sd(X_abund$herb_prp, na.rm = TRUE)
X_abund$open_prp <- (X_abund$open_prp - mean(X_abund$open_prp, na.rm = TRUE)) / sd(X_abund$open_prp, na.rm = TRUE) 
X_abund$woody_mnParea <- (X_abund$woody_mnParea - mean(X_abund$woody_mnParea, na.rm = TRUE)) / sd(X_abund$woody_mnParea, na.rm = TRUE)
X_abund$herb_mnParea <- (X_abund$herb_mnParea - mean(X_abund$herb_mnParea, na.rm = TRUE)) / sd(X_abund$herb_mnParea, na.rm = TRUE)
X_abund$woody_ClmIdx <- (X_abund$woody_ClmIdx - mean(X_abund$woody_ClmIdx, na.rm = TRUE)) / sd(X_abund$woody_ClmIdx, na.rm = TRUE)
X_abund$herb_ClmIdx <- (X_abund$herb_ClmIdx - mean(X_abund$herb_ClmIdx, na.rm = TRUE)) / sd(X_abund$herb_ClmIdx, na.rm = TRUE)
X_abund$woody_ShpInx <- (X_abund$woody_ShpInx - mean(X_abund$woody_ShpInx, na.rm = TRUE)) / sd(X_abund$woody_ShpInx, na.rm = TRUE)
X_abund$herb_ShpInx <- (X_abund$herb_ShpInx - mean(X_abund$herb_ShpInx, na.rm = TRUE)) / sd(X_abund$herb_ShpInx, na.rm = TRUE)
X_abund$woody_lrgPInx <- (X_abund$woody_lrgPInx - mean(X_abund$woody_lrgPInx, na.rm = TRUE)) / sd(X_abund$woody_lrgPInx, na.rm = TRUE)
X_abund$herb_lrgPInx <- (X_abund$herb_lrgPInx - mean(X_abund$herb_lrgPInx, na.rm = TRUE)) / sd(X_abund$herb_lrgPInx, na.rm = TRUE)
X_abund$woody_AggInx <- (X_abund$woody_AggInx - mean(X_abund$woody_AggInx, na.rm = TRUE)) / sd(X_abund$woody_AggInx, na.rm = TRUE)
X_abund$herb_AggInx <- (X_abund$herb_AggInx - mean(X_abund$herb_AggInx, na.rm = TRUE)) / sd(X_abund$herb_AggInx, na.rm = TRUE)
X_abund$woody_EdgDens <- (X_abund$woody_EdgDens - mean(X_abund$woody_EdgDens, na.rm = TRUE)) / sd(X_abund$woody_EdgDens, na.rm = TRUE)
X_abund$herb_EdgDens <- (X_abund$herb_EdgDens - mean(X_abund$herb_EdgDens, na.rm = TRUE)) / sd(X_abund$herb_EdgDens, na.rm = TRUE)
X_abund$woody_Pdens <- (X_abund$woody_Pdens - mean(X_abund$woody_Pdens, na.rm = TRUE)) / sd(X_abund$woody_Pdens, na.rm = TRUE)
X_abund$herb_Pdens <- (X_abund$herb_Pdens - mean(X_abund$herb_Pdens, na.rm = TRUE)) / sd(X_abund$herb_Pdens, na.rm = TRUE)
X_abund$woody_Npatches <- (X_abund$woody_Npatches - mean(X_abund$woody_Npatches, na.rm = TRUE)) / sd(X_abund$woody_Npatches, na.rm = TRUE)
X_abund$herb_Npatches <- (X_abund$herb_Npatches - mean(X_abund$herb_Npatches, na.rm = TRUE)) / sd(X_abund$herb_Npatches, na.rm = TRUE)
X_abund$woody_mnFocal30m <- (X_abund$woody_mnFocal30m - mean(X_abund$woody_mnFocal30m, na.rm = TRUE)) / sd(X_abund$woody_mnFocal30m, na.rm = TRUE)
X_abund$vegDens50m <- (X_abund$vegDens50m - mean(X_abund$vegDens50m, na.rm = TRUE)) / sd(X_abund$vegDens50m, na.rm = TRUE)
X_abund <- as.matrix(X_abund)
print(X_abund)



## Extract and scale detection covariates to matrix ## 
temp_mat <- scale(weather_dat$Temp_degF)
wind_mat <- scale(weather_dat$Wind_mph)
sky_mat <- as.integer(as.factor(weather_dat$Sky_Condition))

# Making a day of year matrix
doy_vec <- yday(as.Date(paste0(formatted_dates, "_2024"), format = "%b_%d_%Y"))
doy_vec <-  as.integer(as.factor(doy_vec))

# Combine into a single dataframe where each row corresponds to a survey
X_det <- as.matrix(data.frame(temp = temp_mat, 
                              wind = wind_mat,
                              sky = sky_mat,
                              doy = doy_vec))

# Area surveyed 
# area <- pi * (200^2) / 4046.86  # in acres
# area <- pi * (200^2) / 10000  # in hectares
# Offset <- rep(area, 27)

# For random effect
n_days = max(J)

# ----------------------
# Bayesian P-value
# ----------------------

S_A <- sum(J_r > 0)
sites_av <- which(J_r > 0)
J_A <- max(J)


# ----------------------
# Bundle Data 
# ----------------------
Wolfe14.data <- list(v = v, 
                     y = y,
                     n = n,
                     k = k, 
                     X_abund = X_abund,
                     X_det = X_det
                     # Offset = area
)


# Check structure
str(Wolfe14.data)

# ----------------------
# Constants 
# ----------------------

constants<- list(S = S, 
                 A_times = A_times,
                 val_times = val_times, 
                 sites_a = sites_a, 
                 sites_av = sites_av,
                 S_val = S_val, 
                 J_val = J_val, 
                 S_A = S_A,
                 J_A = J_A,
                 J_r = J_r,
                 n_days = n_days
                 # days = days,
                 # n.doy = length(unique(doy_vec)),
                 # Sky_Lvls = length(unique(sky_mat))
)

# Check structure
str(constants)
 

# ---------------------------------------------------------- 
# 
#           Acoustic HM with Validated Calls
# 
# ----------------------------------------------------------

# ----------------------
# MCMC Specifications
# ----------------------
niter = 200000
nburnin =  80000
nchains = 3 
nthin = 10

# Posterior samples
(((niter - nburnin) / nthin) * nchains)

# ----------------------
# Model Specifications
# ----------------------

# Parameters monitored
params <- c('lambda', # Abundance
            'N_tot',
            'N',
            'beta0',
            'beta1',
            'beta2',
            'beta3',
            'nu',
            'p_a_tot',# Detection 
            'alpha0', 
            'alpha1', 
            'alpha2',
            'delta_tot', # Vocalization
            # 'gamma0', 
            'mu_j',
            'tau_j',
            'jRE',
            'omega',
            # 'phi',
            # 'mu_phi',
            # 'tau_phi',
            'fit_y',# Posterior Predictive Checks
            'fit_y_pred',
            'fit_v',
            'fit_v_pred',
            'bp_y', # Bayes p-value
            'bp_v'
) # End Params

 

# Initial Values 
inits <- list(
  # Abundance
  N = rep(1, S),
  nu = 2, # >1 is underdispersed        
  beta0 = rnorm(1, 0, 1),
  beta1 = rnorm(1, 0, 1),
  beta2 = rnorm(1, 0, 1),
  beta3 = rnorm(1, 0, 1),
  
  # Detection
  alpha0 = 0,            
  mu_alpha = 0.5,  
  alpha1 = runif(1, 0, 1), 
  alpha2 = rnorm(1, 0, 1),
  alpha3 = rnorm(1, 0, 1),
  
  # Vocalization
  omega = runif(1, 0, 1),  
  gamma0 = log(10),  
  
  # Survey random effect
  jRE = rep(0, n_days),
  mu_j = 1,
  tau_j = 1,
  
  # Overdispersion
  mu_phi = 1,  
  tau_phi = 1,
  phi = matrix(1, nrow = S, ncol = J_A),
  
  # False Positive
  K = k # initializing with observed true vocalizations

)


# ----------------------------- 
# Model Statement 
# ----------------------------- 
acoustic_model <- nimbleCode({
  
  # ----------------------
  # Abundance Priors
  # ----------------------
  
  # Intercept
  beta0 ~ dnorm(0, 10) 
  
  # Covariate effect
  beta1 ~ dnorm(0, 10) # Herbaceous Clumpy Index 
  beta2 ~ dnorm(0, 10) # Woody Aggregation Index 
  beta3 ~ dnorm(0, 10) # Interactive effect
  
  # Underdispersion
  nu ~ dunif(1, 3)
  

  # ------------------------
  # Detection Priors
  # ------------------------
  
  # Intercept
  alpha0 <- logit(mu_alpha) # Constrains alpha0 to be between 0 and 1 on the logit scale (propability)
  mu_alpha ~ dunif(0, 1)
  
  # True individuals
  alpha1 ~ dunif(0, 1000) # Constrained to be positive
  
  # Covariate effect
  alpha2 ~ dnorm(0, 10) # Vegetation Density
  alpha3 ~ dnorm(0, 10) # Wind
  
  # ------------------------
  # Call Rate Priors
  # ------------------------
  
  # False positive rate
  omega  ~ dunif(0, 1000) # From Doser et al.  
  
  # Intercept
  gamma0 ~ dunif(log(1), log(60))
  
  # Survey random effect - Centered
  mu_j ~ dgamma(0.01, 0.01)
  tau_j ~ dgamma(0.01, 0.01)
  for (j in 1:n_days) {
    jRE[j] ~ dnorm(mu_j, tau_j)
  }
  
  # Overdispersion
  mu_phi ~ dgamma(0.01, 0.01)    
  tau_phi ~ dgamma(0.01, 0.01)   
  for (s in 1:S) {
    for (j in 1:J_A) {
      phi[s, j] ~ dgamma(mu_phi, tau_phi)
    }
  }
  
  # ------------------------
  # Likelihood and Process Model
  # ------------------------
  
  # Site
  for (s in 1:S) {
    
    # ---------------------------------
    # Abundance Submodel  
    # ---------------------------------
    
    # # Poisson
    # # Intercept + Herbaceous Clumpy Indx + Woody Aggregation Indx + Interactive effect
    # log(lambda[s]) <- beta0 + beta1 * X_abund[s, 7] + beta2 * X_abund[s, 12] + beta3 * X_abund[s, 7] * X_abund[s, 12]
    # N[s] ~ dpois(lambda[s])
    
    # Conway-Maxwell Poisson
    # Intercept + Herbaceous Clumpy Indx + Woody Aggregation Indx  
    log(lambda[s]) <- beta0 + beta1 * X_abund[s, 7] + beta2 * X_abund[s, 12]
    N[s] ~ dCOMPois_nimble(lambda[s], nu)
    
    # Survey
    for (j in 1:J_A) {
      
    # ---------------------------------
    # Detection Submodel  
    # ---------------------------------
      
    # True Positives + Vegetation Density + Wind
    logit(p_a[s, j]) <- alpha0 + alpha1 * N[s] + alpha2 * X_abund[s, 21] + alpha3 * X_det[j, 2]
      
    # ---------------------------------
    # Call rate Submodel  
    # ---------------------------------
      
    # Intercept + Survey Random Effect
    log(delta[s, j]) <- jRE[j]
      
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
    
    # ---------------------------------
    # Vocalizations  
    # ---------------------------------

    # Surveys with Vocalizations
    for (j in 1:J_r[s]) {

      # Zero Truncated Negative Binomial
      v[s, A_times[s, j]] ~ T(dpois((delta[s, A_times[s, j]] * N[s] + omega) * phi[s, A_times[s, j]] * y[s, A_times[s, j]]), 1, )


      # ---------------------------------
      # PPC calls
      # ---------------------------------
      v_pred[s, j] ~ T(dpois((delta[s, A_times[s, j]] * N[s] + omega) * phi[s, A_times[s, j]] * y[s, A_times[s, j]]), 1, )
      mu_v[s, j] <- ((delta[s, A_times[s, j]] * N[s] + omega) * phi[s, A_times[s, j]]) / (1 - exp(-1 * ((delta[s, A_times[s, j]] * N[s] + omega) * phi[s, A_times[s, j]])))
      resid_v[s, j] <- pow(pow(v[s, A_times[s, j]], 0.5) - pow(mu_v[s, j], 0.5), 2)
      resid_v_pred[s, j] <- pow(pow(v_pred[s, j], 0.5) - pow(mu_v[s, j], 0.5), 2)


      }# End J_R
    
  } # End S
  
  # ------------------------
  # Manual Validation
  # ------------------------

  for (s in 1:S_val) {
    for (j in 1:J_val[s]) {
      K[s, j] ~ dbin(tp[sites_a[s], j], v[sites_a[s], val_times[s, j]])
      k[s, val_times[s, j]] ~ dhyper_nimble(K[s, j], v[sites_a[s], val_times[s, j]] - K[s, j], n[s, val_times[s, j]])
    } # End J
  } # End S

  # -------------------------------------------
  # PPC and Bayesian P-value
  # -------------------------------------------
  for (s in 1:S_A) {
    tmp_v[s] <- sum(resid_v[sites_av[s], 1:J_r[sites_av[s]]])
    tmp_v_pred[s] <- sum(resid_v_pred[sites_av[s], 1:J_r[sites_av[s]]])

  }
  
  for (s in 1:S) {
    tmp_y[s] <- sum(resid_y[sites_a[s], 1:J_A])
    tmp_y_pred[s] <- sum(resid_y_pred[sites_a[s], 1:J_A])
    
  }
  
  fit_y <- sum(tmp_y[1:S])
  fit_y_pred <-sum(tmp_y_pred[1:S])
  fit_v <- sum(tmp_v[1:S_A])
  fit_v_pred <- sum(tmp_v_pred[1:S_A])
  bp_y <- step(fit_y_pred - fit_y)
  bp_v <- step(fit_v_pred - fit_v)
  
  
  
  # -------------------------------------------
  # Derive Parameters
  # -------------------------------------------
  
  # Abundance
  N_tot <- sum(N[1:S])
  
  # Detection 
  p_a_tot <- sum(p_a[1:S, 1:J_A])
  
  # Call Rate
  delta_tot <- sum(delta[1:S, 1:J_A])
  


})
# ---------------------------- End Model ----------------------------

 

# ------------------------
# Fit Model 
# ------------------------

# Fit model with standard NIMBLE MCMC sampler
fm1 <- nimbleMCMC(code = acoustic_model,
                  data = Wolfe14.data,
                  constants = constants,
                  inits = inits,
                  monitors = params,
                  niter = niter,
                  nburnin = nburnin,
                  nchains = nchains,
                  thin = nthin,
                  progressBar = getNimbleOption("MCMCprogressBar"),
                  samplesAsCodaMCMC = TRUE,
                  summary = TRUE)
fm1$summary
# -------------------------------------------------------
# Check Convergence
# -------------------------------------------------------

# Trace plots
mcmcplots::mcmcplot(fm1$samples, 
                    parms = params) 

# Rhat
coda::gelman.diag(fm1$samples)
check_rhat(coda::gelman.diag(fm1$samples), threshold = 1.05) 

# -------------------------------------------------------
# Combine Chains for Posterior inference
# -------------------------------------------------------

# Combine chains
combined_chains <- as.mcmc(do.call(rbind, fm1$samples))

# -------------------------------------------------------
# Posterior Predictive Checks
# -------------------------------------------------------

# ----------------------
# Bayes P-value
# ----------------------

# Extracting samples and calculating mean
bp_y_samples <- combined_chains[, "bp_y"]
bp_v_samples <- combined_chains[, "bp_v"]
mn_bp_y <- mean(bp_y_samples)
mn_bp_v <- mean(bp_v_samples)

# P-value = 0.5 means good fit, = 1 or 0 is a poor fit

# Abundance
cat("Abundance Model Bayesian p-value =", mn_bp_y, "\n")

# Call    
cat("Call Model Bayesian p-value =", mn_bp_v, "\n") 

# ----------------------
# Extract Fits
# ----------------------

# Abundance
fit_y_data <- data.frame(
  Observed = as.vector(combined_chains[, "fit_y"]),       # Observed values
  Predicted = as.vector(combined_chains[, "fit_y_pred"]) # Predicted values
)

# Calls
fit_v_data <- data.frame(
  Observed = as.vector(combined_chains[, "fit_v"]),
  Predicted = as.vector(combined_chains[, "fit_v_pred"])
)

# ----------------------
# Density Plot
# ----------------------

# Abundance
y_PPC_Dens <- ggplot(fit_y_data) +
  geom_density(aes(x = Observed, fill = "Observed"), alpha = 0.5) +   
  geom_density(aes(x = Predicted, fill = "Predicted"), alpha = 0.5) +  
  scale_fill_manual(values = c("Observed" = "blue", "Predicted" = "red")) +  
  labs(title = "Posterior Predictive Check for Abundance", 
       x = "Fit Values", 
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank()) 

# View
print(y_PPC_Dens)

# Export                
ggsave(plot = y_PPC_Dens, "./Figures/PPC/AV_Wolfe_Abund_Density.jpeg", width = 8, height = 5, dpi = 300)
dev.off()


# Call
v_PPC_Dens <- ggplot(fit_v_data) +
  geom_density(aes(x = Observed, fill = "Observed"), alpha = 0.5) +   
  geom_density(aes(x = Predicted, fill = "Predicted"), alpha = 0.5) +  
  scale_fill_manual(values = c("Observed" = "blue", "Predicted" = "red")) +  
  labs(title = "Posterior Predictive Check for Call", 
       x = "Fit Values", 
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank()) 

# View
print(v_PPC_Dens)

# Export                
ggsave(plot = v_PPC_Dens, "./Figures/PPC/AV_Wolfe_Call_Density.jpeg", width = 8, height = 5, dpi = 300)
dev.off()


# -------------------------------------------------------
# Beta Estimates **** interactive effect
# -------------------------------------------------------

# Extract beta estimates
beta0_samples <- combined_chains[, "beta0"]
beta1_samples <- combined_chains[, "beta1"]
beta2_samples <- combined_chains[, "beta2"]

# Means
beta0 <- mean(beta0_samples)
beta1 <- mean(beta1_samples)
beta2 <- mean(beta2_samples)

# Compute 95% CI for each beta
beta_df <- data.frame(
  value = c(beta0_samples, beta1_samples, beta2_samples),  
  parameter = rep(c("beta0", "beta1", "beta2"), each = length(beta0_samples))
) %>%
  group_by(parameter) %>%
  filter(value >= quantile(value, 0.025) & value <= quantile(value, 0.975))  

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

# Clear
dev.off()


# Export beta dataframe
saveRDS(beta_df, "./Data/Model_Data/ARU_WolfeAV_beta_df.rds")

# -------------------------------------------------------
# Covariate Effects
# -------------------------------------------------------

# Covariate names
print(colnames(X_abund))

# Set covariate name 
Cov1_name <- "herb_ClmIdx"
Cov2_name <- "woody_AggInx"

# Create a prediction of covariate values
cov1_scaled <- (X_abund[, Cov1_name] - mean(X_abund[, Cov1_name])) / (max(X_abund[, Cov1_name]) - min(X_abund[, Cov1_name]))
cov2_scaled <- (X_abund[, Cov2_name] - mean(X_abund[, Cov2_name])) / (max(X_abund[, Cov2_name]) - min(X_abund[, Cov2_name]))

# Matrices for storing predictions
cov1_preds <- matrix(NA, nrow = length(beta0_samples), ncol = length(cov1_scaled))
cov2_preds <- matrix(NA, nrow = length(beta0_samples), ncol = length(cov2_scaled))

# Generate predictions
for (i in 1:length(beta0_samples)) {
  cov1_preds[i, ] <- beta0_samples[i] + beta1_samples[i] * cov1_scaled  
  cov2_preds[i, ] <- beta0_samples[i] + beta2_samples[i] * cov2_scaled
}# end loop


# Calculate mean predictions
cov1_preds_mean <- apply(cov1_preds, 2, mean)# mean
cov1_preds_LCI <- apply(cov1_preds, 2, quantile, probs = 0.025) # LCI
cov1_preds_HCI <- apply(cov1_preds, 2, quantile, probs = 0.975) # HCI

cov2_preds_mean <- apply(cov2_preds, 2, mean)
cov2_preds_LCI <- apply(cov2_preds, 2, quantile, probs = 0.025)
cov2_preds_HCI <- apply(cov2_preds, 2, quantile, probs = 0.975)

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

# Clear
dev.off()

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

# Clear
dev.off()

# -------------------------------------------------------
#  Estimating Abundance 
# -------------------------------------------------------

# Extract abundance posterior
Ntot_samples <- combined_chains[ ,"N_tot"]

# Ntotal is the abundance based on 27 acoustic sites at a radius of 200m.
# To correct for density, Ntotal needs to be divided by 27 * area surveyed
area <- pi * (200^2) / 4046.86  # Area in acres
dens_samples <- Ntot_samples / (area * 27)

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
dens_df <- dens_df[dens_df$Density >= dens_summary$Lower_CI & dens_df$Density <= dens_summary$Upper_CI, ]

# Getting total abundance
abund_summary <- dens_summary
abund_summary[,2:4] <- abund_summary[,2:4] * 2710

# Plot Abundance - Violin
abund_df <- dens_df
abund_df$Density <- abund_df$Density * 2710

ggplot(abund_df, aes(x = Model, y = Density, fill = Model)) + 
  geom_violin(trim = FALSE, alpha = 0.6, adjust = 5) +   
  labs(x = "Model", y = "Total Abundance") +
  scale_fill_manual(values = c("AV Wolfe" = "red")) +   
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

# Clear
dev.off()

# Export density dataframe
saveRDS(dens_df, "./Data/Model_Data/ARU_WolfeAV__dens_df.rds")
saveRDS(dens_summary, "./Data/Model_Data/ARU_WolfeAV_dens_summary.rds")
saveRDS(abund_summary, "./Data/Model_Data/ARU_WolfeAV_abund_summary.rds")

# End Script