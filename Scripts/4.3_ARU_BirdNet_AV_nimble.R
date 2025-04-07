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
model_name <- "AV Bnet"

# -------------------------------------------------------
#
#             Variable and Object Definitions
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
# p.a = prob of detecting at least one vocalization in an acoustic recording
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
bnet_dat_all <- read.csv("./Data/Acoustic_Data/NOBO_BirdNETall_2024.csv")

# ARU weather coreadRDS()# ARU weather covariates
weather_dat <- read.csv("./Data/Acoustic_Data/ARU_weathercovs.csv")

# Site covariates
site_covs <- read.csv("./Data/Acoustic_Data/ARU_siteCovs.csv")


# -------------------------------------------------------
#
#    Data Wrangling
#
# -------------------------------------------------------

# Subsetting to a 21 day period
start_date <- as.Date("2024-06-18")
end_date <- start_date + 20  # 21 days total
bnet_dat <- bnet_dat_all[bnet_dat_all$Date >= start_date & bnet_dat_all$Date <= end_date, ]
head(bnet_dat)
nrow(bnet_dat)

# Subsetting to the last 10 mins of the recording
# bnet_dat$Time <- as.POSIXct(bnet_dat$Time, format = "%H:%M:%S")
# start_time <- as.POSIXct("06:45:00", format = "%H:%M:%S")
# end_time   <- as.POSIXct("06:55:00", format = "%H:%M:%S")
# bnet_dat <- bnet_dat[bnet_dat$Time >= start_time & bnet_dat$Time <= end_time, ]
# head(bnet_dat)
# min(bnet_dat$Time)
# max(bnet_dat$Time)
# nrow(bnet_dat)
 

# Creating an occasion column using dates
date_order <- seq(from = start_date, to = end_date, by = "day")
bnet_dat <- bnet_dat %>%  # Adding occasion column
  mutate(Occasion = match(Date, date_order))

# Adding a row counter
bnet_dat$Count <- 1

# Initialize a site by survey matrix
v <- matrix(0, nrow = 27, ncol = 21)        

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
sum(v) # matches nrow(bnet_dat)

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

# Init for K in validation model
K_init <- v[sites_a, ]


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
J <- rep(21, S) # can vary by site

# The max number of surveys
J_A <- max(J)

# J_r contains the number of surveys at each acoustic site that contains at least 1 detected vocalization. 
J_r <- apply(v_mat, 1, function(a) {sum(a != 0)})
J_r <- ifelse(is.na(J_r), 0, J_r)
J_r <- as.numeric(J_r)
print(J_r)

# Sites with at least one acoustic recording
sites_av <- which(J_r > 0)

# Number of sites with at least one acoustic recording
S_A <- sum(J_r > 0)
 
# changing 0 to NA for masking ----------------------------------------------------
# J_r[J_r == 0] <- NA
# 
# J_r_max <- max(J_r, na.rm = TRUE)
# 
# run_j <- matrix(0, nrow = length(J_r), ncol = J_r_max)
# for (s in seq_along(J_r)) {
#   if (!is.na(J_r[s])) {
#     run_j[s, 1:J_r[s]] <- 1
#   }
# }



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
rownames(A_times) <- rownames(v) 
print(A_times)
A_times <- A_times[!apply(is.na(A_times), 1, all), ]

 

# ---------------------------------------
# Manually validated  
# ---------------------------------------

## This is from 14 day every 4 day data set

# # Validated Calls
# # Do not include sites with no calls
# # only include occasions where at least 1 call was validated for a site
# n <- read.csv("./Data/Acoustic_Data/Bnet14day_n.csv", row.names = 1)
# n <- as.matrix(n)
# 
# # True Calls
# # Calls found to be true, same dimension as n
# k <- read.csv("./Data/Acoustic_Data/Bnet14day_k.csv", row.names = 1)
# k <- as.matrix(k)
# 
# # Survey days calls were validated, same dimension as n
# val_times <- read.csv("./Data/Acoustic_Data/Bnet14day_val_times.csv", row.names = 1)
# val_times <- as.matrix(val_times)
# 
# # Total number of sites with manually validated data
# S_val <- nrow(n)
# 
# # How many surveys were validate
# J_val <- apply(!is.na(val_times), 1, sum)



# #### This is fake data ######
# write.csv(v, "bnet_v.csv")

# Validated Calls
# Do not include sites with no calls
# only include occasions where at least 1 call was validated for a site
n <- read.csv("./bnet_n_21days30mins.csv", row.names = 1)
n <- as.matrix(n)
n

# True Calls
# Calls found to be true, same dimension as n
k <- read.csv("./bnet_k_21days30mins.csv", row.names = 1)
k <- as.matrix(k)
k

# # Survey days calls were validated, same dimension as n
val_times <- read.csv("./bnet_valtimes_21days30mins.csv", row.names = 1)
val_times <- as.matrix(val_times)


# Check dimensions
dim(n) 
dim(k)
dim(val_times)


# Total number of sites with manually validated data
S_val <- nrow(n)

# How many surveys were validate
J_val <- apply(!is.na(val_times), 1, sum)
J_val <- J_val[J_val != 0] # removing sites with no validations

# Check dimensions
S_val
length(J_val)


# ---------------------------------------
# Covariates 
# ---------------------------------------

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



# ## Extract and scale detection covariates to matrix ## 
# temp_mat <- scale(weather_dat$Temp_degF)
# wind_mat <- scale(weather_dat$Wind_mph)
# sky_mat <- as.integer(as.factor(weather_dat$Sky_Condition))
# 
# # Making a day of year matrix
# doy_vec <- yday(as.Date(paste0(formatted_dates, "_2024"), format = "%b_%d_%Y"))
# doy_vec <-  as.integer(as.factor(doy_vec))
# 
# # Combine into a single dataframe where each row corresponds to a survey
# X_det <- as.matrix(data.frame(temp = temp_mat, 
#                               wind = wind_mat,
#                               sky = sky_mat,
#                               doy = doy_vec))

# # Area surveyed 
# #area <- pi * (200^2) / 4046.86  # in acres
# area <- pi * (200^2) / 10000  # in hectares
# Offset <- rep(area, 27)

# For random effect
n_days = max(J)
 

# ----------------------
# Bundle Data 
# ----------------------
data <- list(v = v, 
             y = y,
             n = n,
             k = k, 
             X_abund = X_abund
             # X_det = X_det
             # Offset = area
)


# Check structure
str(data)

# ----------------------
# Constants 
# ----------------------

constants<- list(S = S, 
                 # J = J,
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
niter = 1000
nburnin =  0
nchains = 3 
nthin = 1

# ----------------------
# Model Specifications
# ----------------------

# Parameters monitored
params <- c(
            # Abundance
            'lambda', 
            'N_tot',
            'N',
            'beta0',
            'beta1',
            'beta2',
            # 'nu',
            
            # Detection 
            #'p_a',
            'alpha0', 
            'alpha1', 
            'alpha2',
            
            # Vocalization
            #'delta', 
            'tau_j',
            'jRE',
            'omega',
            'phi',
            'mu_phi',
            'tau_phi'
            
            # Posterior Predictive Checks
            # 'fit_y',
            # 'fit_y_pred',
            # 'fit_v',
            # 'fit_v_pred',
            # 'bp_y', # Bayesian p-value
            # 'bp_v'
  
) # End Params



# Initial Values 
inits <- list(
  
      # Abundance
      N = rep(1, S),
      # nu = rep(10, S), # <1 = overdispersed, 1 = Poisson,  >1 = underdispersed        
      beta0 = rnorm(1, 0, 1),
      beta1 = rnorm(1, 0, 1),
      beta2 = rnorm(1, 0, 1),

      # Detection
      alpha0 = 0,            
      mu_alpha = 0.5,  
      alpha1 = runif(1, 0, 1), 
      alpha2 = rnorm(1, 0, 1),
      alpha3 = rnorm(1, 0, 1),
      
      # Vocalization
      omega = 0,  
      # gamma0 = log(10),  
      
      # Survey random effect
      jRE = rep(0, n_days),
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

  # Underdispersion
  # for (s in 1:S) {
  #   #nu[s] ~ T(dgamma(3, 0.5), 1, )
  #   nu[s] ~ T(dgamma(2, 0.002), 1, )
  # }
  
  
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
  
  # Survey random effect  
  tau_j ~ dgamma(0.01, 0.01)
  for (j in 1:n_days) {
    jRE[j] ~ dnorm(0, tau_j)
  }
  
  # Overdispersion
  mu_phi ~ dgamma(0.01, 0.01)
  tau_phi ~ dgamma(0.01, 0.01)
  for (s in 1:S) {
    for (j in 1:J_A) {
      phi[s, j] ~ dgamma(mu_phi, tau_phi)
    }
  }
  
 
  
  # ----------------------------------------
  #
  #     Likelihood and Process Model
  #
  # ----------------------------------------
  
  # Site
  for (s in 1:S) {
    
    # ---------------------------------
    # Abundance   
    # ---------------------------------
    
    # Poisson
    # Intercept + Herbaceous Clumpy Indx + Woody Aggregation Indx
    log(lambda[s]) <- beta0 + beta1 * X_abund[s, 7] + beta2 * X_abund[s, 12]
    N[s] ~ dpois(lambda[s])
    
    # Conway-Maxwell Poisson
    # Intercept + Herbaceous Clumpy Indx + Woody Aggregation Indx  
    # log(lambda[s]) <- beta0 + beta1 * X_abund[s, 7] + beta2 * X_abund[s, 12] 
    # N[s] ~ dCOMPois_nimble(lambda[s], nu[s])
    
    # Survey
    for (j in 1:J_A) {
      
      # ---------------------------------
      # Detection   
      # ---------------------------------
      
      # True Positives + Vegetation Density + Wind
      logit(p_a[s, j]) <- alpha0 + alpha1 * N[s] + alpha2 * X_abund[s, 21] #+ alpha3 * X_det[j, 2]
      
      # Observations
      y[s, j] ~ dbin(p_a[s, j], 1)
      
      # Expected and residuals
      y_pred[s, j] ~ dbin(p_a[s, j], 1)
      resid_y[s, j] <- pow(pow(y[s, j], 0.5) - pow(p_a[s, j], 0.5), 2)
      resid_y_pred[s, j] <- pow(pow(y_pred[s, j], 0.5) - pow(p_a[s, j], 0.5), 2)
      
      # ---------------------------------
      # Call rate   
      # ---------------------------------
      
      # Intercept + Survey Random Effect
      log(delta[s, j]) <- jRE[j]

      # ---------------------------------
      # True Positives 
      # ---------------------------------
      tp[s, j] <- delta[s, j] * N[s] / (delta[s, j] * N[s] + omega)
      
      # ---------------------------------
      # Vocalizations  
      # ---------------------------------
  
      # Zero Truncated Negative Binomial - Implementation as seen in Doser et al. 2021
      v[s,  j] ~ dpois((delta[s, j] * N[s] + omega) * phi[s, j] * y[s, j]) 
      
      # Expected and residuals
      v_pred[s, j] ~ dpois((delta[s, j] * N[s] + omega) * phi[s, j] * y[s, j]) 
      mu_v[s, j] <- ((delta[s, j] * N[s] + omega) * phi[s, j]) / (1 - exp(-1 * ((delta[s, j] * N[s] + omega) * phi[s, j])))
      resid_v[s, j] <- pow(pow(v[s, j], 0.5) - pow(mu_v[s, j], 0.5), 2)
      resid_v_pred[s, j] <- pow(pow(v_pred[s, j], 0.5) - pow(mu_v[s, j], 0.5), 2)

    } # End J
  } # End S

  # ------------------------
  # Manual Validation
  # ------------------------
  
  for (s in 1:S_val) {
    for (j in 1:J_val[s]) {
      K[s, j] ~ dbin(tp[sites_a[s], val_times[s, j]], v[sites_a[s], val_times[s, j]])
      k[s, val_times[s, j]] ~ dhyper_nimble(K[s, j], v[sites_a[s], val_times[s, j]] - K[s, j], n[s, val_times[s, j]])
    } # End J_val
  } # End S_val

  
  # -------------------------------------------
  # PPC and Bayesian P-value
  # -------------------------------------------
  # for (s in 1:S_A) {
  #   tmp_v[s] <- sum(resid_v[sites_av[s], 1:J_r[sites_av[s]]])
  #   tmp_v_pred[s] <- sum(resid_v_pred[sites_av[s], 1:J_r[sites_av[s]]])
  #   
  # }
  # 
  # for (s in 1:S) {
  #   tmp_y[s] <- sum(resid_y[sites_a[s], 1:J_A])
  #   tmp_y_pred[s] <- sum(resid_y_pred[sites_a[s], 1:J_A])
  #   
  # }
  # 
  # fit_y <- sum(tmp_y[1:S])
  # fit_y_pred <-sum(tmp_y_pred[1:S])
  # fit_v <- sum(tmp_v[1:S_A])
  # fit_v_pred <- sum(tmp_v_pred[1:S_A])
  # bp_y <- step(fit_y_pred - fit_y)
  # bp_v <- step(fit_v_pred - fit_v)
 

  
  # -------------------------------------------
  # Derive Parameters
  # -------------------------------------------
  
  # Abundance
  N_tot <- sum(N[1:S])
  
})
# ---------------------------- End Model ----------------------------



# ------------------------
# Fit Model 
# ------------------------

# Fit model with standard NIMBLE MCMC sampler
fm1 <- nimbleMCMC(code = acoustic_model,
                  data = data,
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

# Export model
# saveRDS(fm1, "./Data/Model_Data/ModelFits_AV-BirdNET_fm1.rds")

# fm1 <- readRDS("./Data/Model_Data/ModelFits_AV-BirdNET_fm1.rds")

# -------------------------------------------------------
# Check Convergence
# -------------------------------------------------------

# Extract posterior samples
fm1_samples <- as.mcmc.list(fm1$samples)

# Trace plots
mcmcplots::mcmcplot(fm1_samples,
                    parms = params)

# # Rhat
# coda::gelman.diag(fm1_samples, multivariate = FALSE)
# 
# # Effective sample size
# as.data.frame(coda::effectiveSize(fm1_samples))
# 
# # Summary
# fm1_summary <- summary(fm1_samples)
# View(fm1_summary$statistics)


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
ggsave(plot = y_PPC_Dens, "./Figures/PPC/Abund_Density_AV-BirdNET.jpeg", width = 8, height = 5, dpi = 300)
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
ggsave(plot = v_PPC_Dens, "./Figures/PPC/Call_Density_AV-BirdNET.jpeg", width = 8, height = 5, dpi = 300)
dev.off()
