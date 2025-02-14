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

# Set seed, scientific notation options, and working directory
set.seed(123)
options(scipen = 9999)
setwd(".")

# Setting up cores
Ncores <- parallel::detectCores()
print(Ncores) # Number of available cores
workers <- 2 # Ncores * 0.5 # For low background use 80%, for medium use 50% of Ncores

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
# n.days = number of recording days.
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
# days = variable used to index different recording days. 
# A.times = indexing variable used to determine specific indexes with v > 0.
# K = true number of acosutic vocalizations
# k = true number of manually validated acoustic vocalizations
# n = total number of manually validated acoustic vocalizations 
# sites.a = specific indices of sites where acoustic data were obtained
# R.val = number of validated sites for acoustic data
# Other variables not defined are for computation of Bayesian p-values. 

# -------------------------------------------------------
#
#                    Load Data
#
# -------------------------------------------------------

# BirdNet detections
bnet_dat_all <- read.csv("./Data/Acoustic_Data/NOBO_BirdNETall_2024.csv")

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
bnet_dat <- bnet_dat_all %>% filter(Date %in% date_order)
head(bnet_dat)
nrow(bnet_dat)

# Adding occasion column
bnet_dat <- bnet_dat %>%
  mutate(Occasion = match(Date, date_order))

# Adding a row count
bnet_dat$Count <- 1


# Initialize the matrix
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

# take a look
print(v)
                             
# ARU at site 24 stopped recording on 6/20. Making Columns 8:14 NA
v[24,8:14] <- NA
print(v)  

# Renaming columns to date Month_day
formatted_dates <- format(as.Date(date_order), "%b_%d")
colnames(v) <- formatted_dates
print(v) 

# Adding rownames
rownames(v) <- as.numeric(1:27)

# Get the site numbers with at least one call
v_mat <- as.matrix(v)
sites.a <- which(rowSums(!is.na(v_mat) & v >= 1) > 0)
print(sites.a)



# Creating a binary detection matrix, values >= 1 become 1, and 0 stays 0
v_df <- as.data.frame(v)
y <- v_df %>% mutate(across(where(is.numeric), ~ if_else(. >= 1, 1, 0)))
y <- as.matrix(y)
rownames(y) <- NULL  
colnames(y) <- NULL 
print(y)

# Total number of sites
R <- as.integer(27) 

# Number of repeat visits for each site
J <- rep(14, R)  

# J.r contains the number of surveys at each acoustic site that contains at least 1 detected vocalization. 
J.r <- apply(v_mat, 1, function(a) {sum(a != 0)})
J.r <- ifelse(is.na(J.r), 0, J.r)
J.r <- as.numeric(J.r)
print(J.r)

# A.times is a R x J matrix that contains the indices of the
# specific surveys at each acoustic site that have at least one detected
# vocalization. This determines the specific components of v[i, j] that
# are used in the zero-truncated Poisson.
A.times <- matrix(NA, R, max(J))
tmp <- apply(v, 1, function(a) {which(a != 0)})
for (i in 1:R) {
  if (length(tmp[[i]]) > 0) {
    A.times[i, 1:length(tmp[[i]])] <- tmp[[i]]
  }
}
print(A.times)

# ---------------------------------------
# Manually validated  
# ---------------------------------------


# Validated Calls
# Do not include sites with no calls 
# only include occasions where at least 1 call was validated for a site
n <- read.csv("./Data/Acoustic_Data/Bnet14day_n.csv", row.names = 1)
n <- as.matrix(n)

# True Calls
# Calls found to be true, same dimension as n
k <- read.csv("./Data/Acoustic_Data/Bnet14day_k.csv", row.names = 1)
k <- as.matrix(k)

# Survey days calls were validated, same dimension as n
val.times <- read.csv("./Data/Acoustic_Data/Bnet14day_val.times.csv", row.names = 1)
val.times <- as.matrix(val.times)

# Total number of sites with manually validated data
R.val <- nrow(n)

# How many surveys were validate
J.val <- rep(14, R.val)  
 

# Check
dim(n) # dimensions
dim(k)
dim(val.times)
print(n) # format
print(k)
print(val.times)


# ---------------------------------------
# Covariates 
# ---------------------------------------

# For day random effect on cue production rate
days <- matrix(rep(1:14, each = 27), nrow = 27, ncol = 14, byrow = TRUE)
print(days)

## Extract and scale site covariates for X.abund ##
X.abund <- site_covs[,-c(1:4)] 
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
X.abund <- as.matrix(X.abund)
print(X.abund)


 
## Extract and scale detection covariates to matrix ## 
temp_mat <- scale(weather_dat$Temp_degF)
wind_mat <- scale(weather_dat$Wind_mph)
X.det <- as.matrix(data.frame(temp = temp_mat, wind = wind_mat))
                   
# sky_vec <- as.factor(weather_dat$Sky_Condition)  # Convert to factor
# sky_dum <- model.matrix(~ sky_vec - 1)  # Create dummy variables for sky, without intercept

# Area in acres
area <- pi*(200^2)/4046.86

# ---------------------------------------
# Bayesian P-value
# ---------------------------------------

R.A <- sum(J.r > 0)
sites.a.v <- which(J.r > 0)
J.A <- max(J)


# ---------------------------------------
# Bundle Data 
# ---------------------------------------

Bnet14.data <- list(R = R, 
                    J = J, 
                    v = v, 
                    y = y,
                    n = n,
                    k = k, 
                    val.times = val.times, 
                    sites.a = sites.a, 
                    R.val = R.val, 
                    J.val = J.val, 
                    J.r = J.r, 
                    A.times = A.times, 
                    R.A = R.A, 
                    J.A = J.A, 
                    sites.a.v = sites.a.v, 
                    days = days,
                    n.days = max(J),
                    X.abund = X.abund,
                    X.det = X.det,
                    area = area)

# Check structure
str(Bnet14.data)


# ---------------------------------------------------------- 
# 
#           Acoustic HM with Validated Calls
# 
# ----------------------------------------------------------


# -------------------
# MCMC Specifications
# -------------------
n.iter = 300000
n.burnin = 60000 
n.chains = 3
n.thin = 10


# ----------------------------------------------------------
#                   Model 1 
# Cue Rate =  ran eff of survey/day
# Detection =  Wind
# Abundance = Herbaceous Patch Density + Woody large patch index
# ----------------------------------------------------------



# -------------------------------------------------------
# Model Specifications
# -------------------------------------------------------

# Parameters monitored
params <- c('alpha0', 
            'alpha1', 
            'alpha2',
             'beta0',
             'beta1',
             'tau', 
             'tau.day', 
             'a.phi', 
             'omega', 
             'bp.y', 
             'bp.v',
             'lambda',
             'N',
             'N_tot',
             'D_tot')

# Initial Values 
inits <- function() {
              list(
                N = rep(1, R), 
                beta0 = rnorm(1),
                beta1 = 0,
                omega = runif(1, 0, 10), 
                tau.day = runif(1, 0.1, 1),
                tau = runif(1, 0.1, 1),
                tau.p = runif(1, 0.1, 1),
                tau.day.p = runif(1, 0.1, 1),
                alpha1 = runif(1, 0, 1), 
                alpha2 = 0,
                a.phi = runif(1, 0, 5)
              )
}#end inits
  
  

# -------------------------------------------------------
# Model Statement
# -------------------------------------------------------
cat(" model {
  
  # -------------------------------------------
  # Priors 
  # -------------------------------------------
  beta0 ~ dnorm(0, 1)
  beta1 ~ dnorm(0, 1)
  alpha0 <- logit(mu.alpha)
  mu.alpha ~ dunif(0, 1)
  alpha1 ~ dunif(0, 1000) # Constrained to be positive
  alpha2 ~ dnorm(0, 1)
  omega ~ dunif(0, 1000)
  tau.day ~ dgamma(.01, .01)
  a.phi ~ dunif(0, 100)

  for (i in 1:n.days) {
    gamma.1[i] ~ dnorm(0, tau.day)
  }

  for (i in 1:R) {
    for (j in 1:J.A) {
      phi[i, j] ~ dgamma(a.phi, a.phi)
    }
  }
  
  # -------------------------------------------
  # Likelihood and process model 
  # -------------------------------------------
  for (i in 1:R) {
  
    # Abundance Model
    log(lambda[i]) <- beta0 + beta1 * X.abund[i,17] 
    N[i] ~ dpois(lambda[i])
    

    
    # Acoustic Data 
    for (j in 1:J[i]) {
    
    # Detection Model
    logit(p.a[i, j]) <- alpha0 + alpha1 * N[i] + alpha2 * X.abund[i,21]
    
    # Vocalization Model
      log(delta[i, j]) <- gamma.1[days[i, j]]
      y[i, j] ~ dbin(p.a[i, j], 1)
      tp[i, j] <- delta[i, j] * N[i] / (delta[i, j] * N[i] + omega)
      
      # Posterior predictive checks for Bayesian P-value
      y.pred[i, j] ~ dbin(p.a[i, j], 1)
      resid.y[i, j] <- pow(pow(y[i, j], 0.5) - pow(p.a[i, j], 0.5), 2)
      resid.y.pred[i, j] <- pow(pow(y.pred[i, j], 0.5) - pow(p.a[i, j], 0.5), 2)
    } # j
    
    for (j in 1:J.r[i]) {
      v[i, A.times[i, j]] ~ dpois((delta[i, A.times[i, j]] * N[i] + omega) * phi[i, A.times[i, j]] * y[i, A.times[i, j]]) T(1, )
      v.pred[i, j] ~ dpois((delta[i, A.times[i, j]] * N[i] + omega) * phi[i, A.times[i, j]] * y[i, A.times[i, j]]) T(1, )
      mu.v[i, j] <- ((delta[i, A.times[i, j]] * N[i] + omega) * phi[i, A.times[i, j]]) / (1 - exp(-1 * ((delta[i, A.times[i, j]] * N[i] + omega) * phi[i, A.times[i, j]])))
      resid.v[i, j] <- pow(pow(v[i, A.times[i, j]], 0.5) - pow(mu.v[i, j], 0.5), 2)
      resid.v.pred[i, j] <- pow(pow(v.pred[i, j], 0.5) - pow(mu.v[i, j], 0.5), 2)
    } # j
  } # i
  
  # ------------------------------------------- 
  # Manual validation 
  # -------------------------------------------
  for (i in 1:R.val) {
    for (j in 1:J.val[i]) {
      K[i, j] ~ dbin(tp[sites.a[i], j], v[sites.a[i], val.times[i, j]])
      k[i, val.times[i, j]] ~ dhyper(K[i, j], v[sites.a[i], val.times[i, j]] - K[i, j], n[i, val.times[i, j]], 1)

    } # j
  } # i
  
  # -------------------------------------------
  # Derive Abundance/Density
  # -------------------------------------------
  N_tot <- sum(N[])
  D_tot <- N_tot/area

  
  # -------------------------------------------
  # Bayesian P-value
  # -------------------------------------------
  for (i in 1:R.A) {
    tmp.v[i] <- sum(resid.v[sites.a.v[i], 1:J.r[sites.a.v[i]]])
    tmp.v.pred[i] <- sum(resid.v.pred[sites.a.v[i], 1:J.r[sites.a.v[i]]])
  }
  fit.y <- sum(resid.y[sites.a, 1:J.A])
  fit.y.pred <- sum(resid.y.pred[sites.a, 1:J.A])
  fit.v <- sum(tmp.v[1:R.A])
  fit.v.pred <- sum(tmp.v.pred[1:R.A])
  bp.y <- step(fit.y.pred - fit.y)
  bp.v <- step(fit.v.pred - fit.v)
}
", fill=TRUE, file="./jags_models/Bnet_AV_mod1.txt")
# ------------End Model-------------


# Fit Model
fm.6 <- jags(data = Bnet14.data, 
             inits = inits, 
             parameters.to.save = params, 
             model.file = "./jags_models/Bnet_AV_mod1.txt", 
             n.iter = n.iter,
             n.burnin = n.burnin,
             n.chains = n.chains, 
             n.thin = n.thin,
             parallel = TRUE,
             n.cores = workers,
             DIC = TRUE) 
 
# Check Convergence
mcmcplot(fm.6$samples)# Trace plots
cat("Bayesian p-value =", fm.6$summary["bp.v",1], "\n")# Bayesian P value
fm.6$Rhat# Rhat

# Model summary
print(fm.6, digits = 2)

# Average number of false positives detections
cat("False positives =", fm.6$summary["omega",1], "\n")


# Save Environment
save.image(file = "./ARU_AV_Bnet14day_JAGs.RData")


# -------------------------------------------------------
#
#   Estimating Abundance 
#
# -------------------------------------------------------

# Combine chains
combined_chains <- as.mcmc(do.call(rbind, fm.6$samples))


# Extracting Abundance
Ntot_samples <- combined_chains[ ,"N_tot"]

# Ntotal is the abundance based on 27 acoustic sites at a radius of 200m.
# To correct for density, Ntotal needs to be divided by 27 * area surveyed
area <- pi * (200^2) / 4046.86  # Area in acres
dens_samples <- Ntot_samples / (area * 27)

# Create data frame for density
dens_df <- data.frame(Model = rep("AV Bnet", length(dens_samples)), Density = dens_samples)
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
  scale_y_continuous(limits = c(100,1000),
                     breaks = seq(100, 1000, by = 100),
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
save.image(file = "./ARU_AV_Bnet14day_JAGs.RData")

# Export density dataframe
saveRDS(dens_summary, "./Data/Fitted_Models/ARU_BnetAV_dens_summary.rds")
saveRDS(abund_summary, "./Data/Fitted_Models/ARU_BnetAV_abund_summary.rds")



