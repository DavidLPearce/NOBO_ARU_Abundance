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
workers <- Ncores * 0.7 # For low background use 80%, for medium use 50% of Ncores
print(workers)

# -------------------------------------------------------
#
# Variable and Object Definitions
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
# K = true number of acoustic vocalizations
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


# Initialize the matrix
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

# take a look
print(v)



# Renaming columns to date Month_day
date_order <- c("2024-05-26", "2024-05-30", "2024-06-03", "2024-06-07", 
                "2024-06-11", "2024-06-15", "2024-06-19", "2024-06-23",
                "2024-06-27", "2024-07-01", "2024-07-05", "2024-07-09",
                "2024-07-13", "2024-07-17")
formatted_dates <- format(as.Date(date_order), "%b_%d")
colnames(v) <- formatted_dates
print(v) 
write.csv(v, "v_Wolfe.csv")
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

# ARU at site 24 stopped recording on 6/20. Making Columns 8:14 NA
# v[24,8:14] <- NA
# print(v)

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
n <- read.csv("./Data/Acoustic_Data/Wolfe14day_nTEST.csv", row.names = 1)
n <- as.matrix(n)

# True Calls
# Calls found to be true, same dimension as n
k <- read.csv("./Data/Acoustic_Data/Wolfe14day_kTEST.csv", row.names = 1)
k <- as.matrix(k)

# Survey days calls were validated, same dimension as n
val.times <- read.csv("./Data/Acoustic_Data/Wolfe14day_val.times.csv", row.names = 1)
val.times <- as.matrix(val.times)

# Total number of sites with manually validated data
R.val <- nrow(n)

# How many surveys were validate
J.val <- rep(14, R.val)  


# Check
dim(n) # dimensions
dim(k)
dim(val.times)

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

sky_vec <- as.factor(weather_dat$Sky_Condition)  # Convert to factor
sky_dum <- model.matrix(~ sky_vec - 1)  # Create dummy variables for sky, without intercept


# ---------------------------------------
# Bayesian P-value
# ---------------------------------------

R.A <- sum(J.r > 0)
sites.a.v <- which(J.r > 0)
J.A <- max(J)


# ---------------------------------------
# Bundle Data 
# ---------------------------------------

Wolfe14.data <- list(R = R, 
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
                     sky_dum = sky_dum)

# Check structure
str(Wolfe14.data)



# ---------------------------------------------------------- 
# 
#           Acoustic HM 
# 
# ----------------------------------------------------------

# ---------------------------------------
# MCMC Specifications
# ---------------------------------------
n.iter =  300000
n.burnin =  20000
n.chains = 3
n.thin = 10
n.adapt = 5000

# ----------------------------------------------------------
#                   Model 1
# Abundance = Herb Patch Density 
# Detection = Vegetation Density
# Call Rate = Day Ran Eff + conspecific density
# Vocalization = Ran eff of site & survey?
# ----------------------------------------------------------



# Parameters monitored
params <- c('lambda',
            'N_tot',
            'N',
            'alpha0', 
            'alpha1', 
            'alpha2',
            'beta0',
            'S.raneff',
            'sigma_S',
            'eta',
            'beta1',
            'kappa',
            'gamma1',
            'mu_gamma1',   
            'tau_gamma1',
            'a.phi',
            'omega', 
            'bp.y', 
            'bp.v')

# Initial Values 
inits <- function() {
  list(
    N = rep(1, R), # Abundance
    beta0 = rnorm(1),
    beta1 = rnorm(1), 
    sigma_S = runif(1, 0.1, 2), 
    eta = rnorm(R, 0, 1), 
    alpha0 = rnorm(1, 0, 1), # Detection
    alpha1 = runif(1, 0, 1), 
    alpha2 = rnorm(1, 0, 0.1),
    kappa = runif(1, 0.8, 1.5), # Call Rate
    a.phi = runif(1, 0, 5) # Vocalization
  )
}


# ----------------------------- 
# Model Statement 
# ----------------------------- 
cat(" 
model {

  # ---------------------------------
  # Priors 
  # ---------------------------------
  
  # ----------
  # Abundance
  # ----------
  beta0 ~ dnorm(0, .01)
  beta1 ~ dnorm(0, .01)
  sigma_S ~ dunif(0, 10)  # Prior on standard deviation
  
  # Site random effect on abundance - missing habitat variability
  for (i in 1:R) {
    eta[i] ~ dnorm(0, 1)  # Standard normal for non-centered approach
    S.raneff[i] <- beta0 + sigma_S * eta[i]  
  }

  # ----------
  # Detection
  # ----------
  alpha0 ~ dnorm(0, 1)
  alpha1 ~ dunif(0, 1000) # Constrained to be positive
  alpha2 ~ dnorm(0, 1)
 
  # ----------
  # Call Rate
  # ----------
  mu_gamma1 ~ dlnorm(1, 1)
  tau_gamma1 ~ dgamma(0.01, 0.01)
  kappa ~ dnorm(1, 1) T(0, ) # truncated to always be positive
  
  # Random effect of survey day
  for (j in 1:n.days) {
    gamma1[j] ~ dnorm(mu_gamma1, tau_gamma1) 
  }

  # ----------
  # Vocalizations
  # ----------
  a.phi ~ dunif(0, 100)
  omega ~ dunif(0, 1000)
  
  # random effect of site and survey??
  for (i in 1:R) {
    for (j in 1:J.A) {
      phi[i, j] ~ dgamma(a.phi, a.phi)
    }
  }
  
  # ---------------------------------
  #
  #    Liklihood and process model
  #
  # --------------------------------- 
  
  # Site loop
  for (i in 1:R) {
    
    # ---------------------------------
    # Abundance submodel  
    # ---------------------------------
    log(lambda[i]) <- beta0 + beta1 * X.abund[i,17] + S.raneff[i]
    N[i] ~ dpois(lambda[i])
    

    
    # Survey loop
    for (j in 1:J[i]) {
    
    
      # ---------------------------------
      # Detection submodel
      # ---------------------------------
      logit(p.a[i, j]) <- alpha0 + alpha1*log(N[i]) + alpha2*X.abund[i,21] 
    
    
      # ---------------------------------
      # Call Rate submodel  
      # ---------------------------------
      log(delta[i, j]) <- gamma1[days[i, j]] + kappa*log(N[i])
      y[i, j] ~ dbin(p.a[i, j], 1)
      tp[i, j] <- delta[i, j] * N[i] / (delta[i, j] * N[i] + omega)
      
      # Posterior predictive checks for Bayesian P-value
      y.pred[i, j] ~ dbin(p.a[i, j], 1)
      resid.y[i, j] <- pow(pow(y[i, j], 0.5) - pow(p.a[i, j], 0.5), 2)
      resid.y.pred[i, j] <- pow(pow(y.pred[i, j], 0.5) - pow(p.a[i, j], 0.5), 2)
    } # j
    
    
    for (j in 1:J.r[i]) {
      # ---------------------------------
      # Vocalization submodel  
      # ---------------------------------
      v[i, A.times[i, j]] ~ dpois((delta[i, A.times[i, j]] * N[i] + omega) * phi[i, A.times[i, j]] * y[i, A.times[i, j]]) T(1, )
      v.pred[i, j] ~ dpois((delta[i, A.times[i, j]] * N[i] + omega) * phi[i, A.times[i, j]] * y[i, A.times[i, j]]) T(1, )
      mu.v[i, j] <- ((delta[i, A.times[i, j]] * N[i] + omega) * phi[i, A.times[i, j]]) / (1 - exp(-1 * ((delta[i, A.times[i, j]] * N[i] + omega) * phi[i, A.times[i, j]])))
      resid.v[i, j] <- pow(pow(v[i, A.times[i, j]], 0.5) - pow(mu.v[i, j], 0.5), 2)
      resid.v.pred[i, j] <- pow(pow(v.pred[i, j], 0.5) - pow(mu.v[i, j], 0.5), 2)
    } # j

  } # i
  
  # ---------------------------------
  # Bayesian P-value  
  # ---------------------------------
  for (i in 1:R.A) {
    tmp.v[i] <- sum(resid.v[sites.a.v[i], 1:J.r[sites.a.v[i]]])
    tmp.v.pred[i] <- sum(resid.v.pred[sites.a.v[i], 1:J.r[sites.a.v[i]]])
  }
  
  # ---------------------------------
  # Derived Values  
  # ---------------------------------
  
  # Abundance
  N_tot <- sum(N[])
  
  # Bayes P-value
  fit.y <- sum(resid.y[sites.a, 1:J.A])
  fit.y.pred <- sum(resid.y.pred[sites.a, 1:J.A])
  fit.v <- sum(tmp.v[1:R.A])
  fit.v.pred <- sum(tmp.v.pred[1:R.A])
  bp.y <- step(fit.y.pred - fit.y)
  bp.v <- step(fit.v.pred - fit.v)

}
", fill=TRUE, file="./jags_models/ARU_A_mod1.txt")
# ------------End Model-------------


# Fit Model
out.model.A <- jags(data = Wolfe14.data, 
                   inits = inits, 
                   parameters.to.save = params, 
                   model.file = "./jags_models/ARU_A_mod1.txt", 
                   n.iter = n.iter,
                   n.burnin = n.burnin,
                   n.chains = n.chains, 
                   n.thin = n.thin,
                   parallel = TRUE,
                   n.cores = workers,
                   verbose = TRUE) 

# Trace plots
mcmcplot(out.model.A$samples)

# Rhat
out.model.A$Rhat

# Bayesian P value
cat("Bayesian p-value =", out.model.A$summary["bp.v",1], "\n")


# Model summary
print(out.model.A, digits = 2)



# -------------------------------------------------------
#
#   Estimating Abundance 
#
# -------------------------------------------------------

# load(file = "./Data/Model_Environments/ARU_A_Wolfe14day_JAGs.RData")

# Combine chains
combined_chains <- as.mcmc(do.call(rbind, out.model.A$samples))


# Extracting Abundance
Ntot_samples <- combined_chains[ ,"N_tot"]

# Ntotal is the abundance based on 27 acoustic sites at a radius of 200m.
# To correct for density, Ntotal needs to be divided by 27 * area surveyed
area <- pi * (250^2) / 4046.86  # Area in acres
dens_samples <- Ntot_samples / (area * 27)

# Create data frame for density
dens_df <- data.frame(Model = rep("AV Wolfe", length(dens_samples)), Density = dens_samples)
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

# Plot violin

# Create Abundance dataframe
abund_df <- dens_df
abund_df$Density <- abund_df$Density * 2710
colnames(abund_df)[2] <- 'Abundance'
head(abund_df)

# Plot violin
ggplot(abund_df, aes(x = Model, y = Abundance, fill = Model)) + 
  geom_violin(trim = FALSE, alpha = 0.6, adjust = 5) +  # Adjust bandwidth for smoothing
  labs(x = "Model", y = "Total Abundance") +
  scale_fill_manual(values = c("PC CMR" = "orange", 
                               "PC HDS" = "purple", 
                               "AV Bnet" = "blue",
                               "AV Wolfe" = "red")) +  
  scale_y_continuous(limits = c(0, 900),
                     breaks = seq(0, 900, by = 100),
                     labels = scales::comma) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),  
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
        panel.grid = element_blank(),
        legend.position = "none") 







# Total abundance
mean(dens_df$Density) * 2710





# Export density dataframe
saveRDS(dens_df, "./Data/Fitted_Models/ARU_WolfeA__dens_df.rds")
saveRDS(dens_summary, "./Data/Fitted_Models/ARU_WolfeA_dens_summary.rds")
saveRDS(abund_summary, "./Data/Fitted_Models/ARU_WolfeA_abund_summary.rds")

# Save Environment
save.image(file = "./Data/Model_Environments/ARU_A_Wolfe14day_JAGs.RData")

# End Script