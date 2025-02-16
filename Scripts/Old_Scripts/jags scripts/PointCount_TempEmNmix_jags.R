# -------------------------------------------------------
#
#                    Load libraries
#
# -------------------------------------------------------


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
workers <- Ncores * 0.80 # For low background use 80%, for medium use 50% of Ncores

# -------------------------------------------------------
#
#                    Load Data
#
# -------------------------------------------------------

# Load in capture data
#pc_dat <- read.csv("./Data/Point_Count_Data/NOBO_PC_Summer2024data.csv")
pc_dat <- read.csv("./Data/Point_Count_Data/NOBO_PC_Summer2024data - bin42bin3.csv")

# Load in site covariates
site_covs <- read.csv("./Data/Point_Count_Data/PointCount_siteCovs.csv")

# -------------------------------------------------------
#
#                   Data Wrangling
#
# -------------------------------------------------------

# Remove NAs
pc_dat_NAom <- na.omit(pc_dat)


# creating a matrix that is 4 Surveys  and 10 rows long
det_mat <- matrix(0, nrow = 10, ncol = 4)

# adding a column to state a NOBO was detected, a count column
pc_dat_NAom$count <- 1

# Loop to fill matrix with data
for (i in 1:nrow(pc_dat_NAom)) {
  point_num <- pc_dat_NAom$PointNum[i]
  occasion <- pc_dat_NAom$Survey[i]

  # Fill in the matrix with the number of individuals
  det_mat[point_num, occasion] <- det_mat[point_num, occasion] + pc_dat_NAom$count[i]
  
}#end loop

# Take a look
print(det_mat)


## Observation covariates
# Create matrix for each covariate
temp_mat <- matrix(NA, nrow = 10, ncol = 4)
wind_mat <- matrix(NA, nrow = 10, ncol = 4)
sky_mat <- matrix(NA, nrow = 10, ncol = 4)


# Fill the matrices
for (i in 1:nrow(pc_dat)) {
  # extract site and occasion
  point_num <- pc_dat$PointNum[i]
  occasion <- pc_dat$Survey[i]
  
  # fill mats
  temp_mat[point_num, occasion] <-  pc_dat$Temp.deg.F[i]
  wind_mat[point_num, occasion] <-  pc_dat$Wind.Beau.Code[i]
  sky_mat[point_num, occasion] <-  pc_dat$Sky.Beau.Code[i]

}# end loop

# Take a look
print(temp_mat)
print(wind_mat)
print(sky_mat)


# Extract and scale detection covariates for X.det array
X.det <- array(NA, dim = c(10,  # Number of sites
                           4,  # Number of surveys
                           3),                      # Number of covariates
               dimnames = list(NULL, NULL, c("Temp", "Wind", "Sky")))

# Assign each covariate to the respective slice in the array
X.det[, , "Temp"] <- as.matrix(scale(temp_mat))  # Scaled
X.det[, , "Wind"] <- as.matrix(scale(wind_mat))
X.det[, , "Sky"] <- as.matrix(sky_mat)
print(X.det)
X.det[, , 1]  # Temp
X.det[1, 2, 1] # Row 1, column 2, Observer


# Extract and scale site covariates for X.abund
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
print(X.abund)


# Constances 
K <- 4                          # Number of primary occasions
nsites <- nrow(det_mat)         # Number of sites
area <- pi*(200^2)/4046.86      # Area in acres


# Bundle data
data <- list(y = det_mat, 
             nsites = nsites, 
             K = K, 
             area = area,
             X.det = X.det,
             X.abund = X.abund)

# Look at structure
str(data)


# ---------------------------------------------------------- 
# 
#       Temporary Emigration N-mixture Model
# 
# ----------------------------------------------------------

# -------------------
# MCMC Specifications
# -------------------

n.iter = 300000
n.burnin = 20000
n.chains = 3 
n.thin = 10



# ---------------------------------------------------------- 
#                 Detection Models
# ----------------------------------------------------------


# ---------------------- 
# Model Specifications
# ---------------------- 


# Parameters monitored
params <- c("phi0", 
            "beta0",
            "beta1",
            "beta2",
            "alpha0",
            "alpha1",
            "gamma1", 
            "logit.gamma1",
            "gamma2",
            "lambda",
            "p_Bayes")
            



# Initial values
inits  <- function() {
  list(
    M = apply(det_mat, 1, max) + 5,
    gamma1 = rep(0.5, 4),
    gamma2 = 0,
    beta0 = 0,
    beta1 = 0,
    beta2 = 0,
    alpha0 = 0,
    alpha1 = 0,
    phi0 = 0.5
  )
}






# ---------------------------------------------------------- 
#                 Availability Models
# ----------------------------------------------------------



cat("
model {

  # --------------------
  # Prior
  # --------------------
  beta0 ~ dnorm(0, 0.01)   
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  alpha0 ~ dnorm(0, 0.01) 
  alpha1 ~ dnorm(0, 0.01)

  # Availability parameters
  phi0 ~ dunif(0.1, 0.9)
  logit.phi0 <- log(phi0 / (1 - phi0))
  
  for (k in 1:K) {
    gamma1[k] ~ dunif(0.1, 0.9)  # Availability effects for surveys
    logit.gamma1[k] <- log(gamma1[k] / (1 - gamma1[k]))
  }

  # --------------------
  # Likelihood
  # --------------------
  for (s in 1:nsites) {
    for (k in 1:K) {

      # Availability Model
      logit.phi[s,k] <- logit.gamma1[k]
      phi[s,k] <- exp(logit.phi[s,k]) / (1 + exp(logit.phi[s,k]))

      # Detection Model (logit-link)
      logit(p[s,k]) <- alpha0 + alpha1 * X.det[s,k,2]  
      
      # Probability of detection given availability
      pmarg[s,k] <- p[s,k] * phi[s,k]

      # Observation model (Binomial likelihood)
      y[s,k] ~ dbin(pmarg[s,k], M[s])  
      
      # Posterior Predictive Checks
      y_rep[s,k] ~ dbin(pmarg[s,k], M[s])
      discrepancy_obs[s,k] <- pow(y[s,k] - (pmarg[s,k] * M[s]), 2)
      discrepancy_rep[s,k] <- pow(y_rep[s,k] - (pmarg[s,k] * M[s]), 2)
    }

    # Abundance Model
    log(lambda[s]) <- beta0 + beta1 * X.abund[s,17]  + beta2 * X.abund[s,10]  
    M[s] ~ dpois(lambda[s])  
  }

  # Derived Quantities
  for (k in 1:K){
    Davail[k] <- mean(phi[,k]) * exp(beta0) / area
  }

  # Bayesian p-value Computation
  sum_obs <- sum(discrepancy_obs[,])
  sum_rep <- sum(discrepancy_rep[,])
  p_Bayes <- step(sum_rep - sum_obs)

} # End model
", fill=TRUE, file="./jags_models/Nmix_mod1.txt")
# ------------End Model-------------

# Run JAGs 
fm.1 <- jags(data = data, 
                  inits = inits, 
                  parameters.to.save = params, 
                  model.file = "./jags_models/NMix_mod1.txt",
                  n.iter = n.iter,
                  n.burnin = n.burnin,
                  n.chains = n.chains, 
                  n.thin = n.thin,
                  parallel = TRUE,
                  n.cores = workers,
                  DIC = TRUE)


# Best model fit. P-value = 0.5 means good fit, = 1 or 0 is a poor fit
cat("Bayesian p-value =", fm.1$summary["p_Bayes",1], "\n")

# Check convergence
fm.1$Rhat # Rhat: less than 1.1 means good convergence
mcmcplot(fm.1$samples)# Visually inspect trace plots

# Model summary
summary(fm.1$samples)

# Save Environment
save.image(file = "./NmixTempEm_JAGs.RData")


# -------------------------------------------------------
#
#   Estimating Abundance 
#
# -------------------------------------------------------

# Combine chains
combined_chains <- as.mcmc(do.call(rbind, fm.1$samples))

# Extract lambda estimates
lambda_columns <- grep("^lambda\\[", colnames(combined_chains))
lambda_samples <- combined_chains[ ,lambda_columns]

# mean abundance 
lambda_tot <- rowSums(lambda_samples)

# Area in hectares
#area <- pi*(200^2)/10000

# Area in acres
area <- pi*(200^2)/4046.86

# Getting density
dens_df <- as.data.frame(lambda_tot/area)

# Summarize by row
colnames(dens_df)[1] <- "Density"
dens_df[,2] <- "PC Nmix"
colnames(dens_df)[2] <- "Model"
dens_df <- dens_df[, c("Model", "Density")] # Switch the order of columns
head(dens_df)

# Calculate the 95% Credible Interval
ci_bounds <- quantile(dens_df$Density, probs = c(0.025, 0.975))

# Subset the data frame to 95% CI
dens_df <- subset(dens_df, Density >= ci_bounds[1] & Density <= ci_bounds[2])


# Violin plot
ggplot(dens_df, aes(x = Model, y = Density, fill = Model)) +
  geom_violin() +
  geom_boxplot(aes(x = Model, y = Density),
               width = 0.2, position = position_dodge(width = 0.8)) +
  labs(
    title = "",
    x = "Model",
    y = "Density (N/hectare)") +
  scale_y_continuous(limits = c(0, 20),
                     breaks = seq(0, 5, by = 5),
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


# total density
print(mean(dens_df$Density))
print(min(dens_df$Density))
print(max(dens_df$Density))


# Total abundance
mean(dens_df$Density) * 2710



# Save Environment
save.image(file = "HDS_JAGs.RData")

# Export density dataframe
saveRDS(dens_df, "./Data/Fitted_Models/PC_Nmix_DensityDF.rds")



# End Script

