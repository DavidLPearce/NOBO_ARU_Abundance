# -------------------------------------------------------
#
#                    Load libraries
#
# -------------------------------------------------------

# Install packages (if needed)
# install.packages("nimble")
# install.packages("coda")

# Load library
library(tidyverse)
library(nimble)
library(coda)
library(mcmcplots)

# Set seed, scientific notation options, and working directory
set.seed(123)
options(scipen = 9999)
setwd(".")

# Model name object
model_name <- "PC HDS"

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

# Creating a day of year column
pc_dat$Date <- as.POSIXct(pc_dat$Date, format = "%m/%d/%Y") 
pc_dat$DOY <- yday(pc_dat$Date) 

# ----------------------
# Observation Data  
# ----------------------

# Remove NAs
pc_dat_NAom <- na.omit(pc_dat)

# creating a matrix that is 4 Surveys * 4 Distance bins wide and 10 rows long
det_mat <- matrix(0, nrow = 10, ncol = (4 * 4))

# adding a column to state a NOBO was detected, a count column
pc_dat_NAom$count <- 1

# Loop to fill matrix with data
for (i in 1:nrow(pc_dat_NAom)) {
  point_num <- pc_dat_NAom$PointNum[i]
  occasion <- pc_dat_NAom$Survey[i]
  distance_cat <- as.numeric(pc_dat_NAom$DistBin[i])
  
  # Determine the column in the matrix
  col_index <- (occasion - 1) * 3 + distance_cat
  
  # Fill in the matrix with the number of individuals
  det_mat[point_num, col_index] <- det_mat[point_num, col_index] + pc_dat_NAom$count[i]
  
}#end loop

# Take a look
print(det_mat)

# Create a 3D array
y3d <- array(NA,dim=c(nrow(det_mat), 4, 4) ) # Length of site (10), width of distance bins (4), depth of surveys (4)

# Fill array
y3d[,,1] <- det_mat[,1:4]    
y3d[,,2] <- det_mat[,5:8]  
y3d[,,3] <- det_mat[,9:12]   
y3d[,,4] <- det_mat[,13:16]

# ----------------------
# Observation Covariates  
# ----------------------

# Create matrix for each covariate
obsvr_mat <- matrix(NA, nrow = 10, ncol = 4)
temp_mat <- matrix(NA, nrow = 10, ncol = 4)
wind_mat <- matrix(NA, nrow = 10, ncol = 4)
sky_mat <- matrix(NA, nrow = 10, ncol = 4)
doy_mat <- matrix(NA, nrow = 10, ncol = 4)


# Fill the matrices
for (i in 1:nrow(pc_dat)) {
  point_num <- pc_dat$PointNum[i]
  occasion <- pc_dat$Survey[i]
  
  obsvr_mat[point_num, occasion] <-  pc_dat$Observer[i]
  temp_mat[point_num, occasion] <-  pc_dat$Temp.deg.F[i]
  wind_mat[point_num, occasion] <-  pc_dat$Wind.Beau.Code[i]
  sky_mat[point_num, occasion] <-  pc_dat$Sky.Beau.Code[i]
  doy_mat[point_num, occasion] <-  pc_dat$DOY[i]
  
}# end loop

# Take a look
print(obsvr_mat)
print(temp_mat)
print(wind_mat)
print(sky_mat)
print(doy_mat)

# Categorical covariates
obsvr_det <- apply(obsvr_mat, 2, function(col) as.integer(factor(col))) # Observer

wind_factor <- as.integer(factor(as.vector(wind_mat))) # Wind
wind_det <- matrix(wind_factor, nrow = nrow(wind_mat), ncol = ncol(wind_mat))

sky_factor <- as.integer(factor(as.vector(sky_mat))) # Sky
sky_det <- matrix(sky_factor, nrow = nrow(sky_mat), ncol = ncol(sky_mat))

doy_factor <- as.integer(factor(as.vector(doy_mat))) # Day of Year
doy_det <- matrix(doy_factor, nrow = nrow(doy_mat), ncol = ncol(doy_mat))


# Detection array
X_det <- array(NA, dim = c(10, 4, 5))
X_det[,,1] <- obsvr_det
X_det[,,2] <- temp_mat
X_det[,,3] <- wind_det
X_det[,,4] <- sky_det
X_det[,,5] <- doy_det

print(X_det)

# Categorical covariate levels
Obs_Lvls = length(unique(as.numeric(as.factor(pc_dat$Observer))))
Wind_Lvls = length(unique(as.numeric(as.factor(pc_dat$Wind.Beau.Code))))
Sky_Lvls = length(unique(as.numeric(as.factor(pc_dat$Sky.Beau.Code))))


# ----------------------
# Site Covariates  
# ----------------------

# Remove x, siteid, lat, long
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




# ----------------------
# Constants  
# ----------------------
J <- 4                          # Number of primary occasions
S <- 10                         # Number of sites
nD <- 4                         # Number of distance classes
delta <- 50                     # Class width
B <- 200                        # Maximum distance
midpt <- seq(delta/2, B, delta) # Class midpoint distance
nobs <- apply(y3d, c(1,3), sum) # Total detections per site and occasion
area <- pi*(200^2)/10000        # Area surveyed in hectares

# ----------------------
# Bundle Data  
# ----------------------

# Bundle data for nimble
data <- list(y3d = y3d, 
             nobs = nobs,
             X_abund = X_abund
)


# Take a look
str(data)

# State Constants
constants <- list(S = S, 
                  J = J, 
                  nD = nD, 
                  midpt = midpt, 
                  delta = delta, 
                  B = B,
                  X_det = X_det, # For some reason this has to be here, will give warning
                  Wind_Lvls = Wind_Lvls
                  # Sky_Lvls = Sky_Lvls
)

# Take a look
str(constants)




# ------------------------------------------------------------------------------ 
#
#                             Model Fitting
#
# ------------------------------------------------------------------------------ 

# -------------------
# MCMC Specifications
# -------------------

niter = 1000
nburnin = 100
nchains = 3 
nthin =  1

# ----------------------
# Model Specifications
# ----------------------

# Parameters monitored
params <- c(
            # Abundance
            "N",
            "N_tot",
            "lambda",
            "beta0", 
            "beta1",
            "beta2",

            # Detection
            "pdet",
            "pmarg",
            "sigma0", 
            "alpha1",
            "alpha2",
            
            # Availability
            "gamma0", 
            "jRE",
            "tau_j",
            
            # PPC
            "fit_y",
            "fit_y_pred",
            "bp_y"
)


# Initial Values 
Navail.st <- apply(y3d, c(1,3),sum)
Mst <- apply(Navail.st, c( 1), max,na.rm = TRUE) + 5
inits <- function() {
  list(
    
    # Abundance
    M = Mst,
    Navail = Navail.st,
    beta0 = 0,
    beta1 = 0,
    beta2 = 0,
    
    # Detection
    sigma0 = 100,
    alpha1 = rep(0, Wind_Lvls),
    alpha2 = 0,
    
    # Availability
    gamma0 = 0,
    tau_j = 1,
    jRE = rep(0, J),
    
    # PPC 
    # will get errors if rep is not initialized
    # Initializing at y seems to work
    y3d_rep = y3d  

  )
}


# ----------------------------- 
# Model Statement 
# ----------------------------- 
TEHDS_model <- nimbleCode({
  
  # ---------------------------------
  # Abundance Priors
  # ---------------------------------
  
  # Intercept
  beta0 ~ dnorm(0, 10)
  
  # Covariate
  beta1 ~ dnorm(0, 10)
  beta2 ~ dnorm(0, 10)

  # ---------------------------------
  # Availability Priors
  # ---------------------------------
  
  # Intercept
  gamma0 ~ dnorm(0, 10)

  #Survey random effect
  tau_j ~ dgamma(0.01, 0.01)
  for (j in 1:J) {
    jRE[j] ~ dnorm(gamma0, tau_j)
  }
  
  # ---------------------------------
  # Detection Priors
  # ---------------------------------
  
  # Intercept
  sigma0 ~ dunif(0.1, 500)
  
  # Covariate Effects
  
  
  # Covariate effect
  for (w in 1:Wind_Lvls){ # Wind is a categorical covariate
    alpha1[w] ~ dnorm(0, 10)
  }
  
  alpha2 ~ dnorm(0, 10) # Vegetation Density
  
  # -------------------------------------------
  #
  # Likelihood and Process Model 
  #
  # -------------------------------------------
  
  
  for (s in 1:S) {
    for (j in 1:J) {
      
      # ---------------------------------
      # Availability 
      # ---------------------------------

      logit.phi[s,j] <- jRE[j]
      phi[s,j] <- exp(logit.phi[s,j])/(1+ exp(logit.phi[s,j]))
      
      # ---------------------------------
      # Detection 
      # ---------------------------------
      
      # Distance Sampling = Intercept + Wind + Veg Density 
      log(sigma[s,j]) <- log(sigma0) + alpha1[X_det[s, j, 3]] +  alpha2 * X_abund[s,21] 
      
      # Multinomial cell probability construction
      for(b in 1:nD){
        
        # Half-normal Detection Function
        log(g[s,b,j]) <- -midpt[b]*midpt[b]/(2*sigma[s,j]*sigma[s,j]) 
        
        # Hazard Rate Detection Function
        #cloglog(g[s,b,j]) <- theta*log(sigma[s,j])  - theta*log(midpt[b]) 
        
        # Density function for distance bins
        f[s,b,j] <- (2*midpt[b]*delta)/(B*B)
        cellprobs[s,b,j] <- g[s,b,j]*f[s,b,j]
        cellprobs.cond[s,b,j] <- cellprobs[s,b,j]/sum(cellprobs[s,1:nD,j])
        
      } # End B
      
      # Add probability of undetected individuals
      
      cellprobs[s,nD+1,j]<- 1-sum(cellprobs[s,1:nD,j])
      
      # Detection probabilities
      pdet[s,j] <- sum(cellprobs[s,1:nD,j])
      pmarg[s,j] <- pdet[s,j]*phi[s,j]
      
      # ---------------------------------
      # Observation 
      # ---------------------------------
      
      y3d[s,1:nD,j] ~ dmulti(cellprobs.cond[s,1:nD,j], nobs[s,j])  
      y3d_rep[s,1:nD,j] ~ dmulti(cellprobs.cond[s,1:nD,j], nobs[s,j])
      
      # Number of detected individuals
      nobs[s,j] ~ dbin(pmarg[s,j], M[s])   
      
      # Number of available individuals
      Navail[s,j] ~ dbin(phi[s,j],M[s]) 
      
    } # end J loop
    
    # ---------------------------------
    # Abundance 
    # ---------------------------------
    
    M[s] ~ dpois(lambda[s])
    log(lambda[s]) <- beta0 + beta1 * X_abund[s, 7] +  beta2 * X_abund[s, 12]

  }  # end S loop

  # ---------------------------------
  # Posterior Predictive Checks
  # ---------------------------------
  for (s in 1:S) {
    for (j in 1:J) {
       for (b in 1:nD) {
        discrepancy_obs[s,b,j] <- pow(y3d[s,b,j] - (cellprobs.cond[s,b,j] * nobs[s,j]), 2)
        discrepancy_rep[s,b,j] <- pow(y3d_rep[s,b,j] - (cellprobs.cond[s,b,j] * nobs[s,j]), 2)
       }
      tmp_obs[s, j] <- sum(discrepancy_obs[s, 1:nD, j])
      tmp_rep[s, j] <- sum(discrepancy_rep[s, 1:nD, j])
    }
  }

  # Total fit
  fit_y <- sum(tmp_obs[1:S, 1:J])
  fit_y_pred <- sum(tmp_rep[1:S, 1:J])

  # Bayes p-value
  bp_y <- step(fit_y_pred - fit_y)
  
  # ---------------------------------
  # Derived Quantities
  # ---------------------------------
  
  for (s in 1:S) {
    for (j in 1:J) {
      N_site_j[s,j] <- pdet[s,j] * phi[s,j] * M[s]
    }
    N[s] <- sum(N_site_j[s, 1:J])
  }

  # Abundance
  N_tot <- sum(N[1:S])
  
 
  
})  
# ---------------------------- End Model ----------------------------



# Fit Model
fm1 <- nimbleMCMC(code = TEHDS_model,
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
# saveRDS(fm1, "./Data/Model_Data/ModelFits_PC-TEHDS_fm1.rds")

# -------------------------------------------------------
# Check Convergence
# -------------------------------------------------------

# Extract posterior samples
fm1_samples <- as.mcmc.list(fm1$samples)

# Trace plots
mcmcplots::mcmcplot(fm1_samples, 
                    parms = params) 

# Rhat
coda::gelman.diag(fm1_samples, multivariate = FALSE)

# Effective sample size
as.data.frame(coda::effectiveSize(fm1_samples))


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
mn_bp_y <- mean(bp_y_samples)


# P-value = 0.5 means good fit, = 1 or 0 is a poor fit

# Abundance
cat("Abundance Model Bayesian p-value =", mn_bp_y, "\n")

# ----------------------
# Extract Fits
# ----------------------

# Abundance
fit_y_data <- data.frame(
  Observed = as.vector(combined_chains[, "fit_y"]),       # Observed values
  Predicted = as.vector(combined_chains[, "fit_y_pred"]) # Predicted values
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
ggsave(plot = y_PPC_Dens, "./Figures/PPC/PC-TEHDS_Abund_Density.jpeg", width = 8, height = 5, dpi = 300)
dev.off()

# -------------------------------------------------------
#
#   Beta Estimates and Covariate Effects 
#
# -------------------------------------------------------

# -------------------------------------------------------
# Beta Estimates
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
saveRDS(beta_df, "./Data/Model_Data/Beta_df_PC-TEHDS.rds")


# -------------------------------------------------------
# Covariate Effects
# -------------------------------------------------------

# Covariate names
print(colnames(X_abund))

# Set covariate name 
Cov1_name <- "herb_ClmIdx"
Cov2_name <- "woody_AggInx"

# Create a prediction of covariate values
cov1_pred_vals <- seq(min(X_abund[, Cov1_name]), max(X_abund[, Cov1_name]), length.out = 1000)
cov2_pred_vals <- seq(min(X_abund[, Cov2_name]), max(X_abund[, Cov2_name]), length.out = 1000)

# Mean scaling covariates
cov1_scaled <- (X_abund[, Cov1_name] - mean(X_abund[, Cov1_name])) / (max(X_abund[, Cov1_name]) - min(X_abund[, Cov1_name]))
cov2_scaled <- (X_abund[, Cov2_name] - mean(X_abund[, Cov2_name])) / (max(X_abund[, Cov2_name]) - min(X_abund[, Cov2_name]))

# Matrices for storing predictions
cov1_preds <- matrix(NA, nrow = length(beta0_samples), ncol = length(cov1_scaled))
cov2_preds <- matrix(NA, nrow = length(beta0_samples), ncol = length(cov2_scaled))

# Generate predictions
for (i in 1:length(beta0_samples)) {
  cov1_preds[i, ] <- beta0_samples[i] + beta1_samples[i] * cov1_scaled # Linear
  cov2_preds[i, ] <- beta0_samples[i] + beta2_samples[i] * cov2_scaled
}


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

# # Cov 2
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



# -------------------------------------------------------
#   Estimating Abundance 
# -------------------------------------------------------

# Extract abundance posterior
Ntot_samples <- combined_chains[ ,"N_tot"]

# Ntotal is the abundance based on 10 point counts at a radius of 200m.
# To correct for density, Ntotal needs to be divided by 10 * area surveyed
area <- pi * (200^2) / 10000  # Area in acres
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
dens_df <- dens_df[dens_df$Density >= dens_summary$Lower_CI & dens_df$Density <= dens_summary$Upper_CI, ]

# Getting total abundance
abund_summary <- dens_summary
abund_summary[,2:4] <- abund_summary[,2:4] * 1096

# Plot Abundance - Violin
abund_df <- dens_df
abund_df$Density <- abund_df$Density * 1096

ggplot(abund_df, aes(x = Model, y = Density, fill = Model)) + 
  geom_violin(trim = FALSE, alpha = 0.6, adjust = 5) +   
  labs(x = "Model", y = "Total Abundance ") +
  scale_fill_manual(values = c("PC HDS" = "Purple")) +   
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

# Export Dataframes
saveRDS(dens_df, "./Data/Model_Data/Density_df_PC-TEHDS.rds")
saveRDS(dens_summary, "./Data/Model_Data/Density_summary_PC-TEHDS.rds")
saveRDS(abund_df, "./Data/Model_Data/Abundance_df_PC-TEHDS.rds")
saveRDS(abund_summary, "./Data/Model_Data/Abundance_summary_PC-TEHDS.rds")

# -------------------------- End Script ------------------------------------