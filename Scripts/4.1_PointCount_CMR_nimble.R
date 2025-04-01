
# ****** add header ******

# -------------------------------------------------------
#
#                    Load libraries
#
# -------------------------------------------------------

# Install packages (if needed)
# install.packages("tidyverse")
# install.packages("plotly")
# install.packages("nimble")
# install.packages("coda")
# install.packages("mcmcplots")
 
# Load library
library(tidyverse)
library(plotly)
library(nimble)
library(coda)
library(mcmcplots)


# Set seed, scientific notation options, and working directory
set.seed(123)
options(scipen = 9999)
setwd(".")

# Model name object
model_name <- "PC CMR"


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

# Number of time intervals
J <- 4

# Number of sites
S <- 10

# ----------------------
# Observation Covariates  
# ----------------------

# Creating a day of year column
pc_dat$Date <- mdy(pc_dat$Date)
pc_dat$DOY <- yday(pc_dat$Date)

# Create matrix for each covariate
obsvr_mat <- matrix(NA, nrow = S, ncol = J)
temp_mat <- matrix(NA, nrow = S, ncol = J)
wind_mat <- matrix(NA, nrow = S, ncol = J)
sky_mat <- matrix(NA, nrow = S, ncol = J)
doy_mat <- matrix(NA, nrow = 10, ncol = 4)

# Fill the matrices
for (i in 1:nrow(pc_dat)) {
  # extract site and occasion
  point_num <- pc_dat$PointNum[i]
  occasion <- pc_dat$Survey[i]
  
  # fill mats
  obsvr_mat[point_num, occasion] <-  pc_dat$Observer[i]
  temp_mat[point_num, occasion] <-  pc_dat$Temp.deg.F[i]
  wind_mat[point_num, occasion] <-  pc_dat$Wind.Beau.Code[i]
  sky_mat[point_num, occasion] <-  pc_dat$Sky.Beau.Code[i]
  doy_mat[point_num, occasion] <-  as.integer(pc_dat$DOY[i])
  
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
X_det <- array(NA, dim = c(S, J, 5))
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
# Observation Data  
# ----------------------

# Remove rows with NA values
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


# ----------------------
# Superpopulation  
# ----------------------

# Number of individuals detected
nind <- nrow(y)

# Superpopulation 
print(nind)
M <- nind + 100 # ensure that M is larger than number detected
print(M)

# Data Augmentation
y <- rbind(y, matrix(0, nrow = (M - nind), ncol = 4))

# For unobserved individuals 
site <- as.numeric(factor(pc_CMR$PointNum))

# Add in NA for M
site <- c(site, rep(NA, M-nind))

# Take a look
print(site)



# ----------------------
# Bundle Data  
# ----------------------

# Bundle data for nimble
data <- list(y = y, 
             X_abund = X_abund,
             X_det = X_det,
             group = site)

# Take a look
str(data)

# State Constants
constants <- list(J = J, 
                  M = M,
                  S = S,
                  Wind_Lvls = Wind_Lvls
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

niter = 600000
nburnin = 100000
nchains = 3 
nthin =  15


# ----------------------
# Model Specifications
# ----------------------

# Parameters monitored
params <- c("lambda",
            "N",
            "N_tot",
            "alpha0", 
            "alpha1", 
            "alpha2",
            "beta0", 
            "beta1",
            "beta2",
            "psi",
            "sRE",
            "tau_s",
            "fit_y",
            "fit_y_pred",
            "bp_y"
)

# Initial Values 
inits <- function() {
  list (p0 = runif(1),
        alpha1 = rep(0, Wind_Lvls),
        alpha2 = 0,
        alpha2 = 0,
        lambda = runif(10, min = 0.1, max = 10),
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        tau_s = 1,
        sRE = rep(0, S),
        z = c( rep(1, nind), rep(0, (M - nind))),
        # To avoid nimble warnings, initialize unobserved groups
        group = c(rep(NA,length(pc_dat$PointNum)),
                  sample(1:S, 
                         size = length(site) - length(pc_dat$PointNum), 
                         replace = TRUE)),
        # Will get log prob inf if y_rep is not initialized
        y_rep = matrix(0, nrow = M, ncol = J)
        
)}# end inits

# ----------------------------- 
# Model Statement 
# ----------------------------- 
CMR_model <- nimbleCode({
  
  # ---------------------------------
  # Abundance Priors
  # ---------------------------------
  beta0 ~ dnorm(0, 10)
  beta1 ~ dnorm(0, 10)
  beta2 ~ dnorm(0, 10)

  # ---------------------------------
  # Detection Priors
  # ---------------------------------
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0)) 

  # Covariate effect
  for (w in 1:Wind_Lvls){# Wind is a categorical covariate
    alpha1[w] ~ dnorm(0, 10)
  }
  
  alpha2 ~ dnorm(0, 10)# Vegetation Density
  
  # Site-level random effect for pseudoreplication
  tau_s ~ dgamma(0.01, 0.01) 
  for(s in 1:S){  
    sRE[s] ~ dnorm(alpha0, tau_s)  
  }
  
  # ---------------------------------
  # Individual Encounter/Presence Derived 
  # ---------------------------------
  psi <- sum(lambda[1:S]) / M 

  # -------------------------------------------
  #
  # Likelihood and Process Model 
  #
  # -------------------------------------------
  
  # ---------------------------------
  # Abundance Submodel
  # ---------------------------------

  # Abundance model: Intercept + Herb Clumpy Indx + Woody Agg Indx  
  for(s in 1:S){
    
    log(lambda[s]) <- beta0 + beta1 * X_abund[s, 7] +  beta2 * X_abund[s, 12]
    
    # Individual site probability
    probs[s] <- lambda[s] / sum(lambda[1:S])
    
    # Estimated abundance at site s
    N[s] <- sum(group[] == s)
    
    
  }
  
  # ---------------------------------
  # Presence Submodel 
  # ---------------------------------

  for(i in 1:M){
    
    # Group == site membership
    group[i] ~ dcat(probs[1:S]) 

    # Presence: Data augmentation variables
    z[i] ~ dbern(psi)
    

    # ---------------------------------
    # Detection Submodel
    # ---------------------------------
    
    # Detection model: Site Random Effect[Intercept] + Wind + Veg Density
    for(j in 1:J){

      logit(p[i,j]) <- sRE[group[i]] + alpha1[X_det[group[i], j, 3]]  + alpha2 * X_abund[group[i], 21] 
      pz[i,j] <- p[i,j] * z[i]
      
      # Observation & Generate replicated data
      y[i,j] ~ dbern(pz[i,j])
      y_rep[i,j] ~ dbern(pz[i,j])
      
    } # End J
  } # End M  

      
    # ---------------------------------
    # Posterior Predictive checks
    # ---------------------------------
    
    # Compute discrepancies
    for(i in 1:M){
      for(j in 1:J){
        discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j]) 
        discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])
      }
    }
    
    # Sum by site
    for(i in 1:M){
      tmp_obs[i] <- sum(discrepancy_obs[i, 1:J])
      tmp_rep[i] <- sum(discrepancy_rep[i, 1:J])
    }
    
    # Total fit
    fit_y <- sum(tmp_obs[1:M])
    fit_y_pred <- sum(tmp_rep[1:M])
    
    # Bayes p-value
    bp_y <- step(fit_y_pred - fit_y)
  
  
    # ---------------------------------
    # Derived Metrics
    # ---------------------------------
    
    # Abundance
    N_tot <- sum(z[])

})
# ---------------------------- End Model ----------------------------


# Fit Model
fm1 <- nimbleMCMC(code = CMR_model,
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
saveRDS(fm1, "./Data/Model_Data/ModelFits_PC-CMR_fm1.rds")

# fm1 <- readRDS("./Data/Model_Data/ModelFits_WolfeAV_fm1.rds")

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
ggsave(plot = y_PPC_Dens, "./Figures/PPC/PC-CMR_Abund_Density.jpeg", width = 8, height = 5, dpi = 300)
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
saveRDS(beta_df, "./Data/Model_Data/Beta_df_PC-CMR.rds")

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
  scale_fill_manual(values = c("PC CMR" = "orange")) +   
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
saveRDS(dens_df, "./Data/Model_Data/Density_df_PC-CMR.rds")
saveRDS(dens_summary, "./Data/Model_Data/Density_summary_PC-CMR.rds")
saveRDS(abund_df, "./Data/Model_Data/Abundance_df_PC-CMR.rds")
saveRDS(abund_summary, "./Data/Model_Data/Abundance_summary_PC-CMR.rds")

# -------------------------- End Script ------------------------------------