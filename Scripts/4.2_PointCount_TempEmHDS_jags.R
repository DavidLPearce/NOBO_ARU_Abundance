# -------------------------------------------------------
#
#                    Load libraries
#
# -------------------------------------------------------


# Load library
library(tidyverse)
library(plotly)
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
workers <- Ncores * 0.5 # For low background use 80%, for medium use 50% of Ncores
print(workers)

# Source custom function for checking Rhat values > 1.1
source("./Scripts/Rhat_check_function.R")

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


# Remove NAs
pc_dat_NAom <- na.omit(pc_dat)

# creating a matrix that is 4 Surveys * 4 Distance bins wide and 10 rows long
det_mat <- matrix(0, nrow = 10, ncol = 16)

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


## Observation covariates
# Create matrix for each covariate
obsvr_mat <- matrix(NA, nrow = 10, ncol = 4)
temp_mat <- matrix(NA, nrow = 10, ncol = 4)
wind_mat <- matrix(NA, nrow = 10, ncol = 4)
sky_mat <- matrix(NA, nrow = 10, ncol = 4)
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

# Convert Observer to numeric factor levels
Observer_numeric <- matrix(as.numeric(as.factor(obsvr_mat)), 
                           nrow = nrow(obsvr_mat), 
                           ncol = ncol(obsvr_mat))


# Extract and scale detection covariates for X.det array
X.det <- array(NA, dim = c(10,  # Number of sites
                                 4,  # Number of surveys
                                 5),                      # Number of covariates
                     dimnames = list(NULL, NULL, c("Observer", "Temp", "Wind", "Sky", "DOY")))

# Assign each covariate to the respective slice in the array
X.det[, , "Observer"] <- as.matrix(Observer_numeric)
X.det[, , "Temp"] <- as.matrix(scale(temp_mat))  # Scaled
X.det[, , "Wind"] <- as.matrix(scale(wind_mat))
X.det[, , "Sky"] <- as.matrix(sky_mat)
X.det[, , "DOY"] <- as.matrix(doy_mat)  
print(X.det)
X.det[, , 1]  # Observer
X.det[1, 2, 1] # Row 1, column 2, Observer

# Format X.abund
X.abund <- as.matrix(site_covs[,-c(1:2)]) # Remove X and site id
print(X.abund)


# Create a 3D array
# Length of site (10), width of distance bins (4), depth of surveys (4)
y3d <- array(NA,dim=c(nrow(det_mat), 4, 4) ) 

# Fill array
y3d[,,1] <- det_mat[,1:4]    
y3d[,,2] <- det_mat[,5:8]  
y3d[,,3] <- det_mat[,9:12]   
y3d[,,4] <- det_mat[,13:16]

# Constances 
J <- 4                          # Number of occasions
S <- nrow(det_mat)              # Number of sites
nD <- 4                         # Number of distance classes
delta <- 50                     # Class width
B <- 200                        # Maximum distance
midpt <- seq(delta/2, B, delta) # Class midpoint distance
nobs <- apply(y3d, c(1,3), sum) # Total detections per site and occasion
area <- pi*(200^2)/4046.86      # Area in acres


# Bundle data for nimble
data <- list(y3d = y3d, 
             S = S, 
             J = J, 
             nD = nD, 
             midpt = midpt, 
             delta = delta, 
             B = B,
             nobs = nobs, 
             area = area,
             X.det = X.det,
             X.abund = X.abund)

# Look at structure
str(data)

# ---------------------------------------------------------- 
# 
#       Temporary Emigration Hierarchical Distance Model
# 
# ----------------------------------------------------------

# -------------------------------------------------------
# Model Specifications
# -------------------------------------------------------

## Distributions
# library(ggfortify)
# ggdistribution(dunif, seq(0, 1, 0.001), min = 0, max = 1) # p0
# ggdistribution(dnorm, seq(0, 0.01, 0.0001)) # alpha's and beta's 
# ggdistribution(dgamma, seq(0, 5, 0.01), shape = 0.1, rate = 0.1) # tau


# -------------------
# MCMC Specifications
# -------------------

n.iter = 200000
n.burnin = 60000
n.chains = 3 
n.thin = 10

# Rough idea posterior samples
est_post_samps = (((n.iter - n.burnin) / n.thin) * n.chains)
print(est_post_samps)

# ----------------------------------------------------------
#                   Model 1 
# Availability =  ran eff of survey/day
# Detection =  Wind
# Abundance = Herbaceous Patch Density
# ----------------------------------------------------------

# Parameters monitored
params <- c("beta0", # Abundance
            "beta1",
            "beta2",
            "beta3",
            "lambda",
            "N",
            "N_tot",
            "sigma0", # Detection
            "alpha1",
            "pdet_mean",
            "gamma0", # Availability
            "Jraneff",
            "tau_j",
            "sigma_j",
            "theta",
            "phi_mean",
            "pmarg_mean",
            "log_lik",# Bayes p-value
            "p_Bayes")


# Initial values
inits  <- function() {
  list(
    M = apply(y3d, 1, max) + 5, # Abundance
    beta0 = 0,
    beta1 = 0,
    beta2 = 0,
    beta3 = 0,
    sigma0 = 100, # Detection
    alpha1 = 0,
    Navail = apply(y3d, c(1, 3), sum), # Availability
    gamma0 = 0,
    gamma2 = 0
  )
}


# ----------------------------- 
# Model 1 Statement 
# ----------------------------- 
cat("
model {

  # ---------------------------------
  # Abundance Priors
  # ---------------------------------
  
  # Intercept
  beta0 ~ dnorm(0, 10)
  
  # Covariate
  beta1 ~ dnorm(0, 10)
  beta2 ~ dnorm(0, 10)
  beta3 ~ dnorm(0, 10)
  
  # ---------------------------------
  # Availability Priors
  # ---------------------------------
  
  # Intercept
  gamma0 ~ dnorm(0, 10)
  
  # Survey random effect
  tau_j ~ dgamma(5, 5) # variance
  sigma_j <- sqrt(1 / tau_j) # standard deviation  
  for (j in 1:J) {
    Jraneff[j] ~ dnorm(0, sigma_j) 
  }

  # ---------------------------------
  # Detection Priors
  # ---------------------------------

  # Intercept
  sigma0 ~ dunif(0.1,200)
  
  # Covariate
  alpha1 ~ dnorm(0, 10)

  # Hazard Rate Detection Function
  theta ~ dgamma(1, 1)
  
  # -------------------------------------------
  #
  # Likelihood and Process Model 
  #
  # -------------------------------------------

  for (s in 1:S) {
    for (j in 1:J) {
    
  # ---------------------------------
  # Availability Submodel
  # ---------------------------------
 
      logit.phi[s,j] <- gamma0 + Jraneff[j]
      phi[s,j] <- exp(logit.phi[s,j]) / (1 + exp(logit.phi[s,j]))
      
  # ---------------------------------
  # Detection Submodel
  # ---------------------------------
  
      # Distance Sampling
      log(sigma[s,j]) <- log(sigma0) + alpha1*X.abund[s,21]  
      
      # Multinomial cell probability construction
      for(b in 1:nD){
      
        # Half-normal Detection Function
        log(g[s,b,j]) <- -midpt[b]*midpt[b]/(2*sigma[s,j]*sigma[s,j]) 
        
        # Hazard Rate Detection Function
        #cloglog(g[s,b,j]) <- theta*log(sigma[s,j])  - theta*log(midpt[b]) 
        
        # Density function for distance bins
        f[s,b,j] <- (2 * midpt[b] * delta) / (B * B)
        cellprobs[s,b,j] <- g[s,b,j] * f[s,b,j]
        cellprobs.cond[s,b,j] <- cellprobs[s,b,j] / sum(cellprobs[s,1:nD,j])
      }
      
      # Add probability of undetected individuals
      cellprobs[s,nD+1,j] <- 1 - sum(cellprobs[s,1:nD,j])

      # Detection probabilities
      pdet[s,j] <- sum(cellprobs[s,1:nD,j])
      pmarg[s,j] <- pdet[s,j] * phi[s,j]
      
    # ---------------------------------
    # Observation Submodel
    # ---------------------------------
 
      y3d[s,1:nD,j] ~ dmulti(cellprobs.cond[s,1:nD,j], nobs[s,j])

      # Number of detected individuals
      nobs[s,j] ~ dbin(pmarg[s,j], M[s])

      # Number of available individuals
      Navail[s,j] ~ dbin(phi[s,j], M[s])

    # ---------------------------------
    # Log-Likelihood
    # ---------------------------------

      log_lik[s,j] <- logdensity.multi(y3d[s,1:nD,j], cellprobs.cond[s,1:nD,j], nobs[s,j])
      
    # ---------------------------------
    # Posterior Predictive Checks
    # ---------------------------------

      y3d_rep[s,1:nD,j] ~ dmulti(cellprobs.cond[s,1:nD,j], nobs[s,j])

      for (b in 1:nD) {
        discrepancy_obs[s,b,j] <- pow(y3d[s,b,j] - (cellprobs.cond[s,b,j] * nobs[s,j]), 2)
        discrepancy_rep[s,b,j] <- pow(y3d_rep[s,b,j] - (cellprobs.cond[s,b,j] * nobs[s,j]), 2)
      }
      
    } # End k loop

  # ---------------------------------
  # Abundance Submodel
  # ---------------------------------

    # Poisson
    log(lambda[s]) <- beta0 + beta1 * X.abund[s, 7] +  beta2 * X.abund[s, 12]
    M[s] ~ dpois(lambda[s])

  } # End s loop

  # ---------------------------------
  # Derived Quantities
  # ---------------------------------

  for (s in 1:S) {
      for (j in 1:J) {
        N_site_k[s,j] <- pdet[s,j] * phi[s,j] * M[s]    
      }
    N[s] <- sum(N_site_k[s,])
  }
  
  # Abundance
  N_tot <- sum(N[])
  
  # Detection
  pdet_mean<- mean(pdet[,])
  
  # Availability
  phi_mean<- mean(phi[,])
  
  # Mean marginal detection probability
  pmarg_mean <- mean(pmarg[,])

  # ---------------------------------
  # Derived Quantities
  # ---------------------------------
  # Bayesian p-value Computation
  sum_obs <- sum(discrepancy_obs[, ,])
  sum_rep <- sum(discrepancy_rep[, ,])
  p_Bayes <- step(sum_rep - sum_obs)  # Bayesian p-value

} # End model
", fill=TRUE, file="./jags_models/HDS_abundmod1.txt")
# ------------End Model-------------


# Run JAGs 
fm.1 <- jags(data = data, 
            inits =inits, 
            parameters.to.save = params, 
            model.file = "./jags_models/HDS_abundmod1.txt",
            n.iter = n.iter,
            n.burnin = n.burnin,
            n.chains = n.chains, 
            n.thin = n.thin,
            parallel = TRUE,
            n.cores = workers,
            DIC = TRUE)



# Check convergence
check_rhat(fm.1$Rhat, threshold = 1.1) # Rhat: less than 1.1 means good convergence
mcmcplot(fm.1$samples) # Visually inspect trace plots

# Check autocorrelation
#mcmc_list <- as.mcmc.list(lapply(fm.1$samples, as.mcmc))  # Convert each chain separately
#autocorr.diag(mcmc_list)  # Check autocorrelation per chain

# Check Effective Sample Size
#effectiveSize(mcmc_list)

# Check Fit
# cat("Bayesian p-value =", fm.1$summary["p_Bayes",1], "\n") #  P-value = 0.5 means good fit, = 1 or 0 is a poor fit

# Model summary
#summary(fm.1$samples)
 

# -------------------------------------------------------
#
#   Beta Estimates and Covariate Effects 
#
# -------------------------------------------------------

# Combine chains
combined_chains <- as.mcmc(do.call(rbind, fm.1$samples))

# -------------------------------------------------------
# Beta Estimates
# -------------------------------------------------------


# Extract beta estimates
beta0_samples <- combined_chains[, "beta0"]
beta1_samples <- combined_chains[, "beta1"]
beta2_samples <- combined_chains[, "beta2"]
# beta3_samples <- combined_chains[, "beta3"]

# Means
beta0 <- mean(beta0_samples)
beta1 <- mean(beta1_samples)
# beta2 <- mean(beta2_samples)
# beta3 <- mean(beta3_samples)

# Compute 95% CI for each beta
beta_df <- data.frame(
  value = c(beta0_samples, beta1_samples, beta2_samples), # , beta3_samples 
  parameter = rep(c("beta0", "beta1", "beta2"), each = length(beta0_samples)) # , "beta3" 
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
saveRDS(beta_df, "./Data/Model_Data/PC_HDS_beta_df.rds")


# -------------------------------------------------------
# Covariate Effects
# -------------------------------------------------------

# Covariate names
print(colnames(X.abund))


# Set covariate name 
Cov1_name <- "herb_ClmIdx"
Cov2_name <- "woody_AggInx"

# Create a prediction of covariate values
cov1_pred_vals <- seq(min(X.abund[, Cov1_name]), max(X.abund[, Cov1_name]), length.out = 1000)
cov2_pred_vals <- seq(min(X.abund[, Cov2_name]), max(X.abund[, Cov2_name]), length.out = 1000)

# Mean scaling covariates
cov1_scaled <- (X.abund[, Cov1_name] - mean(X.abund[, Cov1_name])) / (max(X.abund[, Cov1_name]) - min(X.abund[, Cov1_name]))
cov2_scaled <- (X.abund[, Cov2_name] - mean(X.abund[, Cov2_name])) / (max(X.abund[, Cov2_name]) - min(X.abund[, Cov2_name]))

# Matrices for storing predictions
cov1_preds <- matrix(NA, nrow = length(beta0_samples), ncol = length(cov1_scaled))
cov2_preds <- matrix(NA, nrow = length(beta0_samples), ncol = length(cov2_scaled))

# # Create a meshgrid of covariate values for interaction predictions
# interaction_grid <- expand.grid(cov1_scaled = cov1_scaled, cov2_scaled = cov2_scaled)
# interaction_grid$interaction_term <- interaction_grid$cov1_scaled * interaction_grid$cov2_scaled
# 
# # Initialize matrix to store predictions
# interaction_preds <- matrix(NA, nrow = length(beta0_samples), ncol = nrow(interaction_grid))
# 

# Generate predictions
for (i in 1:length(beta0_samples)) {
  cov1_preds[i, ] <- beta0_samples[i] + beta1_samples[i] * cov1_scaled # Linear
  cov2_preds[i, ] <- beta0_samples[i] + beta2_samples[i] * cov2_scaled
  # 
  # interaction_preds[i, ] <- beta0_samples[i] +                       # Interaction
  #   beta1_samples[i] * interaction_grid$cov1_scaled +
  #   beta2_samples[i] * interaction_grid$cov2_scaled +
  #   beta3_samples[i] * interaction_grid$interaction_term

}
  

# Calculate mean predictions
cov1_preds_mean <- apply(cov1_preds, 2, mean)# mean
cov1_preds_LCI <- apply(cov1_preds, 2, quantile, probs = 0.025) # LCI
cov1_preds_HCI <- apply(cov1_preds, 2, quantile, probs = 0.975) # HCI

cov2_preds_mean <- apply(cov2_preds, 2, mean)
cov2_preds_LCI <- apply(cov2_preds, 2, quantile, probs = 0.025)
cov2_preds_HCI <- apply(cov2_preds, 2, quantile, probs = 0.975)

# interaction_preds_mean <- apply(interaction_preds, 2, mean)
# interaction_preds_LCI <- apply(interaction_preds, 2, quantile, probs = 0.025)
# interaction_preds_HCI <- apply(interaction_preds, 2, quantile, probs = 0.975)

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
# 
# interaction_pred_df <- data.frame(
#   interaction_term = interaction_grid$interaction_term,
#   cov1_scaled = interaction_grid$cov1_scaled,
#   cov2_scaled = interaction_grid$cov2_scaled,
#   interaction_preds_mean = interaction_preds_mean,
#   interaction_preds_LCI = interaction_preds_LCI,
#   interaction_preds_HCI = interaction_preds_HCI)
# 

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


# # Interactive effects
# plot_ly(interaction_pred_df, 
#         x = ~cov1_scaled, 
#         y = ~cov2_scaled, 
#         z = ~interaction_preds_mean, 
#         type = "contour",
#         colorscale = "Viridis",
#         ncontours = 80,
#         contours = list(showlabels = FALSE), 
#         colorbar = list(title = "Interaction Effect", titlefont = list(size = 14)) 
#         )%>%
#           layout(
#             title = list(
#               text = paste0(model_name, " | Interaction Effect "),
#               font = list(size = 18),  
#               x = 0.5,   
#               xanchor = "center"   
#             ),
#             xaxis = list(
#               title = list(text = "Herbaceous Clumpy Index (scaled)", 
#                            font = list(size = 14))),
#             yaxis = list(
#               title = list(text = "Woody Number of Patches (scaled)", 
#                            font = list(size = 14))))
# 
#         

# -------------------------------------------------------
#   Estimating Abundance 
# -------------------------------------------------------
 
# Extracting Abundance
Ntot_samples <- combined_chains[ ,"N_tot"]  

# Ntotal is the abundance based on 10 point counts at a radius of 200m.
# To correct for density, Ntotal needs to be divided by 10 * area surveyed
area <- pi * (200^2) / 4046.86  # Area in acres
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
dens_df <- dens_df[dens_df$Density >= dens_summary$Lower_CI 
                   & dens_df$Density <= dens_summary$Upper_CI, ]


# Getting total abundance
abund_summary <- dens_summary
abund_summary[,2:4] <- abund_summary[,2:4] * 2710

# Plot Abundance - Violin
abund_df <- dens_df
abund_df$Density <- abund_df$Density * 2710

ggplot(abund_df, aes(x = Model, y = Density, fill = Model)) + 
  geom_violin(trim = FALSE, alpha = 0.6, adjust = 5) +  
  labs(x = "Model", y = "Density (N/acre)") +
  scale_fill_manual(values = c("PC HDS" = "purple")) +   
  scale_y_continuous(limits = c(0, 1000),
                     breaks = seq(0, 1000, by = 100),
                     labels = scales::comma) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),  
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
        panel.grid = element_blank(),
        legend.position = "none")  # Removes legend


# Total abundance
abund_summary




# Export density dataframe
saveRDS(dens_df, "./Data/Model_Data/PC_HDS_dens_df.rds")
saveRDS(dens_summary, "./Data/Model_Data/PC_HDS_dens_summary.rds")
saveRDS(abund_summary, "./Data/Model_Data/PC_HDS_abund_summary.rds")

# Save Environment
#save.image(file = "./HDS_JAGs.RData")

# End Script