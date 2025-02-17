
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
workers <- 2 #Ncores * 0.60 # For low background use 80%, for medium use 50% of Ncores
print(workers)

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
  doy_mat[point_num, occasion] <-  pc_dat$DOY[i]
  
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
X.det[, , "DOY"] <- as.matrix(scale(doy_mat)) # Scaled
print(X.det)
X.det[, , 1]  # Observer
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

# Create a 3D array
# Length of site (10), width of distance bins (4), depth of surveys (4)
y3d <- array(NA,dim=c(nrow(det_mat), 4, 4) ) 

# Fill array
y3d[,,1] <- det_mat[,1:4]    
y3d[,,2] <- det_mat[,5:8]  
y3d[,,3] <- det_mat[,9:12]   
y3d[,,4] <- det_mat[,13:16]

# Constances 
K <- 4                          # Number of primary occasions
nsites <- nrow(det_mat)         # Number of sites
nD <- 4                         # Number of distance classes
delta <- 50                     # Class width
B <- 200                        # Maximum distance
midpt <- seq(delta/2, B, delta) # Class midpoint distance
nobs <- apply(y3d, c(1,3), sum) # Total detections per site and occasion
area <- pi*(200^2)/4046.86      # Area in acres


# Bundle data for nimble
data <- list(y3d = y3d, 
             nsites = nsites, 
             K = K, 
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

n.iter = 300000
n.burnin = 20000
n.chains = 3 
n.thin = 10


# ----------------------------------------------------------
#                   Model 1 
# Availability =  ran eff of survey/day
# Detection =  Wind
# Abundance = Herbaceous Patch Density
# ----------------------------------------------------------

# Parameters monitored
params <- c("r",
            "sigma0",
            "theta", 
            "phi0", 
            "beta0", 
            "beta1",
            "gamma1", 
            "logit.gamma1",
            "gamma2",
            "alpha1",
            "lambda",
            "N",
            "N_tot",
            "Davail",
            "log_lik",
            "p_Bayes")



# Initial values
inits  <- function() {
  list(
    M = apply(y3d, 1, max) + 5,
    Navail = apply(y3d, c(1, 3), sum),
    sigma0 = 200,
    gamma1 = rep(0.5, 4),
    gamma2 = 0,
    beta0 = 0,
    beta1 = 0,
    beta2 = 0,
    alpha1 = 0,
    phi0 = 0.5,
    theta = 1,
    r = 5
  )
}


# ----------------------------- 
# Model 1 Statement 
# ----------------------------- 
cat("
model {

  # Priors
  beta0 ~ dnorm(0, 1)
  beta1 ~ dnorm(0, 1)

  # Availability parameters
  phi0 ~ dunif(0.1, 0.9)
  logit.phi0 <- log(phi0/(1-phi0))
  gamma2 ~ dnorm(0, 0.01)
  
  for(k in 1:K){
    gamma1[k] ~ dunif(0.1, 0.9) # Availability effects of surveys 1 - 4
    logit.gamma1[k] <- log(gamma1[k]/(1-gamma1[k]))
  }

  # Detection parameters
  sigma0 ~ dunif(0.1,200)   # Intercept
  theta ~ dgamma(0.1, 0.1)
  r ~ dunif(0, 10)
  alpha1 ~ dnorm(0, 0.01)

  for (s in 1:nsites) {
    for (k in 1:K) {

      # Availability Model
      logit.phi[s,k] <- logit.gamma1[k]
      phi[s,k] <- exp(logit.phi[s,k]) / (1 + exp(logit.phi[s,k]))

      # Distance Sampling
      log(sigma[s,k]) <- log(sigma0) + alpha1*X.det[s,k,3]
      
      # Multinomial cell probability construction
      for(b in 1:nD){
        # Half-normal or hazard rate detection functions
        log(g[s,b,k]) <- -midpt[b]*midpt[b]/(2*sigma[s,k]*sigma[s,k]) # half-normal
        #cloglog(g[s,b,k]) <- theta*log(sigma[s,k])  - theta*log(midpt[b])  # hazard
        
        # Density function for distance bins
        f[s,b,k] <- (2 * midpt[b] * delta) / (B * B)
        cellprobs[s,b,k] <- g[s,b,k] * f[s,b,k]
        cellprobs.cond[s,b,k] <- cellprobs[s,b,k] / sum(cellprobs[s,1:nD,k])
      }
      
      # Add probability of undetected individuals
      cellprobs[s,nD+1,k] <- 1 - sum(cellprobs[s,1:nD,k])

      # Detection probabilities
      pdet[s,k] <- sum(cellprobs[s,1:nD,k])
      pmarg[s,k] <- pdet[s,k] * phi[s,k]

      # Observation model (Multinomial likelihood)
      y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[s,1:nD,k], nobs[s,k])

      # Number of detected individuals
      nobs[s,k] ~ dbin(pmarg[s,k], M[s])

      # Number of available individuals
      Navail[s,k] ~ dbin(phi[s,k], M[s])

      # Log-Likelihood Calculation
      log_lik[s,k] <- logdensity.multi(y3d[s,1:nD,k], cellprobs.cond[s,1:nD,k], nobs[s,k])

      # Posterior Predictive Checks (Bayesian p-value)
      y3d_rep[s,1:nD,k] ~ dmulti(cellprobs.cond[s,1:nD,k], nobs[s,k])

      for (b in 1:nD) {
        discrepancy_obs[s,b,k] <- pow(y3d[s,b,k] - (cellprobs.cond[s,b,k] * nobs[s,k]), 2)
        discrepancy_rep[s,b,k] <- pow(y3d_rep[s,b,k] - (cellprobs.cond[s,b,k] * nobs[s,k]), 2)
      }
      
    } # End k loop

    # Abundance Model
    log(lambda[s]) <- beta0 + beta1*X.abund[s,17] 

    # Population size follows a negative binomial distribution
    # M[s] ~ dnegbin(prob[s], r)
    # prob[s] <- r / (r + lambda[s])
    
    # Poisson
    M[s] ~ dpois(lambda[s])
  } # End s loop

  # Derived Quantities
  for (k in 1:K){
    Davail[k] <- mean(phi[,k]) * exp(beta0) / area
  }
  
  for (s in 1:nsites) {
      for (k in 1:K) {
        N_site_k[s,k] <- pdet[s,k] * phi[s,k] * M[s]    
      }
    N[s] <- sum(N_site_k[s,])
  }
  # Abundance and Density
  N_tot <- sum(N[])

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
fm.1$Rhat # Rhat: less than 1.1 means good convergence
mcmcplot(fm.1$samples)# Visually inspect trace plots
cat("Bayesian p-value =", fm.1$summary["p_Bayes",1], "\n")# Best model fit. P-value = 0.5 means good fit, = 1 or 0 is a poor fit


# Model summary
summary(fm.1$samples)

# Save Environment
save.image(file = "./HDS_JAGs.RData")


# -------------------------------------------------------
#
#   Estimating Abundance 
#
# -------------------------------------------------------
 
#load(file = "./Data/Model_Environments/HDS_JAGs.RData")


# Combine chains
combined_chains <- as.mcmc(do.call(rbind, fm.1$samples))


# Extracting Abundance
Ntot_samples <- combined_chains[ ,"N_tot"]  

# Ntotal is the abundance based on 10 point counts at a radius of 200m.
# To correct for density, Ntotal needs to be divided by 10 * area surveyed
area <- pi * (200^2) / 4046.86  # Area in acres
dens_samples <- Ntot_samples / (area * 10)

# Create data frame for density
dens_df <- data.frame(Model = rep("PC HDS", length(dens_samples)), Density = dens_samples)
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
  scale_y_continuous(limits = c(400,800),
                     breaks = seq(400, 800, by = 25),
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




# Export density dataframe
saveRDS(dens_df, "./Data/Fitted_Models/PC_HDS_dens_df.rds")
saveRDS(dens_summary, "./Data/Fitted_Models/PC_HDS_dens_summary.rds")
saveRDS(abund_summary, "./Data/Fitted_Models/PC_HDS_abund_summary.rds")

# Save Environment
save.image(file = "./HDS_JAGs.RData")
