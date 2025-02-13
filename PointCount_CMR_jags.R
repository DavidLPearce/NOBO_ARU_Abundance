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
workers <- 2 # Ncores  * 0.6 # For low background use 80%, for medium use 50% of Ncores

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

# Remove rows with NA values
pc_dat <- na.omit(pc_dat)

# Creating a day of year column
pc_dat$Date <- mdy(pc_dat$Date)
pc_dat$DOY <- yday(pc_dat$Date)
pc_dat$DOYscaled <- ((pc_dat$DOY - 1) / 365)


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

# Number of time intervals
J <- 4

# number of sites
nsites <- 10

# Number of individuals detected
nind <- nrow(y)

# Superpopulation
print(nind)
M <- nind + 100 # ensure that M is larger than number detected
print(M)

# Data Augmentation
y <- rbind(y, matrix(0, nrow=(M-nind), ncol=4))

# Number of "pseudo-groups" added
nz <- M - nind

# Make site (Point_Individual) into an integer
site <- as.numeric(factor(pc_CMR$PointNum))

# Add in NA for M
site <- c(site, rep(NA, M-nind))

# Take a look
print(site)

# Extract and scale site covariates for X.abund
X.abund <- site_covs[,-c(1:4,25)] 
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




# Extract and scale site covariates for X.det
X.det <- pc_dat[,c(2,7:9,18)]
X.det$Observer <- as.numeric(as.factor(X.det$Observer))
X.det$Temp.deg.F <- as.numeric(X.det$Temp.deg.F)
X.det$Wind.Beau.Code <- as.numeric(as.factor(X.det$Wind.Beau.Code)) 
X.det$Sky.Beau.Code <- as.numeric(as.factor(X.det$Sky.Beau.Code))
merged_df <- merge(site_covs, pc_CMR, by = "PointNum") # VegDens to stacked 
X.det$vegDens = as.matrix(merged_df[,c("vegDens50m")]) 
print(X.det)


# Bundle data for JAGs
data <- list(J = J,
             M = M,
             nsites = nsites,
             y = y,
             X.abund = X.abund,
             X.det = X.det,
             group = site)
 
# Check structure 
str(data)

 
# ------------------------------------------------------------------------------
#
#                             Model Fitting
#
# ------------------------------------------------------------------------------

# -------------------------------------------------------
# Model Specifications
# -------------------------------------------------------

## Distributions
# library(ggfortify)
# ggdistribution(dunif, seq(0, 1, 0.001), min = 0, max = 1) # p0
# ggdistribution(dnorm, seq(0, 0.01, 0.0001)) # alpha's and beta's 
# ggdistribution(dgamma, seq(0, 5, 0.01), shape = 0.1, rate = 0.1) # tau

# Parameters monitored
params <- c("lambda",
            "N",
            "N_tot",
            "p0", 
            "alpha0", 
            "alpha1", 
            "beta0", 
            "beta1",
            "beta2",
            "beta3", 
            "psi",
            "S.raneff",
            "tau",
             "p_Bayes",
             "log_lik")
            


# Initial values
 inits  <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        beta3 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
)}

# MCMC 
n.iter = 300000
n.burnin = 20000
n.chains = 3 
n.thin = 10


# -------------------------------------------------------
# Model 33: herb_Pdens + woody_lrgPInx 
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0)) 
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  beta3 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }

  # Abundance model 
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1*X.abund[s,17] + beta2*X.abund[s,10]
    probs[s] <- lambda[s] / sum(lambda[])
    N[s] <- sum(group[] == s)  # Estimated abundance at site s
  }


  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Detection model = intercept + wind + site random effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    
    # Log-Likelihood
    log_lik[i,j] <- logdensity.bern(y[i,j], pz[i,j])
    }
  }
  
  # Total estimated population
  N_tot <- sum(z[])
  
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_mod33.txt")
# ------------End Model-------------

# Fit Model
fm.33 <- jags(data = data, 
              parameters.to.save = params,
              inits = inits, 
              model.file = "./jags_models/CMR_mod33.txt",
              n.iter = n.iter,
              n.burnin = n.burnin,
              n.chains = n.chains, 
              n.thin = n.thin,
              parallel = TRUE,
              n.cores = workers,
              DIC = FALSE)  

# Save Environment
save.image(file = "./CMR_bm_JAGs.RData")

# Trace plots
mcmcplot(fm.33$samples)

# Rhat
fm.33$Rhat

# Model summary
print(fm.33, digits = 2)

# Bayes p-value
cat("Bayesian p-value =", fm.33$summary["p_Bayes",1], "\n") 


# WAIC
log_lik_array <- fm.33$sims.list$log_lik  # Extract log-likelihood samples
log_lik_matrix <- apply(log_lik_array, c(1, 2), sum)  # Summing across J
waic_result <- loo::waic(log_lik_matrix)
print(waic_result)

# # Leave One Out
# loo_result <- loo::loo(log_lik_matrix)
# print(loo_result)



# -------------------------------------------------------
#
#   Estimating Abundance 
#
# -------------------------------------------------------

# Combine chains
combined_chains <- as.mcmc(do.call(rbind, fm.33$samples))

# Extract lambda estimates
N_columns <- grep("^N\\[", colnames(combined_chains))
N_samples <- combined_chains[ ,N_columns]

# mean abundance 
N_tot <- rowSums(N_samples)

# Area in hectares
# area <- pi*(200^2)/10000

# Area in acres
area <- pi*(200^2)/4046.86

# Getting density
dens_df <- as.data.frame(N_tot/area)

# Summarize by row
colnames(dens_df)[1] <- "Density"
dens_df[,2] <- "PC CMR"
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
    y = "Density (N/acre)") +
  scale_y_continuous(limits = c(0, 5),
                     breaks = seq(0, 5, by = 1),
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
save.image(file = "./CMR_bm_JAGs.RData")

# Export density dataframe
saveRDS(dens_df, "./Data/Fitted_Models/PC_CMR_DensityDF.rds")


# End Script