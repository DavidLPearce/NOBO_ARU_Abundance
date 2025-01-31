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

# Number of individuals detected
nind <- nrow(y)

# Superpopulation 
M <- 400

# Data Augmentation
y <- rbind(y, matrix(0, nrow=(M-nind), ncol=4))

# Make site (Point_Individual) into an integer
site <- as.numeric(factor(pc_CMR$PointNum))

# Add in NA for M
site <- c(site, rep(NA, M-nind))

# Take a look
print(site)

# Extract observation covariates
sitecovs = scale(site_covs[,c("woody_prp", 
                              "herb_prp", 
                              "woody_mean_p_Area")])
                     

# Bundle data for nimble
data <- list(y = y, 
             X = sitecovs,
             group = site)

# Take a look
str(data)

# State Constants
constants <- list(J = 4, 
                  M = 400,
                  nsites = 10)

# Take a look
str(constants)


# -------------------------------------------------------
#
#                   Model Fitting
#
# -------------------------------------------------------


# -------------------------------------------------------
#                 Covariate Model
# -------------------------------------------------------
model.cov <- nimbleCode({

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0 / (1-p0))    # same as logit(p0)
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0,0.01)
  beta1 ~ dnorm(0,0.01)
  beta2 ~ dnorm(0,0.01)
  psi <- sum(lambda[1:nsites]) / M   # psi is a derived parameter

  # log-linear model for abundance: Herbaceous Proportion + Woody Mean Patch Area
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1*X[s,2] + beta2*X[s,3]
    probs[s] <- lambda[s] / sum(lambda[1:nsites])
  }

  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[1:nsites])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables

    # Observation model: p depends on woody proportion
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1*X[group[i],1]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
    }
  }
  # Derived quantities
  Lam_mean <- mean(lambda[1:nsites])
}
) # End model statement
# -------------------------------------------------------


# Parameters monitored
params.cov <- c("lambda",
                "p0", 
                "alpha0", 
                "alpha1", 
                "beta0", 
                "beta1",
                "beta2",
                "psi", 
                "lambda",
                "Lam_mean")

# Initial values
inits.cov <- function() {
  list (p0 = runif(1),
        alpha1 = runif(1),
        beta0=runif(1),
        beta1=rnorm(1),
        beta2=rnorm(1),
        z = c( rep(1,nind), rep(0, (M-nind))),
        # To avoid nimble warnings, initialize unobserved groups
        group = c(rep(NA,length(pc_dat$UniqueID)),
                  sample(1:nrow(data$X),
                         size = length(site) - length(pc_dat$UniqueID), 
                         replace = TRUE)),
        lambda = runif(10, min = 0.1, max = 10)
)}# end inits


# Fit Model
fm.covs <- nimbleMCMC(code = model.cov,
                  data = data,
                  constants = constants,
                  inits = inits.cov,
                  monitors = params.cov,
                  niter = 500000,
                  nburnin = 50000,
                  nchains = 3,
                  thin = 10,
                  samplesAsCodaMCMC = TRUE,
                  WAIC = TRUE)

# Save fitted model
saveRDS(fm.covs, "./Data/Fitted_Models/PC_CMR_fmcovs.rds")

# Rhat
#coda::gelman.diag(fm.covs$samples)

# Trace plots
mcmcplot(fm.covs$samples)

# Model Summary
summary(fm.covs$samples)



# -------------------------------------------------------
#               Individual Heterogeneity Model
# -------------------------------------------------------
model.H <- nimbleCode( {
  
  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0 / (1-p0))    # same as logit(p0)
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0,0.01)
  beta1 ~ dnorm(0,0.01)
  beta2 ~ dnorm(0,0.01)
  psi <- sum(lambda[1:nsites])/M
  tau ~ dgamma(0.1,0.1)
  sigma <- 1/sqrt(tau)
  
  # log-linear model for abundance: Herbaceous Proportion + Woody Mean Patch Area
  for(s in 1:nsites){
    log(lambda[s])<- beta0 + beta1*X[s,2] + beta2*X[s,3]
    probs[s]<- lambda[s]/sum(lambda[1:nsites])
  }
  
  # Model for individual encounter histories
  for(i in 1:M){
    eta[i] ~ dnorm(alpha0, tau)  # Individual random effect
    group[i] ~ dcat(probs[1:nsites])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Observation model: Woody Proportion + ind. heterogeneity 
    for(j in 1:J){
      logit(p[i,j]) <-  alpha1*X[group[i],1] + eta[i]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
    }   
  }
  # Derived quantities
  Lam_mean <- mean(lambda[1:nsites])
}
) # End model statement
# -------------------------------------------------------


# Parameters monitored
params.h <- c("lambda",
              "p0", 
              "alpha0", 
              "alpha1",
              "beta0", 
              "beta1",
              "beta2",
              "psi", 
              "sigma",
              "Lam_mean")


# Initial values: add tau
inits.h <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind))),
        group = c(rep(NA,length(pc_dat$UniqueID)),
                  sample(1:nrow(data$X), 
                         size = length(site) - length(pc_dat$UniqueID), 
                         replace = TRUE)),
        lambda = runif(10, min = 0.1, max = 10)  # Adding initial values for lambda
)}# end inits



# Fit Model
fm.h <- nimbleMCMC(code = model.H, 
                   data = data, 
                   constants = constants,
                   inits = inits.h,
                   monitors = params.h,
                   niter = 500000,
                   nburnin = 50000,
                   nchains = 3,
                   thin = 10,
                   samplesAsCodaMCMC = TRUE,
                   WAIC = TRUE)

# Save fitted model
saveRDS(fm.h, "./Data/Fitted_Models/PC_CMR_fmH.rds")

# Rhat
#coda::gelman.diag(fm.h$samples)

# Trace plots
mcmcplot(fm.h$samples)

# Model Summary
summary(fm.h$samples)



# ---------------------------------------------------------- 
# Ranking Detection Models using WAIC
# ----------------------------------------------------------

# Extract the WAIC values for each model
waic_values <- c(fm.covs$WAIC$WAIC,
                 fm.h$WAIC$WAIC)

# Naming models
fitnames <- c("fm.covs", 
              "fm.h")


# Combine model names and WAIC values into a data frame for ranking
waic_df <- data.frame(Model = fitnames, WAIC = waic_values)

# Rank models based on WAIC (lower WAIC is better)
waic_df <- waic_df[order(waic_df$WAIC), ]

# Print the ranked models
print(waic_df)

# ModelDetection model 2 is the best detection model
summary(fm.covs$samples)



#  -------------------------------------------------------
#
#   Saving Data
#
#  -------------------------------------------------------

# Save environment
#save.image(file = "PointCount_CMR_nimble.RData")

#  -------------------------------------------------------
#
#   Estimating Abundance 
#
#  -------------------------------------------------------


# Combine chains
combined_chains <- as.mcmc(do.call(rbind, fm.covs$samples))

# Extract lambda estimates (assuming "lambda[1]", "lambda[2]", etc., are in the output)
lambda_columns <- grep("^lambda\\[", colnames(combined_chains))
lambda_samples <- combined_chains[, lambda_columns]

# mean abundance 
lambda_tot <- rowSums(lambda_samples)

# Area in hectares
area <- pi*(200^2)/10000

# Getting density
dens_df <- as.data.frame(lambda_tot/area)

# Summarize by row
colnames(dens_df)[1] <- "Density"
dens_df[,2] <- "PC CMR"
colnames(dens_df)[2] <- "Model"
dens_df <- dens_df[, c("Model", "Density")]# Switch the order of columns
head(dens_df)

# Calculate the 95% Credible Interval
ci_bounds <- quantile(dens_df$Density, probs = c(0.025, 0.975))


# Subset the data frame to 95% CI
dens_df <- subset(dens_df, Density >= ci_bounds[1] & Density <= ci_bounds[2])


# Plot
ggplot(dens_df, aes(x = Model, y = Density, fill = Model)) +
  geom_violin() +
  geom_boxplot(aes(x = Model, y = Density),
               width = 0.2, position = position_dodge(width = 0.8)) +
  labs(
    title = "Latent Density",
    x = "Model",
    y = "Density") +
  scale_y_continuous(limits = c(0, 10),
                     breaks = seq(0, 10, by = 5),
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
mean_dens <- mean(dens_df$Density)
LCI_dens <- min(dens_df$Density)
HCI_dens <- max(dens_df$Density)

print(mean_dens)
print(LCI_dens)
print(HCI_dens)

# total abundance
mean_dens * 1096.698
LCI_dens * 1096.698
HCI_dens * 1096.698


# Export Density data frame
saveRDS(dens_df, "./Data/Fitted_Models/PC_CMR_Dens.rds")


