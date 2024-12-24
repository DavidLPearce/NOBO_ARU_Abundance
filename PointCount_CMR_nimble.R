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
                "lambda")

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
                  niter = 1000,
                  nburnin = 100,
                  nchains = 3,
                  thin = 5,
                  samplesAsCodaMCMC = TRUE)


# Trace plots
mcmcplot(fm.covs)

# Model Summary
summary(fm.covs)



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
              "sigma")


# Initial values: add tau
inits.h <- function() {
  list (p0 = runif(1),
        alpha1 = runif(1),
        beta0=runif(1),
        beta1=rnorm(1),
        beta2=rnorm(1),
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
                   niter = 1000,
                   nburnin = 100, 
                   nchains = 3,
                   thin = 5,
                   samplesAsCodaMCMC  = TRUE)

# Trace plots
mcmcplot(fm.h)

# Model Summary
summary(fm.h)

# Model structure
str(fm.h)


#  -------------------------------------------------------
#
#   Saving Data
#
#  -------------------------------------------------------

# Save environment
save.image(file = "PointCount_CMR_nimble.RData")

#  -------------------------------------------------------
#
#   Estimating Abundance 
#
#  -------------------------------------------------------


# Combine chains
combined_chains <- as.mcmc(do.call(rbind, fm.h))

# Extract lambda estimates (assuming "lambda[1]", "lambda[2]", etc., are in the output)
lambda_columns <- grep("^lambda\\[", colnames(combined_chains))
lambda_samples <- combined_chains[, lambda_columns]

# Summarize abundance for each site
abundance_summary <- apply(lambda_samples, 2, function(x) {
  c(mean = mean(x), sd = sd(x), quantile(x, probs = c(0.025, 0.5, 0.975)))
})
abundance_summary <- t(abundance_summary)

# Print summary
print(abundance_summary)

# Abundance by site
install.packages("reshape2")

library(ggplot2)
library(reshape2)

# Reshape lambda_samples for ggplot
lambda_df <- as.data.frame(lambda_samples)
colnames(lambda_df) <- paste0("Site", 1:ncol(lambda_df))
lambda_long <- melt(lambda_df, variable.name = "Site", value.name = "Abundance")

# Create violin plot
ggplot(lambda_long, aes(x = Site, y = Abundance)) +
  geom_violin(fill = "skyblue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Posterior Distribution of Abundance by Site",
       x = "Site",
       y = "Abundance")


### Getting total abundance
# ------------------------------


# Assuming fm.h is your MCMC result object
# Check if site-specific lambdas exist
lambda_sites <- sapply(fm.h, function(chain) {
  chain_lambda <- as.matrix(chain[, grep("^lambda\\[", colnames(chain))])
  rowSums(chain_lambda, na.rm = TRUE)
})

# Compute total abundance
total_abundance <- rowMeans(lambda_sites)

# Compute 95% confidence intervals
ci_2.5 <- apply(lambda_sites, 1, quantile, probs = 0.025)
ci_97.5 <- apply(lambda_sites, 1, quantile, probs = 0.975)

# Summary
abundance_summary <- data.frame(
  mean = total_abundance,
  sd = apply(lambda_sites, 1, sd),
  `2.5%` = ci_2.5,
  `50%` = apply(lambda_sites, 1, median),
  `97.5%` = ci_97.5
)

print(abundance_summary)


fm.h_abund <- colMeans(abundance_summary)

# Area = π * r^2
# Point counts were assumed to have a effective sampling radius of 
# 200 meters which is 656.2 feet (1 meter = 3.281 feet).
# So the Area =  pi * 656.2 ^2 = 1352765 feet squared.
# Converting Area from ft^2 to acres is 1352765 ft^2 / 43560 ft^2 per acre
# Area = 31.05521 acres
area = 31.05521


# So the mean density per acre is 
mean_density <- fm.h_abund / area
print(mean_density)

# The study area has a acreage of 2710 acres
study_area = 2710

# The abundance across the study area is 
study_area_abund = mean_density * study_area
print(study_area_abund)


# Creating a matrix of latent density
lat_dens_mat <- (abundance_summary) / area
head(lat_dens_mat)

# Matrix is density per point.
# stacking each point column to one column for plotting
lat_dens_df <- as.data.frame(lat_dens_mat) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Density")

# The point that the estimate came from doesn't matter since comparison is across models
colnames(lat_dens_df)[1] <- "Model"
lat_dens_df[,1] <- "PC ind H"
head(lat_dens_df)

# Plot
ggplot(lat_dens_df, aes(x = Model, y = Density, fill = Model)) +
  geom_violin() + 
  geom_boxplot(aes(x = Model, y = Density), 
               width = 0.2, position = position_dodge(width = 0.8)) +
  labs(
    title = "Latent Density", 
    x = "Model", 
    y = "Density") +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2), labels = scales::comma) + # Customize y-axis
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1), # Tilt x-axis text
    axis.text.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none" 
  )






