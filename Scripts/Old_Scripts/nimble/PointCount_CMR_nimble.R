# -------------------------------------------------------
#
#                    Load libraries
#
# -------------------------------------------------------

# Install packages (if needed)
# install.packages("tidyverse")
# install.packages("nimble")
# install.packages("coda")
# install.packages("mcmcplots")
# install.packages("parallel")

# Load library
library(tidyverse)
library(nimble)
library(coda)
library(mcmcplots)
library(parallel)

# Set seed, scientific notation options, and working directory
set.seed(123)
options(scipen = 9999)
setwd(".")

# Setting up cores
Ncores <- detectCores()
print(Ncores) # number of available cores
workers <- Ncores * 0.5 # may need to change, for low background use 80% num, for medium use 50%


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




                     

# Bundle data for nimble
data <- list(y = y, 
             X.abund = X.abund,
             X.det = X.det,
             group = site)

# Take a look
str(data)

# State Constants
constants <- list(J = 4, 
                  M = 400,
                  nsites = 10)

# Take a look
str(constants)


# ------------------------------------------------------------------------------ 
#
#                             Model Fitting
#
# ------------------------------------------------------------------------------ 


# Parameters monitored
params <- c("lambda",
            "p0", 
            "alpha0", 
            "alpha1", 
            "beta0", 
            "beta1",
            "beta2",
            "beta3",
            "S.raneff",
            "tau",
            "psi")


# Initial values
# inits <- function() {
#   list (p0 = runif(1),
#         alpha1 = runif(1),
#         beta0=runif(1),
#         beta1=rnorm(1),
#         beta2=rnorm(1),
#         beta3=rnorm(1),
#         z = c( rep(1,nind), rep(0, (M-nind))),
#         # To avoid nimble warnings, initialize unobserved groups
#         group = c(rep(NA,length(pc_dat$UniqueID)),
#                   sample(1:nrow(data$X),
#                          size = length(site) - length(pc_dat$UniqueID),
#                          replace = TRUE)),
#         lambda = runif(10, min = 0.1, max = 10)
#   )}# end inits

inits <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        beta1 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits


# -------------------------------------------------------
# Fit Model 0: Null Model
# -------------------------------------------------------
model.0 <- nimbleCode({

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0 / (1-p0)) # same as logit(p0)
  # alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0,0.01)
  # beta1 ~ dnorm(0,0.01)
  # beta2 ~ dnorm(0,0.01)
  # beta3 ~ dnorm(0,0.01)
  psi <- sum(lambda[1:nsites]) / M   # psi is a derived parameter
  
  # Precision for survey-level random effect
  # tau ~ dgamma(0.1, 0.1)  
  # sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # log-linear model for abundance
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 
    probs[s] <- lambda[s] / sum(lambda[1:nsites])
  }
  
  # Site-level random effect
  # for(s in 1:nsites){  
  #   S.raneff[s] ~ dnorm(0, tau)  
  # }

  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[1:nsites])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables

    # Detection model:  Intercept + Wind + Site Random Effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 #+ alpha1*X.det[group[i],3] #+ S.raneff[group[i]] 
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
    }
  }
})
# ------------End Model-------------


# Fit Model
fm.0 <- nimbleMCMC(code = model.0,
                      data = data,
                      constants = constants,
                      inits = inits,
                      monitors = params,
                      niter = 300000,
                      nburnin = 20000,
                      nchains = 3,
                      thin = 10,
                      samplesAsCodaMCMC = TRUE,
                      WAIC = TRUE)



# Rhat
#coda::gelman.diag(fm.covs$samples)

# Trace plots
mcmcplot(fm.covs$samples)

# Model Summary
summary(fm.covs$samples)

 