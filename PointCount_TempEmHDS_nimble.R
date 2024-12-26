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

# All objects created by script
#load("./PointCount_TempEmHDS_nimble.RData")

# -------------------------------------------------------
#
#                   Data Wrangling
#
# -------------------------------------------------------

# Format date
pc_dat$Date <- as.POSIXct(pc_dat$Date, format = "%m/%d/%Y")

# Create a day of year column
pc_dat$DOY <- yday(pc_dat$Date)

# Remove NAs
pc_dat_NAom <- na.omit(pc_dat)

# creating a matrix that is 4 Surveys * 3 Distance bins wide and 10 rows long
det_mat <- matrix(0, nrow = 10, ncol = 12)

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

# Convert Observer to numeric factor levels
Observer_numeric <- matrix(as.numeric(as.factor(obsvr_mat)), 
                           nrow = nrow(obsvr_mat), 
                           ncol = ncol(obsvr_mat))

# Create a 3D array
y3d <- array(NA,dim=c(nrow(det_mat), 3, 4) ) # Length of site (10), width of distance bins (3), depth of surveys (4)

# Fill array
y3d[,,1] <- det_mat[,1:3]    
y3d[,,2] <- det_mat[,4:6]  
y3d[,,3] <- det_mat[,7:9]   
y3d[,,4] <- det_mat[,10:12]

# Constances 
K <- 4                          # Number of primary occasions
nsites <- nrow(det_mat)         # Number of sites
nD <- 3                         # Number of distance classes
midpt <- c(25, 75, 125)         # Class midpoint distance
delta <- 50                     # Class width
B <- 200                        # Maximum distance
nobs <- apply(y3d, c(1,3), sum) # Total detections per site and occasion


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
             HerbPRP = scale(site_covs[,'herb_prp']),
             WoodyPatch = scale(site_covs[,'woody_mean_p_Area']),
             Observer = Observer_numeric,
             Temp = scale(temp_mat),
             Wind = wind_mat,
             Sky = sky_mat,
             DOY = scale(doy_mat)
             )
             
# Look at structure
str(data)



# ---------------------------------------------------------- 
# 
#       Temporary Emigration Hierarcical Distance Model
# 
# ----------------------------------------------------------

# ---------------------------------------------------------- 
#                 Availability
# ----------------------------------------------------------

# ---------------------------------------------------------- 
# Avail Model 0: Null
# ----------------------------------------------------------
availmodel.0 <- nimbleCode({
  
  # Priors
  # Abundance parameters
  beta0 ~ dnorm(0, 0.1)
  beta1 ~ dnorm(0, 0.1)
  beta2 ~ dnorm(0, 0.1)
  
  # Availability parameters
  phi0 ~ dunif(0.1, 0.9)
  logit.phi0 <- log(phi0/(1-phi0))
  for(k in 1:4) {
    gamma1[k] ~ dunif(0.1, 0.9)
    logit.gamma1[k]<- log(gamma1[k]/(1-gamma1[k])) 
  }

  # Detection parameters
  sigma0 ~ dunif(0.1, 50)  
  theta ~ dgamma(0.1, 0.1)
  r ~ dunif(0, 10)
  
  for (s in 1:nsites) {
    for (k in 1:K) {
      # Availability parameter
      logit.phi[s,k] <- logit.gamma1[k]
      phi[s,k] <- exp(logit.phi[s,k])/(1+ exp(logit.phi[s,k]))
      
      # Distance sampling parameter
      log(sigma[s,k]) <- log(sigma0)
      
      # Multinomial cell probability construction
      for(b in 1:nD){
        #log(g[s,b,k]) <- -midpt[b]*midpt[b]/(2*sigma[s,k]*sigma[s,k]) # half-normal
        cloglog(g[s,b,k]) <- theta*log(sigma[s,k])  - theta*log(midpt[b])  # hazard
        f[s,b,k] <- (2*midpt[b]*delta)/(B*B)
        cellprobs[s,b,k] <- g[s,b,k]*f[s,b,k]
        cellprobs.cond[s,b,k] <- cellprobs[s,b,k]/sum(cellprobs[s,1:nD,k])
      }
      cellprobs[s,nD+1,k]<- 1-sum(cellprobs[s,1:nD,k])
      
      #  Conditional 4-part hierarchical model
      pdet[s,k] <- sum(cellprobs[s,1:nD,k])
      pmarg[s,k] <- pdet[s,k]*phi[s,k]
      y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[s,1:nD,k], nobs[s,k]) # Part 4
      nobs[s,k] ~ dbin(pmarg[s,k], M[s])  # Part 3: Number of detected individuals
      Navail[s,k] ~ dbin(phi[s,k],M[s])   # Part 2: Number of available individuals
    } # end k loop
    
    M[s] ~ dnegbin(prob[s], r)
    prob[s] <- r/(r+lambda[s])
    # M[s] ~ dpois(lambda[s])         # Part 1: Abundance model
    log(lambda[s]) <- beta0 + beta1*HerbPRP[s] + beta2*WoodyPatch[s]
  }  # end s loop
  
  # Derived quantities
  Lam_mean <- mean(lambda[1:nsites])
  
}) # End model
# ----------------------------------------------------------

# Inits
availinits.0 <- function() {
  list(
    M = apply(y3d, 1, max) + 2,
    Navail = apply(y3d, c(1, 3), sum),
    sigma0 = 50,
    gamma1 = rep(0.5, 4),
    beta0 = 0,
    beta1 = 0,
    beta2 = 0,
    phi0 = 0.5,
    theta = 1,
    r = 5
  )
}

# Parameters to monitor
availparams.0 <- c("r",
                   "sigma0",
                   "beta0", 
                   "beta1", 
                   "beta2", 
                   "theta", 
                   "phi0", 
                   "gamma1", 
                   "logit.gamma1",
                   "lambda",
                   "Lam_mean")


# Run nimble 
availfm.0 <- nimbleMCMC(code = availmodel.0,
                        constants = data,
                        inits = availinits.0,
                        monitors = availparams.0,
                        thin = 10,
                        niter = 500000,
                        nburnin = 50000,
                        nchains = 3,
                        samplesAsCodaMCMC = TRUE,
                        WAIC = TRUE)


# Rhat
#coda::gelman.diag(availfm.0$samples)

# Inspect
mcmcplot(availfm.0$samples)

# Model Summary
#summary(availfm.0$samples)

# WAIC
#print(availfm.0$WAIC$WAIC)



# ---------------------------------------------------------- 
# Avail Model 1: Temperature
# ----------------------------------------------------------
availmodel.1 <- nimbleCode({
  
  # Priors
  # Abundance parameters
  beta0 ~ dnorm(0, 0.1)
  beta1 ~ dnorm(0, 0.1)
  beta2 ~ dnorm(0, 0.1)

  # Availability parameters
  phi0 ~ dunif(0.1, 0.9)
  logit.phi0 <- log(phi0/(1-phi0))
  for(k in 1:4) {
    gamma1[k] ~ dunif(0.1, 0.9)
    logit.gamma1[k]<- log(gamma1[k]/(1-gamma1[k])) 
  }
  gamma2 ~ dnorm(0, 0.01)
  
  # Detection parameters
  sigma0 ~ dunif(0.1, 50)  
  theta ~ dgamma(0.1, 0.1)
  r ~ dunif(0, 10)
  
  for (s in 1:nsites) {
    for (k in 1:K) {
      # Availability parameter
      logit.phi[s,k] <- logit.gamma1[k] + gamma2*Temp[s,k]
      phi[s,k] <- exp(logit.phi[s,k])/(1+ exp(logit.phi[s,k]))
      
      # Distance sampling parameter
      log(sigma[s,k]) <- log(sigma0)
      
      # Multinomial cell probability construction
      for(b in 1:nD){
        #log(g[s,b,k]) <- -midpt[b]*midpt[b]/(2*sigma[s,k]*sigma[s,k]) # half-normal
        cloglog(g[s,b,k]) <- theta*log(sigma[s,k])  - theta*log(midpt[b])  # hazard
        f[s,b,k] <- (2*midpt[b]*delta)/(B*B)
        cellprobs[s,b,k] <- g[s,b,k]*f[s,b,k]
        cellprobs.cond[s,b,k] <- cellprobs[s,b,k]/sum(cellprobs[s,1:nD,k])
      }
      cellprobs[s,nD+1,k]<- 1-sum(cellprobs[s,1:nD,k])
      
      #  Conditional 4-part hierarchical model
      pdet[s,k] <- sum(cellprobs[s,1:nD,k])
      pmarg[s,k] <- pdet[s,k]*phi[s,k]
      y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[s,1:nD,k], nobs[s,k]) # Part 4
      nobs[s,k] ~ dbin(pmarg[s,k], M[s])  # Part 3: Number of detected individuals
      Navail[s,k] ~ dbin(phi[s,k],M[s])   # Part 2: Number of available individuals
    } # end k loop
    
    M[s] ~ dnegbin(prob[s], r)
    prob[s] <- r/(r+lambda[s])
    # M[s] ~ dpois(lambda[s])         # Part 1: Abundance model
    log(lambda[s]) <- beta0 + beta1*HerbPRP[s] + beta2*WoodyPatch[s]
  }  # end s loop
  
  # Derived quantities
  Lam_mean <- mean(lambda[1:nsites])
  
}) # End model
# ----------------------------------------------------------

# Inits
availinits.1 <- function() {
  list(
    M = apply(y3d, 1, max) + 2,
    Navail = apply(y3d, c(1, 3), sum),
    sigma0 = 50,
    gamma1 = rep(0.5, 4),
    gamma2 = 0,
    beta0 = 0,
    beta1 = 0,
    beta2 = 0,
    phi0 = 0.5,
    theta = 1,
    r = 5
  )
}

# Parameters to monitor
availparams.1 <- c("r",
                   "sigma0",
                   "beta0", 
                   "beta1", 
                   "beta2", 
                   "theta", 
                   "phi0", 
                   "gamma1", 
                   "gamma2", 
                   "logit.gamma1",
                   "lambda",
                   "Lam_mean")


# Run nimble 
availfm.1 <- nimbleMCMC(code = availmodel.1,
                        constants = data,
                        inits = availinits.1,
                        monitors = availparams.1,
                        thin = 10,
                        niter = 500000,
                        nburnin = 50000,
                        nchains = 3,
                        samplesAsCodaMCMC = TRUE,
                        WAIC = TRUE)


# Rhat
#coda::gelman.diag(availfm.1$samples)

# Inspect
mcmcplot(availfm.1$samples)

# Model Summary
#summary(availfm.1$samples)

# WAIC
#print(availfm.1$WAIC$WAIC)
                      


# ---------------------------------------------------------- 
# Avail Model 2: Wind
# ----------------------------------------------------------
availmodel.2 <- nimbleCode({
  
  # Priors
  # Abundance parameters
  beta0 ~ dnorm(0, 0.1)
  beta1 ~ dnorm(0, 0.1)
  beta2 ~ dnorm(0, 0.1)
  
  # Availability parameters
  phi0 ~ dunif(0.1, 0.9)
  logit.phi0 <- log(phi0/(1-phi0))
  for(k in 1:4) {
    gamma1[k] ~ dunif(0.1, 0.9)
    logit.gamma1[k]<- log(gamma1[k]/(1-gamma1[k])) 
  }
  gamma2 ~ dnorm(0, 0.01)
  
  # Detection parameters
  sigma0 ~ dunif(0.1, 50)  
  theta ~ dgamma(0.1, 0.1)
  r ~ dunif(0, 10)
  
  for (s in 1:nsites) {
    for (k in 1:K) {
      # Availability parameter
      logit.phi[s,k] <- logit.gamma1[k] + gamma2*Wind[s,k]
      phi[s,k] <- exp(logit.phi[s,k])/(1+ exp(logit.phi[s,k]))
      
      # Distance sampling parameter
      log(sigma[s,k]) <- log(sigma0)
      
      # Multinomial cell probability construction
      for(b in 1:nD){
        #log(g[s,b,k]) <- -midpt[b]*midpt[b]/(2*sigma[s,k]*sigma[s,k]) # half-normal
        cloglog(g[s,b,k]) <- theta*log(sigma[s,k])  - theta*log(midpt[b])  # hazard
        f[s,b,k] <- (2*midpt[b]*delta)/(B*B)
        cellprobs[s,b,k] <- g[s,b,k]*f[s,b,k]
        cellprobs.cond[s,b,k] <- cellprobs[s,b,k]/sum(cellprobs[s,1:nD,k])
      }
      cellprobs[s,nD+1,k]<- 1-sum(cellprobs[s,1:nD,k])
      
      #  Conditional 4-part hierarchical model
      pdet[s,k] <- sum(cellprobs[s,1:nD,k])
      pmarg[s,k] <- pdet[s,k]*phi[s,k]
      y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[s,1:nD,k], nobs[s,k]) # Part 4
      nobs[s,k] ~ dbin(pmarg[s,k], M[s])  # Part 3: Number of detected individuals
      Navail[s,k] ~ dbin(phi[s,k],M[s])   # Part 2: Number of available individuals
    } # end k loop
    
    M[s] ~ dnegbin(prob[s], r)
    prob[s] <- r/(r+lambda[s])
    # M[s] ~ dpois(lambda[s])         # Part 1: Abundance model
    log(lambda[s]) <- beta0 + beta1*HerbPRP[s] + beta2*WoodyPatch[s]
  }  # end s loop
  
  # Derived quantities
  Lam_mean <- mean(lambda[1:nsites])
  
}) # End model
# ----------------------------------------------------------

# Inits
availinits.2 <- function() {
  list(
    M = apply(y3d, 1, max) + 2,
    Navail = apply(y3d, c(1, 3), sum),
    sigma0 = 50,
    gamma1 = rep(0.5, 4),
    gamma2 = 0,
    beta0 = 0,
    beta1 = 0,
    beta2 = 0,
    phi0 = 0.5,
    theta = 1,
    r = 5
  )
}

# Parameters to monitor
availparams.2 <- c("r",
                   "sigma0",
                   "beta0", 
                   "beta1", 
                   "beta2", 
                   "theta", 
                   "phi0", 
                   "gamma1", 
                   "gamma2", 
                   "logit.gamma1",
                   "lambda",
                   "Lam_mean")


# Run nimble 
availfm.2 <- nimbleMCMC(code = availmodel.2,
                        constants = data,
                        inits = availinits.2,
                        monitors = availparams.2,
                        thin = 10,
                        niter = 500000,
                        nburnin = 50000,
                        nchains = 3,
                        samplesAsCodaMCMC = TRUE,
                        WAIC = TRUE)


# Rhat
#coda::gelman.diag(availfm.2$samples)

# Inspect
mcmcplot(availfm.2$samples)

# Model Summary
#summary(availfm.2$samples)

# WAIC
#print(availfm.2$WAIC$WAIC)




# ---------------------------------------------------------- 
# Avail Model 3: Sky
# ----------------------------------------------------------
availmodel.3 <- nimbleCode({
  
  # Priors
  # Abundance parameters
  beta0 ~ dnorm(0, 0.1)
  beta1 ~ dnorm(0, 0.1)
  beta2 ~ dnorm(0, 0.1)
  
  # Availability parameters
  phi0 ~ dunif(0.1, 0.9)
  logit.phi0 <- log(phi0/(1-phi0))
  for(k in 1:4) {
    gamma1[k] ~ dunif(0.1, 0.9)
    logit.gamma1[k]<- log(gamma1[k]/(1-gamma1[k])) 
  }
  gamma2 ~ dnorm(0, 0.01)
  
  # Detection parameters
  sigma0 ~ dunif(0.1, 50)  
  theta ~ dgamma(0.1, 0.1)
  r ~ dunif(0, 10)
  
  for (s in 1:nsites) {
    for (k in 1:K) {
      # Availability parameter
      logit.phi[s,k] <- logit.gamma1[k] + gamma2*Sky[s,k]
      phi[s,k] <- exp(logit.phi[s,k])/(1+ exp(logit.phi[s,k]))
      
      # Distance sampling parameter
      log(sigma[s,k]) <- log(sigma0)
      
      # Multinomial cell probability construction
      for(b in 1:nD){
        #log(g[s,b,k]) <- -midpt[b]*midpt[b]/(2*sigma[s,k]*sigma[s,k]) # half-normal
        cloglog(g[s,b,k]) <- theta*log(sigma[s,k])  - theta*log(midpt[b])  # hazard
        f[s,b,k] <- (2*midpt[b]*delta)/(B*B)
        cellprobs[s,b,k] <- g[s,b,k]*f[s,b,k]
        cellprobs.cond[s,b,k] <- cellprobs[s,b,k]/sum(cellprobs[s,1:nD,k])
      }
      cellprobs[s,nD+1,k]<- 1-sum(cellprobs[s,1:nD,k])
      
      #  Conditional 4-part hierarchical model
      pdet[s,k] <- sum(cellprobs[s,1:nD,k])
      pmarg[s,k] <- pdet[s,k]*phi[s,k]
      y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[s,1:nD,k], nobs[s,k]) # Part 4
      nobs[s,k] ~ dbin(pmarg[s,k], M[s])  # Part 3: Number of detected individuals
      Navail[s,k] ~ dbin(phi[s,k],M[s])   # Part 2: Number of available individuals
    } # end k loop
    
    M[s] ~ dnegbin(prob[s], r)
    prob[s] <- r/(r+lambda[s])
    # M[s] ~ dpois(lambda[s])         # Part 1: Abundance model
    log(lambda[s]) <- beta0 + beta1*HerbPRP[s] + beta2*WoodyPatch[s]
  }  # end s loop
  
  # Derived quantities
  Lam_mean <- mean(lambda[1:nsites])
  
}) # End model
# ----------------------------------------------------------

# Inits
availinits.3 <- function() {
  list(
    M = apply(y3d, 1, max) + 2,
    Navail = apply(y3d, c(1, 3), sum),
    sigma0 = 50,
    gamma1 = rep(0.5, 4),
    gamma2 = 0,
    beta0 = 0,
    beta1 = 0,
    beta2 = 0,
    phi0 = 0.5,
    theta = 1,
    r = 5
  )
}

# Parameters to monitor
availparams.3 <- c("r",
                   "sigma0",
                   "beta0", 
                   "beta1", 
                   "beta2", 
                   "theta", 
                   "phi0", 
                   "gamma1", 
                   "gamma2", 
                   "logit.gamma1",
                   "lambda",
                   "Lam_mean")


# Run nimble 
availfm.3 <- nimbleMCMC(code = availmodel.3,
                        constants = data,
                        inits = availinits.3,
                        monitors = availparams.3,
                        thin = 10,
                        niter = 500000,
                        nburnin = 50000,
                        nchains = 3,
                        samplesAsCodaMCMC = TRUE,
                        WAIC = TRUE)


# Rhat
#coda::gelman.diag(availfm.3$samples)

# Inspect
mcmcplot(availfm.3$samples)

# Model Summary
#summary(availfm.3$samples)

# WAIC
#print(availfm.3$WAIC$WAIC)




# ---------------------------------------------------------- 
# Avail Model 4: DOY
# ----------------------------------------------------------
availmodel.4 <- nimbleCode({
  
  # Priors
  # Abundance parameters
  beta0 ~ dnorm(0, 0.1)
  beta1 ~ dnorm(0, 0.1)
  beta2 ~ dnorm(0, 0.1)
  
  # Availability parameters
  phi0 ~ dunif(0.1, 0.9)
  logit.phi0 <- log(phi0/(1-phi0))
  for(k in 1:4) {
    gamma1[k] ~ dunif(0.1, 0.9)
    logit.gamma1[k]<- log(gamma1[k]/(1-gamma1[k])) 
  }
  gamma2 ~ dnorm(0, 0.01)
  
  # Detection parameters
  sigma0 ~ dunif(0.1, 50)  
  theta ~ dgamma(0.1, 0.1)
  r ~ dunif(0, 10)
  
  for (s in 1:nsites) {
    for (k in 1:K) {
      # Availability parameter
      logit.phi[s,k] <- logit.gamma1[k] + gamma2*DOY[s,k]
      phi[s,k] <- exp(logit.phi[s,k])/(1+ exp(logit.phi[s,k]))
      
      # Distance sampling parameter
      log(sigma[s,k]) <- log(sigma0)
      
      # Multinomial cell probability construction
      for(b in 1:nD){
        #log(g[s,b,k]) <- -midpt[b]*midpt[b]/(2*sigma[s,k]*sigma[s,k]) # half-normal
        cloglog(g[s,b,k]) <- theta*log(sigma[s,k])  - theta*log(midpt[b])  # hazard
        f[s,b,k] <- (2*midpt[b]*delta)/(B*B)
        cellprobs[s,b,k] <- g[s,b,k]*f[s,b,k]
        cellprobs.cond[s,b,k] <- cellprobs[s,b,k]/sum(cellprobs[s,1:nD,k])
      }
      cellprobs[s,nD+1,k]<- 1-sum(cellprobs[s,1:nD,k])
      
      #  Conditional 4-part hierarchical model
      pdet[s,k] <- sum(cellprobs[s,1:nD,k])
      pmarg[s,k] <- pdet[s,k]*phi[s,k]
      y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[s,1:nD,k], nobs[s,k]) # Part 4
      nobs[s,k] ~ dbin(pmarg[s,k], M[s])  # Part 3: Number of detected individuals
      Navail[s,k] ~ dbin(phi[s,k],M[s])   # Part 2: Number of available individuals
    } # end k loop
    
    M[s] ~ dnegbin(prob[s], r)
    prob[s] <- r/(r+lambda[s])
    # M[s] ~ dpois(lambda[s])         # Part 1: Abundance model
    log(lambda[s]) <- beta0 + beta1*HerbPRP[s] + beta2*WoodyPatch[s]
  }  # end s loop
  
  # Derived quantities
  Lam_mean <- mean(lambda[1:nsites])
  

}) # End model
# ----------------------------------------------------------

# Inits
availinits.4 <- function() {
  list(
    M = apply(y3d, 1, max) + 2,
    Navail = apply(y3d, c(1, 3), sum),
    sigma0 = 50,
    gamma1 = rep(0.5, 4),
    gamma2 = 0,
    beta0 = 0,
    beta1 = 0,
    beta2 = 0,
    phi0 = 0.5,
    theta = 1,
    r = 5
  )
}

# Parameters to monitor
availparams.4 <- c("r",
                   "sigma0",
                   "beta0", 
                   "beta1", 
                   "beta2", 
                   "theta", 
                   "phi0", 
                   "gamma1", 
                   "gamma2", 
                   "logit.gamma1",
                   "lambda",
                   "Lam_mean")


# Run nimble 
availfm.4 <- nimbleMCMC(code = availmodel.4,
                        constants = data,
                        inits = availinits.4,
                        monitors = availparams.4,
                        thin = 10,
                        niter = 500000,
                        nburnin = 50000,
                        nchains = 3,
                        samplesAsCodaMCMC = TRUE,
                        WAIC = TRUE)


# Rhat
#coda::gelman.diag(availfm.4$samples)

# Inspect
mcmcplot(availfm.4$samples)

# Model Summary
#summary(availfm.4$samples)

# WAIC
#print(availfm.4$WAIC$WAIC)




# ---------------------------------------------------------- 
# Ranking Availability Models using WAIC
# ----------------------------------------------------------

# Extract the WAIC values for each model
avail_waic_values <- c(availfm.0$WAIC$WAIC,
                       availfm.1$WAIC$WAIC,
                       availfm.2$WAIC$WAIC,
                       availfm.3$WAIC$WAIC,
                       availfm.4$WAIC$WAIC)

# Naming models
avail_fitnames <- c("availfm.0", 
                    "availfm.1", 
                    "availfm.2",
                    "availfm.3", 
                    "availfm.4")


# Combine model names and WAIC values into a data frame for ranking
avail_waic_df <- data.frame(Model = avail_fitnames, WAIC = avail_waic_values)

# Rank models based on WAIC (lower WAIC is better)
avail_waic_df <- avail_waic_df[order(avail_waic_df$WAIC), ]

# Print the ranked models
print(avail_waic_df)


# Model 4: DOY is the most supported





# ----------------------------------------------------------
#                         Detection
# ----------------------------------------------------------


# ---------------------------------------------------------- 
# Det Model 0: Null
# ----------------------------------------------------------
detmod.0 <- nimbleCode({
  
  # Priors
  # Abundance parameters
  beta0 ~ dnorm(0, 0.1)
  beta1 ~ dnorm(0, 0.1)
  beta2 ~ dnorm(0, 0.1)
  
  # Availability parameters
  phi0 ~ dunif(0.1, 0.9)
  logit.phi0 <- log(phi0/(1-phi0))
  for(k in 1:4) {
    gamma1[k] ~ dunif(0.1, 0.9)
    logit.gamma1[k]<- log(gamma1[k]/(1-gamma1[k])) 
  }
  gamma2 ~ dnorm(0, 0.01)
  
  # Detection parameters
  sigma0 ~ dunif(0.1, 50)  
  theta ~ dgamma(0.1, 0.1)
  r ~ dunif(0, 10)
  
  for (s in 1:nsites) {
    for (k in 1:K) {
      # Availability parameter
      logit.phi[s,k] <- logit.gamma1[k] + gamma2*DOY[s,k]
      phi[s,k] <- exp(logit.phi[s,k])/(1+ exp(logit.phi[s,k]))
      
      # Distance sampling parameter
      log(sigma[s,k]) <- log(sigma0)
      
      # Multinomial cell probability construction
      for(b in 1:nD){
        #log(g[s,b,k]) <- -midpt[b]*midpt[b]/(2*sigma[s,k]*sigma[s,k]) # half-normal
        cloglog(g[s,b,k]) <- theta*log(sigma[s,k])  - theta*log(midpt[b])  # hazard
        f[s,b,k] <- (2*midpt[b]*delta)/(B*B)
        cellprobs[s,b,k] <- g[s,b,k]*f[s,b,k]
        cellprobs.cond[s,b,k] <- cellprobs[s,b,k]/sum(cellprobs[s,1:nD,k])
      }
      cellprobs[s,nD+1,k]<- 1-sum(cellprobs[s,1:nD,k])
      
      #  Conditional 4-part hierarchical model
      pdet[s,k] <- sum(cellprobs[s,1:nD,k])
      pmarg[s,k] <- pdet[s,k]*phi[s,k]
      y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[s,1:nD,k], nobs[s,k]) # Part 4
      nobs[s,k] ~ dbin(pmarg[s,k], M[s])  # Part 3: Number of detected individuals
      Navail[s,k] ~ dbin(phi[s,k],M[s])   # Part 2: Number of available individuals
    } # end k loop
    
    M[s] ~ dnegbin(prob[s], r)
    prob[s] <- r/(r+lambda[s])
    # M[s] ~ dpois(lambda[s])         # Part 1: Abundance model
    log(lambda[s]) <- beta0 + beta1*HerbPRP[s] + beta2*WoodyPatch[s]
  }  # end s loop
  
  # Derived quantities
  Lam_mean <- mean(lambda[1:nsites])
  
}) # End model
# ----------------------------------------------------------

# Inits
detinits.0 <- function() {
  list(
    M = apply(y3d, 1, max) + 2,
    Navail = apply(y3d, c(1, 3), sum),
    sigma0 = 50,
    gamma1 = rep(0.5, 4),
    gamma2 = 0,
    beta0 = 0,
    beta1 = 0,
    beta2 = 0,
    phi0 = 0.5,
    theta = 1,
    r = 5
  )
}

# Parameters to monitor
detparams.0 <- c("r",
                   "sigma0",
                   "beta0", 
                   "beta1", 
                   "beta2", 
                   "theta", 
                   "phi0", 
                   "gamma1", 
                   "gamma2", 
                   "logit.gamma1",
                   "lambda",
                 "Lam_mean")


# Run nimble 
detfm.0 <- nimbleMCMC(code = detmod.0,
                        constants = data,
                        inits = detinits.0,
                        monitors = detparams.0,
                        thin = 10,
                        niter = 500000,
                        nburnin = 50000,
                        nchains = 3,
                        samplesAsCodaMCMC = TRUE,
                        WAIC = TRUE)


# Rhat
#coda::gelman.diag(detfm.0$samples)

# Inspect
mcmcplot(detfm.0$samples)

# Model Summary
#summary(detfm.0$samples)

# WAIC
#print(detfm.0$WAIC$WAIC)



# ----------------------------------------------------------
# Det Model 1: Observer
# ----------------------------------------------------------
detmod.1 <- nimbleCode({
  
  # Priors
  # Abundance parameters
  beta0 ~ dnorm(0, 0.1)
  beta1 ~ dnorm(0, 0.1)
  beta2 ~ dnorm(0, 0.1)
  
  # Availability parameters
  phi0 ~ dunif(0.1, 0.9)
  logit.phi0 <- log(phi0/(1-phi0))
  for(k in 1:4) {
    gamma1[k] ~ dunif(0.1, 0.9)
    logit.gamma1[k]<- log(gamma1[k]/(1-gamma1[k]))
  }
  gamma2 ~ dnorm(0, 0.01)
  
  # Detection parameters
  sigma0 ~ dunif(0.1, 50)
  alpha2 ~ dnorm(0, 0.01)
  theta ~ dgamma(0.1, 0.1)
  r ~ dunif(0, 10)
  
  for (s in 1:nsites) {
    for (k in 1:K) {
      # Availability parameter
      logit.phi[s,k] <- logit.gamma1[k] + gamma2*DOY[s,k]
      phi[s,k] <- exp(logit.phi[s,k])/(1+ exp(logit.phi[s,k]))
      
      # Distance sampling parameter
      log(sigma[s,k]) <- log(sigma0) + alpha2*Observer[s,k]
      
      # Multinomial cell probability construction
      for(b in 1:nD){
        #log(g[s,b,k]) <- -midpt[b]*midpt[b]/(2*sigma[s,k]*sigma[s,k]) # half-normal
        cloglog(g[s,b,k]) <- theta*log(sigma[s,k])  - theta*log(midpt[b])  # hazard
        f[s,b,k] <- (2*midpt[b]*delta)/(B*B)
        cellprobs[s,b,k] <- g[s,b,k]*f[s,b,k]
        cellprobs.cond[s,b,k] <- cellprobs[s,b,k]/sum(cellprobs[s,1:nD,k])
      }
      cellprobs[s,nD+1,k]<- 1-sum(cellprobs[s,1:nD,k])
      
      #  Conditional 4-part hierarchical model
      pdet[s,k] <- sum(cellprobs[s,1:nD,k])
      pmarg[s,k] <- pdet[s,k]*phi[s,k]
      y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[s,1:nD,k], nobs[s,k]) # Part 4
      nobs[s,k] ~ dbin(pmarg[s,k], M[s])  # Part 3: Number of detected individuals
      Navail[s,k] ~ dbin(phi[s,k],M[s])   # Part 2: Number of available individuals
    } # end k loop
    
    M[s] ~ dnegbin(prob[s], r)
    prob[s] <- r/(r+lambda[s])
    # M[s] ~ dpois(lambda[s])         # Part 1: Abundance model
    log(lambda[s]) <- beta0 + beta1*HerbPRP[s] + beta2*WoodyPatch[s]
  }  # end s loop
  
  # Derived quantities
  Lam_mean <- mean(lambda[1:nsites])
  
}) # End model
# ----------------------------------------------------------

# Inits
detinits.1 <- function() {
  list(
    M = apply(y3d, 1, max) + 2,
    Navail = apply(y3d, c(1, 3), sum),
    sigma0 = 50,
    alpha2 = 0,
    gamma1 = rep(0.5, 4),
    gamma2 = 0,
    beta0 = 0,
    beta1 = 0,
    beta2 = 0,
    phi0 = 0.5,
    theta = 1,
    r = 5
  )
}

# Parameters to monitor
detparams.1 <- c("r",
                 "sigma0",
                 "beta0",
                 "beta1",
                 "beta2",
                 "alpha2",
                 "theta",
                 "phi0",
                 "gamma1",
                 "gamma2",
                 "logit.gamma1",
                 "lambda",
                 "Lam_mean")


# Run nimble
detfm.1 <- nimbleMCMC(code = detmod.1,
                      constants = data,
                      inits = detinits.1,
                      monitors = detparams.1,
                      thin = 10,
                      niter = 500000,
                      nburnin = 50000,
                      nchains = 3,
                      samplesAsCodaMCMC = TRUE,
                      WAIC = TRUE)
                    

# Rhat
#coda::gelman.diag(detfm.1$samples)

# Inspect
mcmcplot(detfm.1$samples)

# Model Summary
#summary(detfm.1$samples)

# WAIC
#print(detfm.1$WAIC$WAIC)



# ----------------------------------------------------------
# Det Model 2: Temp
# ----------------------------------------------------------
detmod.2 <- nimbleCode({

  # Priors
  # Abundance parameters
  beta0 ~ dnorm(0, 0.1)
  beta1 ~ dnorm(0, 1)
  beta2 ~ dnorm(0, 1)

  # Availability parameters
  phi0 ~ dunif(0.1, 0.9)
  logit.phi0 <- log(phi0/(1-phi0))
  for(k in 1:4) {
    gamma1[k] ~ dunif(0.1, 0.9)
    logit.gamma1[k]<- log(gamma1[k]/(1-gamma1[k]))
  }
  gamma2 ~ dnorm(0, 0.01)

  # Detection parameters
  sigma0 ~ dunif(0.1, 50)
  alpha2 ~ dnorm(0, 0.01)
  theta ~ dgamma(0.1, 0.1)
  r ~ dunif(0, 10)

  for (s in 1:nsites) {
    for (k in 1:K) {
      # Availability parameter
      logit.phi[s,k] <- logit.gamma1[k] + gamma2*DOY[s,k]
      phi[s,k] <- exp(logit.phi[s,k])/(1+ exp(logit.phi[s,k]))

      # Distance sampling parameter
      log(sigma[s,k]) <- log(sigma0) + alpha2*Temp[s,k]

      # Multinomial cell probability construction
      for(b in 1:nD){
        #log(g[s,b,k]) <- -midpt[b]*midpt[b]/(2*sigma[s,k]*sigma[s,k]) # half-normal
        cloglog(g[s,b,k]) <- theta*log(sigma[s,k])  - theta*log(midpt[b])  # hazard
        f[s,b,k] <- (2*midpt[b]*delta)/(B*B)
        cellprobs[s,b,k] <- g[s,b,k]*f[s,b,k]
        cellprobs.cond[s,b,k] <- cellprobs[s,b,k]/sum(cellprobs[s,1:nD,k])
      }
      cellprobs[s,nD+1,k]<- 1-sum(cellprobs[s,1:nD,k])

      #  Conditional 4-part hierarchical model
      pdet[s,k] <- sum(cellprobs[s,1:nD,k])
      pmarg[s,k] <- pdet[s,k]*phi[s,k]
      y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[s,1:nD,k], nobs[s,k]) # Part 4
      nobs[s,k] ~ dbin(pmarg[s,k], M[s])  # Part 3: Number of detected individuals
      Navail[s,k] ~ dbin(phi[s,k],M[s])   # Part 2: Number of available individuals
    } # end k loop

    M[s] ~ dnegbin(prob[s], r)
    prob[s] <- r/(r+lambda[s])
    # M[s] ~ dpois(lambda[s])         # Part 1: Abundance model
    log(lambda[s]) <- beta0 + beta1*HerbPRP[s,1] + beta2*WoodyPatch[s,1]
  }  # end s loop
  
  # Derived quantities
  Lam_mean <- mean(lambda[1:nsites])
  
}) # End model
# ----------------------------------------------------------

# Inits
detinits.2 <- function() {
  list(
    M = apply(y3d, 1, max) + 2,
    Navail = apply(y3d, c(1, 3), sum),
    sigma0 = 50,
    alpha2 = 0,
    gamma1 = rep(0.5, 4),
    gamma2 = 0,
    beta0 = 0,
    beta1 = 0,
    beta2 = 0,
    phi0 = 0.5,
    theta = 1,
    r = 5
  )
}

# Parameters to monitor
detparams.2 <- c("r",
            "sigma0",
            "beta0",
            "beta1",
            "beta2",
            "alpha2",
            "theta",
            "phi0",
            "gamma1",
            "gamma2",
            "logit.gamma1",
            "lambda",
            "Lam_mean")


# Run nimble
detfm.2 <- nimbleMCMC(
                  code = detmod.2,
                  constants = data,
                  inits = detinits.2,
                  monitors = detparams.2,
                  thin = 10,
                  niter = 1000000,
                  nburnin = 50000,
                  nchains = 3,
                  samplesAsCodaMCMC = TRUE,
                  WAIC = TRUE)

# Save model
saveRDS(detfm.2, "./Data/Fitted_Models/PC_TempEmHDS_detfm2.rds")

# Rhat
coda::gelman.diag(detfm.2$samples)

# Inspect
mcmcplot(detfm.2$samples)

# Model Summary
summary(detfm.2$samples)

# WAIC
print(detfm.2$WAIC$WAIC)




# ----------------------------------------------------------
# Det Model 3: Wind
# ----------------------------------------------------------
detmod.3 <- nimbleCode({
  
  # Priors
  # Abundance parameters
  beta0 ~ dnorm(0, 0.1)
  beta1 ~ dnorm(0, 0.1)
  beta2 ~ dnorm(0, 0.1)
  
  # Availability parameters
  phi0 ~ dunif(0.1, 0.9)
  logit.phi0 <- log(phi0/(1-phi0))
  for(k in 1:4) {
    gamma1[k] ~ dunif(0.1, 0.9)
    logit.gamma1[k]<- log(gamma1[k]/(1-gamma1[k]))
  }
  gamma2 ~ dnorm(0, 0.01)
  
  # Detection parameters
  sigma0 ~ dunif(0.1, 50)
  alpha2 ~ dnorm(0, 0.01)
  theta ~ dgamma(0.1, 0.1)
  r ~ dunif(0, 10)
  
  for (s in 1:nsites) {
    for (k in 1:K) {
      # Availability parameter
      logit.phi[s,k] <- logit.gamma1[k] + gamma2*DOY[s,k]
      phi[s,k] <- exp(logit.phi[s,k])/(1+ exp(logit.phi[s,k]))
      
      # Distance sampling parameter
      log(sigma[s,k]) <- log(sigma0) + alpha2*Wind[s,k]
      
      # Multinomial cell probability construction
      for(b in 1:nD){
        #log(g[s,b,k]) <- -midpt[b]*midpt[b]/(2*sigma[s,k]*sigma[s,k]) # half-normal
        cloglog(g[s,b,k]) <- theta*log(sigma[s,k])  - theta*log(midpt[b])  # hazard
        f[s,b,k] <- (2*midpt[b]*delta)/(B*B)
        cellprobs[s,b,k] <- g[s,b,k]*f[s,b,k]
        cellprobs.cond[s,b,k] <- cellprobs[s,b,k]/sum(cellprobs[s,1:nD,k])
      }
      cellprobs[s,nD+1,k]<- 1-sum(cellprobs[s,1:nD,k])
      
      #  Conditional 4-part hierarchical model
      pdet[s,k] <- sum(cellprobs[s,1:nD,k])
      pmarg[s,k] <- pdet[s,k]*phi[s,k]
      y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[s,1:nD,k], nobs[s,k]) # Part 4
      nobs[s,k] ~ dbin(pmarg[s,k], M[s])  # Part 3: Number of detected individuals
      Navail[s,k] ~ dbin(phi[s,k],M[s])   # Part 2: Number of available individuals
    } # end k loop
    
    M[s] ~ dnegbin(prob[s], r)
    prob[s] <- r/(r+lambda[s])
    # M[s] ~ dpois(lambda[s])         # Part 1: Abundance model
    log(lambda[s]) <- beta0 + beta1*HerbPRP[s] + beta2*WoodyPatch[s]
  }  # end s loop
  
  # Derived quantities
  Lam_mean <- mean(lambda[1:nsites])
  
}) # End model
# ----------------------------------------------------------

# Inits
detinits.3 <- function() {
  list(
    M = apply(y3d, 1, max) + 2,
    Navail = apply(y3d, c(1, 3), sum),
    sigma0 = 50,
    alpha2 = 0,
    gamma1 = rep(0.5, 4),
    gamma2 = 0,
    beta0 = 0,
    beta1 = 0,
    beta2 = 0,
    phi0 = 0.5,
    theta = 1,
    r = 5
  )
}

# Parameters to monitor
detparams.3 <- c("r",
                 "sigma0",
                 "beta0",
                 "beta1",
                 "beta2",
                 "alpha2",
                 "theta",
                 "phi0",
                 "gamma1",
                 "gamma2",
                 "logit.gamma1",
                 "lambda",
                 "Lam_mean")


# Run nimble
detfm.3 <- nimbleMCMC(
                      code = detmod.3,
                      constants = data,
                      inits = detinits.3,
                      monitors = detparams.3,
                      thin = 10,
                      niter = 500000,
                      nburnin = 50000,
                      nchains = 3,
                      samplesAsCodaMCMC = TRUE,
                      WAIC = TRUE)

# Rhat
#coda::gelman.diag(detfm.3$samples)

# Inspect
mcmcplot(detfm.3$samples)

# Model Summary
#summary(detfm.3$samples)

# WAIC
#print(detfm.3$WAIC$WAIC)


# ----------------------------------------------------------
# Det Model 4: Sky
# ----------------------------------------------------------
detmod.4 <- nimbleCode({
  
  # Priors
  # Abundance parameters
  beta0 ~ dnorm(0, 0.1)
  beta1 ~ dnorm(0, 0.1)
  beta2 ~ dnorm(0, 0.1)
  
  # Availability parameters
  phi0 ~ dunif(0.1, 0.9)
  logit.phi0 <- log(phi0/(1-phi0))
  for(k in 1:4) {
    gamma1[k] ~ dunif(0.1, 0.9)
    logit.gamma1[k]<- log(gamma1[k]/(1-gamma1[k]))
  }
  gamma2 ~ dnorm(0, 0.01)
  
  # Detection parameters
  sigma0 ~ dunif(0.1, 50)
  alpha2 ~ dnorm(0, 0.01)
  theta ~ dgamma(0.1, 0.1)
  r ~ dunif(0, 10)
  
  for (s in 1:nsites) {
    for (k in 1:K) {
      # Availability parameter
      logit.phi[s,k] <- logit.gamma1[k] + gamma2*DOY[s,k]
      phi[s,k] <- exp(logit.phi[s,k])/(1+ exp(logit.phi[s,k]))
      
      # Distance sampling parameter
      log(sigma[s,k]) <- log(sigma0) + alpha2*Sky[s,k]
      
      # Multinomial cell probability construction
      for(b in 1:nD){
        #log(g[s,b,k]) <- -midpt[b]*midpt[b]/(2*sigma[s,k]*sigma[s,k]) # half-normal
        cloglog(g[s,b,k]) <- theta*log(sigma[s,k])  - theta*log(midpt[b])  # hazard
        f[s,b,k] <- (2*midpt[b]*delta)/(B*B)
        cellprobs[s,b,k] <- g[s,b,k]*f[s,b,k]
        cellprobs.cond[s,b,k] <- cellprobs[s,b,k]/sum(cellprobs[s,1:nD,k])
      }
      cellprobs[s,nD+1,k]<- 1-sum(cellprobs[s,1:nD,k])
      
      #  Conditional 4-part hierarchical model
      pdet[s,k] <- sum(cellprobs[s,1:nD,k])
      pmarg[s,k] <- pdet[s,k]*phi[s,k]
      y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[s,1:nD,k], nobs[s,k]) # Part 4
      nobs[s,k] ~ dbin(pmarg[s,k], M[s])  # Part 3: Number of detected individuals
      Navail[s,k] ~ dbin(phi[s,k],M[s])   # Part 2: Number of available individuals
    } # end k loop
    
    M[s] ~ dnegbin(prob[s], r)
    prob[s] <- r/(r+lambda[s])
    # M[s] ~ dpois(lambda[s])         # Part 1: Abundance model
    log(lambda[s]) <- beta0 + beta1*HerbPRP[s] + beta2*WoodyPatch[s]
  }  # end s loop
  
  # Derived quantities
  Lam_mean <- mean(lambda[1:nsites])
  
}) # End model
# ----------------------------------------------------------

# Inits
detinits.4 <- function() {
  list(
    M = apply(y3d, 1, max) + 2,
    Navail = apply(y3d, c(1, 3), sum),
    sigma0 = 50,
    alpha2 = 0,
    gamma1 = rep(0.5, 4),
    gamma2 = 0,
    beta0 = 0,
    beta1 = 0,
    beta2 = 0,
    phi0 = 0.5,
    theta = 1,
    r = 5
  )
}

# Parameters to monitor
detparams.4 <- c("r",
                 "sigma0",
                 "beta0",
                 "beta1",
                 "beta2",
                 "alpha2",
                 "theta",
                 "phi0",
                 "gamma1",
                 "gamma2",
                 "logit.gamma1",
                 "lambda",
                 "Lam_mean")


# Run nimble
detfm.4 <- nimbleMCMC(code = detmod.4,
                      constants = data,
                      inits = detinits.4,
                      monitors = detparams.4,
                      thin = 10,
                      niter = 500000,
                      nburnin = 50000,
                      nchains = 3,
                      samplesAsCodaMCMC = TRUE,
                      WAIC = TRUE)
                      

# Rhat
#coda::gelman.diag(detfm.4$samples)

# Inspect
mcmcplot(detfm.4$samples)

# Model Summary
#summary(detfm.4$samples)

# WAIC
#print(detfm.4$WAIC$WAIC)



# ---------------------------------------------------------- 
# Ranking Detection Models using WAIC
# ----------------------------------------------------------

# Extract the WAIC values for each model
det_waic_values <- c(detfm.0$WAIC$WAIC,
                     detfm.1$WAIC$WAIC,
                     detfm.2$WAIC$WAIC,
                     detfm.3$WAIC$WAIC,
                     detfm.4$WAIC$WAIC)

# Naming models
det_fitnames <- c("detfm.0", 
                  "detfm.1", 
                  "detfm.2",
                  "detfm.3", 
                  "detfm.4")


# Combine model names and WAIC values into a data frame for ranking
det_waic_df <- data.frame(Model = det_fitnames, WAIC = det_waic_values)

# Rank models based on WAIC (lower WAIC is better)
det_waic_df <- det_waic_df[order(det_waic_df$WAIC), ]

# Print the ranked models
print(det_waic_df)

# Detection model 2 is the best detection model
summary(detfm.2$samples)


#  -------------------------------------------------------
#
#   Estimating Abundance 
#
#  -------------------------------------------------------

# Combine chains
combined_chains <- as.mcmc(do.call(rbind, detfm.2$samples))

# Extract Lam_mean column
Lam_mean_column <- combined_chains[, "Lam_mean"]

# Area in hectares
area <- pi*(200^2)/10000

# Getting density
dens_df <- as.data.frame(Lam_mean_column/area)

# Summarize by row
colnames(dens_df)[1] <- "Density"
dens_df[,2] <- "PC Temp Em HDS"
colnames(dens_df)[2] <- "Model"
# Switch the order of columns
dens_df <- dens_df[, c("Model", "Density")]
head(dens_df)

# Calculate the 95% Credible Interval
# ci_bounds <- quantile(dens_df$Density, probs = c(0.025, 0.975))

# Calculate the 85% Credible Interval
ci_bounds <- quantile(dens_df$Density, probs = c(0.025, 0.875))

# Subset the data frame to 95% CI
dens_df <- subset(dens_df, Density >= ci_bounds[1] & Density <= ci_bounds[2])

print(paste("min =", min(dens_df$Density), 
            "mean =", mean(dens_df$Density),
            "max =", max(dens_df$Density)))




# Export Density data frame
saveRDS(dens_df, "./Data/Fitted_Models/PC_TempEmHDS_Dens.rds")

# Plot
ggplot(dens_df, aes(x = Model, y = Density, fill = Model)) +
  geom_violin() +
  geom_boxplot(aes(x = Model, y = Density),
               width = 0.2, position = position_dodge(width = 0.8)) +
  labs(
    title = "Latent Density",
    x = "Model",
    y = "Density") +
  scale_y_continuous(limits = c(0, 25),
                     breaks = seq(0, 25, by = 5),
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

#  -------------------------------------------------------
#
#   Saving Data
#
#  -------------------------------------------------------



# Save environment
#save.image(file = "PointCount_TempEmHDS_nimble.RData")

# Clear environment
#rm(list = ls())
