# -------------------------------------------------------
#
#                    Load libraries
#
# -------------------------------------------------------

# Load packages
library(tidyverse)
library(jagsUI)
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

# Format date
pc_dat$Date <- as.POSIXct(pc_dat$Date, format = "%m/%d/%Y")

# Create a day of year column
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

# Create a 3D array
y3d <- array(NA,dim=c(nrow(det_mat), 4, 4) ) # Length of site (10), width of distance bins (4), depth of surveys (4)

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
area <- pi*(200^2)/10000        # Area surveyed in hectares

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
             HerbPRP = as.vector(site_covs[,'herb_prp']),
             WoodyPatch = as.vector(scale(site_covs[,'woody_Parea'])),
             Observer = Observer_numeric,
             Temp = scale(temp_mat),
             Wind = wind_mat,
             Sky = sky_mat,
             DOY = as.matrix(scale(doy_mat))
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
# Avail Model Hazad Det Fun
# ----------------------------------------------------------
cat("
model {

  # Priors
  # Abundance parameters
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)

  # Availability parameters
  phi0 ~ dunif(0.1, 0.9)
  logit.phi0 <- log(phi0/(1-phi0))
  for(k in 1:4){
    gamma1[k] ~ dunif(0.1, 0.9) # Availability effects of surveys 1 - 4
    logit.gamma1[k]<- log(gamma1[k]/(1-gamma1[k]))
  }
  gamma2 ~ dnorm(0, 0.01)

  # Detection parameters
  sigma0 ~ dunif(0.1,50)   # Intercept
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

      # Conditional 4-part hierarchical model
      pdet[s,k] <- sum(cellprobs[s,1:nD,k])
      pmarg[s,k] <- pdet[s,k]*phi[s,k]
      
      # Part 4
      y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[s,1:nD,k], nobs[s,k])
      
      # Part 3: Number of detected individuals
      nobs[s,k] ~ dbin(pmarg[s,k], M[s])  
      
      # Part 2: Number of available individuals
      Navail[s,k] ~ dbin(phi[s,k],M[s])   
    } # end k loop

     M[s] ~ dnegbin(prob[s], r)
     prob[s] <- r/(r+lambda[s])
     # M[s] ~ dpois(lambda[s]) 
     
     # Part 1: Abundance model
     log(lambda[s]) <- beta0 + beta1*HerbPRP[s] + beta2*WoodyPatch[s]
     
  }  # end s loop

  # Derived quantities
  for(k in 1:K){
    Davail[k] <- mean(phi[,k])*exp(beta0)/area
  }
  Mtotal <- sum(M[])
  Dtotal <- exp(beta0)/area
} # End model
" ,fill=TRUE, file="availmodel0.txt")
# ----------------------------------------------------------
 
# Inits
availinits.0 <- function() {
  list(
    M = apply(y3d, 1, max) + 5,
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
                   "gamma2"
                   # "lambda",
                   # "Mtotal",
                   # "Dtotal"
                   )


# Run JAGs 
availfm.hz <- jags(data = data, 
                  inits = availinits.0, 
                  parameters.to.save = availparams.0, 
                  model.file = "availmodel0.txt", 
                  n.thin = 10,
                  n.iter = 100000,
                  n.burnin = 10000,
                  n.chains = 3, 
                  parallel = TRUE)

# DIC
print(availfm.hz$DIC)

# # model summary
# summary(availfm.hz)
# 
# # Convergence
# print(availfm.hz$Rhat) # Rhat
# mcmcplot(availfm.0$samples) # trace plots






