#---(begin AHMnimble header)---
# This file contains code adapted from the file R_BUGS_code_AHM_Vol_1_20170519.R
# available at https://www.mbr-pwrc.usgs.gov/pubanalysis/keryroylebook/, with
# permission of the authors.  It may have been modified so that examples originally
# written to work with JAGS or WinBUGS will work with NIMBLE.  More information
# about NIMBLE can be found at https://r-nimble.org, https://github.com/nimble-dev/nimble,
# and https://cran.r-project.org/web/packages/nimble/index.html.
#
# The file  R_BUGS_code_AHM_Vol_1_20170519.R contains the following header:
# -----(begin R_BUGS_code_AHM_Vol_1_20170519.R header)-----
# =========================================================================
#
#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#
#   Marc Kéry & J. Andy Royle
#
#   *** This is the text file with all R and BUGS code from the book ***
#
#   Created 2 Dec 2015 based on draft from 21 Oct 2015
#
# =========================================================================
### Last change: 19 May 2017 by Mike Meredith
# Incorporated errata up to 19 May 2017
# Code updated to implement new names and changes to functions in AHMbook package 0.1.4
# Use built-in data sets instead of .csv files
# In Chapter 11, replaced 'Y', 'Ysum' and 'Yaug' with lower case 'y', 'ysum' and 'yaug'
#  to match the code in the printed book.
#
# -----(end R_BUGS_code_AHM_Vol_1_20170519.R header)-----

# This file was created by:
#
# Jacob Levine and Perry de Valpine
#
#---(end AHMnimble header)---
# DOES NOT WORK
# 9.5.4.2 Bayesian analysis of the Wagtail data
# ------------------------------------------------------------------------

##########################################################################
## From Mike Meredith AHM1 Ch 9 ##
library(AHMbook)
library(unmarked)

# Load the wagtail data, investigate NA patterns in detection data
data("wagtail")
str(wagtail)
Y <- wagtail$Y
table(n.missing <- rowSums(is.na(Y))) # Frequency distribution of number of NAs per site
n.missing
keep <- which(n.missing == 0)     # Sites with complete distance data
Y <- Y[keep,]                     # restrict analysis to those

# Harvest other data for sites with complete distance data
potato <- wagtail$potato[keep]   ;   grass <- wagtail$grass[keep]
lscale <- wagtail$lscale[keep]   ;   hour <- wagtail$hour[keep,]
date <- wagtail$date[keep,]   ;   rep <- wagtail$rep[keep,]
breaks <- wagtail$breaks

# Look at the distance data
str(Y)
tmp <- apply(Y, 2, sum, na.rm = TRUE)
matplot(1:6, t(matrix(tmp, nrow = 4, byrow= TRUE)), type = "b", ylim = c(0, 90),
        xlab = "Distance class", ylab = "Number of wagtails", frame = FALSE,
        lwd = 3, lty = 1)

# Standardize all continuous covariates
mn.potato <- mean(potato)   ;   sd.potato <- sd(potato)
mn.grass <- mean(grass)   ;   sd.grass <- sd(grass)
mn.lscale <- mean(lscale)   ;   sd.lscale <- sd(lscale)
mn.date <- mean(date)    ;   sd.date <- sd(c(date))
mn.hour <- mean(hour)   ;   sd.hour <- sd(c(hour))
POTATO <- (potato - mn.potato) / sd.potato
GRASS <- (grass - mn.grass) / sd.grass
LSCALE <- (lscale - mn.lscale) / sd.lscale
DATE <- (date - mn.date) / sd.date
HOUR <- (hour - mn.hour) / sd.hour

##########################################################################

y3d <- array(NA,dim=c(nrow(Y), 6, 4) )          # Create 3d array
y3d[,,1] <- Y[,1:6]    
y3d[,,2] <- Y[,7:12]  
y3d[,,3] <- Y[,13:18]   
y3d[,,4] <- Y[,19:24]

K <- 4                          # Number of primary occasions
nsites <- nrow(Y)               # Number of sites
nD <- 6                         # Number of distance classes
midpt <- seq(25,275,50)         # Class midpoint distance
delta <- 50                     # Class width
B <- 300                        # Maximum distance
nobs <- apply(y3d, c(1,3), sum) # Total detections per site and occasion

# Area in hectares
area <- pi*(300^2)/10000

# Bundle data
data <- list(y3d = y3d, 
             nsites = nsites, 
             K = K, 
             nD = nD, 
             midpt = midpt, 
             delta = delta, 
             B = B,
             nobs = nobs, 
             area = area,
             POTATO = POTATO,
             GRASS = GRASS,
             LSCALE = LSCALE,
             DATE = DATE,
             HOUR = HOUR)

str(data)

# ---------------------------------------------------------- 
# 
# No Covariates
# 
# ----------------------------------------------------------
Section9p5p4p2_code <- nimbleCode({
    
    # Priors
    # Abundance parameters
    beta0 ~ dnorm(0, 0.01)
    beta1 ~ dnorm(0, 0.01)
    
    # Availability parameter
    phi ~ dunif(0,1)
    
    # Detection parameter
    sigma ~ dunif(0,500)  
    
    # Multinomial cell probabilities
    for(b in 1:nD){
      log(g[b]) <- -midpt[b]*midpt[b]/(2*sigma*sigma)  # Half-normal model
      f[b] <- (2*midpt[b]*delta)/(B*B) # Scaled radial density function
      cellprobs[b] <- g[b]*f[b]
      cellprobs.cond[b] <- cellprobs[b]/sum(cellprobs[1:nD])
    }
    cellprobs[nD+1] <- 1-sum(cellprobs[1:nD])
    
    for (s in 1:nsites) {
      for (k in 1:K) {
      # Conditional 4-part version of the model
        pdet[s,k] <- sum(cellprobs[1:nD])
        pmarg[s,k] <- pdet[s,k]*phi
        y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[1:nD], nobs[s,k]) # Part 4: distance  
        nobs[s,k] ~ dbin(pmarg[s,k], M[s])  # Part 3: number of detected individuals
        Navail[s,k] ~ dbin(phi,M[s])        # Part 2: Number of available individuals
      }  # end k loop
    
      M[s] ~ dpois(lambda[s])    #  Part 1: Abundance model
      log(lambda[s]) <- beta0    #  Habitat variables would go here
    }  # end s loop
    
    # Derived quantities
    for(k in 1:K){
      Davail[k] <- phi*exp(beta0)/area
    }
    Mtotal <- sum(M[1:nsites])
    Dtotal<- exp(beta0)/area
  } # end model
)
# ----------------------------------------------------------
# # Inits
# Navail.st <- apply(y3d, c(1,3),sum)  
# Mst <- apply(Navail.st, c( 1), max,na.rm=TRUE) + 2
# inits <- function() list(M=Mst, sigma = 100.0)
# initsValues <- inits
# 
# # Parameters to save
# params <- c("sigma", "phi", "beta0", "beta1", "Mtotal", "Davail", "Dtotal")
# 
# # MCMC settings
# ni <- 12000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3
# 
# ## run nimble from R
# 
# wag1 <- nimbleMCMC(
#   code = Section9p5p4p2_code,
#   constants = data,
#   inits = initsValues,
#   monitors = params,
#   nburnin = nb,
#   niter = ni,
#   samplesAsCodaMCMC = TRUE
# )
# 
# library(coda)
# library(mcmcplots)
# colnames(wag1)
# params_of_interest <- c("Davail[1]", "Dtotal", "Mtotal", "beta0", "beta1", "phi", "sigma")
# mcmcplot(wag1[, params_of_interest])
# ## mixing looks good


# ---------------------------------------------------------- 
# 
# Covariates
# 
# ----------------------------------------------------------
Section9p5p4p2_code.1 <- nimbleCode({
    
    # Priors
    # Abundance parameters
    beta0 ~ dnorm(0, 0.01)
    beta1 ~ dnorm(0, 0.01)
    beta2 ~ dnorm(0, 0.01)
    beta3 ~ dnorm(0, 0.01)
    
    # Availability parameters
    phi0 ~ dunif(0,1)
    logit.phi0 <- log(phi0/(1-phi0))
    for(k in 1:4){
      gamma1[k] ~ dunif(0, 1) # Availability effects of surveys 1 - 4
      logit.gamma1[k]<- log(gamma1[k]/(1-gamma1[k])) 
    }
    gamma2 ~ dnorm(0, 0.01)  
    gamma3 ~ dnorm(0, 0.01)
    gamma4 ~ dnorm(0, 0.01)
    
    # Detection parameters
    sigma0 ~ dunif(0.1,500)   # Intercept  
    alpha2 ~ dnorm(0, 0.01)   # effect of DATE (linear)
    alpha3 ~ dnorm(0, 0.01)   # effect of DATE (squared)
    alpha4 ~ dnorm(0, 0.01)   # effect of HOUR
    theta ~ dgamma(0.1, 0.1)
    r ~ dunif(0, 10)
    
    for (s in 1:nsites) {
      for (k in 1:K) {
        # Availability parameter
        logit.phi[s,k] <- logit.gamma1[k] + gamma2*POTATO[s] + gamma3*GRASS[s] + gamma4*LSCALE[s]
        phi[s,k] <- exp(logit.phi[s,k])/(1+ exp(logit.phi[s,k]))
        
        # Distance sampling parameter
        log(sigma[s,k]) <- log(sigma0) + alpha2*DATE[s,k] + alpha3*pow(HOUR[s,k],2) + alpha4*HOUR[s,k]
        
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
      log(lambda[s]) <- beta0 + beta1*POTATO[s] + beta2*GRASS[s] + beta3*LSCALE[s]
    }  # end s loop
    
    # Derived quantities
    for(k in 1:K){
    Davail[k] <- mean(phi[1:nsites,k])*exp(beta0)/area
    }
    Mtotal <- sum(M[1:nsites])
    Dtotal <- exp(beta0)/area
  } # End model
)
# ----------------------------------------------------------
# Inits
# Navail.st <- apply(y3d, c(1,3),sum)
# Mst <- apply(Navail.st, c( 1), max,na.rm=TRUE) + 2

Navail = apply(y3d, c(1, 3), sum, na.rm = TRUE)
M = apply(Navail.st, 1, max, na.rm = TRUE) + 2

inits <- function() {
  list(
    M = Mst,
    Navail = Navail.st,
    sigma0 = 100,
    alpha2 = 0,
    alpha3 = 0,
    alpha4 = 0,
    gamma1 = rep(0.5, 4),
    gamma2 = 0,
    gamma3 = 0,
    gamma4 = 0,
    beta0 = 0,
    beta1 = 0,
    beta2 = 0,
    beta3 = 0,
    phi0 = 0.5,
    theta = 1,
    r = 5
  )
}

# Parameters to save
params <- c("r",
            "sigma0",
            "beta0", 
            "beta1", 
            "beta2", 
            "beta3", 
            "Mtotal", 
            "alpha2", 
            "alpha3",  
            "alpha4", 
            "theta", 
            "Dtotal", 
            "Davail", 
            "phi0", 
            "gamma1", 
            "gamma2", 
            "gamma3", 
            "gamma4" ,
            "logit.gamma1")


## Run nimble 
wag2 <- nimbleMCMC(
  code = Section9p5p4p2_code.1,
  constants = data,
  inits = inits,
  monitors = params,
  thin = 5,
  niter = 10000,
  nburnin = 1000,
  nchains = 3,
  samplesAsCodaMCMC = TRUE
)


mcmcplot(wag2)

# ## compare efficiency between nimble and jags
# Section9p5p4p2_compare <- compareMCMCs(
#   modelInfo = list(
#     code = Section9p5p4p2_code.1,
#     data = data,
#     inits = initsValues
#   ),
#   monitors = params,
#   MCMCs = c('nimble', 'jags'),
#   summary = FALSE,
#   burnin = nb,
#   niter = ni
# )
# 
# outputDir = "~/Desktop/NIMBLE_EXAMPLES/AHMnimble/Chapter_9/comparison_pages"
# make_MCMC_comparison_pages(Section9p5p4p2_compare, modelNames = "Section9p5p4p2", dir = outputDir)
# 
# browseURL(file.path(outputDir, "Section9p5p4p2.html"))
