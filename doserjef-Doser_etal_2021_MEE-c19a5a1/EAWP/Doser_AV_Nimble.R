# Author: Jeffrey W. Doser
# Description: Integrated model combining acoustic and point-count data to 
#              estimate abundance of a single species, assuming constant
#              abundance and detection probability across sites used to 
#              fit data from Marsh-Billings-Rockefeller National Historic
#              Park to estimate Eastern Wood-pewee abundance. 

# This code extends code from the following sources: 
#     1. Chambert, T., Waddle, J. H., Miller, D. A., Walls, S. C., 
#        and Nichols, J. D. (2018b). A new framework for analysing 
#        automated acoustic species detection data: Occupancy estimation 
#        and optimization of recordings post-processing. 
#        Methods in Ecology and Evolution, 9(3):560–570.
#     2. Kery, M. and Royle, J. A. (2020). Applied hierarchical modeling 
#        in ecology: Analysis of distribution, abundance, and species 
#        richness in R and BUGS: Volume 2: Dynamic and advanced models. 
#        (in press). Academic Press.

# Citation: 
#     Doser, J.W., Finley, A. O., Weed, A. S., and Zipkin, E. F. (2020).
#     Integrating automated acoustic vocalization data and point count 
#     surveys for efficient estimation of bird abundance. 
#     arXiv preprint arXiv:2011.11047.

# -------------------------------------------------------
# Load libraries
# -------------------------------------------------------
library(tidyverse)
library(coda)
library(nimble)
library(mcmcplots)

# Set seed, scientific notation options, and working directory
set.seed(123)
options(scipen = 9999)
setwd(".")

# -------------------------------------------------------
# Variable Definitions
# -------------------------------------------------------
# beta.0 = abundance intercept 
# beta.1 = abundance trend estimate
# alpha.0 = prob (on logit scale) of detecting at least one vocalization at a site
#           that is not occupied.
# alpha.1 = additional prob (on logit scale) of detecting at least one vocalization at 
#           a site that is not occupied. 
# omega = mean # of false positive acoustic detections
# p = detection probability of an individual in point count data
# tau.day = precision for random day effect on true vocalization detection rate. 
# a.phi = overdispersion parameter for zero-truncated negative binomial. 
# gamma.1 = random day effect on true vocalization detection rate
# n.days = number of recording days.
# N = latent abundance process
# tp = true positive rate
# p.a = prob of detecting at least one vocalization in an acoustic recording
# v = acoustic vocalization data from clustering algorithm
# y = binary summary of v, takes value 0 if v = 0 and value 1 if v > 0. 
# c = point count data
# R = number of total sites
# J = number of repeat visits for acoustic data at each site
# J.A = max number of repeat visits at each acoustic data site. 
# n.count = number of repeat visits for count data
# R.val = number of sites where validation of acoustic data occurred
# days = variable used to index different recording days. 
# A.times = indexing variable used to determine specific indexes
#           with v > 0. 
# K = true number of acosutic vocalizations
# k = true number of manually validated acoustic vocalizations
# n = total number of manually validated acoustic vocalizations 
# sites.a = specific indices of sites where acoustic data were obtained
# R.val = number of validated sites for acoustic data
# Other variables not defined are for computation of Bayesian p-values. 
# ------------------------------------------------------------------------------

# -------------------------------------------------------
# Load and Bundle Data
# -------------------------------------------------------

# data from Doser et al 2021
load("./doserjef-Doser_etal_2021_MEE-c19a5a1/EAWP/eawp-data.R")

c          # c = point count data
dates.c
J          # J = number of repeat visits for acoustic data at each site
k          # k = true number of manually validated acoustic vocalizations
n          # n = total number of manually validated acoustic vocalizations 
n.count    # n.count = number of repeat visits for count data
R          # R = number of total sites
R.val      # R.val = number of validated sites for acoustic data
sites.a    # sites.a = specific indices of sites where acoustic data were obtained
v          # v = acoustic vocalization data from clustering algorithm
X.lambda   # Covatiate matrix for abundance
y          # y = binary summary of v, takes value 0 if v = 0 and value 1 if v > 0. 

# J.r contains the number of surveys at each acoustic site that contains
# at least 1 detected vocalization. 
J.r <- apply(v, 1, function(a) {sum(a != 0)})
J.r <- ifelse(is.na(J.r), 0, J.r)

# A.times is a n.sites x J matrix that contains the indices of the
# specific surveys at each acoustic site that have at least one detected
# vocalization. This determines the specific components of v[i, j] that
# are used in the zero-truncated Poisson.
A.times <- matrix(NA, R, max(J))
tmp <- apply(v, 1, function(a) {which(a != 0)})
for (i in 1:R) {
  if (length(tmp[[i]]) > 0) {
    A.times[i, 1:length(tmp[[i]])] <- tmp[[i]]
  }
}

J.val <- apply(v[sites.a, ], 1, function(a){sum(!is.na(a))})
val.times  <- matrix(NA, R.val, max(J))
tmp <- apply(v[sites.a, ], 1, function(a) {which(!is.na(a))})
for (i in 1:R.val) {
  val.times[i, 1:length(tmp[[i]])] <- tmp[[i]]
}

# For day random effect on cue production rate
n.days <- 14
days <- matrix(NA, R, max(J))
days[76, ] <- rep(1:n.days, each = 3)
days[77, ] <- rep((n.days + 1):(n.days * 2), each = 3)
days[78, ] <- rep((n.days * 2 +  1):(n.days * 3), each = 3)
days[79, ] <- rep(1:n.days, each = 3)
days[80, ] <- rep((n.days + 1):(n.days * 2), each = 3)
days[81, ] <- rep((n.days * 2 +  1):(n.days * 3), each = 3)
days[82, ] <- rep(1:n.days, each = 3)
days[83, ] <- rep((n.days + 1):(n.days * 2), each = 3)
days[84, ] <- rep((n.days * 2 +  1):(n.days * 3), each = 3)
days[85, ] <- rep(1:n.days, each = 3)
days[86, ] <- rep((n.days + 1):(n.days * 2), each = 3)
days[87, ] <- rep((n.days * 2 +  1):(n.days * 3), each = 3)

# For Bayesian P-value
sites.c <- which(n.count != 0)
n.count.c <- n.count[which(n.count != 0)]
R.A <- sum(J.r > 0)
sites.a.v <- which(J.r > 0)
J.A <- max(J)

# Bundle data  
data <- list(R = R, 
             n.count = n.count, 
             J = J, 
             X.lambda = X.lambda,
             v = v, 
             y = y, 
             k = k, 
             n = n, 
             sites.a = sites.a, 
             R.val = R.val, 
             J.val = J.val, 
             J.r = J.r, 
             A.times = A.times, 
             val.times = val.times, 
             R.A = R.A, 
             J.A = J.A, 
             sites.a.v = sites.a.v, 
             n.sites = R / 3, 
             days = days,
             n.days = max(days, na.rm = TRUE))

# Constants
constants <- list(n.days = max(days, na.rm = TRUE),
                  R = R,
                  J = J,
                  J.A = J.A,
                  J.r = J.r,
                  R.val = R.val, 
                  J.val = J.val)

# Initial Values 
N.init <- rep(1, R)
c.max <- apply(c, 1, max, na.rm = TRUE)
N.init <- ifelse(c.max > 1, c.max + 1, 1)
inits <- function() {
  list(
    N = N.init, 
    beta.0 = rnorm(1),
    beta.1 = rnorm(1), 
    omega = runif(1, 0, 10), 
    tau.day = runif(1, 0.1, 1),
    tau = runif(1, 0.1, 1),
    tau.p = runif(1, 0.1, 1),
    tau.day.p = runif(1, 0.1, 1),
    alpha.1 = runif(1, 0, 1), 
    a.phi = runif(1, 0, 5)
  )
}

# -------------------------------------------------------
# Hypergeometric Function
# -------------------------------------------------------

# Nimble does not have a hypergeometric function
# need to create one


# Custom binomial coefficient function for nimble
binom_coeff <- nimbleFunction(
  run = function(n = integer(0), k = integer(0)) {
    if (k > n) return(0)
    result <- 1
    for (i in 1:k) {
      result <- result * (n - i + 1) / i
    }
    returnType(double(0))
    return(result)
  }
)

# Define custom distribution function
dhyper_nimble <- nimbleFunction(
  run = function(x = double(0),   # Number of successes in sample (use 'x' here)
                 m = double(0),   # Number of successes in population
                 nf = double(0),  # Number of failures in population
                 k = double(0),   # Sample size
                 log = logical(0)) {  # Log-transform flag
    # Hypergeometric probability density function:
    num <- binom_coeff(m, x) * binom_coeff(nf, k - x)
    denom <- binom_coeff(m + nf, k)
    result <- num / denom
    if (log) {
      result <- log(result)
    }
    returnType(double(0))  # Set the return type to double(0)
    return(result)
  }
)

# Define custom random generation function
rhyper_nimble <- nimbleFunction(
  run = function(n = integer(0),   # Number of random samples to generate
                 m = double(0),    # Number of successes in population
                 nf = double(0),   # Number of failures in population
                 k = double(0)) {   # Sample size
    # Initialize output vector for random samples
    result <- double(n)  # Change to double(0) to match the expected return type
    
    # Generate n samples
    for (i in 1:n) {
      successes <- 0
      failures <- 0
      # Create a pool of successes and failures
      pool <- c(rep(1, m), rep(0, nf))  # 1s for successes, 0s for failures
      
      # Shuffle the pool manually using Nimble's random numbers
      for (j in 1:k) {
        # Generate a random index from the pool size
        idx <- floor(runif(1, 1, m + nf + 1))  # Random index between 1 and m+nf
        
        # Increment success or failure based on the chosen index
        if (pool[idx] == 1) {
          successes <- successes + 1
        } else {
          failures <- failures + 1
        }
        
        # Remove the selected element by setting it to a special value (e.g., -1)
        pool[idx] <- -1
      }
      
      result[i] <- successes
    }
    
    returnType(double(0))  # Set the return type to double(0)
    return(result)
  }
)

# Register the custom distribution
registerDistributions(list(
  dhyper_nimble = list(
    BUGSdist = "dhyper_nimble(x, m, nf, k)",
    discrete = TRUE,
    range = c(0, Inf),
    types = c(
      "value = double(0)",  # The result of the distribution (output)
      "x = double(0)",      # Number of successes in sample (the variable)
      "m = double(0)",      # Parameter for successes in population
      "nf = double(0)",     # Parameter for failures in population
      "k = double(0)"       # Sample size
    )
  )
))


# -------------------------------------------------------
# Model Statement
# -------------------------------------------------------

model_AV <- nimbleCode({
  # Priors
  beta.0 ~ dnorm(0, .01)
  beta.1 ~ dnorm(0, .01)
  alpha.0 <- logit(mu.alpha)
  mu.alpha ~ dunif(0, 1)
  alpha.1 ~ dunif(0, 1000)
  omega ~ dunif(0, 1000)
  tau.day ~ dgamma(.01, .01)
  a.phi ~ dunif(0, 100)
  
  for (i in 1:n.days) {
    gamma.1[i] ~ dnorm(0, tau.day)
  }
  
  for (i in 1:R) {
    for (j in 1:J.A) {
      phi[i, j] ~ dgamma(a.phi, a.phi)
    }
  }
  
  # Likelihood and process model
  for (i in 1:R) {
    log(lambda[i]) <- beta.0 + beta.1 * X.lambda[i, 2]
    N[i] ~ dpois(lambda[i])
    logit(p.a[i]) <- alpha.0 + alpha.1 * N[i]
    
    # Acoustic Data
    for (j in 1:J[i]) {  # Use precomputed J_i
      log(delta[i, j]) <- gamma.1[days[i, j]]
      y[i, j] ~ dbin(p.a[i], 1)
      tp[i, j] <- delta[i, j] * N[i] / (delta[i, j] * N[i] + omega)
      
      # Posterior predictive checks for Bayesian P-value
      y.pred[i, j] ~ dbin(p.a[i], 1)
      resid.y[i, j] <- pow(pow(y[i, j], 0.5) - pow(p.a[i], 0.5), 2)
      resid.y.pred[i, j] <- pow(pow(y.pred[i, j], 0.5) - pow(p.a[i], 0.5), 2)
    }
    
    for (j in 1:J.r[i]) {   
      idx <- A.times[i, j]  # Ensure A.times[i, j] is treated as a scalar index
      v[i, idx] ~ dpois((delta[i, idx] * N[i] + omega) * phi[i, idx] * y[i, idx]) 
      v[i, idx] <- max(v[i, idx], 1)
      mu.v[i, j] <- ((delta[i, idx] * N[i] + omega) * phi[i, idx]) / (1 - exp(-1 * ((delta[i, idx] * N[i] + omega) * phi[i, idx])))
      resid.v[i, j] <- pow(pow(v[i, idx], 0.5) - pow(mu.v[i, j], 0.5), 2)
      resid.v.pred[i, j] <- pow(pow(v.pred[i, j], 0.5) - pow(mu.v[i, j], 0.5), 2)
    }
  }
  
  # Manual validation
  for (i in 1:R.val) {
    for (j in 1:J.val[i]) {
      # Directly assign idx_val from val.times, ensuring it's a valid index before passing it into the model
      idx_val <- val.times[i, j]  # Assuming val.times[i, j] produces valid indices
      K[i, j] ~ dbin(tp[sites.a[i], j], v[sites.a[i], idx_val])
      k[i, idx_val] ~ dhyper_nimble(K[i, j], v[sites.a[i], idx_val] - K[i, j], n[i, idx_val], 1)
    }
  }
})

  
  # # Bayesian P-value
  # for (i in 1:R.A) {
  #   tmp.v[i] <- 0  # Initialize the sum for tmp.v
  #   tmp.v.pred[i] <- 0  # Initialize the sum for tmp.v.pred
  #   
  #   # Manually sum over the indices instead of using dynamic indexing
  #   n_jr_i_v <- J.r[sites.a.v[i]]  # Get the number of valid indices
  #   for (j in 1:n_jr_i_v) {
  #     tmp.v[i] <- tmp.v[i] + resid.v[sites.a.v[i], j]  # Sum residuals for each site
  #     tmp.v.pred[i] <- tmp.v.pred[i] + resid.v.pred[sites.a.v[i], j]  # Sum predicted residuals for each site
  #   }
  # }
  # 
  # fit.y <- 0
  # fit.y.pred <- 0
  # for (i in 1:R) {
  #   for (j in 1:J.A) {
  #     fit.y <- fit.y + resid.y[i, j]
  #     fit.y.pred <- fit.y.pred + resid.y.pred[i, j]
  #   }
  # }
  # 
  # fit.v <- 0
  # fit.v.pred <- 0
  # for (i in 1:R.A) {
  #   fit.v <- fit.v + tmp.v[i]
  #   fit.v.pred <- fit.v.pred + tmp.v.pred[i]
  # }
  # 
  # bp.y <- step(fit.y.pred - fit.y)
  # bp.v <- step(fit.v.pred - fit.v)
#})

    

# --------------------------------------------------------------------------------------------------------


# Parameters monitored ----------------------------------------------------
params.A <- c('alpha.0', 
              'alpha.1', 
              'beta.0',
              'beta.1', 
              'tau', 
              'tau.day', 
              'a.phi', 
              'omega', 
              'bp.y', 
              'bp.v')



# Fit Model ----------------------------------------------------
fm.AV <- nimbleMCMC(code = model_AV, 
                   data = data,
                   constants = constants,
                   inits = inits,
                   monitors = params.A,
                   niter = 200000,
                   nburnin = 60000,
                   nchains = 3,
                   thin = 50,
                   samplesAsCodaMCMC = TRUE,
                   WAIC = TRUE)



