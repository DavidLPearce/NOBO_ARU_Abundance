# Time-of-Detection Removal Model with Correct Data Structure
# y[site, interval, occasion]

library(R2jags)

# JAGS Model Code - CORRECTED (no loop over intervals)
cat("
model {
  
  # Priors for abundance regression
  alpha.lam ~ dnorm(0, 0.01)
  beta.lam1 ~ dnorm(0, 0.01)
  beta.lam2 ~ dnorm(0, 0.01)
  
  # Priors for detection probability regression
  alpha.p ~ dnorm(0, 0.01)
  beta.p1 ~ dnorm(0, 0.01)
  beta.p2 ~ dnorm(0, 0.01)
  
  # SURVEYED SITES (with detection data)
  for (i in 1:n_surveyed) {
    
    # Ecological model
    log(lambda_surveyed[i]) <- alpha.lam + beta.lam1*site_cov1_surveyed[i] + beta.lam2*site_cov2_surveyed[i]
    N[i] ~ dpois(lambda_surveyed[i])
    
    # Loop over occasions (j) - NO loop over intervals!
    for (j in 1:noccasions) {
      
      # Detection model varies by occasion (not by interval)
      logit(p[i,j]) <- alpha.p + beta.p1*obs_cov1[i,j] + beta.p2*obs_cov2[i,j]
      
      # Removal probabilities for the 4 intervals
      # pi[interval, site, occasion]
      pi[1,i,j] <- p[i,j]
      pi[2,i,j] <- (1 - p[i,j]) * p[i,j]
      pi[3,i,j] <- (1 - p[i,j]) * (1 - p[i,j]) * p[i,j]
      pi[4,i,j] <- 1 - p[i,j] - (1-p[i,j])*p[i,j] - (1-p[i,j])*(1-p[i,j])*p[i,j]
      
      # Multinomial observation model
      # y[site, intervals 1:4, occasion] ~ multinomial
      y[i,1:4,j] ~ dmulti(pi[1:4,i,j], N[i])
    }
  }
  
  # UNSURVEYED SITES (prediction only)
  for (u in 1:n_unsurveyed) {
    log(lambda_unsurveyed[u]) <- alpha.lam + beta.lam1*site_cov1_unsurveyed[u] + beta.lam2*site_cov2_unsurveyed[u]
    N_un[u] ~ dpois(lambda_unsurveyed[u])
  }
  
  # Derived quantities
  total_surveyed <- sum(N[])
  total_unsurveyed <- sum(N_un[])
  total_property <- total_surveyed + total_unsurveyed
  
}
", fill = TRUE, file = "./jags_models/Model_TOD_Claude.txt")

# Simulated data with CORRECT structure: y[site, interval, occasion]
set.seed(789)
n_surveyed <- 10
n_unsurveyed <- 30
noccasions <- 3
nintervals <- 4  # 3 time intervals + never detected

# Site-level covariates
site_cov1_surveyed <- rnorm(n_surveyed)
site_cov1_unsurveyed <- rnorm(n_unsurveyed)
site_cov2_surveyed <- rnorm(n_surveyed)
site_cov2_unsurveyed <- rnorm(n_unsurveyed)

# Observation-level covariates [site, occasion]
obs_cov1 <- matrix(rnorm(n_surveyed * noccasions), n_surveyed, noccasions)
obs_cov2 <- matrix(sample(0:1, n_surveyed*noccasions, replace=TRUE), n_surveyed, noccasions)

# Simulate true abundance
alpha.lam <- 2
beta.lam1 <- 0.5
beta.lam2 <- -0.3

lambda_surveyed_true <- exp(alpha.lam + beta.lam1*site_cov1_surveyed + beta.lam2*site_cov2_surveyed)
lambda_unsurveyed_true <- exp(alpha.lam + beta.lam1*site_cov1_unsurveyed + beta.lam2*site_cov2_unsurveyed)

true_N_surveyed <- rpois(n_surveyed, lambda_surveyed_true)
true_N_unsurveyed <- rpois(n_unsurveyed, lambda_unsurveyed_true)

# Simulate detections: y[site, interval, occasion]
alpha.p <- -1
beta.p1 <- -0.5
beta.p2 <- 0.8

y <- array(0, dim = c(n_surveyed, nintervals, noccasions))

for (i in 1:n_surveyed) {
  for (k in 1:noccasions) {
    p_ik <- plogis(alpha.p + beta.p1*obs_cov1[i,k] + beta.p2*obs_cov2[i,k])
    probs <- c(p_ik, 
               (1-p_ik)*p_ik, 
               (1-p_ik)^2*p_ik,
               1 - p_ik - (1-p_ik)*p_ik - (1-p_ik)^2*p_ik)
    y[i,,k] <- rmultinom(1, true_N_surveyed[i], probs)
  }
}

# Check data structure
cat("\n=== DATA STRUCTURE ===\n")
cat("y dimensions [sites, intervals, occasions]:", dim(y), "\n")
cat("Example - Site 1, all intervals, occasion 1:", y[1,,1], "\n")
cat("Example - Site 1, all intervals, occasion 2:", y[1,,2], "\n\n")

# Check that N is large enough for each occasion
cat("=== CHECKING N CONSTRAINTS ===\n")
for(i in 1:n_surveyed) {
  max_detections <- max(apply(y[i,1:3,], 2, sum))  # Max across occasions
  cat("Site", i, ": true N =", true_N_surveyed[i], 
      ", max detections =", max_detections, "\n")
}

# Data for JAGS
jags_data <- list(
  y = y,
  n_surveyed = n_surveyed,
  n_unsurveyed = n_unsurveyed,
  noccasions = noccasions,
  site_cov1_surveyed = site_cov1_surveyed,
  site_cov1_unsurveyed = site_cov1_unsurveyed,
  site_cov2_surveyed = site_cov2_surveyed,
  site_cov2_unsurveyed = site_cov2_unsurveyed,
  obs_cov1 = obs_cov1,
  obs_cov2 = obs_cov2
)

# Parameters to monitor
params <- c("alpha.lam", "beta.lam1", "beta.lam2",
            "alpha.p", "beta.p1", "beta.p2",
            "N", "N_un",
            "total_surveyed", "total_unsurveyed", "total_property")

# Initial values
inits <- function() {
  # For each site, N must be >= max detections across any occasion
  N_surveyed_init <- numeric(n_surveyed)
  for(i in 1:n_surveyed) {
    max_at_site <- max(apply(y[i,1:3,], 2, sum))  # Max detections in any occasion
    N_surveyed_init[i] <- max_at_site + rpois(1, 3)
  }
  
  N_unsurveyed_init <- rpois(n_unsurveyed, 10)
  
  list(alpha.lam = rnorm(1, 2, 0.5),
       beta.lam1 = rnorm(1, 0, 0.5),
       beta.lam2 = rnorm(1, 0, 0.5),
       alpha.p = rnorm(1, -1, 0.5),
       beta.p1 = rnorm(1, 0, 0.5),
       beta.p2 = rnorm(1, 0, 0.5),
       N = N_surveyed_init,
       N_un = N_unsurveyed_init)
}

# Run JAGS
fm1 <- jags(data = jags_data,
            inits = inits,
            parameters.to.save = params,
            model.file = "./jags_models/Model_TOD_Claude.txt",
            n.chains = 3,
            n.iter = 20000,
            n.burnin = 10000,
            n.thin = 5)

print(fm1)

# Compare estimates to true values
cat("\n=== PARAMETER ESTIMATES ===\n")
cat("alpha.lam:", alpha.lam, "vs", fm1$BUGSoutput$mean$alpha.lam, "\n")
cat("beta.lam1:", beta.lam1, "vs", fm1$BUGSoutput$mean$beta.lam1, "\n")
cat("beta.lam2:", beta.lam2, "vs", fm1$BUGSoutput$mean$beta.lam2, "\n")
cat("alpha.p:", alpha.p, "vs", fm1$BUGSoutput$mean$alpha.p, "\n")
cat("beta.p1:", beta.p1, "vs", fm1$BUGSoutput$mean$beta.p1, "\n")
cat("beta.p2:", beta.p2, "vs", fm1$BUGSoutput$mean$beta.p2, "\n")

cat("\n=== ABUNDANCE ESTIMATES ===\n")
cat("True total (surveyed):", sum(true_N_surveyed), "\n")
cat("Estimated total (surveyed):", fm1$BUGSoutput$mean$total_surveyed, "\n")
cat("True total (unsurveyed):", sum(true_N_unsurveyed), "\n")
cat("Estimated total (unsurveyed):", fm1$BUGSoutput$mean$total_unsurveyed, "\n")
cat("True total (property):", sum(true_N_surveyed) + sum(true_N_unsurveyed), "\n")
cat("Estimated total (property):", fm1$BUGSoutput$mean$total_property, 
    "±", fm1$BUGSoutput$sd$total_property, "\n")