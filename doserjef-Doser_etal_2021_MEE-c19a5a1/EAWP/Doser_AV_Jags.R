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


library(jagsUI)
library(coda)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(dplyr)
library(mcmcplots)

# Variable Definitions --------------------------------------------------------- 

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
# A.times = indexing variable used to determine specific indexes with v > 0.
# K = true number of acosutic vocalizations
# k = true number of manually validated acoustic vocalizations
# n = total number of manually validated acoustic vocalizations 
# sites.a = specific indices of sites where acoustic data were obtained
# R.val = number of validated sites for acoustic data
# Other variables not defined are for computation of Bayesian p-values. 
# ------------------------------------------------------------------------------

#load("eawp-data.R")
load("./doserjef-Doser_etal_2021_MEE-c19a5a1/EAWP/eawp-data.R")

# c          # c = point count data
# dates.c
J          # J = number of repeat visits for acoustic data at each site
k          # k = true number of manually validated acoustic vocalizations
n          # n = total number of manually validated acoustic vocalizations 
# n.count    # n.count = number of repeat visits for count data
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
i = 87
for (i in 1:R) {
  if (length(tmp[[i]]) > 0) {
    A.times[i, 1:length(tmp[[i]])] <- tmp[[i]]
  }
}

# write.csv(A.times, "A.times.csv")
# write.csv(v, "v.csv")

J.val <- apply(v[sites.a, ], 1, function(a){sum(!is.na(a))})



# survey dates v was validated
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

# Bundle data ------------------------------------------------------------
bugs.data.A <- list(R = R, 
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

str(bugs.data.A)




# Initial Values ----------------------------------------------------------
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

# MCMC settings -----------------------------------------------------------
n.iter <- 200000
n.thin <- 50
n.burn <- 60000
n.chain <- 3
n.adapt <- 5000



# Model Statement -----------------------------------------------------------
model_AV <- "
model {
  
  # -------------------------------------------
  # Priors 
  # -------------------------------------------
  beta.0 ~ dnorm(0, .01)
  beta.1 ~ dnorm(0, .01)
  alpha.0 <- logit(mu.alpha)
  mu.alpha ~ dunif(0, 1)
  alpha.1 ~ dunif(0, 1000) # Constrained to be positive
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
  
  # -------------------------------------------
  # Likelihood and process model 
  # -------------------------------------------
  
  for (i in 1:R) {
    log(lambda[i]) <- beta.0 + beta.1 * X.lambda[i, 2]
    N[i] ~ dpois(lambda[i])
    logit(p.a[i]) <- alpha.0 + alpha.1 * N[i]
    # Acoustic Data -------------------
    for (j in 1:J[i]) {
      log(delta[i, j]) <- gamma.1[days[i, j]]
      y[i, j] ~ dbin(p.a[i], 1)
      tp[i, j] <- delta[i, j] * N[i] / (delta[i, j] * N[i] + omega)
      # Posterior predictive checks for Bayesian P-value
      y.pred[i, j] ~ dbin(p.a[i], 1)
      resid.y[i, j] <- pow(pow(y[i, j], 0.5) - pow(p.a[i], 0.5), 2)
      resid.y.pred[i, j] <- pow(pow(y.pred[i, j], 0.5) - pow(p.a[i], 0.5), 2)
    } # j
    
    for (j in 1:J.r[i]) {
      v[i, A.times[i, j]] ~ dpois((delta[i, A.times[i, j]] * N[i] + omega) * phi[i, A.times[i, j]] * y[i, A.times[i, j]]) T(1, )
      v.pred[i, j] ~ dpois((delta[i, A.times[i, j]] * N[i] + omega) * phi[i, A.times[i, j]] * y[i, A.times[i, j]]) T(1, )
      mu.v[i, j] <- ((delta[i, A.times[i, j]] * N[i] + omega) * phi[i, A.times[i, j]]) / (1 - exp(-1 * ((delta[i, A.times[i, j]] * N[i] + omega) * phi[i, A.times[i, j]])))
      resid.v[i, j] <- pow(pow(v[i, A.times[i, j]], 0.5) - pow(mu.v[i, j], 0.5), 2)
      resid.v.pred[i, j] <- pow(pow(v.pred[i, j], 0.5) - pow(mu.v[i, j], 0.5), 2)
    } # j
  } # i
  
  # ------------------------------------------- 
  # Manual validation 
  # -------------------------------------------
  
  for (i in 1:R.val) {
    for (j in 1:J.val[i]) {
      K[i, j] ~ dbin(tp[sites.a[i], j], v[sites.a[i], val.times[i, j]])
      k[i, val.times[i, j]] ~ dhyper(K[i, j], v[sites.a[i], val.times[i, j]] - K[i, j], n[i, val.times[i, j]], 1)
    } # j
  } # i
  
  # -------------------------------------------
  # Bayesian P-value
  # -------------------------------------------
  
  for (i in 1:R.A) {
    tmp.v[i] <- sum(resid.v[sites.a.v[i], 1:J.r[sites.a.v[i]]])
    tmp.v.pred[i] <- sum(resid.v.pred[sites.a.v[i], 1:J.r[sites.a.v[i]]])
  }
  fit.y <- sum(resid.y[sites.a, 1:J.A])
  fit.y.pred <- sum(resid.y.pred[sites.a, 1:J.A])
  fit.v <- sum(tmp.v[1:R.A])
  fit.v.pred <- sum(tmp.v.pred[1:R.A])
  bp.y <- step(fit.y.pred - fit.y)
  bp.v <- step(fit.v.pred - fit.v)
}
"
# Write the model to a temporary file
model_file <- tempfile(pattern = "model_", fileext = ".txt")
writeLines(model_AV, con = model_file)

# Run Model -----------------------------------------------------------
out.model.A <- jags(
  data = bugs.data.A, 
  inits = inits, 
  parameters.to.save = params.A, 
  model.file = model_file, 
  n.iter = n.iter, 
  n.thin = n.thin, 
  n.burnin = n.burn, 
  n.chains = n.chain, 
  n.adapt = n.adapt, 
  parallel = TRUE
)

# # Clean up the temporary file if no longer needed
# unlink(model_file)

# Rhat
coda::gelman.diag(out.model.A$samples)

# Trace plots
mcmcplot(out.model.A$samples)

# Model Summary
summary(out.model.A$samples)

# Abundance Intercept -----------------
summary(mcmc(out.model.A$sims.list$beta.0))$quantiles

# Year Effect -------------------------
summary(mcmc(out.model.A$sims.list$beta.1))$quantiles
#sum(out.model.A$sims.list$beta.1 < 0) / length(out.model.A.C$sims.list$beta.1)

# Bayesian p-values -------------------------------------------------------
mean(out.model.A$sims.list$bp.v)

# Figure of trend estimates -----------------------------------------------

# Extract the summary statistics
a.v.trends <- data.frame(
  Parameter = "beta.1",
  Mean = mean(out.model.A$sims.list$beta.1),
  Lower = quantile(out.model.A$sims.list$beta.1, probs = 0.025),
  Upper = quantile(out.model.A$sims.list$beta.1, probs = 0.975)
)

ggplot(a.v.trends, aes(x = Parameter, y = Mean)) +
  geom_pointrange(aes(ymin = Lower, ymax = Upper), size = 1) +
  geom_point(size = 3, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Posterior Estimate of beta.1",
    x = "Parameter",
    y = "Estimate"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
