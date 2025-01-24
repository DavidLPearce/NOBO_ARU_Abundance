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
#             Variable and Object Definitions
#
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
# A.times = indexing variable used to determine specific indexes with v > 0.
# K = true number of acosutic vocalizations
# k = true number of manually validated acoustic vocalizations
# n = total number of manually validated acoustic vocalizations 
# sites.a = specific indices of sites where acoustic data were obtained
# R.val = number of validated sites for acoustic data
# Other variables not defined are for computation of Bayesian p-values. 

# -------------------------------------------------------
#
#                    Load Data
#
# -------------------------------------------------------

# BirdNet detections
bnet_dat_all <- read.csv("./Data/Acoustic_Data/ARU_BirdNET_alldates.csv")


# -------------------------------------------------------
#
#                                   Data Wrangling
#
# -------------------------------------------------------


# Subset to same dates as point count data
bnet_dat <- bnet_dat_all %>%
  filter(Date %in% c("2024-04-28", "2024-05-19", "2024-06-18", "2024-07-12"))

## Organize into a site x occasion matrix of counts of calls
# Define the full list of site numbers
all_sites <- tibble(Site_Number = 1:27)

# Count observations by Site_Number and Date ***************** add in NA for missing survey dates
countobs <- bnet_dat %>%
  group_by(Site_Number, Date) %>%
  summarise(Count = n(), .groups = "drop")

# Ensure all site numbers and dates are represented
v <- countobs %>%
  complete(Site_Number = all_sites$Site_Number, Date, fill = list(Count = 0)) %>%
  pivot_wider(names_from = Date, values_from = Count, values_fill = 0)

# Formatting to numeric
v <- v %>%
  mutate(across(-1, as.numeric))

# View the result
print(v, n = 27)  

# Get the site numbers with at least one call
sites.a <- v %>%
  rowwise() %>%
  mutate(Any_Call = sum(c_across(where(is.numeric))) > 0) %>%
  filter(Any_Call) %>%
  pull(Site_Number)

# Format the output as a vector of site numbers
print(sites.a) 

# Creating a binary detection matrix, values >= 1 become 1, and 0 stays 0
y <- v %>%
  mutate(across(where(is.numeric), ~ if_else(. >= 1, 1, 0)))

# Removing first row and defining as a matrix
y <- y[,-1]
y <- as.matrix(y)
rownames(y) <- NULL  
colnames(y) <- NULL 

# View the result
print(y)

# Sites are organized correctly, removing site_number columnm defining as a matrix and setting to default row/col names
v <- v[,-1]
v <- as.matrix(v)
rownames(v) <- NULL  
colnames(v) <- NULL  
write.csv(v, "v.csv")

# View the result
print(v) # used as v = acoustic vocalization data from clustering algorithm

# Total number of sites
R <- as.integer(27) 

# Number of repeat visits
J <- rep(4, R)  



# ---------------------------------------
# Manually validated info ****** made up
# ---------------------------------------

#sum(det_mat, na.rm = TRUE)

#calls2eval<- sample(177, size = (108*.3), replace = FALSE)
#write.csv(det_mat, "det_mat.csv")

# Need to redo this - i think not all sites were validated 
n <- read.csv("./totalval_fake.csv") # n = total number of manually validated acoustic vocalizations 
k <- read.csv("./totaltrue_fake.csv") # k = true number of manually validated acoustic vocalizations

# Convert to a matrix and numeric
n <- as.matrix(n)
k <- as.matrix(k)

# set row and col names to default ([,1], [,2])
rownames(n) <- NULL  
colnames(n) <- NULL  
rownames(k) <- NULL  
colnames(k) <- NULL 

# Total number of sites with manually validated data
# R.val <- sum(apply(n, 1, function(row) any(row >= 1)))
# print(R.val) # used as R.val = number of validated sites for acoustic data
R.val = R

# ---------------------------------------
# Covariates ****** made up
# ---------------------------------------
set.seed(123)

# Generate a random vector for 27 sites (e.g., between 0 and 100)
random_covariate <- runif(R, min = 0, max = 100)

# Scale the covariate between 0 and 1 using min-max scaling
X.lambda <- (random_covariate - min(random_covariate)) / (max(random_covariate) - min(random_covariate))

# convert to matrix
X.lambda <- as.matrix(X.lambda)

# Display the scaled covariate
print(X.lambda)


# ---------------------------------------

# J.r contains the number of surveys at each acoustic site that contains
# at least 1 detected vocalization. 
J.r <- apply(v, 1, function(a) {sum(a != 0)})
J.r <- ifelse(is.na(J.r), 0, J.r)
J.r <- as.numeric(J.r)


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

# How many valid observations across visits
J.val <- apply(v[sites.a, ], 1, function(a){sum(!is.na(a))})


# Times validation occured for sites that were validated
# val.times  <- matrix(NA, R.val, max(J))
# tmp <- apply(v[sites.a, ], 1, function(a) list(which(!is.na(a))))
# i = 1
# for (i in 1:R.val) {
#   val.times[i, 1:length(tmp[[i]])] <- tmp[[i]]
# }

# All times are valid (no NAs from survey failure)
val.times <- matrix(rep(1:4, times = nrow(v)), nrow = nrow(v), byrow = TRUE)



# i=1
# j=1
# 
# for (i in 1:R.val) {
#   for (j in 1:J.val[i]) {
#     K[i, j] ~ dbin(tp[sites.a[i], j], v[sites.a[i], val.times[i, j]])
#     k[i, val.times[i, j]] ~ dhyper(K[i, j], v[sites.a[i], val.times[i, j]] - K[i, j], n[i, val.times[i, j]], 1)
#   } # j
# } # i

# For day random effect on cue production rate
days <- matrix(rep(1:4, times = nrow(v)), nrow = nrow(v), byrow = TRUE)

# For Bayesian P-value
R.A <- sum(J.r > 0)
sites.a.v <- which(J.r > 0)
J.A <- max(J)


# Bundle data 
bugs.data.A <- list(R = R, 
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


# Initial Values 
N.init <- rep(1, R)
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



# Parameters monitored
params.A <- c('alpha.0', 
              'alpha.1', 
              'beta.0',
              'beta.1', 
              'tau', 
              'tau.day', 
              'a.phi', 
              'omega', 
              'bp.y', 
              'bp.v',
              'lambda')

# MCMC settings 
n.iter <- 200000
n.thin <- 50
n.burn <- 60000
n.chain <- 3
n.adapt <- 5000



# Model Statement 
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
    log(lambda[i]) <- beta.0 + beta.1 * X.lambda[i, 1]
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

# Run Model
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


# Rhat
coda::gelman.diag(out.model.A$samples)

# Trace plots
mcmcplot(out.model.A$samples)

# Model Summary
summary(out.model.A$samples)




