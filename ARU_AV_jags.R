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

# Site covariates
site_covs <- read.csv("./Data/Acoustic_Data/ARU_siteCovs.csv")


# -------------------------------------------------------
#
#    Data Wrangling
#
# -------------------------------------------------------


# Subset to 14 days starting at May 26 and ending on July 17. Dates are every 4 days.
bnet_dat <- bnet_dat_all %>%
  filter(Date %in% c("2024-05-26", "2024-05-30", "2024-06-03",  "2024-06-07",
                     "2024-06-11", "2024-06-15", "2024-06-19", "2024-06-23",
                     "2024-06-27", "2024-07-01", "2024-07-05", "2024-07-09",
                     "2024-07-13",  "2024-07-17"))
                     
# Dates and their corresponding occasion numbers
date_order <- c("2024-05-26", "2024-05-30", "2024-06-03", "2024-06-07",
                "2024-06-11", "2024-06-15", "2024-06-19", "2024-06-23",
                "2024-06-27", "2024-07-01", "2024-07-05", "2024-07-09",
                "2024-07-13", "2024-07-17")

# Adding occasion column
bnet_dat <- bnet_dat %>%
  mutate(Occasion = match(Date, date_order))

# Adding a row count
bnet_dat$Count <- 1


# Initialize the matrix
v <- matrix(0, nrow = 27, ncol = 14)        

# Extract count data
for (i in 1:nrow(bnet_dat)) {
  
  # Extracting plot ID
  site <- bnet_dat$Site_Number[i]
  
  # Extracting Occasion
  occasion <- bnet_dat$Occasion[i]
  
  # Fill in the matrix with the number of individuals
  v[site, occasion] =  v[site, occasion] + bnet_dat$Count[i]
  
} # end loop 

# take a look
print(v)
                             
# ARU at site 24 stopped recording on 6/20. Making Columns 8:14 NA
v[24,8:14] <- NA
print(v)  

# Renaming columns to date Month_day
formatted_dates <- format(as.Date(date_order), "%b_%d")
colnames(v) <- formatted_dates
print(v) 

# Adding rownames
rownames(v) <- as.numeric(1:27)

# Get the site numbers with at least one call
v_mat <- as.matrix(v)
rownames(v_mat) <- NULL  
sites.a <- which(rowSums(!is.na(v_mat) & v >= 1) > 0)
print(sites.a)



# Creating a binary detection matrix, values >= 1 become 1, and 0 stays 0
v_df <- as.data.frame(v)
y <- v_df %>%
  mutate(across(where(is.numeric), ~ if_else(. >= 1, 1, 0)))
y <- as.matrix(y)
rownames(y) <- NULL  
colnames(y) <- NULL 
print(y)

# Total number of sites
R <- as.integer(27) 

# Number of repeat visits
J <- rep(14, R)  



# ---------------------------------------
# Manually validated info ****** made up ----- Fix this 
# ---------------------------------------

#write.csv(v, "v.csv")

# calls validated, do not include sites with no calls
# only include occasions where at least 1 call was validated for a site
n <- as.matrix(read.csv("./n.csv", row.names = 1))

# calls found to be true, same dimension as n
k <- as.matrix(read.csv("./k.csv", row.names = 1))

# Survey dates call was validated, same dimension as n
# should be the column number 
val.times  <- as.matrix(read.csv("val.times.csv", row.names = 1))

# Total number of sites with manually validated data
R.val <- nrow(n)

# How many surveys were validate
J.val <- rep(14, R.val)  


# ---------------------------------------
# Covariates *** site covs made up ***
# ---------------------------------------

# Extract site covariates
herbPdens <- as.vector(scale(site_covs[,c("herb_Pdens")])) 
woodyParea <- as.vector(site_covs[,c("woody_Parea")]) 

# For day random effect on cue production rate
days <- matrix(NA, nrow = 27, ncol = 14)  
for (col in 1:14) {
  days[, col] <- col
}
print(days)



# J.r contains the number of surveys at each acoustic site that contains at least 1 detected vocalization. 
J.r <- apply(v_mat, 1, function(a) {sum(a != 0)})
J.r <- ifelse(is.na(J.r), 0, J.r)
J.r <- as.numeric(J.r)
print(J.r)

# A.times is a R x J matrix that contains the indices of the
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
print(A.times)



# For Bayesian P-value
R.A <- sum(J.r > 0)
sites.a.v <- which(J.r > 0)
J.A <- max(J)


# Bundle data 
bugs.data.A <- list(R = R, 
                    J = J, 
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
                    days = days,
                    n.days = max(J),
                    herbPdens = herbPdens,
                    woodyParea = woodyParea)

str(bugs.data.A)




# -------------------------------------------------------
# Acoustic Model
# -------------------------------------------------------
cat(" model {
  
  # -------------------------------------------
  # Priors 
  # -------------------------------------------
  beta.0 ~ dnorm(0, .01)
  beta.1 ~ dnorm(0, .01)
  beta.2 ~ dnorm(0, .01)
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
  
    # Abundance Model
    log(lambda[i]) <- beta.0 + beta.1 * herbPdens[i] + beta.2 * woodyParea[i]
    N[i] ~ dpois(lambda[i])
    
    # Detection Model
    logit(p.a[i]) <- alpha.0 + alpha.1 * N[i]
    
    # Acoustic Data 
    for (j in 1:J[i]) {
    
    # Vocalization Model
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
", fill=TRUE, file="./jags_models/ARU_mod0.txt")
# ------------End Model-------------



# Parameters monitored
params.A <- c('alpha.0', 
              'alpha.1', 
              'beta.0',
              'beta.1',
              'beta.2',
              'tau', 
              'tau.day', 
              'a.phi', 
              'omega', 
              'bp.y', 
              'bp.v',
              'lambda')

# Initial Values 
N.init <- rep(1, R)
inits <- function() {
  list(
    N = N.init, 
    beta.0 = rnorm(1),
    beta.1 = rnorm(1), 
    beta.2 = rnorm(1),
    omega = runif(1, 0, 10), 
    tau.day = runif(1, 0.1, 1),
    tau = runif(1, 0.1, 1),
    tau.p = runif(1, 0.1, 1),
    tau.day.p = runif(1, 0.1, 1),
    alpha.1 = runif(1, 0, 1), 
    a.phi = runif(1, 0, 5)
  )
}#end inits




# Fit Model
fm.0 <- jags(data = bugs.data.A, 
                    inits = inits, 
                    parameters.to.save = params.A, 
                    model.file = "./jags_models/ARU_mod.txt", 
                    n.iter = 200000, 
                    n.burnin = 60000, 
                    n.chains = 3,
                    n.thin = 50,
                    n.adapt = 5000, 
                    parallel = TRUE,
                    n.cores = 8,
                    DIC = TRUE) 
                    

# Rhat
fm.0$Rhat

# Model summary
print(fm.0, digits = 3)

# Trace plots
mcmcplot(fm.0$samples)

# Bayesian P value
cat("Bayesian p-value =", fm.0$summary["bp.v",1], "\n")

# Average number of false positives detections
cat("False positives =", fm.0$summary["omega",1], "\n")

#  -------------------------------------------------------
#
#   Estimating Abundance 
#
#  -------------------------------------------------------


# Combine chains
combined_chains <- as.mcmc(do.call(rbind, out.model.A$samples))

# Extract lambda estimates
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
dens_df[,2] <- "ARU Bnet"
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
  scale_y_continuous(limits = c(0, 20),
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
#saveRDS(dens_df, "./Data/Fitted_Models/ARU_Bnet.rds")


