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

# Number of time intervals
J <- 4

# number of sites
nsites <- 10

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

# Extract site covariates
herbPdens = as.vector(scale(site_covs[,c("herb_Pdens")])) 
woodyParea = as.vector(site_covs[,c("woody_Parea")]) 


# Bundle data for JAGs
data <- list(J = J,
             M = M,
             nsites = nsites,
             y = y,
             herbPdens = herbPdens,
             woodyParea = woodyParea,
             group = site)
 
# Take a look
str(data)

 


# -------------------------------------------------------
#
#                   Model Fitting
#
# -------------------------------------------------------


# -------------------------------------------------------
# Covariate Model
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0 / (1-p0))    # same as logit(p0)
  # alpha1 ~ dnorm(0, 0.01)
  # alpha2 ~ dnorm(0,0.01)
  beta0 ~ dnorm(0,0.01)
  beta1 ~ dnorm(0,0.01)
  beta2 ~ dnorm(0,0.01)
  psi <- sum(lambda[]) / M   # psi is a derived parameter

  # log-linear model for abundance
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1 * herbPdens[s] + beta2 * woodyParea[s]
    probs[s] <- lambda[s] / sum(lambda[])
  }

  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables

    # Observation model: p 
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
    }
  }
}
",fill=TRUE,file="cov_mod.txt")
# ------------------------- End Model ------------------------------

# Parameters monitored
params.cov <- c("lambda",
                "p0", 
                "alpha0", 
                #"alpha1", 
                "beta0", 
                "beta1",
                "beta2",
                "psi")

# Initial values
inits.cov <- function(){
      list (p0 = runif(1), 
            alpha1 = runif(1), 
            alpha2 = rnorm(1),
            beta0 = runif(1),
            beta1 = rnorm(1),
            beta2 = rnorm(1),
            z = c( rep(1,nind), rep(0, (M-nind)))
)}# end inits

fm.cov <- jags(data = data, 
               parameters.to.save = params.cov, 
               inits = inits.cov, 
               model.file = "cov_mod.txt",
               n.iter = 500000,
               n.burnin = 5000,
               n.chains = 3, 
               n.thin = 10,
               parallel = TRUE,
               n.cores = 8,
               DIC = TRUE)  


# Rhat
fm.cov$Rhat

# Trace plots
mcmcplot(fm.cov$samples)

# Model summary
print(fm.cov, digits = 3)





# -------------------------------------------------------
# Covariate Model with Survey Random Effect
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  #alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau_survey ~ dgamma(0.1, 0.1)  
  sigma_survey <- 1/sqrt(tau_survey)

  # Log-linear model for abundance
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1 * herbPdens[s] + beta2 * woodyParea[s]
    probs[s] <- lambda[s] / sum(lambda[])
  }

  # Survey-level random effect
  for(j in 1:J){  
    survey_effect[j] ~ dnorm(0, tau_survey)  
  }

  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Observation model: p depends on survey-specific effect
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + survey_effect[j] 
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
    }
  }
}
", fill=TRUE, file="raneff_cov_mod.txt")
# ------------------------- End Model ------------------------------

# Parameters monitored
params.rcov <- c("lambda",
              "p0", 
              "alpha0", 
              #"alpha1",
              "beta0", 
              "beta1",
              "beta2",
              "psi", 
              "tau_survey",
              "sigma_survey")


# Initial values: add tau
inits.rcov <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        tau_survey = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
  )}# end inits




# Fit Model
fm.rcov <- jags(data = data, 
             parameters.to.save = params.rcov,
             inits = inits.rcov, 
             model.file = "raneff_cov_mod.txt",
             n.iter = 500000,
             n.burnin = 50000,
             n.chains = 3, 
             n.thin = 10,
             parallel = TRUE,
             n.cores = 8,
             DIC = TRUE)  



# Rhat
fm.rcov$Rhat

# Model summary
print(fm.rcov, digits = 3)

# Trace plots
mcmcplot(fm.rcov$samples)








# -------------------------------------------------------
#               Individual Heterogeneity Model
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  alpha1 ~ dnorm(0, 0.01)
  #alpha2 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0,0.01)
  psi <- sum(lambda[])/M
  tau ~ dgamma(0.1,0.1) # precision of ind. random effects
  sigma <- 1/sqrt(tau)

  # log-linear model for abundance: lambda depends on WOODY
  for(s in 1:nsites){
    log(lambda[s])<- beta0 + beta1 * herbPdens[s] + beta2 * woodyParea[s]
    probs[s]<- lambda[s]/sum(lambda[])
  }

  # Model for individual encounter histories
  for(i in 1:M){
    eta[i] ~ dnorm(alpha0, tau)  # Individual random effect
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Observation model: p depends ind. heterogeneity
    for(j in 1:J){
      logit(p[i,j]) <- alpha1 * X[group[i],2] + eta[i] #+ alpha2*X[group[i],3] 
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
    }
  }
}
",fill=TRUE,file="het_mod.txt")
# ------------------------- End Model ------------------------------


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
        alpha1 = 0,
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
)}# end inits




# Fit Model
fm.h <- jags(data = data, 
             parameters.to.save = params.h,
             inits = inits.h, 
             model.file = "het_mod.txt",
             n.iter = 500000,
             n.burnin = 50000,
             n.chains = 3, 
             n.thin = 10,
             parallel = TRUE,
             n.cores = 8,
             DIC = TRUE)  

# Model summary
print(fm.h, digits = 3)

# Rhat
coda::gelman.diag(fm.h$samples)

# Trace plots
mcmcplot(fm.h$samples)






# -------------------------------------------------------
#
#   Compare DIC
#
# -------------------------------------------------------

print(fm.cov$DIC)
print(fm.h$DIC)


# -------------------------------------------------------
#
#   Estimating Abundance 
#
# -------------------------------------------------------


# Combine chains
combined_chains <- as.mcmc(do.call(rbind, fm.covs$samples))

# Extract lambda estimates (assuming "lambda[1]", "lambda[2]", etc., are in the output)
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
dens_df[,2] <- "PC CMR"
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
  scale_y_continuous(limits = c(0, 10),
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
saveRDS(dens_df, "./Data/Fitted_Models/PC_CMR_Dens.rds")

















