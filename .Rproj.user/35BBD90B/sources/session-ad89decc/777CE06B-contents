# Author: David L. Pearce
# Description:
#             TBD

# This code extends code from the following sources: 
#     1. Kery, M. and Royle, J. A. (2020). Applied hierarchical modeling 
#           in ecology: Analysis of distribution, abundance, and species 
#           richness in R and BUGS: Volume 2: Dynamic and advanced models. 
#           Academic Press.


# Citation: 
#      TBD

# -------------------------------------------------------
#
#                    Load libraries
#
# -------------------------------------------------------

# Install packages (if needed)
# install.packages("tidyverse")
# install.packages("jagsUI")
# install.packages("coda")
# install.packages("MCMCvis")

# Load library
library(tidyverse)
library(jagsUI)
library(coda)
library(MCMCvis)

# Check JAGs Version
# Download here: https://sourceforge.net/projects/mcmc-jags/
rjags::jags.version() # Code written under JAGS 4.3.1.  

# Set seed, scientific notation options, and working directory
set.seed(123)
options(scipen = 9999)
setwd(".")

# Setting up cores
Ncores <- parallel::detectCores()
print(Ncores) # Number of available cores
workers <- Ncores  * 0.3 # For low background use 80%, for medium use 50% of Ncores
print(workers)

# Source custom function for checking Rhat values > 1.1
source("./Scripts/Rhat_check_function.R")

# Model name object
model_name <- "PC MCR"

# -------------------------------------------------------
#
# Variable and Object Definitions ******* Finish this
# 
# -------------------------------------------------------

# beta0 = abundance intercept 

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

# Rename SiteID to PointNum for matching
colnames(site_covs)[1] <- "PointNum"


# Remove rows with NA values (no birds detected)
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

# Number of sites
S <- 10

# Number of individuals detected
nind <- nrow(y)

# Superpopulation
print(nind)
M <- nind + 100 # ensure that M is larger than number detected
print(M)

# Data Augmentation
y <- rbind(y, matrix(0, nrow=(M-nind), ncol=4))

# Number of "pseudo-groups" added
nz <- M - nind

# Make site (Point_Individual) into an integer
site <- as.numeric(factor(pc_CMR$PointNum))

# Add in NA for M
site <- c(site, rep(NA, M-nind))

# Take a look
print(site)


# Extract and scale site covariates for detection
Wind <- pc_dat[, 'Wind.Beau.Code']
Wind <- as.matrix(as.numeric(as.factor(Wind)))
Wind_Lvls = nrow(unique(Wind))
merged_df <- merge(site_covs, pc_CMR, by = "PointNum") # VegDens to stacked 
VegDens <- merged_df[,c("vegDens50m")]
VegDens <- scale(VegDens)

# Inspect
head(Wind)
head(VegDens)

# Format abundance covariates
Herb_COH <- as.matrix(scale(site_covs[,'herb_COH'])) #  herb_ClmIdx
Woody_SPLIT <- as.matrix(scale(site_covs[, 'woody_SPLIT']))#  woody_AggInx

# Inspect
head(Herb_COH)
head(Woody_SPLIT)

# Bundle data for JAGs
data <- list(J = J,
             M = M,
             S = S,
             y = y,
             Wind = Wind,
             Wind_Lvls = Wind_Lvls,
             VegDens = VegDens,
             Herb_COH = Herb_COH,
             Woody_SPLIT = Woody_SPLIT,
             group = site
)

# Check structure 
str(data)


# ------------------------------------------------------------------------------
#
#                             Model Fitting
#
# ------------------------------------------------------------------------------

# -------------------
# MCMC Specifications
# -------------------

n_iter = 200000
n_burnin = 60000
n_chains = 3 
n_thin = 10
n_adapt = 5000

# Posterior samples
n_post_samps = ((n_iter * n_chains) - n_burnin) / n_thin
print(n_post_samps)


# Parameters monitored
params <- c("lambda",
            "N",
            "N_tot",
            "alpha0", 
            "alpha1", 
            "alpha2", 
            "p",
            "beta0", 
            "beta1",
            "beta2",
            "psi",
            "S_RE_lam",
            'tau_lam',
            "S_RE_det",
            "tau_det",
            "p_Bayes",
            "sum_obs",
            "sum_rep")



# Initial values
inits  <- function() {
  list (p0 = runif(1),
        alpha1 = rep(0, Wind_Lvls),
        alpha2 = 0,
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
)}


# -------------------------------------------------------
# Model Statement
# -------------------------------------------------------
cat("
model {

  # ---------------------------------
  # Abundance Priors
  # ---------------------------------
  beta0 ~ dnorm(0, 0.001) # Intercept
  beta1 ~ dnorm(0, 0.001) # Herbaceous cohesion index
  beta2 ~ dnorm(0, 0.001) # Woody splitting index
  
  # # Site random effect
  # tau_lam ~ dgamma(0.001, 0.001)
  # for (s in 1:S) {
  #    S_RE_lam[s] ~ dnorm(beta0, tau_lam)
  # }

  # ---------------------------------
  # Detection Priors
  # ---------------------------------
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0)) 
  
  # Wind is a categorical covariate
  for (w in 1:Wind_Lvls){
  alpha1[w] ~ dnorm(0, 0.001)
  }
  
  # Vegetation Density
  alpha2 ~ dnorm(0, 0.001)
  
  # Site-level random effect for pseudoreplication
  tau_det ~ dgamma(0.001, 0.001) 
  for(s in 1:S){  
    S_RE_det[s] ~ dnorm(0, tau_det)  
  }
  
  # ---------------------------------
  # Individual Encounter/Presence Derived 
  # ---------------------------------
  psi <- sum(lambda[])/M
  
  # -------------------------------------------
  #
  # Likelihood and Process Model 
  #
  # -------------------------------------------

  # ---------------------------------
  # Abundance Submodel
  # ---------------------------------

  for(s in 1:S){
    log(lambda[s]) <- beta0 + beta1 * Herb_COH[s, 1] +  beta2 * Woody_SPLIT[s, 1]
    
    # Individual site probability
    probs[s] <- lambda[s] / sum(lambda[])
    
    # Estimated abundance at site s
    N[s] <- sum(group[] == s)  
  }

  # ---------------------------------
  # Presence Submodel 
  # ---------------------------------

  for(i in 1:M){
  
    # Group == site membership
    group[i] ~ dcat(probs[])  
    
    # Presence: Data augmentation variables
    z[i] ~ dbern(psi)         
    
  # ---------------------------------
  # Detection Submodel
  # ---------------------------------
 
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1[Wind[group[i], 1]] + alpha2 * VegDens[group[i],1] + S_RE_det[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
  # ---------------------------------
  # Posterior Predictive checks
  # ---------------------------------

    # Generate replicated data
    y_rep[i,j] ~ dbern(pz[i,j]) 
    
    # Discrepancy for observed data
    # discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j]) 
    discrepancy_obs[i,j] <- pow(y[i,j] - p[i,j], 2)
    
    # Discrepancy for replicated data  
    # discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  
    discrepancy_rep[i,j] <- pow(y_rep[i,j] - p[i,j], 2)


    } # End J
  } # End M
  
  # ---------------------------------
  # Derived Metrics
  # ---------------------------------
 
  # Abundance
  N_tot <- sum(z[])
  
  # ---------------------------------
  # Bayesian p-value
  # ---------------------------------
  
  # Sum of discrepancies for observed data
  sum_obs <- sum(discrepancy_obs[,])  
  
  # Sum of discrepancies for replicated data
  sum_rep <- sum(discrepancy_rep[,])  
  
  # Proportion of times replicated > observed
  p_Bayes <- step(sum_rep - sum_obs)  
  
}

", fill=TRUE, file = "./jags_models/PC_MCR_Model.txt")
# ------------End Model-------------

# Fit Model
fm1 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/PC_MCR_Model.txt",
             n.iter = n_iter,
             n.burnin = n_burnin,
             n.chains = n_chains, 
             n.thin = n_thin,
             n.adapt = n_adapt,
             parallel = TRUE,
             n.cores = workers,
             verbose = TRUE,
             DIC = FALSE)

# -------------------------------------------------------
# Check Convergence
# -------------------------------------------------------

# Rhat
check_rhat(fm1$Rhat, threshold = 1.1)

# Trace plots
MCMCvis::MCMCtrace(fm1, 
                   params = c("lambda",
                              "N",
                              "N_tot",
                              "alpha0", 
                              "alpha1", 
                              "alpha2", 
                              "beta0", 
                              "beta1",
                              "beta2",
                              "psi",
                              # "S_RE_lam",
                              # 'tau_lam',
                              "S_RE_det",
                              "tau_det"
                   ),
                   pdf = T,
                   filename = "PC-MCR_TracePlots.pdf",
                   wd = "./Figures"
)

# -------------------------------------------------------
# Combine Chains for Posterior inference
# -------------------------------------------------------

# Combine chains
combined_chains <- as.mcmc(do.call(rbind, fm1$samples))

# ----------------------
# Extract Fit
# ----------------------

# Abundance
fit_y_data <- data.frame(
  Observed = as.vector(combined_chains[, "sum_obs"]), # Observed values
  Predicted = as.vector(combined_chains[, "sum_rep"]) # Predicted values
)

# ----------------------
# Density Plot
# ----------------------

# Bayes P-value
# P-value = 0.5 means good fit, = 1 or 0 is a poor fit
mn_bpy <- round(mean(fm1$summary["p_Bayes",1]), 2) 

# y
y_PPC_Dens <- ggplot(fit_y_data) +
  geom_density(aes(x = Observed, fill = "Observed"), alpha = 0.5) +   
  geom_density(aes(x = Predicted, fill = "Predicted"), alpha = 0.5) +  
  scale_fill_manual(values = c("Observed" = "blue", "Predicted" = "red")) +  
  labs(title = "", 
       x = "Fit Values", 
       y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") + 
  annotate("text", x = 100, y = 0.03, label = paste0("Bayes p-value = ", mn_bpy), hjust = 0) 

# View
print(y_PPC_Dens)

# Export                
ggsave(plot = y_PPC_Dens, "./Figures/PPC_BP_PC-MCR_Yobs.jpeg", width = 8, height = 5, dpi = 300)
dev.off()

# -------------------------------------------------------
#
#   Beta Estimates and Covariate Effects 
#
# -------------------------------------------------------


# -------------------------------------------------------
# Beta Estimates
# -------------------------------------------------------

# Extract beta estimates
beta0_samples <- combined_chains[, "beta0"]
beta1_samples <- combined_chains[, "beta1"]
beta2_samples <- combined_chains[, "beta2"]

# Compute 95% CI for each beta
beta_df <- data.frame(
  value = c(beta0_samples, beta1_samples, beta2_samples ),
  parameter = rep(c("beta0", "beta1", "beta2" ), each = length(beta0_samples))
) %>%
  group_by(parameter) %>%
  filter(value >= quantile(value, 0.025) & value <= quantile(value, 0.975))  # Keep only values within 95% CI

# Add model
beta_df$Model <- model_name

# Create Summary
beta_summary <- beta_df %>%
  group_by(parameter, Model) %>%
  summarise(
    mean = mean(value),
    LCI = min(value),
    UCI = max(value),
    .groups = "drop"  
  )
# view
print(beta_summary)

# Export beta dataframe and summary
saveRDS(beta_df, "./Data/Model_Data/PC_MCR_beta_df.rds")
write.csv(beta_summary, "./Figures/PC_MCR_BetaSummary.csv")

# -------------------------------------------------------
# Alpha Estimates
# -------------------------------------------------------

# Extract alpha estimates
alpha0_samples <- combined_chains[, "alpha0"]
alpha1_samples <- combined_chains[, grepl("^alpha1\\[", colnames(combined_chains))]
alpha1_samples_1 <- alpha1_samples[,1]
alpha1_samples_2 <- alpha1_samples[,2]
alpha1_samples_3 <- alpha1_samples[,3] 
alpha1_samples_4 <- alpha1_samples[,4]
alpha1_samples_5 <- alpha1_samples[,5]
alpha2_samples <- combined_chains[, "alpha2"]

# Extract site random effect
sRE_samples <- combined_chains[, c("S_RE_det[1]", "S_RE_det[2]","S_RE_det[3]","S_RE_det[4]","S_RE_det[5]",
                                   "S_RE_det[6]", "S_RE_det[7]", "S_RE_det[8]", "S_RE_det[9]", "S_RE_det[10]")]
sRE_samples <- rowMeans(sRE_samples) # Row means

# Compute 95% CI  
alpha_df <- data.frame(
  value = c(alpha0_samples, 
            alpha1_samples_1, 
            alpha1_samples_2,
            alpha1_samples_3,
            alpha1_samples_4,
            alpha1_samples_5,
            alpha2_samples,
            sRE_samples),  
  parameter = rep(c("alpha0", 
                    "alpha1_1", 
                    'alpha1_2',
                    'alpha1_3',
                    'alpha1_4',
                    'alpha1_5',
                    "alpha2",
                    "sRE_det"), each = length(alpha0_samples))
) %>%
  group_by(parameter) %>%
  filter(value >= quantile(value, 0.025) & value <= quantile(value, 0.975))  

# Add model
alpha_df$Model <- model_name

# Create summary
alpha_summary <- alpha_df %>%
  group_by(parameter, Model) %>%
  summarise(
    mean = mean(value),
    LCI = min(value),
    UCI = max(value),
    .groups = "drop"  
  )

# view
print(alpha_summary)

# Export alpha summary
saveRDS(alpha_summary, "./Data/Model_Data/PC_MCR_AlphaSummary.rds")


# -------------------------------------------------------
# Detection probability 
# -------------------------------------------------------

# Extract detection posterior samples
p_samples <- combined_chains[, grepl("^p\\[", colnames(combined_chains))]

# Combine all samples into a single long vector
all_p_samples <- as.vector(p_samples)

# Create summary
p_summary <- data.frame(
  mean = mean(all_p_samples),
  LCI = quantile(all_p_samples, 0.025),
  UCI = quantile(all_p_samples, 0.975)
)

# Add model
p_summary$Model <- model_name

# Add parameter name
p_summary$Parameter <- "detection"

# view
print(p_summary)

# -------------------------------------------------------
#   Estimating Abundance 
# -------------------------------------------------------

# Extracting Abundance
Ntot_samples <- combined_chains[ ,"N_tot"]

# Ntotal is the abundance based on 10 point counts at a radius of 200m.

# Area in hectares of a 200m radius circle
area <- pi * (200^2) / 10000  # Area in hectares

# Calculate density (individuals per hectare)
dens_samples <- Ntot_samples / (area * 10)

# Create data frame for density
dens_df <- data.frame(Model = rep(model_name, length(dens_samples)), Density = dens_samples)
colnames(dens_df)[2] <- "Density"
head(dens_df)

# Abundance estimates
abund_df <- dens_df # posterior estimates
abund_df$Density <- abund_df$Density * 1096  

abund_summary <- abund_df %>% # summary
  group_by(Model) %>%
  summarise(
    mean = mean(Density),
    LCI = quantile(Density, 0.025),
    UCI = quantile(Density, 0.975)
  )

# Add model
abund_summary$Model <- model_name

# Add parameter name
abund_summary$Parameter <- "abundance"

# View
print(abund_summary)

# Combine with detection
param_summary <- rbind(p_summary, abund_summary)

# View
print(param_summary)

# Trim abund df to 95% CI
abund_95df <- abund_df %>%
  left_join(abund_summary, by = "Model") %>%
  filter(Density >= LCI & Density <= UCI) %>%
  select(-LCI, -UCI)


  
# Export alpha summary
saveRDS(abund_95df, "./Data/Model_Data/PC_MCR_abund_df.rds")
saveRDS(param_summary, "./Data/Model_Data/PC_MCR_param_summary.rds")

# ------------ End Script -----------------