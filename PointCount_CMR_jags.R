# -------------------------------------------------------
#
#                    Load libraries
#
# -------------------------------------------------------

# Install packages (if needed)
# install.packages("tidyverse")
# install.packages("jagsUI")
# install.packages("coda")
# install.packages("mcmcplots")

# Load library
library(tidyverse)
library(jagsUI)
library(coda)
library(mcmcplots)


# Check JAGs Version
# Latest version as of (5 Feb 2025) JAGS 4.3.1.  
# Download here: https://sourceforge.net/projects/mcmc-jags/
rjags::jags.version() 

# Set seed, scientific notation options, and working directory
set.seed(123)
options(scipen = 9999)
setwd(".")

# Setting up cores
Ncores <- parallel::detectCores()
print(Ncores) # number of available cores
workers <- Ncores * 0.5 # may need to change, for low background use 80% num, for medium use 50%

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

# Creating a day of year column
pc_dat$Date <- mdy(pc_dat$Date)
pc_dat$DOY <- yday(pc_dat$Date)
pc_dat$DOYscaled <- ((pc_dat$DOY - 1) / 365)


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

# Number of "pseudo-groups" added
nz <- M - nind

# Make site (Point_Individual) into an integer
site <- as.numeric(factor(pc_CMR$PointNum))

# Add in NA for M
site <- c(site, rep(NA, M-nind))

# Take a look
print(site)

# Extract and scale site covariates for X.abund
X.abund <- site_covs[,-c(1:4,25)] 
X.abund$woody_lrgPInx <- scale(X.abund$woody_lrgPInx)
X.abund$herb_lrgPInx  <- scale(X.abund$herb_lrgPInx)
X.abund$woody_AggInx <- scale(X.abund$woody_AggInx)
X.abund$herb_AggInx  <- scale(X.abund$herb_AggInx)
X.abund$woody_EdgDens <- scale(X.abund$woody_EdgDens)
X.abund$herb_EdgDens  <- scale(X.abund$herb_EdgDens)
X.abund$woody_Pdens <- scale(X.abund$woody_Pdens)
X.abund$herb_Pdens  <- scale(X.abund$herb_Pdens)
X.abund$woody_Npatches <- scale(X.abund$woody_Npatches)
X.abund$herb_Npatches  <- scale(X.abund$herb_Npatches)
X.abund$mnElev  <- scale(X.abund$mnElev)
print(X.abund)



# Extract and scale site covariates for X.det
X.det <- pc_dat[,c(2,7:9,18)]
X.det$Observer <- as.numeric(as.factor(X.det$Observer))
X.det$Temp.deg.F <- as.numeric(X.det$Temp.deg.F)
X.det$Wind.Beau.Code <- as.numeric(as.factor(X.det$Wind.Beau.Code)) 
X.det$Sky.Beau.Code <- as.numeric(as.factor(X.det$Sky.Beau.Code))
merged_df <- merge(site_covs, pc_CMR, by = "PointNum") # VegDens to stacked 
X.det$vegDens = as.matrix(merged_df[,c("vegDens50m")]) 
print(X.det)


# Bundle data for JAGs
data <- list(J = J,
             M = M,
             nsites = nsites,
             y = y,
             X.abund = X.abund,
             X.det = X.det,
             group = site)
 
# Take a look
str(data)

 
# ------------------------------------------------------------------------------
#
#                             Model Fitting
#
# ------------------------------------------------------------------------------

## Distributions
# library(ggfortify)
# ggdistribution(dunif, seq(0, 1, 0.001), min = 0, max = 1) # p0
# ggdistribution(dnorm, seq(0, 0.01, 0.0001)) # alpha's and beta's 
# ggdistribution(dgamma, seq(0, 5, 0.01), shape = 0.1, rate = 0.1) # tau

# Parameters monitored
params <- c("lambda",
            "p0", 
            "alpha0", 
            "alpha1", 
            "beta0", 
            "beta1",
            "beta2",
            "beta3", 
            "psi",
            "S.raneff",
            "tau",
            "sigma",
            "p_Bayes")


# Initial values
 inits  <- function() {
  list (p0 = runif(1),
        alpha1 = 0,
        beta0 = 0,
        beta1 = 0,
        beta2 = 0,
        beta3 = 0,
        tau = 1,
        z = c( rep(1,nind), rep(0, (M-nind)))
)}





# -------------------------------------------------------
# Model 0: Null Model
# -------------------------------------------------------
cat("
model {

  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0/(1-p0))
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  beta3 ~ dnorm(0, 0.01)
  psi <- sum(lambda[])/M

  # Precision for survey-level random effect
  tau ~ dgamma(0.1, 0.1)  
  sigma <- 1/sqrt(tau) # Standard deviation of survey effects

  # Abundance model
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 
    probs[s] <- lambda[s] / sum(lambda[])
  }

  # Site-level random effect
  for(s in 1:nsites){  
    S.raneff[s] ~ dnorm(0, tau)  
  }
  
  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Detection model
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 *  X.det[group[i],3] + S.raneff[group[i]]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
      
    # Posterior predictive checks
    y_rep[i,j] ~ dbern(pz[i,j])  # Generate replicated data
    discrepancy_obs[i,j] <- abs(y[i,j] - p[i,j])  # Discrepancy for observed data
    discrepancy_rep[i,j] <- abs(y_rep[i,j] - p[i,j])  # Discrepancy for replicated data  
    }
  }
  # Bayesian p-value
  sum_obs <- sum(discrepancy_obs[,])  # Sum of discrepancies for observed data
  sum_rep <- sum(discrepancy_rep[,])  # Sum of discrepancies for replicated data
  p_Bayes <- step(sum_rep - sum_obs)  # Proportion of times replicated > observed
}
", fill=TRUE, file = "./jags_models/CMR_fm0.txt")
# ------------End Model-------------

# Fit Model
fm.0 <- jags(data = data, 
             parameters.to.save = params,
             inits = inits, 
             model.file = "./jags_models/CMR_fm0.txt",
             n.iter = 300000,
             n.burnin = 20000,
             n.chains = 3, 
             n.thin = 10,
             parallel = TRUE,
             n.cores = workers,
             DIC = TRUE)  

# Save Environment
save.image(file = "CMR_JAGs.RData")
 




# -------------------------------------------------------
# Ranking Models
# -------------------------------------------------------

# Total number of models? fm.0 > fm.? = 
total_fits <- 20 # change as needed
dic_values <- numeric(totmods)  # For DIC values
fitnames <- character(totmods)   # For model names


for (i in 0:(total_fits - 1)) {
  model_name <- paste0("fm.", i)
  dic_var_name <- paste0("dic_fm.", i)
  
  assign(dic_var_name, get(model_name)$DIC)  # Extract DIC
  dic_values[i + 1] <- get(dic_var_name)  # Store DIC value
  fitnames[i + 1] <- model_name  # Store model name
}
                     
# Combine into data frame
DIC_df <- data.frame(Model = fitnames, DIC = dic_values)
                               
# Order by DIC (ascending order)
DIC_df <- DIC_df[order(DIC_df$DIC),]

# Best DIC? Lowest is better
print(DIC_df)

# Best model?
bm <- fm.?  
  
# Best model fit. P-value = 0.5 means good fit, = 1 or 0 is a poor fit
cat("Bayesian p-value =", bm$summary["p_Bayes",1], "\n")

# Check convergence
bm$Rhat # Rhat: less than 1.1 means good convergence
mcmcplot(bm$samples)# Visually inspect trace plots

# Model summary
summary(bm$samples)




# -------------------------------------------------------
#
#   Estimating Abundance 
#
# -------------------------------------------------------

# Combine chains
combined_chains <- as.mcmc(do.call(rbind, bm$samples))

# Extract lambda estimates
lambda_columns <- grep("^lambda\\[", colnames(combined_chains))
lambda_samples <- combined_chains[ ,lambda_columns]

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
dens_df <- dens_df[, c("Model", "Density")] # Switch the order of columns
head(dens_df)

# Calculate the 95% Credible Interval
ci_bounds <- quantile(dens_df$Density, probs = c(0.025, 0.975))

# Subset the data frame to 95% CI
dens_df <- subset(dens_df, Density >= ci_bounds[1] & Density <= ci_bounds[2])


# Violin plot
ggplot(dens_df, aes(x = Model, y = Density, fill = Model)) +
  geom_violin() +
  geom_boxplot(aes(x = Model, y = Density),
               width = 0.2, position = position_dodge(width = 0.8)) +
  labs(
    title = "",
    x = "Model",
    y = "Density (N/hectare)") +
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

# Total abundance
mean_dens * 1096.698
LCI_dens * 1096.698
HCI_dens * 1096.698


# Save Environment
save.image(file = "CMR_JAGs.RData")


# End Script