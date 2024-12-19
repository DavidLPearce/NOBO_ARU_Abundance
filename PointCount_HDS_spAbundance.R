# -------------------------------------------------------
#
#                    Load libraries
#
# -------------------------------------------------------

# Install packages (if needed)
# install.packages("tidyverse")
# install.packages("spAbundance")

# Load library
library(tidyverse)
library(spAbundance)

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

# Removing NAs, an NA means that there were no detections
pc_dat <- na.omit(pc_dat)

# -------------------------------------------------------
# Detection Matrix
# ------------------------------------------------------- 

# Organize data into counts by distance bins for each row
det_mat <- pc_dat %>%
  group_by(PointNum, DistBin) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = DistBin, values_from = Count, values_fill = 0)%>%
  rename_with(~paste0("DistBin.", .), -PointNum) %>% # colnames to DistBin..1 ... DistBin.3
  column_to_rownames(var = "PointNum") # rownames are point count ID 

# Take a look 
print(det_mat)


# -------------------------------------------------------
# Covariates Matrix
# ------------------------------------------------------- 

# Since counts are summarized across (pooled) four survey occasions
# the observation covatiates would be a average across survey days.
# So, will only use the distance detection function to estimate detection


# -------------------------------------------------------
# Format long
# ------------------------------------------------------- 

# Formatting long
pc_spA_HDS <- list(y = det_mat,
                   covs = site_covs[,c('herb_prp', 
                                       'woody_mean_p_Area', 
                                       'woody_c_clumpy')],
                   dist.breaks = c(0, 0.05, 0.1, 0.2),
                   offset = 31.05521) # in acres
                                                            

# Check formatting
str(pc_spA_HDS)

# -------------------------------------------------------
#
#               Hierarchical Distance Models
#
# -------------------------------------------------------

# MCMC Specifications
batch.length <- 20
n.batch <- 10000 
batch.length * n.batch # Total number of MCMC samples per chain
n.burn <- 20000 
n.thin <- 10
n.chains <- 3

# Initial values for model start
inits <- list(
              alpha = 0.1,                                   # Detection covariates
              beta = 0.1,                                    # Single abundance intercept
              N = apply(pc_spA_HDS$y, 1, sum)) # Abundance


# Set Priors
priors <- list(
              alpha.normal = list(mean = 0, var = 10), # Prior for detection - Narrower variance for faster convergence
              beta.normal = list(mean = 0, var = 10))  # Prior for abundance - Narrower variance for faster convergence


# Tuning
tuning <- list(
              alpha = 0.25, # Tuning for detection 
              beta = 0.25)  # Tuning for abundance 



# -------------------------------------------------------
#                    Detection Function
# ------------------------------------------------------- 

# Half-normal detection function
fn_hn <- DS(abund.formula = ~ 1,
            det.formula = ~ 1,
            data = pc_spA_HDS,
            family = 'Poisson',
            det.func = 'halfnormal',
            transect = 'point',
            inits = inits,
            priors = priors,
            tuning = tuning,
            accept.rate = 0.43,
            n.batch = n.batch,
            batch.length = batch.length,
            n.burn = n.burn,
            n.thin = n.thin,
            n.chains = n.chains,
            n.omp.threads = 1,
            n.report = 1000,
            verbose = TRUE)


# Negative exponential detection function
fn_ne <- DS(abund.formula = ~ 1,
            det.formula = ~ 1,
            data = pc_spA_HDS,
            family = 'Poisson',
            det.func = 'negexp',
            transect = 'point',
            inits = inits,
            priors = priors,
            tuning = tuning,
            accept.rate = 0.43,
            n.batch = n.batch,
            batch.length = batch.length,
            n.burn = n.burn,
            n.thin = n.thin,
            n.chains = n.chains,
            n.omp.threads = 1,
            n.report = 1000,
            verbose = TRUE) 



## Checking convergence ##

# Half-normal detection function
plot(fn_hn, 'beta', density = FALSE) # Abundance parameters
plot(fn_hn, 'alpha', density = FALSE) # Detection parameters
fn_hn$rhat # Rhat values of 1.0 to 1.1 indicate good mixing
dev.off() # clear plots

# Negative exponential detection function
plot(fn_ne, 'beta', density = FALSE) # Abundance parameters
plot(fn_ne, 'alpha', density = FALSE) # Detection parameters
fn_ne$rhat # Rhat values of 1.0 to 1.1 indicate good mixing
dev.off() # clear plots



## Ranking Detection functions Models ##

# Calculating WAIC
waic_fn_hn <- waicAbund(fn_hn)
waic_fn_ne<- waicAbund(fn_ne)
 
# Extract the WAIC values for each model
waic_values <- c(waic_fn_hn["WAIC"],
                 waic_fn_ne["WAIC"])

# Create a named vector with model names
fitnames <- c("fn_hn", 
              "fn_ne")

# Combine model names and WAIC values into a data frame for ranking
model_waic_df <- data.frame(Model = fitnames, WAIC = waic_values)

# Rank models based on WAIC (lower WAIC is better)
model_waic_df <- model_waic_df[order(model_waic_df$WAIC), ]

# Print the ranked models
print(model_waic_df)

## Negative exponential is the better detection function


# -------------------------------------------------------
#                    Abundance Models 
# ------------------------------------------------------- 

# -------------------------------------------------------
# Abundance Fit 1: Null
# ------------------------------------------------------- 
fm.0 <- DS(abund.formula = ~ 1,
           det.formula = ~ 1,
           data = pc_spA_HDS,
           family = 'Poisson',
           det.func = 'negexp',
           transect = 'point',
           inits = inits,
           priors = priors,
           tuning = tuning,
           accept.rate = 0.43,
           n.batch = n.batch,
           batch.length = batch.length,
           n.burn = n.burn,
           n.thin = n.thin,
           n.chains = n.chains,
           n.omp.threads = 1,
           n.report = 1000,
           verbose = TRUE) 


# -------------------------------------------------------
# Abundance Fit 1: Herbaceous Proportion 
# ------------------------------------------------------- 
fm.1 <- DS(abund.formula = ~ herb_prp,
           det.formula = ~ 1,
           data = pc_spA_HDS,
           family = 'Poisson',
           det.func = 'negexp',
           transect = 'point',
           inits = inits,
           priors = priors,
           tuning = tuning,
           accept.rate = 0.43,
           n.batch = n.batch,
           batch.length = batch.length,
           n.burn = n.burn,
           n.thin = n.thin,
           n.chains = n.chains,
           n.omp.threads = 1,
           n.report = 1000,
           verbose = TRUE) 


# -------------------------------------------------------
# Abundance Fit 2: Mean Woody Patch Area 
# ------------------------------------------------------- 
fm.2 <- DS(abund.formula = ~ woody_mean_p_Area,
           det.formula = ~ 1,
           data = pc_spA_HDS,
           family = 'Poisson',
           det.func = 'negexp',
           transect = 'point',
           inits = inits,
           priors = priors,
           tuning = tuning,
           accept.rate = 0.43,
           n.batch = n.batch,
           batch.length = batch.length,
           n.burn = n.burn,
           n.thin = n.thin,
           n.chains = n.chains,
           n.omp.threads = 1,
           n.report = 1000,
           verbose = TRUE) 


# -------------------------------------------------------
# Abundance Fit 3: Woody Clumpy Index 
# ------------------------------------------------------- 
fm.3 <- DS(abund.formula = ~ woody_c_clumpy,
           det.formula = ~ 1,
           data = pc_spA_HDS,
           family = 'Poisson',
           det.func = 'negexp',
           transect = 'point',
           inits = inits,
           priors = priors,
           tuning = tuning,
           accept.rate = 0.43,
           n.batch = n.batch,
           batch.length = batch.length,
           n.burn = n.burn,
           n.thin = n.thin,
           n.chains = n.chains,
           n.omp.threads = 1,
           n.report = 1000,
           verbose = TRUE) 


# -------------------------------------------------------
# Abundance Fit 4: Herbaceous Proportion + Woody Patch
# ------------------------------------------------------- 
fm.4 <- DS(abund.formula = ~ herb_prp + woody_mean_p_Area, 
           det.formula = ~ 1,
           data = pc_spA_HDS,
           family = 'Poisson',
           det.func = 'negexp',
           transect = 'point',
           inits = inits,
           priors = priors,
           tuning = tuning,
           accept.rate = 0.43,
           n.batch = n.batch,
           batch.length = batch.length,
           n.burn = n.burn,
           n.thin = n.thin,
           n.chains = n.chains,
           n.omp.threads = 1,
           n.report = 1000,
           verbose = TRUE) 


# -------------------------------------------------------
# Abundance Fit 5: Herbaceous Proportion + Woody Clumpy
# ------------------------------------------------------- 
fm.5 <- DS(abund.formula = ~ herb_prp + woody_c_clumpy, 
           det.formula = ~ 1,
           data = pc_spA_HDS,
           family = 'Poisson',
           det.func = 'negexp',
           transect = 'point',
           inits = inits,
           priors = priors,
           tuning = tuning,
           accept.rate = 0.43,
           n.batch = n.batch,
           batch.length = batch.length,
           n.burn = n.burn,
           n.thin = n.thin,
           n.chains = n.chains,
           n.omp.threads = 1,
           n.report = 1000,
           verbose = TRUE) 


# -------------------------------------------------------
# Checking Convergence
# -------------------------------------------------------

# Inspect Trace plots
# If convergence is good and the MCMC chains are mixing well
# the plots will look like tall grass

# Abundance Fit 0: Null Model
plot(fm.0, 'beta', density = FALSE) # Abundance parameters
plot(fm.0, 'alpha', density = FALSE) # Detection parameters
fm.0$rhat # Rhat values of 1.0 to 1.1 indicate good mixing
dev.off() # clear plots

# Abundance Fit 1: Herbaceous Proportion
plot(fm.1, 'beta', density = FALSE)  
plot(fm.1, 'alpha', density = FALSE)  
fm.1$rhat
dev.off()

# Abundance Fit 2: Mean Woody Patch Area
plot(fm.2, 'beta', density = FALSE)  
plot(fm.2, 'alpha', density = FALSE)  
fm.2$rhat
dev.off()

# Abundance Fit 3: Woody Clumpy Index
plot(fm.3, 'beta', density = FALSE)  
plot(fm.3, 'alpha', density = FALSE)  
fm.3$rhat
dev.off()

# Abundance Fit 4: Herbaceous Proportion + Woody Patch
plot(fm.4, 'beta', density = FALSE)  
plot(fm.4, 'alpha', density = FALSE)  
fm.4$rhat
dev.off()

# Abundance Fit 5: Herbaceous Proportion + Woody Clumpy
plot(fm.5, 'beta', density = FALSE)  
plot(fm.5, 'alpha', density = FALSE)  
fm.5$rhat
dev.off()


# -------------------------------------------------------
# Ranking Abundance Models
# -------------------------------------------------------

# Calculating WAIC
waic_fm.0 <- waicAbund(fm.0)
waic_fm.1 <- waicAbund(fm.1)
waic_fm.2 <- waicAbund(fm.2)
waic_fm.3 <- waicAbund(fm.3)  
waic_fm.4 <- waicAbund(fm.4)
waic_fm.5 <- waicAbund(fm.5)  


# Extract the WAIC values for each model
waic_values <- c(waic_fm.0["WAIC"],
                 waic_fm.1["WAIC"],
                 waic_fm.2["WAIC"],
                 waic_fm.3["WAIC"],
                 waic_fm.4["WAIC"],
                 waic_fm.5["WAIC"])

# Create a named vector with model names
fitnames <- c("fm.0", 
              "fm.1", 
              "fm.2",
              "fm.3", 
              "fm.4",
              "fm.5")


# Combine model names and WAIC values into a data frame for ranking
model_waic_df <- data.frame(Model = fitnames, WAIC = waic_values)

# Rank models based on WAIC (lower WAIC is better)
model_waic_df <- model_waic_df[order(model_waic_df$WAIC), ]

# Print the ranked models
print(model_waic_df)

# Model 4: Herbaceous Proportion + Woody Patch is the best abundance model
summary(fm.4)

# Check model posterior predictive check
ppc1 <- ppcAbund(fm.4, fit.stat = 'freeman-tukey', group = 1)
summary(ppc1)

# Bayesian p-value = 0.6683, using group of site
# using group of 2, replicate, PPC Bayesian p-value = 0.5288
# ppc2 <- ppcAbund(fm.4, fit.stat = 'freeman-tukey', group = 2)
# summary(ppc2)

















