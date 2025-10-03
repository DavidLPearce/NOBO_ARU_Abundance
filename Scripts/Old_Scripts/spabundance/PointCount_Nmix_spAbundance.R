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

# -------------------------------------------------------
# Detection Matrix
# ------------------------------------------------------- 

# Summarize detections by site and survey
det_mat <- pc_dat %>%
  group_by(PointNum, Survey) %>% # grouping detections by site and survey
  summarise(Detections = n(), .groups = 'drop') %>%
  pivot_wider(names_from = Survey, values_from = Detections, values_fill = 0) %>% # site x survey
  rename_with(~ paste0("y.", .), -PointNum) %>% # colnames to y.1 ... y.4
  column_to_rownames(var = "PointNum") # rownames are point count ID 

# View the result
print(det_mat)

# -------------------------------------------------------
# Observation Matrix
# ------------------------------------------------------- 

# Format date
pc_dat$Date <- as.POSIXct(pc_dat$Date, format = "%m/%d/%Y")

# Create a day of year column
pc_dat$DOY <- yday(pc_dat$Date)

# Function for organizing detection matrix
cov_matFUN <- function(data, covariate_name) {
  data %>%
    group_by(PointNum, Survey) %>%
    summarise(Value = first(.data[[covariate_name]]), .groups = "drop") %>%
    pivot_wider(names_from = Survey, values_from = Value, names_prefix = paste0(covariate_name, ".")) %>%
    column_to_rownames(var = "PointNum")
} # end function

# Create matrices for each covariate
obs_mat <- cov_matFUN(pc_dat, "Observer")
temp_mat <- cov_matFUN(pc_dat, "Temp.deg.F")
wind_mat <- cov_matFUN(pc_dat, "Wind.Beau.Code")
sky_mat <- cov_matFUN(pc_dat, "Sky.Beau.Code")
doy_mat <- cov_matFUN(pc_dat, "DOY")

# View the Observer matrix as an example
print(obs_mat)
print(temp_mat)
print(wind_mat)
print(sky_mat)
print(doy_mat)



# -------------------------------------------------------
# Format long
# ------------------------------------------------------- 

# Formatting long
pc_spA_Nmix <- list(y = det_mat,
                    abund.covs = site_covs[,-c(1:4)],
                    det.covs = list(obsvr = obs_mat,
                                    temp = temp_mat,
                                    wind = wind_mat,
                                    sky = sky_mat,
                                    doy= doy_mat))

# Check structure
str(pc_spA_Nmix)

# Observor, Wind, Sky covariates need to be a factor
pc_spA_Nmix$det.covs$obsvr <- as.data.frame(lapply(pc_spA_Nmix$det.covs$obsvr, as.factor))
pc_spA_Nmix$det.covs$wind <- as.data.frame(lapply(pc_spA_Nmix$det.covs$wind, as.factor))
pc_spA_Nmix$det.covs$doy <- as.data.frame(lapply(pc_spA_Nmix$det.covs$doy, as.factor))

# Recheck structure
str(pc_spA_Nmix)



# -------------------------------------------------------
#
#                   N-mixture Models
#
# -------------------------------------------------------

# -------------------------------------------------------
#                    Detection Model
# ------------------------------------------------------- 


##  Determining best detection model ##

# MCMC Specifications
batch.length <- 20
n.batch <- 25000
batch.length * n.batch # Total number of MCMC samples per chain
n.burn <- 50000
n.thin <- 10
n.chains <- 3

# Initial values for model start
inits <- list(
              alpha = 0.1,                                   # Detection covariates
              beta = 0.1,                                    # Single abundance intercept
              N = apply(pc_spA_Nmix$y, 1, max, na.rm = TRUE)) # Abundance


# Set Priors
priors <- list(
               alpha.normal = list(mean = 0, var = 10),   # Prior for detection
               beta.normal = list(mean = 0, var = 10))   # Prior for abundance
         

# Tuning
tuning <- list(
                 alpha = 0.25,        # Tuning for detection 
                 beta = 0.25)         # Tuning for abundance 

# -------------------------------------------------------
# Detection Fit 0: Null model 
# ------------------------------------------------------- 
detfm.0 <- NMix(abund.formula = ~ 1, 
              det.formula = ~ 1, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Detection Fit 1: Observer
# ------------------------------------------------------- 
detfm.1 <- NMix(abund.formula = ~ 1, 
                det.formula = ~ obsvr, 
                data = pc_spA_Nmix,
                family = 'Poisson',
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
                n.report = 5000,
                verbose = TRUE)

# -------------------------------------------------------
# Detection Fit 2: Temperature
# ------------------------------------------------------- 
detfm.2 <- NMix(abund.formula = ~ 1, 
                det.formula = ~ temp, 
                data = pc_spA_Nmix,
                family = 'Poisson',
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
                n.report = 5000,
                verbose = TRUE)

# -------------------------------------------------------
# Detection Fit 3: Wind
# ------------------------------------------------------- 
detfm.3 <- NMix(abund.formula = ~ 1, 
                det.formula = ~ wind, 
                data = pc_spA_Nmix,
                family = 'Poisson',
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
                n.report = 5000,
                verbose = TRUE)

# -------------------------------------------------------
# Detection Fit 4: Sky
# ------------------------------------------------------- 
detfm.4 <- NMix(abund.formula = ~ 1, 
                det.formula = ~ sky, 
                data = pc_spA_Nmix,
                family = 'Poisson',
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
                n.report = 5000,
                verbose = TRUE)

# -------------------------------------------------------
# Detection Fit 5: Day of Year
# ------------------------------------------------------- 
detfm.5 <- NMix(abund.formula = ~ 1, 
                det.formula = ~ doy, 
                data = pc_spA_Nmix,
                family = 'Poisson',
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
                n.report = 5000,
                verbose = TRUE)

# -------------------------------------------------------
# Detection Fit 6: Day of Year + Observer
# ------------------------------------------------------- 
detfm.6 <- NMix(abund.formula = ~ 1, 
                det.formula = ~ doy + obsvr, 
                data = pc_spA_Nmix,
                family = 'Poisson',
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
                n.report = 5000,
                verbose = TRUE)

# -------------------------------------------------------
# Checking Convergence
# -------------------------------------------------------

# Inspect Trace plots
# If convergence is good and the MCMC chains are mixing well
# the plots will look like tall grass

# Detection Fit 0: Null Model
plot(detfm.0, 'beta', density = FALSE) # Abundance parameters
plot(detfm.0, 'alpha', density = FALSE) # Detection parameters
detfm.0$rhat # Rhat values of 1.0 to 1.1 indicate good mixing
dev.off() # clear plots

# Detection Fit 1: Observer
plot(detfm.1, 'beta', density = FALSE)  
plot(detfm.1, 'alpha', density = FALSE)  
detfm.1$rhat
dev.off()

# Detection Fit 2: Temperature
plot(detfm.2, 'beta', density = FALSE)  
plot(detfm.2, 'alpha', density = FALSE) 
detfm.2$rhat
dev.off()

# Detection Fit 3: Wind
plot(detfm.3, 'beta', density = FALSE)  
plot(detfm.3, 'alpha', density = FALSE)  
detfm.3$rhat
dev.off()

# Detection Fit 4: Sky
plot(detfm.4, 'beta', density = FALSE)  
plot(detfm.4, 'alpha', density = FALSE) 
detfm.4$rhat
dev.off()

# Detection Fit 5: Day of Year
plot(detfm.5, 'beta', density = FALSE)  
plot(detfm.5, 'alpha', density = FALSE)  
detfm.5$rhat
dev.off()

# Detection Fit 6: Day of Year + Observer
plot(detfm.6, 'beta', density = FALSE)  
plot(detfm.6, 'alpha', density = FALSE)  
detfm.6$rhat
dev.off()


# -------------------------------------------------------
# Ranking Detection Models
# -------------------------------------------------------


# Calculating WAIC
waic_detfm.0 <- waicAbund(detfm.0) # Message is for setting upper limit for integration 
waic_detfm.1 <- waicAbund(detfm.1)
waic_detfm.2 <- waicAbund(detfm.2)
waic_detfm.3 <- waicAbund(detfm.3)  
waic_detfm.4 <- waicAbund(detfm.4)
waic_detfm.5 <- waicAbund(detfm.5)  
waic_detfm.6 <- waicAbund(detfm.6)

# Extract the WAIC values for each model
waic_detvalues <- c(waic_detfm.0["WAIC"],
                    waic_detfm.1["WAIC"],
                    waic_detfm.2["WAIC"],
                    waic_detfm.3["WAIC"],
                    waic_detfm.4["WAIC"],
                    waic_detfm.5["WAIC"],
                    waic_detfm.6["WAIC"])

# Create a named vector with model names
detfitnames <- c("detfm.0", 
                 "detfm.1", 
                 "detfm.2",
                 "detfm.3", 
                 "detfm.4",
                 "detfm.5",
                 "detfm.6")


# Combine model names and WAIC values into a data frame for ranking
detmodel_waic_df <- data.frame(Model = detfitnames, WAIC = waic_detvalues)

# Rank models based on WAIC (lower WAIC is better)
detmodel_waic_df <- detmodel_waic_df[order(detmodel_waic_df$WAIC), ]

# Print the ranked models
print(detmodel_waic_df)

## Detection model 1 observer, is the best detection model
summary(detfm.1)

# -------------------------------------------------------
#                    Abundance Models
# ------------------------------------------------------- 

## Getting woody and herb cov model combinations ##

# Select woody and herbaceous covariates
woody_covs <- site_covs %>% select(starts_with("woody"))
herb_covs <- site_covs %>% select(starts_with("herb"))

# Get column names
woody_vars <- colnames(woody_covs)
herb_vars <- colnames(herb_covs)

# Generate all combinations
model_combinations <- expand.grid(Woody = woody_vars, Herbaceous = herb_vars, stringsAsFactors = FALSE)

# View the result
print(model_combinations)

# -------------------------------------------------------
# Abundance Fit 0: Null
# ------------------------------------------------------- 
fm.0 <- NMix(abund.formula = ~ 1, 
             det.formula = ~ obsvr, 
             data = pc_spA_Nmix,
             family = 'Poisson',
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
             n.report = 5000,
             verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 1: Woody Proportion 
# ------------------------------------------------------- 
fm.1 <- NMix(abund.formula = ~ woody_prp, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 2: Mean Woody Patch Area 
# ------------------------------------------------------- 
fm.2 <- NMix(abund.formula = ~ woody_Parea, 
             det.formula = ~ obsvr, 
             data = pc_spA_Nmix,
             family = 'Poisson',
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
             n.report = 5000,
             verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 3: Woody Clumpy Index 
# ------------------------------------------------------- 
fm.3 <- NMix(abund.formula = ~ woody_ClmIdx, 
             det.formula = ~ obsvr, 
             data = pc_spA_Nmix,
             family = 'Poisson',
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
             n.report = 5000,
             verbose = TRUE)


# -------------------------------------------------------
# Abundance Fit 4: Woody Patch Density
# ------------------------------------------------------- 
fm.4 <- NMix(abund.formula = ~ scale(woody_Pdens), 
             det.formula = ~ obsvr, 
             data = pc_spA_Nmix,
             family = 'Poisson',
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
             n.report = 5000,
             verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 5: Woody Shape Index
# ------------------------------------------------------- 
fm.5 <- NMix(abund.formula = ~ woody_ShpInx, 
             det.formula = ~ obsvr, 
             data = pc_spA_Nmix,
             family = 'Poisson',
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
             n.report = 5000,
             verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 6: Woody Largest Patch Index
# ------------------------------------------------------- 
fm.6 <- NMix(abund.formula = ~ scale(woody_lrgPInx), 
             det.formula = ~ obsvr, 
             data = pc_spA_Nmix,
             family = 'Poisson',
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
             n.report = 5000,
             verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 7: Woody Aggregation Index
# ------------------------------------------------------- 
fm.7 <- NMix(abund.formula = ~ scale(woody_AggInx), 
             det.formula = ~ obsvr, 
             data = pc_spA_Nmix,
             family = 'Poisson',
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
             n.report = 5000,
             verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 8: Woody edge Density
# ------------------------------------------------------- 
fm.8 <- NMix(abund.formula = ~ scale(woody_EdgDens), 
             det.formula = ~ obsvr, 
             data = pc_spA_Nmix,
             family = 'Poisson',
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
             n.report = 5000,
             verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 9: herb Proportion 
# ------------------------------------------------------- 
fm.9 <- NMix(abund.formula = ~ herb_prp, 
             det.formula = ~ obsvr, 
             data = pc_spA_Nmix,
             family = 'Poisson',
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
             n.report = 5000,
             verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 10: Mean herb Patch Area 
# ------------------------------------------------------- 
fm.10 <- NMix(abund.formula = ~ herb_Parea, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 11: herb Clumpy Index 
# ------------------------------------------------------- 
fm.11 <- NMix(abund.formula = ~ herb_ClmIdx, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)


# -------------------------------------------------------
# Abundance Fit 12: herb Patch Density
# ------------------------------------------------------- 
fm.12 <- NMix(abund.formula = ~ scale(herb_Pdens), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 13: herb Shape Index
# ------------------------------------------------------- 
fm.13 <- NMix(abund.formula = ~ herb_ShpInx, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 14: herb Largest Patch Index
# ------------------------------------------------------- 
fm.14 <- NMix(abund.formula = ~ scale(herb_lrgPInx), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 15: herb Aggregation Index
# ------------------------------------------------------- 
fm.15 <- NMix(abund.formula = ~ scale(herb_AggInx), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 16: herb Edge Density
# ------------------------------------------------------- 
fm.16 <- NMix(abund.formula = ~ scale(herb_EdgDens), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 17:  woody_Parea + herb_prp
# ------------------------------------------------------- 
fm.17 <- NMix(abund.formula = ~ woody_Parea + herb_prp, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 18: woody_ClmIdx + herb_prp
# ------------------------------------------------------- 
fm.18 <- NMix(abund.formula = ~ woody_ClmIdx + herb_prp, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 19: scale(woody_Pdens) + herb_prp
# ------------------------------------------------------- 
fm.19 <- NMix(abund.formula = ~ scale(woody_Pdens) + herb_prp, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 20: woody_ShpInx +  herb_prp
# ------------------------------------------------------- 
fm.20 <- NMix(abund.formula = ~ woody_ShpInx + herb_prp, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)


# -------------------------------------------------------
# Abundance Fit 21: scale(woody_lrgPInx) + herb_prp
# ------------------------------------------------------- 
fm.21 <- NMix(abund.formula = ~ scale(woody_lrgPInx) + herb_prp, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 22: scale(woody_AggInx) + herb_prp
# ------------------------------------------------------- 
fm.22 <- NMix(abund.formula = ~ scale(woody_AggInx) + herb_prp, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 23: scale(woody_EdgDens) + herb_prp
# ------------------------------------------------------- 
fm.23 <- NMix(abund.formula = ~ scale(woody_EdgDens) + herb_prp, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 24: woody_prp +	herb_Parea
# ------------------------------------------------------- 
fm.24 <- NMix(abund.formula = ~ woody_prp +	herb_Parea, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 25: 
# ------------------------------------------------------- 
fm.25 <- NMix(abund.formula = ~ woody_Parea +	herb_Parea, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 26: woody_ClmIdx	+	herb_Parea
# ------------------------------------------------------- 
fm.26 <- NMix(abund.formula = ~ woody_ClmIdx	+	herb_Parea, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 27: scale(woody_Pdens)	+	herb_Parea
# ------------------------------------------------------- 
fm.27 <- NMix(abund.formula = ~ scale(woody_Pdens)	+	herb_Parea, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 28: woody_ShpInx	+	herb_Parea
# ------------------------------------------------------- 
fm.28 <- NMix(abund.formula = ~ woody_ShpInx	+	herb_Parea, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)


# -------------------------------------------------------
# Abundance Fit 29: scale(woody_lrgPInx)	+	herb_Parea
# ------------------------------------------------------- 
fm.29 <- NMix(abund.formula = ~ scale(woody_lrgPInx)	+	herb_Parea, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 30: scale(woody_AggInx)	+	herb_Parea
# ------------------------------------------------------- 
fm.30 <- NMix(abund.formula = ~ scale(woody_AggInx)	+	herb_Parea, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 31: scale(woody_EdgDens)	+	herb_Parea
# ------------------------------------------------------- 
fm.31 <- NMix(abund.formula = ~ scale(woody_EdgDens)	+	herb_Parea, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 32: woody_prp	+	herb_ClmIdx
# ------------------------------------------------------- 
fm.32 <- NMix(abund.formula = ~ woody_prp	+	herb_ClmIdx, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 33: woody_Parea	+	herb_ClmIdx
# ------------------------------------------------------- 
fm.33 <- NMix(abund.formula = ~ woody_Parea	+	herb_ClmIdx, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 34: woody_ClmIdx	+	herb_ClmIdx
# ------------------------------------------------------- 
fm.34 <- NMix(abund.formula = ~ woody_ClmIdx	+	herb_ClmIdx, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)


# -------------------------------------------------------
# Abundance Fit 35: scale(woody_Pdens)	+	herb_ClmIdx
# ------------------------------------------------------- 
fm.35 <- NMix(abund.formula = ~ scale(woody_Pdens)	+	herb_ClmIdx, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 36: woody_ShpInx	+	herb_ClmIdx
# ------------------------------------------------------- 
fm.36 <- NMix(abund.formula = ~ woody_ShpInx	+	herb_ClmIdx, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)


# -------------------------------------------------------
# Abundance Fit 37: scale(woody_lrgPInx)	+	herb_ClmIdx
# ------------------------------------------------------- 
fm.37 <- NMix(abund.formula = ~ scale(woody_lrgPInx)	+	herb_ClmIdx, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 38: scale(woody_AggInx)	+	herb_ClmIdx
# ------------------------------------------------------- 
fm.38 <- NMix(abund.formula = ~ scale(woody_AggInx)	+	herb_ClmIdx, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 39: scale(woody_EdgDens)	+	herb_ClmIdx
# ------------------------------------------------------- 
fm.39 <- NMix(abund.formula = ~ scale(woody_EdgDens)	+	herb_ClmIdx, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 40: woody_prp	+	scale(herb_Pdens)
# ------------------------------------------------------- 
fm.40 <- NMix(abund.formula = ~ woody_prp	+	scale(herb_Pdens), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 41: woody_Parea	+	scale(herb_Pdens)
# ------------------------------------------------------- 
fm.41 <- NMix(abund.formula = ~ woody_Parea	+	scale(herb_Pdens), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 42: woody_ClmIdx	+	scale(herb_Pdens)
# ------------------------------------------------------- 
fm.42 <- NMix(abund.formula = ~ woody_ClmIdx	+	scale(herb_Pdens), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 43: scale(woody_Pdens)	+	scale(herb_Pdens)
# ------------------------------------------------------- 
fm.43 <- NMix(abund.formula = ~ scale(woody_Pdens)	+	scale(herb_Pdens), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 44: woody_ShpInx	+	scale(herb_Pdens)
# ------------------------------------------------------- 
fm.44 <- NMix(abund.formula = ~ woody_ShpInx	+	scale(herb_Pdens), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)


# -------------------------------------------------------
# Abundance Fit 45: scale(woody_lrgPInx)	+	scale(herb_Pdens)
# ------------------------------------------------------- 
fm.45 <- NMix(abund.formula = ~ scale(woody_lrgPInx)	+	scale(herb_Pdens), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 46: scale(woody_AggInx)	+	scale(herb_Pdens)
# ------------------------------------------------------- 
fm.46 <- NMix(abund.formula = ~ scale(woody_AggInx)	+	scale(herb_Pdens), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 47: scale(woody_EdgDens)	+	scale(herb_Pdens)
# ------------------------------------------------------- 
fm.47 <- NMix(abund.formula = ~ scale(woody_EdgDens)	+	scale(herb_Pdens), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 48: woody_prp	+	herb_ShpInx
# ------------------------------------------------------- 
fm.48 <- NMix(abund.formula = ~ woody_prp	+	herb_ShpInx, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 49: woody_Parea	+	herb_ShpInx
# ------------------------------------------------------- 
fm.49 <- NMix(abund.formula = ~ woody_Parea	+	herb_ShpInx, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 50:  woody_ClmIdx	+	herb_ShpInx
# ------------------------------------------------------- 
fm.50 <- NMix(abund.formula = ~ woody_ClmIdx	+	herb_ShpInx, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 51:  scale(woody_Pdens)	+	herb_ShpInx
# ------------------------------------------------------- 
fm.51 <- NMix(abund.formula = ~ scale(woody_Pdens)	+	herb_ShpInx, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 52:  woody_ShpInx	+	herb_ShpInx

# ------------------------------------------------------- 
fm.52 <- NMix(abund.formula = ~  woody_ShpInx	+	herb_ShpInx, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 53:  scale(woody_lrgPInx)	+	herb_ShpInx
# ------------------------------------------------------- 
fm.53 <- NMix(abund.formula = ~  scale(woody_lrgPInx)	+	herb_ShpInx, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 54:  scale(woody_AggInx)	+	herb_ShpInx
# ------------------------------------------------------- 
fm.54 <- NMix(abund.formula = ~ scale(woody_AggInx)	+	herb_ShpInx, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 55: scale(woody_EdgDens)	+	herb_ShpInx
# ------------------------------------------------------- 
fm.55 <- NMix(abund.formula = ~  scale(woody_EdgDens)	+	herb_ShpInx, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 56:  woody_prp	+	scale(herb_lrgPInx)
# ------------------------------------------------------- 
fm.56 <- NMix(abund.formula = ~  woody_prp	+	scale(herb_lrgPInx), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 57:  woody_Parea	+	scale(herb_lrgPInx)
# ------------------------------------------------------- 
fm.57 <- NMix(abund.formula = ~  woody_Parea	+	scale(herb_lrgPInx), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 58: woody_ClmIdx	+	scale(herb_lrgPInx)
# ------------------------------------------------------- 
fm.58 <- NMix(abund.formula = ~  woody_ClmIdx	+	scale(herb_lrgPInx), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)


# -------------------------------------------------------
# Abundance Fit 59:  scale(woody_Pdens)	+	scale(herb_lrgPInx)
# ------------------------------------------------------- 
fm.59 <- NMix(abund.formula = ~  scale(woody_Pdens)	+	scale(herb_lrgPInx), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 60: woody_ShpInx	+	scale(herb_lrgPInx)
# ------------------------------------------------------- 
fm.60 <- NMix(abund.formula = ~  woody_ShpInx	+	scale(herb_lrgPInx), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 61: scale(woody_lrgPInx)	+	scale(herb_lrgPInx)
# ------------------------------------------------------- 
fm.61 <- NMix(abund.formula = ~  scale(woody_lrgPInx)	+	scale(herb_lrgPInx),
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 62:  scale(woody_AggInx)	+	scale(herb_lrgPInx)
# ------------------------------------------------------- 
fm.62 <- NMix(abund.formula = ~ scale(woody_AggInx)	+	scale(herb_lrgPInx), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 63:  scale(woody_EdgDens)	+	scale(herb_lrgPInx)
# ------------------------------------------------------- 
fm.63 <- NMix(abund.formula = ~  scale(woody_EdgDens)	+	scale(herb_lrgPInx), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 64:  woody_prp	+	scale(herb_AggInx)
# ------------------------------------------------------- 
fm.64 <- NMix(abund.formula = ~  woody_prp	+	scale(herb_AggInx), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)


# -------------------------------------------------------
# Abundance Fit 65: woody_Parea	+	scale(herb_AggInx)
# ------------------------------------------------------- 
fm.65 <- NMix(abund.formula = ~  woody_Parea	+	scale(herb_AggInx), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 66:  woody_ClmIdx	+	scale(herb_AggInx)
# ------------------------------------------------------- 
fm.66 <- NMix(abund.formula = ~ woody_ClmIdx	+	scale(herb_AggInx), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)


# -------------------------------------------------------
# Abundance Fit 67:  scale(woody_Pdens)	+	scale(herb_AggInx)
# ------------------------------------------------------- 
fm.67 <- NMix(abund.formula = ~  scale(woody_Pdens)	+	scale(herb_AggInx),
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 68:  woody_ShpInx	+	scale(herb_AggInx)
# ------------------------------------------------------- 
fm.68 <- NMix(abund.formula = ~  woody_ShpInx	+	scale(herb_AggInx), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 69:  scale(woody_lrgPInx)	+	scale(herb_AggInx)
# ------------------------------------------------------- 
fm.69 <- NMix(abund.formula = ~  scale(woody_lrgPInx)	+	scale(herb_AggInx), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 70:  scale(woody_AggInx)	+	scale(herb_AggInx)
# ------------------------------------------------------- 
fm.70 <- NMix(abund.formula = ~ scale(woody_AggInx)	+	scale(herb_AggInx), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 71:  scale(woody_EdgDens)	+	scale(herb_AggInx)
# ------------------------------------------------------- 
fm.71 <- NMix(abund.formula = ~  scale(woody_EdgDens)	+	scale(herb_AggInx), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 72:  woody_prp	+	scale(herb_EdgDens)	
# ------------------------------------------------------- 
fm.72 <- NMix(abund.formula = ~  woody_prp	+	scale(herb_EdgDens)	, 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 73:  woody_Parea	+	scale(herb_EdgDens)
# ------------------------------------------------------- 
fm.73 <- NMix(abund.formula = ~  woody_Parea	+	scale(herb_EdgDens), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 74:  woody_ClmIdx	+	scale(herb_EdgDens)
# ------------------------------------------------------- 
fm.74 <- NMix(abund.formula = ~  woody_ClmIdx	+	scale(herb_EdgDens), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)


# -------------------------------------------------------
# Abundance Fit 75:  scale(woody_Pdens)	+	scale(herb_EdgDens)
# ------------------------------------------------------- 
fm.75 <- NMix(abund.formula = ~  scale(woody_Pdens)	+	scale(herb_EdgDens), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 76:  woody_ShpInx	+	scale(herb_EdgDens)
# ------------------------------------------------------- 
fm.76 <- NMix(abund.formula = ~ woody_ShpInx	+	scale(herb_EdgDens), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 77:  scale(woody_lrgPInx)	+	scale(herb_EdgDens)
# ------------------------------------------------------- 
fm.77 <- NMix(abund.formula = ~  scale(woody_lrgPInx)	+	scale(herb_EdgDens), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 78: scale(woody_AggInx)	+	scale(herb_EdgDens)
# ------------------------------------------------------- 
fm.78 <- NMix(abund.formula = ~  scale(woody_AggInx)	+	scale(herb_EdgDens), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 79:  scale(woody_EdgDens)	+	scale(herb_EdgDens)
# ------------------------------------------------------- 
fm.79 <- NMix(abund.formula = ~  scale(woody_EdgDens)	+	scale(herb_EdgDens), 
              det.formula = ~ obsvr, 
              data = pc_spA_Nmix,
              family = 'Poisson',
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
              n.report = 5000,
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

# Abundance Fit 1 
plot(fm.1, 'beta', density = FALSE)  
plot(fm.1, 'alpha', density = FALSE)  
fm.1$rhat
dev.off()

# Abundance Fit 2 
plot(fm.2, 'beta', density = FALSE)  
plot(fm.2, 'alpha', density = FALSE)  
fm.2$rhat
dev.off()

# Abundance Fit 3 
plot(fm.3, 'beta', density = FALSE)  
plot(fm.3, 'alpha', density = FALSE)  
fm.3$rhat
dev.off()

# Abundance Fit 4 
plot(fm.4, 'beta', density = FALSE)  
plot(fm.4, 'alpha', density = FALSE)  
fm.4$rhat
dev.off()

# Abundance Fit 5 
plot(fm.5, 'beta', density = FALSE)  
plot(fm.5, 'alpha', density = FALSE)  
fm.5$rhat
dev.off()

# Abundance Fit 6 
plot(fm.6, 'beta', density = FALSE)  
plot(fm.6, 'alpha', density = FALSE)  
fm.6$rhat
dev.off()

# Abundance Fit 7 
plot(fm.7, 'beta', density = FALSE)  
plot(fm.7, 'alpha', density = FALSE)  
fm.7$rhat
dev.off()

# Abundance Fit 8 
plot(fm.8, 'beta', density = FALSE)  
plot(fm.8, 'alpha', density = FALSE)  
fm.8$rhat
dev.off()

# -------------------------------------------------------
# Ranking Abundance Models
# -------------------------------------------------------

# Calculating WAIC for models fm.0 to fm.79
waic_values <- numeric(80)  # Create an empty numeric vector for WAIC values
fitnames <- character(80)   # Create an empty character vector for model names

for (i in 0:79) {
  model_name <- paste0("fm.", i)
  waic_var_name <- paste0("waic_fm.", i)
  
  assign(waic_var_name, waicAbund(get(model_name)))  # Calculate WAIC
  waic_values[i + 1] <- get(waic_var_name)["WAIC"]  # Store WAIC value
  fitnames[i + 1] <- model_name  # Store model name
}


# Combine model names and WAIC values into a data frame for ranking
model_waic_df <- data.frame(Model = fitnames, WAIC = waic_values)

# Rank models based on WAIC (lower WAIC is better)
model_waic_df <- model_waic_df[order(model_waic_df$WAIC), ]

# Print the ranked models
print(model_waic_df)

# Best Model 
summary(fm.3)

# Check model posterior predictive check
ppc_fm.3 <- ppcAbund(fm.3, fit.stat = 'freeman-tukey', group = 1)
summary(ppc_fm.3)

# Bayesian p-value = 0.8226, using group of site



# -------------------------------------------------------
#                    Abundance
# -------------------------------------------------------

# Mean abundance per point
print(mean(fm.3$N.samples)) # Latent Abundance
print(mean(fm.3$mu.samples)) # Expected abundance

# Area in hectares
area <- pi*(200^2)/10000

# So the mean density per acre is 
nmix_density <- mean(fm.4$N.samples) / area
print(nmix_density)


# The abundance across the study area is 
study_area_abund = nmix_density * 1096.698
print(study_area_abund)

# Creating a matrix of latent density
lat_dens_mat <- (fm.4$N.samples) / area
head(lat_dens_mat)

# Matrix is density per point.
# stacking each point column to one column for plotting
lat_dens_df <- as.data.frame(lat_dens_mat) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Density")

# The point that the estimate came from doesn't matter since comparison is across models
colnames(lat_dens_df)[1] <- "Model"
lat_dens_df[,1] <- "PC Nmix"
head(lat_dens_df)

# Export density dataframe
saveRDS(lat_dens_df, "./Data/Fitted_Models/PC_Nmix_Dens.rds")

# In case of crash
# lat_dens_df <- readRDS("./Data/Fitted_Models/PC_Nmix_Dens.rds")

# Plot
ggplot(lat_dens_df, aes(x = Model, y = Density, fill = Model)) +
  geom_violin() + 
  geom_boxplot(aes(x = Model, y = Density), 
               width = 0.2, position = position_dodge(width = 0.8)) +
  labs(
    title = "Latent Density", 
    x = "Model", 
    y = "Density") +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2), labels = scales::comma) + # Customize y-axis
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1), # Tilt x-axis text
    axis.text.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none" 
  )

# Export Point Count N-mixture Latent Density
saveRDS(lat_dens_df, "./Data/Fitted_Models/PC_Nmix_Dens.rds")


