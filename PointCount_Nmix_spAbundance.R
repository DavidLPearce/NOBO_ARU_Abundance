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

# Site 3 is missing covariate information for survey 1 
# adding in information based on other survey 1 info
obs_mat[3,1] <- "Dave"
temp_mat[3,1] <- 76
wind_mat[3,1] <- 4
sky_mat[3,1] <- 2
doy_mat[3,1] <- 118


# -------------------------------------------------------
# Format long
# ------------------------------------------------------- 

# Formatting long
pc_spA_Nmix <- list(y = det_mat,
                    abund.covs = site_covs[,4:8],
                    det.covs = list(obsvr = obs_mat,
                                    temp = temp_mat,
                                    wind = wind_mat,
                                    sky = sky_mat,
                                    doy= doy_mat))

# Check structure
str(pc_spA_Nmix)

# Detection covariates need to be a factor
pc_spA_Nmix$det.covs$obsvr <- as.data.frame(lapply(pc_spA_Nmix$det.covs$obsvr, as.factor))
pc_spA_Nmix$det.covs$temp <- as.data.frame(lapply(pc_spA_Nmix$det.covs$temp, as.factor))
pc_spA_Nmix$det.covs$wind <- as.data.frame(lapply(pc_spA_Nmix$det.covs$wind, as.factor))
pc_spA_Nmix$det.covs$sky <- as.data.frame(lapply(pc_spA_Nmix$det.covs$sky, as.factor))
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
n.batch <- 10000 
batch.length * n.batch # Total number of MCMC samples per chain
n.burn <- 20000 
n.thin <- 10
n.chains <- 3

# Initial values for model start
inits <- list(
              alpha = 0.1,                                   # Detection covariates
              beta = 0.1,                                    # Single abundance intercept
              N = apply(pc_spA_Nmix$y, 1, max, na.rm = TRUE)) # Abundance


# Set Priors
priors <- list(
                 alpha.normal = list(mean = 0, var = 10),   # Prior for detection - Narrower variance for faster convergence
                 beta.normal = list(mean = 0, var = 10))   # Prior for abundance - Narrower variance for faster convergence
         

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
              n.report = 1000,
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
                n.report = 1000,
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
                n.report = 1000,
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
                n.report = 1000,
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
                n.report = 1000,
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
                n.report = 1000,
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
                n.report = 1000,
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

## Detection model 5 (day of year), is the best detection model
summary(detfm.5)

# -------------------------------------------------------
#                    Abundance Models
# ------------------------------------------------------- 

# -------------------------------------------------------
# Abundance Fit 1: Null
# ------------------------------------------------------- 
fm.0 <- NMix(abund.formula = ~ 1, 
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
             n.report = 1000,
             verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 1: Herbaceous Proportion 
# ------------------------------------------------------- 
fm.1 <- NMix(abund.formula = ~ herb_prp, 
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
              n.report = 1000,
              verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 2: Mean Woody Patch Area 
# ------------------------------------------------------- 
fm.2 <- NMix(abund.formula = ~ woody_mean_p_Area, 
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
             n.report = 1000,
             verbose = TRUE)

# -------------------------------------------------------
# Abundance Fit 3: Woody Clumpy Index 
# ------------------------------------------------------- 
fm.3 <- NMix(abund.formula = ~ woody_c_clumpy, 
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
             n.report = 1000,
             verbose = TRUE)


# -------------------------------------------------------
# Abundance Fit 4: Herbaceous Proportion + Woody Patch
# ------------------------------------------------------- 
fm.4 <- NMix(abund.formula = ~ herb_prp + woody_mean_p_Area, 
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
             n.report = 1000,
             verbose = TRUE)


# -------------------------------------------------------
# Abundance Fit 5: Herbaceous Proportion + Woody Clumpy
# ------------------------------------------------------- 
fm.5 <- NMix(abund.formula = ~ herb_prp + woody_c_clumpy, 
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


# -------------------------------------------------------
#                    Abundance
# -------------------------------------------------------

# Mean abundance per point
print(mean(fm.4$N.samples)) # Latent Abundance
print(mean(fm.4$mu.samples)) # Expected abundance

# Area = π * r^2
# Point counts were assumed to have a effective sampling radius of 
# 200 meters which is 656.2 feet (1 meter = 3.281 feet).
# So the Area =  pi * 656.2 ^2 = 1352765 feet squared.
# Converting Area from ft^2 to acres is 1352765 ft^2 / 43560 ft^2 per acre
# Area = 31.05521 acres
area = 31.05521

# So the mean density per acre is 
nmix_density <- mean(fm.4$N.samples) / area
print(nmix_density)

# The study area has a acreage of 2710 acres
study_area = 2710

# The abundance across the study area is 
study_area_abund = nmix_density * study_area
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

# Plot
ggplot(lat_dens_df, aes(x = Model, y = Density, fill = Model)) +
  geom_violin() + 
  geom_boxplot(aes(x = Model, y = Density), 
               width = 0.2, position = position_dodge(width = 0.8)) +
  labs(
    title = "Latent Density", 
    x = "Model", 
    y = "Density") +
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, by = 0.25), labels = scales::comma) + # Customize y-axis
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


