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

# Take a look at Point Count data structure
str(pc_dat)

# Format date
pc_dat$Date <- as.POSIXct(pc_dat$Date, format = "%m/%d/%Y")

# Create a day of year column
pc_dat$DOY <- yday(pc_dat$Date)

# Stack data by PointNum_survey 
pcDist_dat <- pc_dat %>%
      group_by(PointNum, Survey, DistBin) %>% # Organize data by Point_Survey and distance bin
      summarise(Count = n(), .groups = 'drop') %>% # Summarize counts in dist bins 
      pivot_wider(names_from = DistBin, values_from = Count, values_fill = 0) %>% # formating wide
      mutate(PointNum_Survey = paste(PointNum, Survey, sep = "_")) %>% # Point_Survey
      left_join(                                             # Adding covariates
        pc_dat %>%
          select(PointNum, Survey, Observer, Temp.deg.F, Wind.Beau.Code, Sky.Beau.Code, DOY) %>%
          distinct(),
        by = c("PointNum", "Survey")
      ) %>%
      rename_with(~paste0("DistBin.", .), matches("^[0-9]+$")) %>% # Naming distance bins
      select(PointNum_Survey, PointNum, Observer, Survey, Temp = Temp.deg.F, # Renaming covariate columns
             Wind = Wind.Beau.Code, Sky = Sky.Beau.Code, DOY, everything()) %>%
      replace(is.na(.), 0) %>% # Replace NA with 0 for counts
      select(-any_of("NA"))   # Remove the NA column

# Converting to a dataframe
pcDist_dat <- as.data.frame(pcDist_dat)

# Take a look
print(pcDist_dat)
str(pcDist_dat)


# -------------------------------------------------------
# Detection Matrix
# ------------------------------------------------------- 

# Subset distance data
det_mat <- pcDist_dat[,9:11]

# name row names using PointNum_survey
row.names(det_mat) <- pcDist_dat$PointNum_Survey

# Take a look
head(det_mat)

# -------------------------------------------------------
# Covariates Matrix
# ------------------------------------------------------- 

# Adding site covariates to pcDist_dat using PointNum
pcDist_dat <- pcDist_dat %>%
  left_join(site_covs, by = "PointNum")

# Take a look
head(pcDist_dat)

# Subsetting covariate matrix
cov_mat <- pcDist_dat[,c(2:8, 14:18)]

# Take a look
head(cov_mat)


# -------------------------------------------------------
# Format long
# ------------------------------------------------------- 

# Formatting long
pc_spA_HDS <- list(y = det_mat,
                   covs = cov_mat,
                   dist.breaks = c(0, 0.05, 0.1, 0.2),
                   offset = 31.05521) # in acres

# Detection covariates need to be a factor
pc_spA_HDS$covs$Observer <- factor(pc_spA_HDS$covs$Observer)
pc_spA_HDS$covs$Survey <- factor(pc_spA_HDS$covs$Survey)
pc_spA_HDS$covs$Temp <- factor(pc_spA_HDS$covs$Temp)
pc_spA_HDS$covs$Wind <- factor(pc_spA_HDS$covs$Wind)
pc_spA_HDS$covs$Sky <- factor(pc_spA_HDS$covs$Sky)
pc_spA_HDS$covs$DOY <- factor(pc_spA_HDS$covs$DOY)

# Check formatting
str(pc_spA_HDS)

# -------------------------------------------------------
#
#               Hierarchical Distance Models
#
# -------------------------------------------------------

# MCMC Specifications
batch.length <- 20
n.batch <- 25000
batch.length * n.batch # Total number of MCMC samples per chain
n.burn <- 50000
n.thin <- 10
n.chains <- 3

# Initial values for model start
inits <- list(
              alpha = 0.1,               
              beta = 0.1,                
              sigma.sq.mu = 0.1,
              N = apply(pc_spA_HDS$y, 1, sum)) 


# Set Priors
priors <- list(
              alpha.normal = list(mean = 0, var = 10),  
              beta.normal = list(mean = 0, var = 10),
              sigma.sq.mu.ig = list(mean = 0, var = 2.72))


# Tuning
tuning <- list(
              alpha = 0.25,  
              beta = 0.25,
              beta.star = 0.25)   



# -------------------------------------------------------
#                    Detection Function
# ------------------------------------------------------- 

# PointNum is included as a random effect in the abundance model
# to account for pseudoreplication in the stacked data

# Half-normal detection function
fn_hn <- DS(abund.formula = ~ (1|PointNum),
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
            n.report = 5000,
            verbose = TRUE)


# Negative exponential detection function
fn_ne <- DS(abund.formula = ~ (1|PointNum),
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
            n.report = 5000,
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
waic_fn_values <- c(waic_fn_hn["WAIC"],
                 waic_fn_ne["WAIC"])

# Create a named vector with model names
fnfitnames <- c("fn_hn", 
                "fn_ne")

# Combine model names and WAIC values into a data frame for ranking
model_fnwaic_df <- data.frame(Model = fnfitnames, WAIC = waic_fn_values)

# Rank models based on WAIC (lower WAIC is better)
model_fnwaic_df <- model_fnwaic_df[order(model_fnwaic_df$WAIC), ]

# Print the ranked models
print(model_fnwaic_df)

## Negative exponential is the better detection function



# -------------------------------------------------------
#                    Detection Models 
# ------------------------------------------------------- 


# -------------------------------------------------------
# Detection Fit 0: Null model 
# ------------------------------------------------------- 
detfm.0 <- DS(abund.formula = ~ (1|PointNum),
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
              n.report = 5000,
              verbose = TRUE) 

# -------------------------------------------------------
# Detection Fit 1: Observer
# ------------------------------------------------------- 
detfm.1 <- DS(abund.formula = ~ (1|PointNum),
              det.formula = ~ Observer,
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
              n.report = 5000,
              verbose = TRUE) 

# -------------------------------------------------------
# Detection Fit 2: Temperature
# ------------------------------------------------------- 
detfm.2 <- DS(abund.formula = ~ (1|PointNum),
              det.formula = ~ Temp,
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
              n.report = 5000,
              verbose = TRUE) 

# -------------------------------------------------------
# Detection Fit 3: Wind
# ------------------------------------------------------- 
detfm.3 <- DS(abund.formula = ~ (1|PointNum),
              det.formula = ~ Wind,
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
              n.report = 5000,
              verbose = TRUE) 

# -------------------------------------------------------
# Detection Fit 4: Sky
# ------------------------------------------------------- 
detfm.4 <- DS(abund.formula = ~ (1|PointNum),
              det.formula = ~ Sky,
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
              n.report = 5000,
              verbose = TRUE) 

# -------------------------------------------------------
# Detection Fit 5: Day of Year
# ------------------------------------------------------- 
detfm.5 <- DS(abund.formula = ~ (1|PointNum),
              det.formula = ~ DOY,
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
              n.report = 5000,
              verbose = TRUE) 

# -------------------------------------------------------
# Detection Fit 6: Day of Year + Observer
# ------------------------------------------------------- 
detfm.6 <- DS(abund.formula = ~ (1|PointNum),
              det.formula = ~ Observer + DOY,
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
              n.report = 5000,
              verbose = TRUE) 

# -------------------------------------------------------
# Detection Fit 7: Woody Proportion
# ------------------------------------------------------- 
detfm.7 <- DS(abund.formula = ~ (1|PointNum),
              det.formula = ~ woody_prp,
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

# Detection Fit 7: Woody Proportion
plot(detfm.7, 'beta', density = FALSE)  
plot(detfm.7, 'alpha', density = FALSE)  
detfm.7$rhat
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
waic_detfm.7 <- waicAbund(detfm.7)

# Extract the WAIC values for each model
waic_detvalues <- c(waic_detfm.0["WAIC"],
                    waic_detfm.1["WAIC"],
                    waic_detfm.2["WAIC"],
                    waic_detfm.3["WAIC"],
                    waic_detfm.4["WAIC"],
                    waic_detfm.5["WAIC"],
                    waic_detfm.6["WAIC"],
                    waic_detfm.7["WAIC"])

# Create a named vector with model names
detfitnames <- c("detfm.0", 
                 "detfm.1", 
                 "detfm.2",
                 "detfm.3", 
                 "detfm.4",
                 "detfm.5",
                 "detfm.6",
                 "detfm.7")


# Combine model names and WAIC values into a data frame for ranking
detmodel_waic_df <- data.frame(Model = detfitnames, WAIC = waic_detvalues)

# Rank models based on WAIC (lower WAIC is better)
detmodel_waic_df <- detmodel_waic_df[order(detmodel_waic_df$WAIC), ]

# Print the ranked models
print(detmodel_waic_df)

## Detection model 1 (Observer), is the best detection model
summary(detfm.1)


# -------------------------------------------------------
#                    Abundance Models 
# ------------------------------------------------------- 

# -------------------------------------------------------
# Abundance Fit 1: Null
# ------------------------------------------------------- 
fm.0 <- DS(abund.formula = ~ (1|PointNum),
           det.formula = ~ Observer,
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
           n.report = 5000,
           verbose = TRUE) 


# -------------------------------------------------------
# Abundance Fit 1: Herbaceous Proportion 
# ------------------------------------------------------- 
fm.1 <- DS(abund.formula = ~ herb_prp + (1|PointNum),
           det.formula = ~ Observer,
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
           n.report = 5000,
           verbose = TRUE) 


# -------------------------------------------------------
# Abundance Fit 2: Mean Woody Patch Area 
# ------------------------------------------------------- 
fm.2 <- DS(abund.formula = ~ woody_mean_p_Area + (1|PointNum),
           det.formula = ~ Observer,
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
           n.report = 5000,
           verbose = TRUE) 


# -------------------------------------------------------
# Abundance Fit 3: Woody Clumpy Index 
# ------------------------------------------------------- 
fm.3 <- DS(abund.formula = ~ woody_c_clumpy + (1|PointNum),
           det.formula = ~ Observer,
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
           n.report = 5000,
           verbose = TRUE) 


# -------------------------------------------------------
# Abundance Fit 4: Herbaceous Proportion + Woody Patch
# ------------------------------------------------------- 
fm.4 <- DS(abund.formula = ~ herb_prp + woody_mean_p_Area  + (1|PointNum), 
           det.formula = ~ Observer,
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
           n.report = 5000,
           verbose = TRUE) 

# Save fitted model
saveRDS(fm.4, "./Data/Fitted_Models/PC_HDS_spAbund_fm4.rds")


# -------------------------------------------------------
# Abundance Fit 5: Herbaceous Proportion + Woody Clumpy
# ------------------------------------------------------- 
fm.5 <- DS(abund.formula = ~ herb_prp + woody_c_clumpy + (1|PointNum), 
           det.formula = ~ Observer,
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

# Model 3: Woody Clumpy Index is the best abundance model
# Model 4: Herbaceous Proportion + Woody Patch also shows support
summary(fm.3)
summary(fm.4)

# Check model posterior predictive check
ppc_fm.3 <- ppcAbund(fm.3, fit.stat = 'freeman-tukey', group = 1)
summary(ppc_fm.3)

ppc_fm.4 <- ppcAbund(fm.4, fit.stat = 'freeman-tukey', group = 1)
summary(ppc_fm.4)

# Model 2 Bayesian p-value = 0.4905 
# Model 4 Bayesian p-value = 0.4882

# Model shows a slightly better model fit

# -------------------------------------------------------
#                    Abundance
# -------------------------------------------------------

# Mean abundance per point
print(mean(fm.4$N.samples)) # Latent Abundance
print(mean(fm.4$mu.samples)) # Expected abundance


# Area in hectares
area <- pi*(200^2)/10000

# So the mean density per acre is 
hds_density <- mean(fm.4$N.samples) / area
print(hds_density)


# The abundance across the study area is 
study_area_abund = hds_density * 1096.698
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
lat_dens_df[,1] <- "PC HDS"
head(lat_dens_df)

# Export density dataframe
saveRDS(lat_dens_df, "./Data/Fitted_Models/PC_HDS_Dens.rds")

# In case of crash
# lat_dens_df <- readRDS("./Data/Fitted_Models/PC_HDS_spAbund_fm4.rds")

# Plot
ggplot(lat_dens_df, aes(x = Model, y = Density, fill = Model)) +
  geom_violin() + 
  geom_boxplot(aes(x = Model, y = Density), 
               width = 0.2, position = position_dodge(width = 0.8)) +
  labs(
    title = "Latent Density", 
    x = "Model", 
    y = "Density") +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 1), labels = scales::comma) + # Customize y-axis
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1), # Tilt x-axis text
    axis.text.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none" 
  )










