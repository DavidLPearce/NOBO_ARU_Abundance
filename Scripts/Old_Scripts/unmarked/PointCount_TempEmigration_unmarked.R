# -------------------------------------------------------
#
#                    Load libraries
#
# -------------------------------------------------------

# Install packages (if needed)
# install.packages("tidyverse")
# install.packages("unmarked")
# install.packages("AICcmodavg")

# Load library
library(tidyverse)
library(unmarked)
library(AICcmodavg)
library(ggplot2)

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

# Remove NAs
pc_dat_NAom <- na.omit(pc_dat)

# creating a matrix that is 4 Surveys * 3 Distance bins wide
# and 10 rows long
det_mat <- matrix(0, nrow = 10, ncol = 12)

# adding a column to state a NOBO was detected, a count column
pc_dat_NAom$count <- 1

# Loop to fill matrix with data
for (i in 1:nrow(pc_dat_NAom)) {
  point_num <- pc_dat_NAom$PointNum[i]
  occasion <- pc_dat_NAom$Survey[i]
  distance_cat <- as.numeric(pc_dat_NAom$DistBin[i])
  
  # Determine the column in the matrix
  col_index <- (occasion - 1) * 3 + distance_cat
  
  # Fill in the matrix with the number of individuals
  det_mat[point_num, col_index] <- det_mat[point_num, col_index] + pc_dat_NAom$count[i]
  
}#end loop

# Take a look
print(det_mat)


## Observation covariates
# Create matrix for each covariate
obsvr_mat <- matrix(NA, nrow = 10, ncol = 4)
temp_mat <- matrix(NA, nrow = 10, ncol = 4)
wind_mat <- matrix(NA, nrow = 10, ncol = 4)
sky_mat <- matrix(NA, nrow = 10, ncol = 4)
doy_mat <- matrix(NA, nrow = 10, ncol = 4)


# Fill the matrices
for (i in 1:nrow(pc_dat)) {
  point_num <- pc_dat$PointNum[i]
  occasion <- pc_dat$Survey[i]
  
  obsvr_mat[point_num, occasion] <-  pc_dat$Observer[i]
  temp_mat[point_num, occasion] <-  pc_dat$Temp.deg.F[i]
  wind_mat[point_num, occasion] <-  pc_dat$Wind.Beau.Code[i]
  sky_mat[point_num, occasion] <-  pc_dat$Sky.Beau.Code[i]
  doy_mat[point_num, occasion] <-  pc_dat$DOY[i]

}# end loop

# take a look
print(obsvr_mat)
print(temp_mat)
print(wind_mat)
print(sky_mat)
print(doy_mat)


# -------------------------------------------------------
#
#             Temporary Emigration HDS Model
#
# -------------------------------------------------------

## format data for unmarked


# creating breaks for unmarked frame
breaks <- c(0, 50, 100, 200)

# unmarked fram
pc_umf <- unmarkedFrameGDS(y = det_mat, 
                           survey = "point", 
                           unitsIn = "m", 
                           dist.breaks = breaks, 
                           numPrimary = 4,  # number of occasions
                           siteCovs = data.frame(herb_prp = site_covs$herb_prp,
                                                 woody_mean_p_Area = site_covs$woody_mean_p_Area),
                          yearlySiteCovs=list(
                                              observer = obsvr_mat,
                                              temp = temp_mat,
                                              wind = wind_mat,
                                              sky = sky_mat,
                                              doy = doy_mat ))# end function

# Take a look at the structure
str(pc_umf)



# -------------------------------------------------------
#
#            Determining Detection Function
#
# -------------------------------------------------------

# exponential detection function 
fm0.exp <- gdistsamp(lambdaformula = ~1, 
                     phiformula = ~1, 
                     pformula = ~1, 
                     keyfun = "exp", 
                     output = "abund",
                     mixture = "P", 
                     K = 100, 
                     se = TRUE, 
                     data = pc_umf)
# Model Sumamry
summary(fm0.exp)




# hazard detection function 
fm0.haz <- gdistsamp(lambdaformula = ~1, 
                     phiformula = ~1, 
                     pformula = ~1, 
                     keyfun = "haz", 
                     output = "abund",
                     mixture = "P", 
                     K = 100, 
                     se = TRUE, 
                     data = pc_umf)

# Model Sumamry
summary(fm0.haz)


# half-normal detection function 
fm0.hn <- gdistsamp(lambdaformula = ~1, 
                    phiformula = ~1, 
                    pformula = ~1, 
                    keyfun = "halfnorm", 
                    output = "abund",
                    mixture = "P", 
                    K = 100, 
                    se = TRUE, 
                    data = pc_umf,
                    control=list(trace=TRUE, REPORT=1))

# Model Sumamry
summary(fm0.haz)


## Model ranking to determine best detection function

# List of models
det_func_models <- list(
  exp = fm0.exp,
  haz = fm0.haz,
  hn = fm0.hn
)

# Saving model names
det_func_model_names <- c("Exponential", 
                          "Hazard", 
                          "Half-normal")


# Comparing models
det_func_aic <- aictab(cand.set = det_func_models, 
                       modnames = det_func_model_names)

# Print comparison
print(det_func_aic)


# Exponential function is the better model


# -------------------------------------------------------
#
#                 Detection Models
#
# -------------------------------------------------------

# Fit 0: Null Model
detfm.0 <- gdistsamp(lambdaformula = ~1, 
                     phiformula = ~1, 
                     pformula = ~1, 
                     keyfun = "exp", 
                     output = "abund", 
                     mixture = "P", 
                     K = 100, 
                     se = TRUE, 
                     data = pc_umf)
# Fit 1: observer
detfm.1 <- gdistsamp(lambdaformula = ~1, 
                     phiformula = ~1, 
                     pformula = ~observer, 
                     keyfun = "exp", 
                     output = "abund", 
                     mixture = "P", 
                     K = 100, 
                     se = TRUE, 
                     data = pc_umf)

# Fit 2: temp
detfm.2 <- gdistsamp(lambdaformula = ~1, 
                     phiformula = ~1, 
                     pformula = ~temp, 
                     keyfun = "exp", 
                     output = "abund", 
                     mixture = "P", 
                     K = 100, 
                     se = TRUE, 
                     data = pc_umf)
# Fit 3: wind
detfm.3 <- gdistsamp(lambdaformula = ~1, 
                     phiformula = ~1, 
                     pformula = ~wind, 
                     keyfun = "exp", 
                     output = "abund", 
                     mixture = "P", 
                     K = 100, 
                     se = TRUE, 
                     data = pc_umf)

# Fit 4: sky
detfm.4 <- gdistsamp(lambdaformula = ~1, 
                     phiformula = ~1, 
                     pformula = ~sky, 
                     keyfun = "exp", 
                     output = "abund", 
                     mixture = "P", 
                     K = 100, 
                     se = TRUE, 
                     data = pc_umf)
# Fit 5: doy
detfm.5 <- gdistsamp(lambdaformula = ~1, 
                     phiformula = ~1, 
                     pformula = ~doy, 
                     keyfun = "exp", 
                     output = "abund", 
                     mixture = "P", 
                     K = 100, 
                     se = TRUE, 
                     data = pc_umf)

## Model ranking to determine which covariates influence detection
 
# List of models
det_cov_models <- list(
                      null = detfm.0,
                      observer = detfm.1,
                      temp = detfm.2,
                      wind = detfm.3,
                      sky = detfm.4,
                      doy = detfm.5)

# Saving model names
det_cov_model_names <- c("null",
                         "observer", 
                         "temp", 
                         "wind",
                         "sky",
                         "doy")


# Comparing models
det_cov_aic_comp <- aictab(cand.set = det_cov_models, 
                           modnames = det_cov_model_names)

# Print comparison
print(det_cov_aic_comp)

# Sky is the best detection model

# -------------------------------------------------------
#
#                 Detection Models
#
# -------------------------------------------------------

# Fit 0: Null
availfm.0 <- gdistsamp(lambdaformula = ~1, 
                       phiformula = ~1, 
                       pformula = ~sky, 
                       keyfun = "exp", 
                       output = "abund", 
                       mixture = "P", 
                       K = 100, 
                       se = TRUE, 
                       data = pc_umf)


# Fit 1: temp
availfm.1 <- gdistsamp(lambdaformula = ~1, 
                       phiformula = ~temp, 
                       pformula = ~sky, 
                       keyfun = "exp", 
                       output = "abund", 
                       mixture = "P", 
                       K = 100, 
                       se = TRUE, 
                       data = pc_umf)


# Fit 2: wind
availfm.2 <- gdistsamp(lambdaformula = ~1, 
                       phiformula = ~wind, 
                       pformula = ~sky, 
                       keyfun = "exp", 
                       output = "abund", 
                       mixture = "P", 
                       K = 100, 
                       se = TRUE, 
                       data = pc_umf)


# Fit 3: sky
availfm.3 <- gdistsamp(lambdaformula = ~1, 
                       phiformula = ~sky, 
                       pformula = ~sky, 
                       keyfun = "exp", 
                       output = "abund", 
                       mixture = "P", 
                       K = 100, 
                       se = TRUE, 
                       data = pc_umf)


# Fit 4: doy
availfm.4 <- gdistsamp(lambdaformula = ~1, 
                       phiformula = ~doy, 
                       pformula = ~sky, 
                       keyfun = "exp", 
                       output = "abund", 
                       mixture = "P", 
                       K = 100, 
                       se = TRUE, 
                       data = pc_umf)

# List of models
avail_cov_models <- list(
                       null = detfm.0,
                       temp = detfm.1,
                       wind = detfm.2,
                       sky = detfm.3,
                       doy = detfm.4)

# Saving model names
avail_cov_model_names <- c("null",
                         "temp", 
                         "wind",
                         "sky",
                         "doy")


# Comparing models
avail_cov_aic_comp <- aictab(cand.set = avail_cov_models, 
                           modnames = avail_cov_model_names)

# Print comparison
print(avail_cov_aic_comp)

# null model is the best model

# -------------------------------------------------------
#
#                 Abundance Models
#
# -------------------------------------------------------

# Fit 0: Null
abundfm.0 <- gdistsamp(lambdaformula = ~1, 
                       phiformula = ~1, 
                       pformula = ~sky, 
                       keyfun = "exp", 
                       output = "abund", 
                       mixture = "P", 
                       K = 100, 
                       se = TRUE, 
                       data = pc_umf)

# Fit 1: herb_prp
abundfm.1 <- gdistsamp(lambdaformula = ~herb_prp, 
                       phiformula = ~1, 
                       pformula = ~sky, 
                       keyfun = "exp", 
                       output = "abund", 
                       mixture = "P", 
                       K = 100, 
                       se = TRUE, 
                       data = pc_umf)

# Fit 2: herb_prp
abundfm.2 <- gdistsamp(lambdaformula = ~woody_mean_p_Area, 
                       phiformula = ~1, 
                       pformula = ~sky, 
                       keyfun = "exp", 
                       output = "abund", 
                       mixture = "P", 
                       K = 100, 
                       se = TRUE, 
                       data = pc_umf)

# Fit 3: herb_prp + woody_mean_p_Area
abundfm.3 <- gdistsamp(lambdaformula = ~herb_prp + woody_mean_p_Area, 
                       phiformula = ~1, 
                       pformula = ~sky, 
                       keyfun = "exp", 
                       output = "abund", 
                       mixture = "P", 
                       K = 100, 
                       se = TRUE, 
                       data = pc_umf)


# List of models
abund_cov_models <- list(
                         null = abundfm.0,
                         herb_prp = abundfm.1,
                         woody_mean_p_Area = abundfm.2,
                         herb_prp_woody_mean_p_Area = abundfm.3)

# Saving model names
abund_cov_model_names <- c("null",
                           "herb_prp", 
                           "woody_mean_p_Area",
                           "herb_prp_woody_mean_p_Area")


# Comparing models
abund_cov_aic_comp <- aictab(cand.set = abund_cov_models, 
                             modnames = abund_cov_model_names)

# Print comparison
print(abund_cov_aic_comp)


# herb proportion is the best model

 
str(abundfm.1)
summary(abundfm.1)

# check fit 
yhat <- fitted(abundfm.1)
sum((pc_umf@y - yhat)^2 / yhat, na.rm=TRUE)

Nmix.gof.test(abundfm.1, nsim = 1000)


# -------------------------------------------------------
#
#              Estimatates of Abundance
#
# -------------------------------------------------------


# Compute total abundance and confidence intervals
abundance_pred <- predict(abundfm.1, type = "lambda")
total_abundance <- sum(abundance_pred$Predicted)
total_abundance_se <- sqrt(sum(abundance_pred$SE^2))  # Combined SE using propagation

# Confidence interval
ci_lower <- total_abundance - 1.96 * total_abundance_se
ci_upper <- total_abundance + 1.96 * total_abundance_se

# Create a data frame for plotting
total_abundance_df <- data.frame(
  TotalAbundance = total_abundance,
  CI_Lower = ci_lower,
  CI_Upper = ci_upper
)

# Plot
ggplot(total_abundance_df, aes(x = "Total Abundance", y = TotalAbundance)) +
  geom_point(size = 5, color = "blue") +  # Use size to adjust point size
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, color = "black") +
  labs(title = "Temp E Nmix",
       x = "",
       y = "Abundance") +
  theme_minimal() +
  ylim(0, 300) +  # Set y-axis limits
  theme(axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank())




#### Abundance by site



# Load required library
library(ggplot2)

# Predict abundance and confidence intervals
abundance_pred <- predict(abundfm.1, type = "lambda")

# Add site identifiers if not present
abundance_pred$Site <- 1:nrow(abundance_pred)

# Plot predicted abundance with confidence intervals
ggplot(abundance_pred, aes(x = Site, y = Predicted)) +
  geom_point(size = 3, color = "blue") +  # Predicted abundance as points
  geom_errorbar(aes(ymin = Predicted - 1.96 * SE, ymax = Predicted + 1.96 * SE),
                width = 0.2, color = "black") +  # CI bars
  labs(title = "Site-Specific Abundance Predictions",
       x = "Site",
       y = "Predicted Abundance") +
  theme_minimal()
