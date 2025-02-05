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
#                    Load Data
#
# -------------------------------------------------------

# BirdNet detections
wolfe_dat <- read.csv("C:/Users/davep/OneDrive - Texas A&M University/Rprojects/CH2_NOBO_Abundance/Data/Acoustic_Data/NOBO_Wolfe_2024.csv")

# Site covariates
site_covs <- read.csv("C:/Users/davep/OneDrive - Texas A&M University/Rprojects/CH2_NOBO_Abundance/Data/Acoustic_Data/ARU_siteCovs.csv")


# -------------------------------------------------------
#
#    Data Wrangling
#
# -------------------------------------------------------

# Subset to 14 days starting at May 26 and ending on July 17. Dates are every 4 days.
date_order <- c("2024-05-26", "2024-05-30", "2024-06-03", "2024-06-07", # in ascending order
                "2024-06-11", "2024-06-15", "2024-06-19", "2024-06-23",
                "2024-06-27", "2024-07-01", "2024-07-05", "2024-07-09",
                "2024-07-13", "2024-07-17")

# Adding occasion column
wolfe_dat <- wolfe_dat %>%
  mutate(Occasion = match(Date, date_order))

# Adding a row count
wolfe_dat$Count <- 1

# Initialize the matrix
v <- matrix(0, nrow = 27, ncol = 14)        

# Extract count data
for (i in 1:nrow(wolfe_dat)) {
  
  # Extracting plot ID
  site <- wolfe_dat$Site_Number[i]
  
  # Extracting Occasion
  occasion <- wolfe_dat$Occasion[i]
  
  # Fill in the matrix with the number of individuals
  v[site, occasion] =  v[site, occasion] + wolfe_dat$Count[i]
  
} # end loop 

# take a look
print(v)

# ARU at site 24 stopped recording on 6/20. Making Columns 8:14 NA
v[24,8:14] <- NA
 
# Adding rownames
rownames(v) <- as.numeric(1:27)

# Take a look
print(v)

# Get the site numbers with at least one call
v_mat <- as.matrix(v)
rownames(v_mat) <- NULL  
sites.a <- which(rowSums(!is.na(v_mat) & v >= 1) > 0)
print(sites.a)

# Creating a binary detection matrix, values >= 1 become 1, and 0 stays 0
v_df <- as.data.frame(v)
y <- v_df %>% mutate(across(where(is.numeric), ~ if_else(. >= 1, 1, 0)))
y <- as.matrix(y)
rownames(y) <- NULL  
colnames(y) <- NULL 
print(y)

# Total number of sites
R <- as.integer(27) 

# Number of repeat visits
J <- rep(14, R)  


