# ------------------------------------------------------------------------------
#
#                               Load Packages
#
# ------------------------------------------------------------------------------

install.packages("sf")
install.packages("raster")
install.packages("ggplot2")
install.packages("terra")
install.packages("landscapemetrics")
install.packages("progress")

library(sf)
library(raster)
library(ggplot2)
library(terra)
library(landscapemetrics)
library(progress)

set.seed(123)
options(scipen = 9999)
setwd(".")

# ------------------------------------------------------------------------------
#
#                                 Load Data
#
# ------------------------------------------------------------------------------

# Read in raster data
lulc_rast <- stack("D:/LaCopita_GIS_Data/LaCopitaLULC/LaCopitaLULC_60cm_SVM/LULC_60cm_SVM_Raster/LaCopitaLULC_60cm_SVM.tif")

# Read in site locations
site_dat <- read.csv("./Data/Point_Count_Data/NOBO_PointCount_Locations.csv")


# Extracting each class to a new raster
woody_class <- lulc_rast == 0   # Woody
herb_class <- lulc_rast == 1    # Herbaceous
baregnd_class <- lulc_rast == 2 # Bare Ground
water_class <- lulc_rast == 3   # Water
dev_class <- lulc_rast == 4     # Developed

# Convert to SpatRaster
lulc_rast <- rast(lulc_rast)
woody_class <- rast(woody_class) 
herb_class <- rast(herb_class)
baregnd_class <- rast(baregnd_class)
water_class <- rast(water_class)
dev_class <- rast(dev_class)

# Plots
plot(woody_class, main = "Woody")
plot(herb_class, main = "Herbaceous")
plot(baregnd_class, main = "Bare Ground")
plot(water_class, main = "Water")
plot(dev_class, main = "Developed")


# ------------------------------------------------------------------------------
#
#                                 Extracting
#
# ------------------------------------------------------------------------------


# Initialize the progress bar
pb <- progress_bar$new(
  format = "  Processing [:bar] :percent in :elapsed, remaining: :eta",
  total = NROW(site_dat),
  clear = FALSE,
  width = 100
)


# ------------------------
#    Loop to extract
# ------------------------
for (index in 1:NROW(site_dat[,'Number'])) {
  

  
  # Subset the site
  site_sub <- site_dat[which(site_dat[,'Number'] == index),]
  
  # Setting projection
  site_sub_coords <- SpatialPoints(coords = site_sub[, c("Long", "Lat")],
                                   proj4string = CRS("+proj=longlat +datum=WGS84"))

  # Create a buffer around the site
  site_buffer <- terra::buffer(site_sub_coords, width = 200)
  site_buffer_terra <- vect(site_buffer) # SpatVector
  site_buffer_terra <- terra::project(site_buffer_terra, crs(lulc_rast))# CRS of lulc
  
  
  ## Proportion of Class ##
  
  # Extract raster cover type within buffer
  woody_prp <- terra::extract(woody_class, site_buffer_terra, fun = mean, na.rm = TRUE)
  herb_prp <- terra::extract(herb_class, site_buffer_terra, fun = mean, na.rm = TRUE)
  baregnd_prp <- terra::extract(baregnd_class, site_buffer_terra, fun = mean, na.rm = TRUE) 
  
  
  ## Landscape metrics ##
  
  # list_lsm() list of available metrics.
  
  # Extract and crop the raster for the buffer
  lulc_clip <- terra::crop(lulc_rast, site_buffer_terra)
  lulc_clip <- terra::mask(lulc_clip, site_buffer_terra)
  #plot(lulc_clip)
  
  # Woody
  woody_clip <- terra::crop(woody_class, site_buffer_terra)
  woody_clip <- terra::mask(woody_clip, site_buffer_terra)
  
  # Herbaceous
  herb_clip <- terra::crop(herb_class, site_buffer_terra)
  herb_clip <- terra::mask(herb_clip, site_buffer_terra)
  
  # Bareground
  baregnd_clip <- terra::crop(baregnd_class, site_buffer_terra)
  baregnd_clip <- terra::mask(baregnd_clip, site_buffer_terra)
  
  
  ## Patch-level metrics ##
  
  # Mean Patch Area
  # lsm_p_area: Calculates the area of each patch.(m^2)
  p_area <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_p_area")
  
  # Woody
  woody_p_area <- p_area[which(p_area$class == 0),]
  woody_p_area_mean <- mean(woody_p_area$value, na.rm = TRUE) # mean patch area
  
  # Clumpy Index
  # lsm_c_clumpy: Quantifies the clumpiness (spatial aggregation), 1 indicates highly aggregated, -1 disaggregation
  c_clumpy <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_c_clumpy")
  
  # Woody
  woody_c_clumpy <- c_clumpy[which(c_clumpy$class == 0),]

  
  
  ## Add extracted values ##
  
  site_dat[site_dat$Number == index, "woody_prp"] <- woody_prp[1, 'layer']           # Proportion woody
  site_dat[site_dat$Number == index, "herb_prp"] <- herb_prp[1, 'layer']             # Proportion herbaceous
  site_dat[site_dat$Number == index, "baregnd_prp"] <- baregnd_prp[1, 'layer']       # Proportion bare ground
  site_dat[site_dat$Number == index, "woody_mean_p_Area"] <- woody_p_area_mean       # Woody mean patch area
  site_dat[site_dat$Number == index, "woody_c_clumpy"] <- woody_c_clumpy[1, 'value'] # Woody clumpy index
  
  # Update the progress bar
  pb$tick()

} # end extraction loop  
  
# Take a look
print(site_dat)


# Export data
write.csv(site_dat, "./Data/Point_Count_Data/PointCount_siteCovs.csv")# .csv


  
  
  
  