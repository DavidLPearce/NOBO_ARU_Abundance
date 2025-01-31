# ------------------------------------------------------------------------------
#
#                               Load Packages
#
# ------------------------------------------------------------------------------

# Install packages (if needed)
# install.packages("sf")
# install.packages("raster")
# install.packages("ggplot2")
# install.packages("terra")
# install.packages("landscapemetrics")
# install.packages("progress")

# Load packages
library(sf)
library(raster)
library(ggplot2)
library(terra)
library(landscapemetrics)
library(progress)

# Set seed, scientific notation, and workplace
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
site_dat <- read.csv("./Data/Acoustic_Data/ARU_sites.csv")

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
for (row in 1:NROW(site_dat)) {
  
    # Subset the site
  site_sub <- site_dat[row, ]
  
  # Setting projection
  site_sub_coords <- SpatialPoints(coords = site_sub[, c("Long", "Lat")],
                                   proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  # Create a buffer around the site
  site_buffer <- terra::buffer(site_sub_coords, width = 200)
  site_buffer_terra <- vect(site_buffer) # SpatVector
  site_buffer_terra <- terra::project(site_buffer_terra, crs(lulc_rast))# CRS of lulc
  
  

  # Extract and crop the raster for the buffer
  lulc_clip <- terra::crop(lulc_rast, site_buffer_terra)
  lulc_clip <- terra::mask(lulc_clip, site_buffer_terra)
  # plot(lulc_clip)
  woody_clip <- lulc_clip == 0
  herb_clip <- lulc_clip == 1 
  baregnd_clip <- lulc_clip == 2
 
  ## list_lsm() list of available metrics ##
  
  ## Proportion of Class ##
  # Extract raster cover type within buffer
  woody_prp <- terra::extract(woody_class, site_buffer_terra, fun = mean, na.rm = TRUE)
  herb_prp <- terra::extract(herb_class, site_buffer_terra, fun = mean, na.rm = TRUE)
  baregnd_prp <- terra::extract(baregnd_class, site_buffer_terra, fun = mean, na.rm = TRUE) 
  
  ## Mean Patch Area ##
  # lsm_p_area: Calculates the area of each patch.(m^2)
  p_area <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_p_area")
  woody_p_area <- p_area[which(p_area$class == 0),] # Woody
  woody_p_area_mean <- mean(woody_p_area$value, na.rm = TRUE)
  herb_p_area <- p_area[which(p_area$class == 1),] # Herbaceous
  herb_p_area_mean <- mean(herb_p_area$value, na.rm = TRUE) 
  baregnd_p_area <- p_area[which(p_area$class == 2),] # Bareground
  baregnd_p_area_mean <- mean(baregnd_p_area$value, na.rm = TRUE)
  
  ## Clumpy Index ##
  # lsm_c_clumpy: Quantifies the clumpiness (spatial aggregation)
  # Equals -1 for maximally disaggregated, 0 for randomly distributed and 1 for maximally aggregated classes.
  c_clumpy <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_c_clumpy")
  woody_c_clumpy <- c_clumpy[which(c_clumpy$class == 0),] # Woody
  herb_c_clumpy <- c_clumpy[which(c_clumpy$class == 1),] # Herbaceous
  baregnd_c_clumpy <- c_clumpy[which(c_clumpy$class == 2),] # Bareground
  
  ## Patch Density ##
  # lsm_c_pd: Number per 100 hectares
  # Increases as the landscape gets more patchy. Reaches its maximum if every cell is a different patch.
  c_pd <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_c_pd")
  woody_c_pd <- c_pd[which(c_pd$class == 0),] # Woody
  herb_c_pd <- c_pd[which(c_pd$class == 1),] # Herbaceous
  baregnd_c_pd <- c_pd[which(c_pd$class == 2),] # Bareground
  
  ## Mean Shape Index ##
  # lsm_l_shape_mn: summarised as the mean of all patches in the landscape
  # Equals SHAPE_MN = 1 if all patches are squares. Increases, without limit, as the shapes of patches become more complex.
  woody_shape_mn <- landscapemetrics::calculate_lsm(woody_clip, what = "lsm_l_shape_mn") # Woody
  herb_shape_mn <- landscapemetrics::calculate_lsm(herb_clip, what = "lsm_l_shape_mn") # Herbaceous
  baregnd_shape_mn <- landscapemetrics::calculate_lsm(baregnd_clip, what = "lsm_l_shape_mn") # Bareground
  
  ## Largest Patch Index ##
  # lsm_c_lpi: an 'Area and edge metric'. Percentage of landscape covered by the corresponding largest patch of each class i
  # Approaches LPI = 0 when the largest patch is becoming small and equals LPI = 100 when only one patch is present
  c_lpi <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_c_lpi")
  woody_c_lpi <- c_lpi[which(c_lpi$class == 0),] # Woody
  herb_c_lpi <- c_lpi[which(c_lpi$class == 1),] # Herbaceous
  baregnd_c_lpi <- c_lpi[which(c_lpi$class == 2),] # Bareground
  
  
  ## Aggregation Index ##
  # lsm_c_ai: equals the number of like adjacencies divided by the theoretical maximum possible number of like adjacencies for that class
  # Equals 0 for maximally disaggregated and 100 for maximally aggregated classes.
  c_ai <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_c_ai")
  woody_c_ai <- c_ai[which(c_ai$class == 0),] # Woody
  herb_c_ai <- c_ai[which(c_ai$class == 1),] # Herbaceous
  baregnd_c_ai <- c_ai[which(c_ai$class == 2),] # Bareground
  
  ## Edge Density ##
  # lsm_c_ed: equals the sum of all edges of class i in relation to the landscape area
  # = 0 if only one patch is present (and the landscape boundary is not included) and increases, without limit, as the landscapes becomes more patchy
  c_ed <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_c_ed")
  woody_c_ed <- c_ed[which(c_ed$class == 0),] # Woody
  herb_c_ed <- c_ed[which(c_ed$class == 1),] # Herbaceous
  baregnd_c_ed <- c_ed[which(c_ed$class == 2),] # Bareground
  
  
  ## Add extracted values ##
  # Proportion
  site_dat[row, "woody_prp"] <- woody_prp[1, 'layer']          
  site_dat[row, "herb_prp"] <- herb_prp[1, 'layer']              
  site_dat[row, "baregnd_prp"] <- baregnd_prp[1, 'layer']
  # Mean Patch Area
  site_dat[row, "woody_Parea"] <- woody_p_area_mean           
  site_dat[row, "herb_Parea"] <- herb_p_area_mean
  site_dat[row, "baregnd_Parea"] <- baregnd_p_area_mean
  # Woody Clumpy Index
  site_dat[row, "woody_ClmIdx"] <- woody_c_clumpy[1, 'value'] 
  site_dat[row, "herb_ClmIdx"] <- herb_c_clumpy[1, 'value']
  site_dat[row, "baregnd_ClmIdx"] <- baregnd_c_clumpy[1, 'value']  
  # Patch Density
  site_dat[row, "woody_Pdens"] <- woody_c_pd[1, 'value']      
  site_dat[row, "herb_Pdens"] <- herb_c_pd[1, 'value']  
  site_dat[row, "baregnd_Pdens"] <-baregnd_c_pd[1, 'value'] 
  # Mean Shape Index
  site_dat[row, "woody_ShpInx"] <- woody_shape_mn[1, 'value']   
  site_dat[row, "herb_ShpInx"] <- herb_shape_mn[1, 'value']  
  site_dat[row, "baregnd_ShpInx"] <- baregnd_shape_mn[1, 'value']
  # Largest Patch Index
  site_dat[row, "woody_lrgPInx"] <- woody_c_lpi[1, 'value']   
  site_dat[row, "herb_lrgPInx"] <- herb_c_lpi[1, 'value']  
  site_dat[row, "herb_lrgPInx"] <- baregnd_c_lpi[1, 'value'] 
  # Aggregation Index
  site_dat[row, "woody_AggInx"] <- woody_c_ai[1, 'value']     
  site_dat[row, "herb_AggInx"] <- herb_c_ai[1, 'value']  
  site_dat[row, "baregnd_AggInx"] <- baregnd_c_ai[1, 'value']  
  # Edge Density
  site_dat[row, "woody_EdgDens"] <- woody_c_ed[1, 'value']   
  site_dat[row, "herb_EdgDens"] <- herb_c_ed[1, 'value']  
  site_dat[row, "baregnd_EdgDens"] <- baregnd_c_ed[1, 'value']
  
  # Update the progress bar
  pb$tick()
  
} # end extraction loop  

# Take a look
str(site_dat)
print(site_dat)


# Export data
write.csv(site_dat, "./Data/Acoustic_Data/ARU_siteCovs.csv")# .csv



