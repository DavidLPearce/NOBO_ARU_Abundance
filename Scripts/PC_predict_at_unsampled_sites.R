# ------------------------------------------------------------------------------
#
#                               Load Packages
#
# ------------------------------------------------------------------------------

# Install packages (if needed)
# install.packages("sf")
# install.packages("raster")
# install.packages("ggplot2")
# install.packages("viridis")
# install.packages("terra")
# install.packages(igraph)
# install.packages("lidR")
# install.packages("RCSF")
# install.packages("landscapemetrics")
# install.packages("progress")
# install.packages("psych")

# Load packages
library(sf)
library(sp)
library(raster)
library(ggplot2)
library(viridis)
library(terra)
library(igraph)
library(lidR)
library(RCSF)
library(landscapemetrics)
library(progress)
library(psych)

# Set seed, scientific notation, and workplace
set.seed(123)
options(scipen = 9999)
setwd(".")


# ------------------------------------------------------------------------------
#
#                                 Load Data
#
# ------------------------------------------------------------------------------

# La Copita Extent
lcr_extent <- st_read("E:/Projects/Current_Projects/NOBO Abund S TX/GIS/ShapeFiles/LaCopitaExtent/LaCopita.shp")
print(lcr_extent)

# Read in raster data
lulc_rast <- stack("E:/Projects/Current_Projects/NOBO Abund S TX/GIS/LaCopita_LULC_2022&2024/LaCopita_LULC_2024/LaCopita_LULC_2024.tif")
print(lulc_rast)

# Read in site locations
site_dat <- read.csv("./Data/Point_Count_Data/NOBO_PointCount_Locations.csv")
site_dat_sf <- st_as_sf(site_dat, coords = c("Long", "Lat"), crs = 4326)  
site_dat_sf <- st_transform(site_dat_sf, crs = st_crs(lcr_extent))   

# Directory for plots
output_dir <- "./Figures/LULC/PC_Predict/"

# Extracting each class to a new raster
woody_class <- lulc_rast == 0   # Woody
herb_class <- lulc_rast == 1    # Herbaceous
brgnd_class <- lulc_rast == 2   # Bare Ground 
# Water = 3
# Developed = 4

# Convert to SpatRaster
lulc_rast <- rast(lulc_rast)
woody_class <- rast(woody_class) 
herb_class <- rast(herb_class)
brgnd_class <- rast(brgnd_class)

# # Plots
# plot(woody_class, main = "Woody")
# plot(herb_class, main = "Herbaceous")
# plot(brgnd_class, main = "Bare Ground")
# plot(water_class, main = "Water")
# plot(dev_class, main = "Developed")

# -----------------------------------
# Create Prediction Polygons
# -----------------------------------

# Set seed for reproducibility
set.seed(123)

# Surveyed site buffers
site_buffers <- st_buffer(site_dat_sf, dist = 227)

# Remove sampled areas from the extent
unsampled_area <- st_difference(lcr_extent, st_union(site_buffers))

# Plot
plot(st_geometry(lcr_extent), border = "black", lwd = 2)
plot(st_geometry(site_buffers), add = TRUE, col = "lightgrey", border = "grey")
plot(st_geometry(site_dat_sf), add = TRUE, pch = 16, col = "black")
plot(st_geometry(unsampled_area), add = TRUE, col = "lightgreen", border = "darkgreen")


# Calculate area of unsampled region 
unsampled_area_m2 <- st_area(unsampled_area)

# Using a radius of 227m which is the typical homerange size of a quail
# and what the sampled sites covariates were extracted at
site_area <- pi * 227^2

# Number of prediction sites at rougly same area
N_pred_sites <- as.numeric(unsampled_area_m2) / site_area
print(paste("Number of Prediction Sites =", round(N_pred_sites)))


# Generate regularly-spaced points in the unsampled area
pred_points <- st_sample(unsampled_area, 
                         size = round(N_pred_sites), 
                         type = "regular")

plot(st_geometry(lcr_extent), border = "black", lwd = 2)
plot(st_geometry(site_buffers), add = TRUE, col = "lightgrey", border = "grey")
plot(st_geometry(site_dat_sf), add = TRUE, pch = 16, col = "black")
plot(st_geometry(unsampled_area), add = TRUE, col = "lightgreen", border = "darkgreen")
plot(st_geometry(pred_points), add = TRUE, pch = 16, col = "red")


# Create Voronoi polygons (Thiessen polygons) around each point
# These will divide the area into roughly equal segments
pred_voronoi <- st_voronoi(st_union(pred_points))
pred_voronoi <- st_collection_extract(pred_voronoi, "POLYGON")
pred_voronoi <- st_sf(geometry = pred_voronoi)

# Clip to unsampled area boundary
pred_segments <- st_intersection(pred_voronoi, unsampled_area)
pred_segments$site_id <- 1:nrow(pred_segments)

# Calculate actual areas to verify
pred_segments$area_m2 <- st_area(pred_segments)
print(paste("Mean segment area:", mean(pred_segments$area_m2), "m²"))
print(paste("Target area:", site_area, "m²"))

# Visualize
plot(st_geometry(lcr_extent), border = "black", lwd = 2)
plot(st_geometry(site_buffers), add = TRUE, col = "lightgrey", border = "grey")
plot(st_geometry(pred_segments), add = TRUE, col = "lightgreen", border = "darkgreen")
plot(st_geometry(site_dat_sf), add = TRUE, pch = 16, col = "black", cex = 1)
plot(st_geometry(pred_points), add = TRUE, pch = 16, col = "orange", cex = 1)

# just prediction areas
plot(st_geometry(lcr_extent), border = "black", lwd = 2)
plot(st_geometry(pred_segments), add = TRUE, col = "lightgreen", border = "darkgreen")
plot(st_geometry(pred_points), add = TRUE, pch = 16, col = "black", cex = 1)

# -----------------------------------
# Landscape Metrics Extraction
# -----------------------------------

# Remove geometry
predict_covs <- st_drop_geometry(pred_segments)
predict_covs <- predict_covs[,44:45]
predict_covs$area_m2 <- as.numeric(predict_covs$area_m2) # Remove [m^2]
head(predict_covs)

# Progress bar
pb <- progress_bar$new(
  format = "  Processing [:bar] :percent | :site | Elapsed: :elapsed| ETA: :eta",
  total = NROW(pred_segments),
  clear = FALSE,
  width = 120
)


# row = 1

# Loop
for (row in 1:NROW(pred_segments)) {
  
  ## list_lsm() list of available metrics 
  ## lsm_tbl <- list_lsm()
  ## View(lsm_tbl)
  
  # Subset the site
  site_sub <- pred_segments[row, ]
  
  # Get site ID
  SiteID <- site_sub[,'site_id']
  SiteID <- st_drop_geometry(SiteID)
  
  # Update the progress bar
  pb$tick(tokens = list(site = SiteID),)
  
  
  # Extract and crop the raster for the buffer
  lulc_clip <- terra::crop(lulc_rast, site_sub)
  lulc_clip <- terra::mask(lulc_clip, site_sub)
  plot(lulc_clip, col = c("forestgreen", "lightgreen", "saddlebrown"), main = SiteID)
  
  # Subset to class
  woody_clip <- lulc_clip == 0
  herb_clip <- lulc_clip == 1
  brgnd_clip <- lulc_clip == 2
  
  # Plot & Export
  lulc_Plotfile <- paste0(output_dir, SiteID, "_LULC.png")
  lulc_df <- as.data.frame(lulc_clip, xy = TRUE)
  lulc_plot <- ggplot(lulc_df) +
    geom_raster(aes(x = x, y = y, fill = factor(Class_name))) +
    scale_fill_manual(name = "Land Cover Type",
                      values =  c("Shrubland" = "forestgreen",
                                  "Grassland" = "lightgreen",
                                  "Bare Ground" = "saddlebrown")) +
    coord_fixed() +
    labs(title = paste("LULC Site", SiteID)) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA))
  ggsave(filename = lulc_Plotfile, plot = lulc_plot, width = 8, height = 6, dpi = 300)
  
  # ------------------------------
  # Proportion of Class - Class
  # ------------------------------
  # Extract raster cover type within buffer
  woody_prp <- terra::extract(woody_class, site_sub, fun = mean, na.rm = TRUE)
  herb_prp <- terra::extract(herb_class, site_sub, fun = mean, na.rm = TRUE)
  predict_covs[row, "woody_prp"] <- woody_prp[1, 'layer']
  predict_covs[row, "herb_prp"] <- herb_prp[1, 'layer']
  
  
  
  # ------------------------------
  # Mean Patch Area - Patch
  # ------------------------------
  # lsm_p_area: Calculates the area of each patch.(m^2)
  p_area <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_p_area")
  woody_p_area <- p_area[which(p_area$class == 0),] # Woody
  woody_p_area_mean <- mean(woody_p_area$value, na.rm = TRUE)
  herb_p_area <- p_area[which(p_area$class == 1),] # Herbaceous
  herb_p_area_mean <- mean(herb_p_area$value, na.rm = TRUE)
  predict_covs[row, "woody_mnParea"] <- woody_p_area_mean
  predict_covs[row, "herb_mnParea"] <- herb_p_area_mean
  
  
  # ------------------------------
  # Clumpy Index - Class
  # ------------------------------
  # lsm_c_clumpy: Quantifies the clumpiness (spatial aggregation)
  # Equals -1 for maximally disaggregated, 0 for randomly distributed and 1 for maximally aggregated classes.
  c_clumpy <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_c_clumpy")
  woody_c_clumpy <- c_clumpy[which(c_clumpy$class == 0),] # Woody
  herb_c_clumpy <- c_clumpy[which(c_clumpy$class == 1),] # Herbaceous
  predict_covs[row, "woody_ClmIdx"] <- woody_c_clumpy[1, 'value']
  predict_covs[row, "herb_ClmIdx"] <- herb_c_clumpy[1, 'value']
  
  # ------------------------------
  # Largest Patch Index - Class
  # ------------------------------
  # lsm_c_lpi: an 'Area and edge metric'. Percentage of landscape covered by the corresponding largest patch of each class i
  # Approaches LPI = 0 when the largest patch is becoming small and equals LPI = 100 when only one patch is present
  c_lpi <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_c_lpi")
  woody_c_lpi <- c_lpi[which(c_lpi$class == 0),] # Woody
  herb_c_lpi <- c_lpi[which(c_lpi$class == 1),] # Herbaceous
  predict_covs[row, "woody_lrgPInx"] <- woody_c_lpi[1, 'value']
  predict_covs[row, "herb_lrgPInx"] <- herb_c_lpi[1, 'value']
  
  
  # ------------------------------
  # Aggregation Index - Class
  # ------------------------------
  # lsm_c_ai: equals the number of like adjacencies divided by the theoretical maximum possible number of like adjacencies for that class
  # Equals 0 for maximally disaggregated and 100 for maximally aggregated classes.
  c_ai <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_c_ai")
  woody_c_ai <- c_ai[which(c_ai$class == 0),] # Woody
  herb_c_ai <- c_ai[which(c_ai$class == 1),] # Herbaceous
  predict_covs[row, "woody_AggInx"] <- woody_c_ai[1, 'value']
  predict_covs[row, "herb_AggInx"] <- herb_c_ai[1, 'value']
  
  
  # ------------------------------
  # Edge Density - Class
  # ------------------------------
  # lsm_c_ed: equals the sum of all edges of class i in relation to the landscape area
  # = 0 if only one patch is present (and the landscape boundary is not included) and increases, without limit, as the landscapes becomes more patchy
  c_ed <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_c_ed")
  woody_c_ed <- c_ed[which(c_ed$class == 0),] # Woody
  herb_c_ed <- c_ed[which(c_ed$class == 1),] # Herbaceous
  predict_covs[row, "woody_EdgDens"] <- woody_c_ed[1, 'value']
  predict_covs[row, "herb_EdgDens"] <- herb_c_ed[1, 'value']
  
  
  # ------------------------------
  # Patch Density- Class
  # ------------------------------
  # lsm_c_pd: Number per 100 hectares
  # Increases as the landscape gets more patchy. Reaches its maximum if every cell is a different patch.
  c_pd <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_c_pd")
  woody_c_pd <- c_pd[which(c_pd$class == 0),] # Woody
  herb_c_pd <- c_pd[which(c_pd$class == 1),] # Herbaceous
  predict_covs[row, "woody_Pdens"] <- woody_c_pd[1, 'value']
  predict_covs[row, "herb_Pdens"] <- herb_c_pd[1, 'value']
  
  # ------------------------------
  # Interspersion and Juxtaposition Index (IJI) - Class
  # ------------------------------
  c_iji <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_c_iji")
  woody_c_iji <- c_iji[which(c_iji$class == 0),] # Woody
  herb_c_iji <- c_iji[which(c_iji$class == 1),] # Herbaceous
  predict_covs[row, "woody_IJI"] <- woody_c_iji[1, 'value']
  predict_covs[row, "herb_IJI"] <- herb_c_iji[1, 'value']
  
  # ------------------------------
  # Cohesion Index - Class
  # ------------------------------
  c_coh <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_c_cohesion")
  woody_c_coh <- c_coh[which(c_coh$class == 0),] # Woody
  herb_c_coh <- c_coh[which(c_coh$class == 1),] # Herbaceous
  predict_covs[row, "woody_COH"] <- woody_c_coh[1, 'value']
  predict_covs[row, "herb_COH"] <- herb_c_coh[1, 'value']
  
  # ------------------------------
  # Division - Class
  # ------------------------------
  c_div <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_c_division")
  woody_c_div <- c_div[which(c_div$class == 0),] # Woody
  herb_c_div <- c_div[which(c_div$class == 1),] # Herbaceous
  predict_covs[row, "woody_DIV"] <- woody_c_div[1, 'value']
  predict_covs[row, "herb_DIV"] <- herb_c_div[1, 'value']
  
  # ------------------------------
  # Split - Class
  # ------------------------------
  c_split <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_c_split")
  woody_c_split <- c_split[which(c_split$class == 0),] # Woody
  herb_c_split <- c_split[which(c_split$class == 1),] # Herbaceous
  predict_covs[row, "woody_SPLIT"] <- woody_c_split[1, 'value']
  predict_covs[row, "herb_SPLIT"] <- herb_c_split[1, 'value']
  
} # -------------------- End Extraction Loop -----------------------------


# Take a look
str(predict_covs)
print(predict_covs)

# Export data
write.csv(predict_covs, "./Data/Point_Count_Data/PC_PredictsiteCovs.csv")
