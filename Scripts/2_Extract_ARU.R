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
# install.packages("lidR")
# install.packages("RCSF")
# install.packages("landscapemetrics")
# install.packages("progress")

# Load packages
library(sf)
library(sp)
library(raster)
library(ggplot2)
library(viridis)
library(terra)
library(lidR)
library(RCSF)
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
print(lulc_rast)

# Read in LiDAR data
las_folder <- "D:/LaCopita_GIS_Data/LIDAR2018_70cm/LaCopita_LiDAR_tiles/LAS"# Set directory to las files
lasFiles <- readLAScatalog(las_folder, filter = "-keep_class 2 3 4 5")# Read in las files
summary(lasFiles)

# Read in site locations
site_dat <- read.csv("./Data/Acoustic_Data/ARU_sites.csv")
print(site_dat)

# Directory for plots
output_dir <- "./Figures/LULC/ARU/"

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
# plot(baregnd_class, main = "Bare Ground")
# plot(water_class, main = "Water")
# plot(dev_class, main = "Developed")


# ------------------------------------------------------------------------------
#
#                                 Extracting
#
# ------------------------------------------------------------------------------

# ------------------------
#    Extract Covs
# ------------------------

# Initialize the progress bar
pb <- progress_bar$new(
  format = "  Processing [:bar] :percent in :elapsed, remaining: :eta",
  total = NROW(site_dat),
  clear = FALSE,
  width = 100
)

# row = 1

# Loop
for (row in 1:NROW(site_dat)) {
  
  ## list_lsm() list of available metrics ##
  
  # Subset the site
  site_sub <- site_dat[row, ]
  
  SiteID <- site_sub[,1]
  
  # Setting projection
  site_sub_coords <- SpatialPoints(coords = site_sub[, c("Long", "Lat")],
                                   proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  # ------------------------
  #
  #         LULC
  #
  # ------------------------
  
  # Create a buffer around the site
  site_buffer <- terra::buffer(site_sub_coords, width = 200)
  site_buffer_terra <- vect(site_buffer) # SpatVector
  site_buffer_terra <- terra::project(site_buffer_terra, terra::crs(lulc_rast))# CRS of lulc
  
  # Extract and crop the raster for the buffer
  lulc_clip <- terra::crop(lulc_rast, site_buffer_terra)
  lulc_clip <- terra::mask(lulc_clip, site_buffer_terra)
  
  # Subset to class
  woody_clip <- lulc_clip == 0
  herb_clip <- lulc_clip == 1
  brgnd_clip <- lulc_clip == 2
  
  # Plot & Export
  lulc_Plotfile <- paste0(output_dir, SiteID, "_LULC.png")
  lulc_df <- as.data.frame(lulc_clip, xy = TRUE)
  lulc_plot <- ggplot(lulc_df) +
    geom_raster(aes(x = x, y = y, fill = factor(Classvalue))) +
    scale_fill_manual(name = "Land Cover Type",
                      values =  c("1" = "forestgreen",
                                  "2" = "lightgreen",
                                  "3" = "saddlebrown"),
                      labels = c("1" = "Woody",
                                 "2" = "Herbaceous",
                                 "3" = "Bareground")) +
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
  # Proportion of Class
  # ------------------------------
  # Extract raster cover type within buffer
  woody_prp <- terra::extract(woody_class, site_buffer_terra, fun = mean, na.rm = TRUE)
  herb_prp <- terra::extract(herb_class, site_buffer_terra, fun = mean, na.rm = TRUE)
  brgnd_prp <- terra::extract(brgnd_class, site_buffer_terra, fun = mean, na.rm = TRUE)
  site_dat[row, "woody_prp"] <- woody_prp[1, 'layer']          
  site_dat[row, "herb_prp"] <- herb_prp[1, 'layer']              
  site_dat[row, "open_prp"] <- herb_prp[1, 'layer'] + brgnd_prp[1, 'layer'] 
  
  
  # ------------------------------
  # Mean Patch Area
  # ------------------------------
  # lsm_p_area: Calculates the area of each patch.(m^2)
  p_area <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_p_area")
  woody_p_area <- p_area[which(p_area$class == 0),] # Woody
  woody_p_area_mean <- mean(woody_p_area$value, na.rm = TRUE)
  herb_p_area <- p_area[which(p_area$class == 1),] # Herbaceous
  herb_p_area_mean <- mean(herb_p_area$value, na.rm = TRUE) 
  site_dat[row, "woody_mnParea"] <- woody_p_area_mean    
  site_dat[row, "herb_mnParea"] <- herb_p_area_mean
  
  
  # ------------------------------
  # Clumpy Index
  # ------------------------------
  # lsm_c_clumpy: Quantifies the clumpiness (spatial aggregation)
  # Equals -1 for maximally disaggregated, 0 for randomly distributed and 1 for maximally aggregated classes.
  c_clumpy <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_c_clumpy")
  woody_c_clumpy <- c_clumpy[which(c_clumpy$class == 0),] # Woody
  herb_c_clumpy <- c_clumpy[which(c_clumpy$class == 1),] # Herbaceous
  site_dat[row, "woody_ClmIdx"] <- woody_c_clumpy[1, 'value'] 
  site_dat[row, "herb_ClmIdx"] <- herb_c_clumpy[1, 'value']
  
  
  # ------------------------------
  # Mean Shape Index
  # ------------------------------
  # lsm_l_shape_mn: summarised as the mean of all patches in the landscape
  # Equals SHAPE_MN = 1 if all patches are squares. Increases, without limit, as the shapes of patches become more complex.
  woody_shape_mn <- landscapemetrics::calculate_lsm(woody_clip, what = "lsm_l_shape_mn") # Woody
  herb_shape_mn <- landscapemetrics::calculate_lsm(herb_clip, what = "lsm_l_shape_mn") # Herbaceous
  site_dat[row, "woody_ShpInx"] <- woody_shape_mn[1, 'value']   
  site_dat[row, "herb_ShpInx"] <- herb_shape_mn[1, 'value']  
  
  
  # ------------------------------
  # Largest Patch Index
  # ------------------------------
  # lsm_c_lpi: an 'Area and edge metric'. Percentage of landscape covered by the corresponding largest patch of each class i
  # Approaches LPI = 0 when the largest patch is becoming small and equals LPI = 100 when only one patch is present
  c_lpi <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_c_lpi")
  woody_c_lpi <- c_lpi[which(c_lpi$class == 0),] # Woody
  herb_c_lpi <- c_lpi[which(c_lpi$class == 1),] # Herbaceous
  site_dat[row, "woody_lrgPInx"] <- woody_c_lpi[1, 'value']   
  site_dat[row, "herb_lrgPInx"] <- herb_c_lpi[1, 'value']  
  
  
  # ------------------------------
  # Aggregation Index
  # ------------------------------
  # lsm_c_ai: equals the number of like adjacencies divided by the theoretical maximum possible number of like adjacencies for that class
  # Equals 0 for maximally disaggregated and 100 for maximally aggregated classes.
  c_ai <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_c_ai")
  woody_c_ai <- c_ai[which(c_ai$class == 0),] # Woody
  herb_c_ai <- c_ai[which(c_ai$class == 1),] # Herbaceous
  site_dat[row, "woody_AggInx"] <- woody_c_ai[1, 'value']     
  site_dat[row, "herb_AggInx"] <- herb_c_ai[1, 'value']  
  
  
  # ------------------------------
  # Edge Density 
  # ------------------------------
  # lsm_c_ed: equals the sum of all edges of class i in relation to the landscape area
  # = 0 if only one patch is present (and the landscape boundary is not included) and increases, without limit, as the landscapes becomes more patchy
  c_ed <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_c_ed")
  woody_c_ed <- c_ed[which(c_ed$class == 0),] # Woody
  herb_c_ed <- c_ed[which(c_ed$class == 1),] # Herbaceous
  site_dat[row, "woody_EdgDens"] <- woody_c_ed[1, 'value']   
  site_dat[row, "herb_EdgDens"] <- herb_c_ed[1, 'value'] 
  
  
  # ------------------------------
  # Patch Density
  # ------------------------------
  # lsm_c_pd: Number per 100 hectares
  # Increases as the landscape gets more patchy. Reaches its maximum if every cell is a different patch.
  c_pd <- landscapemetrics::calculate_lsm(lulc_clip, what = "lsm_c_pd")
  woody_c_pd <- c_pd[which(c_pd$class == 0),] # Woody
  herb_c_pd <- c_pd[which(c_pd$class == 1),] # Herbaceous
  site_dat[row, "woody_Pdens"] <- woody_c_pd[1, 'value']
  site_dat[row, "herb_Pdens"] <- herb_c_pd[1, 'value']
  
  
  # ------------------------------
  # Clumpy Patches
  # ------------------------------
  
  # Classify patches
  woody_rast <- raster(woody_clip)
  herb_rast <- raster(herb_clip)
  woody_clumps <- clump(woody_rast, directions = 8) # Id clumps
  herb_clumps <- clump(herb_rast, directions = 8)  
  
  # Plot and Export Woody Patches
  woodyclumps_Plotfile <- paste0(output_dir, SiteID, "_WoodyClumps.png")
  png(woodyclumps_Plotfile, width = 8, height = 6, units = "in", res = 300)  # Set resolution to 300 DPI
  plot(woody_clumps, main = "Woody Clumps", col = turbo(cellStats(woody_clumps, stat = "max"), direction = 1))
  dev.off()
  
  # Number of patches
  site_dat[row, "woody_Npatches"] <- cellStats(woody_clumps, stat = "max")
  site_dat[row, "herb_Npatches"] <- cellStats(herb_clumps, stat = "max")
  
  # ------------------------------
  # Focal Statistics 
  # ------------------------------
  
  # Create classification matrix: from, to, new value
  woody_FocClass <- (lulc_clip == 0)  # Create a mask where woody = 1, others = 0
  woody_FocClass <- as.numeric(woody_FocClass)  # Convert TRUE/FALSE to 1/0
  # plot(woody_FocClass) # check
  
  # 30m circle
  kernel30m <- focalMat(lulc_clip, d = 15, type = "circle")  
  woody_focal30m <- focal(woody_FocClass, w = kernel30m, fun = sum, na.rm = TRUE)
  mean_woody_focal30m <- global(woody_focal30m$focal_sum, fun = "mean", na.rm = TRUE)
  site_dat[row, "woody_mnFocal30m"] <- mean_woody_focal30m[1,1]
  
  
  # ------------------------
  #
  #         LiDAR
  #
  # ------------------------
  
  # Convert to SpatVector
  site_sub_vect <- vect(site_sub_coords)
  
  # Reproject to match `lasFiles` (NAD83 UTM Zone 14N)
  site_sub_vect_utm <- project(site_sub_vect, "+proj=utm +zone=14 +datum=NAD83 +units=m +no_defs")
  
  # Convert back to sp  
  site_sub_coords_utm <- as(site_sub_vect_utm, "Spatial")
  
  # Clipping LiDAR data
  sub_las_50m <- clip_roi(lasFiles, site_sub_coords_utm, radius = 50)
  sub_las_200m <- clip_roi(lasFiles, site_sub_coords_utm, radius = 200)
  #plot(sub_las_50m, color = "Classification", bg = "black", size = 5)
  #plot(sub_las_200m, color = "Classification", bg = "black", size = 5)
  
  # ------------------------------
  # Classifying Ground
  # ------------------------------
  
  # # Custom Cloth Simulation Function
  sub_las_mycsf <- classify_ground(sub_las_200m,
                                   csf(sloop_smooth = TRUE,
                                       class_threshold = 1,
                                       cloth_resolution = 1,
                                       time_step = 1))
  
  # Take a look
  #plot(sub_las_mycsf, color = "Classification", bg = "black", size = 5)
  
  # ------------------------------
  # Proportion Vegetation
  # ------------------------------
  
  # Putting extracted data into a table
  LiDAR_table <- table(sub_las_50m$Classification)
  
  # Calculate the total count for each habitat type
  total_count <- sum(LiDAR_table)
  
  # transform table to dataframe
  LiDAR_table <- as.data.frame(LiDAR_table)
  
  # Subset rows where Var1 is 3 (med veg) or 4 (high veg)
  LiDAR_subset <- LiDAR_table[LiDAR_table$Var1 %in% c(3, 4), ]
  
  # Calculate the proportions of med + high vegetation
  vegdens50m <- (sum(LiDAR_subset$Freq) / total_count)
  
  # Veg Density (prp LiDAR Veg points)
  site_dat[row, "vegDens50m"] <- vegdens50m 
  
  # ------------------------------
  # Digital Elevation/Terrain Model
  # ------------------------------
  
  # Create model
  dtm <- lidR::rasterize_terrain(sub_las_200m, res = 0.7, algorithm = knnidw(k = 10))
  mn_elev <- terra::global(dtm, fun = "mean", na.rm = TRUE) # Calculate mean
  site_dat[row, "mnElev"] <-  mn_elev[,1]
  # plot(dtm)
  
  
  # ------------------------------
  # Normalize Point Cloud
  # ------------------------------
  
  # Normalize the point cloud (subtracts ground elevation)
  sub_las_norm <- normalize_height(sub_las_200m, dtm)
  
  
  # ------------------------------
  # Canopy Height Model
  # ------------------------------
  
  # Create model using normalized point cloud
  chm <- rasterize_canopy(sub_las_norm, res = 0.7, pitfree(subcircle = 0.2))
  
  # Plot and Export
  chm_df <- as.data.frame(chm, xy = TRUE)
  chm_Plotfile <- paste0(output_dir, SiteID, "_CHM.png")
  chm_plot <- ggplot(chm_df) +
    geom_raster(aes(x = x, y = y, fill = Z)) +  # 'value' is the raster values column
    scale_fill_viridis_c(option = "H") +  # You can adjust the color scale here
    coord_fixed() +
    labs(title = paste("CHM Site", SiteID)) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA))
  ggsave(filename = chm_Plotfile, plot = chm_plot, width = 8, height = 6, dpi = 300)
  
  # ------------------------------
  # Focal Statistics
  # ------------------------------
  
  # Brush between 0.5 and 2 meters, heights above and below are in separate classes
  optim_chm <- classify(chm, c(0, 0.5, 2, Inf), c(0, 1, 0))
  optim_chm <- (optim_chm == 1)  # This returns a boolean mask
  optim_chm <- as.numeric(optim_chm)  # Convert TRUE/FALSE to 1/0
  #plot(optim_chm)
  
  # Convert CHM raster to a dataframe
  optim_chm_df <- as.data.frame(optim_chm, xy = TRUE)
  colnames(optim_chm_df) <- c("x", "y", "Class")
  optim_chm_df <- optim_chm_df[which(optim_chm_df[,'Class'] == 1),]
  
  # Reproject lulc_clip to the CRS of optim_chm
  lulc_clip_utm <- project(lulc_clip, terra::crs(optim_chm))
  
  # Resample lulc_clip to match the resolution and extent of optim_chm
  lulc_clip__utm <- resample(lulc_clip_utm, optim_chm, method = "near")
  lulc_utm_df <- as.data.frame(lulc_clip__utm, xy = TRUE)
  
  # Plot LULC with optimal CHM overlay
  overlayPlotfile <- paste0(output_dir, SiteID, "_Overlay.png")
  overlayPlot <- ggplot() +
    geom_raster(data = lulc_utm_df, aes(x = x, y = y, fill = factor(Classvalue))) +
    scale_fill_manual(name = "Land Cover Type",
                      values = c("1" = "forestgreen",
                                 "2" = "lightgreen",
                                 "3" = "saddlebrown"),
                      labels = c("1" = "Woody",
                                 "2" = "Herbaceous",
                                 "3" = "Bareground")) +
    geom_tile(data = optim_chm_df, aes(x = x, y = y), fill = "yellow", alpha = 0.50) + 
    coord_fixed() +
    labs(title = paste("LULC and CHM (0.5-2m [grey]) Overlay for", SiteID )) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA))
  # Export
  ggsave(filename = overlayPlotfile, plot = overlayPlot, width = 8, height = 6, dpi = 300)
  
  
  # 30m at Focal Statistics
  kernel30m <- focalMat(chm, d = 15, type = "circle")  
  optim_focal30m <- focal(optim_chm, w = kernel30m, fun = sum, na.rm = TRUE)
  mean_optim_focal30m <- global(optim_focal30m$focal_sum, fun = "mean", na.rm = TRUE)
  site_dat[row, "optim_mnFocal30m"] <- mean_optim_focal30m[1,1]
  
  # ------------------------------
  # Update the progress bar
  # ------------------------------
  pb$tick()
  
} # -------------------- End Extraction Loop -----------------------------


# Take a look
str(site_dat)
print(site_dat)


# Export data
write.csv(site_dat, "./Data/Acoustic_Data/ARU_siteCovs.csv") # .csv




# --------------------End Script -----------------------------

site_dat<- read.csv("./Data/Acoustic_Data/ARU_siteCovs.csv")

cor_mat <- cor(site_dat[-c((1:4))])
 

# Heat Map
library(reshape2)
library(psych)

# Assuming cor_mat is your correlation matrix
cor_mat_melted <- melt(cor_mat)

# Create the heatmap
ggplot(cor_mat_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Correlation Heatmap", x = "", y = "")

# Checking correlation
pairs.panels(site_dat[-c((1:4))],
             gap = 0,
             bg = c("blue", "red"),
             pch = 21, main = "ARU Site Covs Correlations")

# Filter to less than 30% correlated
low_corr <- which(abs(cor_mat) < 0.3, arr.ind = TRUE)

# Remove redundant and diagonal elements
low_corr <- low_corr[low_corr[, 1] > low_corr[, 2], ]

# Filter to keep only pairs where one variable starts with "woody_" and the other with "herb_"
woody_herb_pairs <- low_corr[
  (grepl("^woody_", rownames(cor_mat)[low_corr[, 1]]) & grepl("^herb_", colnames(cor_mat)[low_corr[, 2]])) |
    (grepl("^herb_", rownames(cor_mat)[low_corr[, 1]]) & grepl("^woody_", colnames(cor_mat)[low_corr[, 2]])), 
]

# Create formatted output of low correlation pairs
low_corr_pairs <- apply(woody_herb_pairs, 1, function(idx) 
  paste(rownames(cor_mat)[idx[1]], colnames(cor_mat)[idx[2]], sep = " - ")
)

print(low_corr_pairs)
str(low_corr_pairs)

# Now convert the character vector into a dataframe
ARU_combo <- data.frame(pairs = low_corr_pairs)


PC_combo <-read.csv("PC_model_combinations.csv")

# Find common pairs
common_pairs <- intersect(ARU_combo$pairs, PC_combo$pairs)

# Print common pairs
print(common_pairs)
# "woody_EdgDens - herb_prp"  "woody_lrgPInx - herb_ClmIdx"  "woody_AggInx - herb_ShpInx"



