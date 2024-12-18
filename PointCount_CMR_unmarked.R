# Load Library
library(tidyverse)
library(unmarked)

# Load in capture data
pc_dat <- read.csv("./Data/Point_Count_Data/NOBO_PC_Summer2024data.csv")
print(pc_dat)
pc_dat <- na.omit(pc_dat)# Remove rows with NA values

# Load in site covariates
site_covs <- read.csv("./Data/Point_Count_Data/PointCount_siteCovs.csv")

# Match Match site covs with point count data on the Point count number
head(pc_dat)
head(site_covs)
pc_dat <- merge(pc_dat, site_covs, by = "PointNum", all.x = TRUE)

# Format date
pc_dat$Date <- as.POSIXct(pc_dat$Date, format = "%m/%d/%Y")

# Create a day of year column
pc_dat$DOY <- yday(pc_dat$Date)

# Add a unique id column based on Point Count number, NOBO number, and Day of Year
pc_dat$UniqueID <- paste(pc_dat$PointNum, pc_dat$NOBOnum, pc_dat$DOY, sep = "_")

# Subset the capture-mark data
pc_CMR <- pc_dat[, c( "UniqueID", "Survey","int1", "int2", "int3", "int4",
                      "Observor", "Temp.deg.F", "Wind.Beau.Code", "Sky.Beau.Code",
                      "DOY", "woody_prp", "herb_prp", "baregnd_prp", 
                      "woody_mean_p_Area", "woody_c_clumpy"),
                 drop = FALSE]
head(pc_CMR)


# Adding a capture history
pc_CMR$captureHistory <- paste(pc_CMR$int1, pc_CMR$int2, pc_CMR$int3, pc_CMR$int4, sep="")

# Capture history as a factor
pc_CMR$captureHistory <- factor(pc_CMR$captureHistory,
                              levels=c("0001","0010","0011","0100",
                                       "0101","0110","0111","1000",
                                       "1001","1010","1011","1100",
                                       "1101","1110","1111"))

# Setting detected NOBOs to a factor
pc_CMR$UniqueID <- factor(pc_CMR$UniqueID)

# Tabulating the detentions
pc_CMRtab <- table(pc_CMR$UniqueID, pc_CMR$captureHistory)
class(pc_CMRtab) <- "matrix"
nrow(pc_CMRtab)

# Obs cov for time interval
inter_mat <- matrix(c('1','2','3', '4'), nrow(pc_CMRtab), 4, byrow=TRUE)

# Obs cov for observer
obser_mat <- matrix(NA, nrow(pc_CMRtab), 4, byrow=TRUE)
obser_mat[,1] <- pc_CMR$Observor
obser_mat[,2] <- pc_CMR$Observor
obser_mat[,3] <- pc_CMR$Observor
obser_mat[,4] <- pc_CMR$Observor

# Obs cov for day of year
doy_mat <- matrix(NA, nrow(pc_CMRtab), 4, byrow=TRUE)
doy_mat[,1] <- pc_CMR$DOY
doy_mat[,2] <- pc_CMR$DOY
doy_mat[,3] <- pc_CMR$DOY
doy_mat[,4] <- pc_CMR$DOY

# Temp.deg.F
temp_mat <- matrix(NA, nrow(pc_CMRtab), 4, byrow=TRUE)
temp_mat[,1] <- pc_CMR$Temp.deg.F
temp_mat[,2] <- pc_CMR$Temp.deg.F
temp_mat[,3] <- pc_CMR$Temp.deg.F
temp_mat[,4] <- pc_CMR$Temp.deg.F

# Wind.Beau.Code
wind_mat <- matrix(NA, nrow(pc_CMRtab), 4, byrow=TRUE)
wind_mat[,1] <- pc_CMR$Wind.Beau.Code
wind_mat[,2] <- pc_CMR$Wind.Beau.Code
wind_mat[,3] <- pc_CMR$Wind.Beau.Code
wind_mat[,4] <- pc_CMR$Wind.Beau.Code

# Sky.Beau.Code
sky_mat <- matrix(NA, nrow(pc_CMRtab), 4, byrow=TRUE)
sky_mat[,1] <- pc_CMR$Sky.Beau.Code
sky_mat[,2] <- pc_CMR$Sky.Beau.Code
sky_mat[,3] <- pc_CMR$Sky.Beau.Code
sky_mat[,4] <- pc_CMR$Sky.Beau.Code

# Matrix describing the relationship between obs covs and y
o2y <- matrix(1, 4, 15)

# Capture Probability statement
crPiFun <- function(p) {
  p1 <- p[,1]
  p2 <- p[,2]
  p3 <- p[,3]
  p4 <- p[,4]
  cbind(
        "0001" = (1 - p1) * (1 - p2) * (1 - p3) *      p4,
        "0010" = (1 - p1) * (1 - p2) *      p3  * (1 - p4),
        "0011" = (1 - p1) * (1 - p2) *      p3  *      p4,
        "0100" = (1 - p1) *      p2  * (1 - p3) * (1 - p4),
        "0101" = (1 - p1) *      p2  * (1 - p3) *      p4,
        "0110" = (1 - p1) *      p2  *      p3  *      p4,
        "0111" = (1 - p1) *      p2 *       p3  *      p4,
        "1000" =      p1  * (1 - p2) * (1 - p3) * (1 - p4),
        "1001" =      p1  * (1 - p2) * (1 - p3) *      p4,
        "1010" =      p1  * (1 - p2) *      p3  * (1 - p4),
        "1011" =      p1  * (1 - p2) *      p3  *      p4,
        "1100" =      p1 *       p2  * (1 - p3) * (1 - p4), 
        "1101" =      p1 *       p2  * (1 - p3) *      p4,
        "1110" =      p1 *       p2  *      p3  * (1 - p4),    
        "1111" =      p1 *       p2  *      p3  *      p4
  )
} # end statement

head(pc_CMR)
head(pc_CMRtab)

# Reorder pc_CMR to match the rownames of pc_CMRtab
pc_CMR <- pc_CMR[match(rownames(pc_CMRtab), pc_CMR$UniqueID), ]

# Bundle into unmarked object
cmr_umf <- unmarkedFrameMPois(y = pc_CMRtab,
                              siteCovs = pc_CMR[,12:16],
                              obsCovs = list(interval = inter_mat,
                                             observor = obser_mat,
                                             doy = doy_mat,
                                             temp = temp_mat,
                                             wind = wind_mat,
                                             sky = sky_mat), 
                              obsToY=o2y, 
                              piFun="crPiFun")

# Take a look
str(cmr_umf)


# Null
fm.0 <- multinomPois(~ 1 ~ 1, cmr_umf)
print(fm.0)

# Interval
detfm.1 <- multinomPois(~ interval - 1 ~ 1, cmr_umf)
print(detfm.1)

# observor
detfm.2 <- multinomPois(~ observor ~ 1, cmr_umf)
print(detfm.2)

# doy
detfm.3 <- multinomPois(~ doy ~ 1, cmr_umf)
print(detfm.3)

# temp
detfm.4 <- multinomPois(~ temp ~ 1, cmr_umf)
print(detfm.4)

# wind
detfm.5 <- multinomPois(~ wind ~ 1, cmr_umf)
print(detfm.5)

# sky
detfm.6 <- multinomPois(~ sky ~ 1, cmr_umf)
print(detfm.6)

# Model selection
model_list <- modSel(fitList(fm.0, 
                             detfm.1,
                             detfm.2,
                             detfm.3,
                             detfm.4,
                             detfm.5,
                             detfm.6))


print(model_list)

# sky shows most support, so does doy
# sky + doy
detfm.7 <- multinomPois(~ sky + doy ~ 1, cmr_umf)
print(detfm.7)

# Model selection
model_list <- modSel(fitList(fm.0, 
                             detfm.1,
                             detfm.2,
                             detfm.3,
                             detfm.4,
                             detfm.5,
                             detfm.6,
                             detfm.7))


print(model_list)

# Sky is the most informative detection covariate

# woody_prp
fm.1 <- multinomPois(~ sky 
                     ~ woody_prp, 
                     cmr_umf)


# herb_prp
fm.2 <- multinomPois(~ sky 
                     ~ herb_prp, 
                     cmr_umf)

# baregnd_prp
fm.3 <- multinomPois(~ sky 
                     ~ baregnd_prp, 
                     cmr_umf)

# woody_mean_p_Area
fm.4 <- multinomPois(~ sky 
                     ~ baregnd_prp, 
                     cmr_umf)

# woody_c_clumpy
fm.5 <- multinomPois(~ sky 
                     ~ woody_c_clumpy, 
                     cmr_umf)

# herb_prp + woody_mean_p_Area
fm.6 <- multinomPois(~ sky 
                     ~ herb_prp + woody_mean_p_Area, 
                     cmr_umf)

# herb_prp + woody_c_clumpy
fm.7 <- multinomPois(~ sky 
                     ~ herb_prp + woody_c_clumpy, 
                     cmr_umf)

# herb_prp + woody_mean_p_Area + woody_c_clumpy
fm.8 <- multinomPois(~ sky 
                     ~ herb_prp + woody_mean_p_Area + woody_c_clumpy, 
                     cmr_umf)

# Model selection
model_list <- modSel(fitList(fm.0, 
                             #fm.1,
                             fm.2,
                             #fm.3,
                             #fm.4,
                             fm.5,
                             fm.6,
                             fm.7,
                             fm.8))


print(model_list)




# 7.8.3 Models with individual heterogeneity model
# ------------------------------------------------------------------------
parID <- matrix(c('p','sig','sig', 'sig'), nrow(pc_CMRtab), 4, byrow=TRUE)

umf.cr2 <- unmarkedFrameMPois(y = alfl.H1,
                              siteCovs = NULL,
                              obsCovs=list(parID=parID), obsToY=o2y, piFun="MhPiFun")


Mtx <- multinomPois(~ parID-1 ~ woody, umf.cr2)











