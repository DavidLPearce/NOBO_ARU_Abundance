
library(tidyverse)

# Load in data
pc_dat <- read.csv("./NOBO_PC_Summer2024data.csv")

print(pc_dat)

# Remove rows with NA values
pc_dat <- na.omit(pc_dat)

# Format date
pc_dat$Date <- as.POSIXct(pc_dat$Date, format = "%m/%d/%Y")

# Creating a survey number column
pc_dat$Survey <- as.numeric(as.factor(pc_dat$Date))

# Create a day of year column
pc_dat$DOY <- yday(pc_dat$Date)

# Add a unique id column based on Point Count number, NOBO number, and Day of Year
pc_dat$UniqueID <- paste(pc_dat$PointNum, pc_dat$NOBOnum, pc_dat$DOY, sep = "_")

# Columns to subset
cols_sub <- c( "UniqueID", "Survey","int1", "int2", "int3", "int4")

# Subset the capture-mark data
pc_CMR <- pc_dat[, cols_sub, drop = FALSE]
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

# Covariate for time interval
intervalMat <- matrix(c('1','2','3', '4'), nrow(pc_CMRtab), 4, byrow=TRUE)

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


# Bundle into unmarked object
cmr_umf <- unmarkedFrameMPois(y = pc_CMRtab,
                              siteCovs = NULL,
                              obsCovs = list(interval=intervalMat), obsToY=o2y, piFun="crPiFun")



M0 <- multinomPois(~ 1 ~ 1, cmr_umf)
print(M0)

Mt <- multinomPois(~ interval - 1 ~ 1, cmr_umf)
print(Mt)


fl <- modSel(fitList(M0, Mt))
print(fl)


# 7.8.3 Models with individual heterogeneity model
# ------------------------------------------------------------------------
parID <- matrix(c('p','sig','sig', 'sig'), nrow(pc_CMRtab), 4, byrow=TRUE)

umf.cr2 <- unmarkedFrameMPois(y = alfl.H1,
                              siteCovs = NULL,
                              obsCovs=list(parID=parID), obsToY=o2y, piFun="MhPiFun")


Mtx <- multinomPois(~ parID-1 ~ woody, umf.cr2)











