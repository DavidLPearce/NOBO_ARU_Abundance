# ------------------------------------------------------------------------------
#
#                               Load Packages
#
# ------------------------------------------------------------------------------

library(reshape2)


# ------------------------------------------------------------------------------
#
#                               Read in Data
#
# ------------------------------------------------------------------------------

pc_covs <- read.csv("./Data/Point_Count_Data/PointCount_siteCovs.csv")
aru_dat<- read.csv("./Data/Acoustic_Data/ARU_siteCovs.csv")

# ------------------------------------------------------------------------------
#
#                               Check Correlations
#
# ------------------------------------------------------------------------------


# Remove site, lat, long, veg density 
pc_covs <- aru_dat[-c( )]
aru_dat <- aru_dat[-c( )]
 

# Create correlation matrix
pc_cor_mat <- cor(site_dat[-c((1:4))])
aru_cor_mat <- cor(site_dat[-c((1:4))])

# Format for heat map
pc_cor_mat_melted <- melt(pc_cor_mat)
aru_cor_mat_melted <- melt(aru_cor_mat)
 
# Plot Heatmap

# PC
ggplot(pc_cor_mat_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(hjust = 1, size = 14)) +
  labs(title = "Point Count Correlation Heatmap", x = "", y = "")

# ARU
ggplot(aru_cor_mat_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(hjust = 1, size = 14)) +
  labs(title = "ARU Correlation Heatmap", x = "", y = "")


# Filter covarite combinations less than 45% correlated
pc_low_corr <- which(abs(pc_cor_mat) < 0.45, arr.ind = TRUE)
aru_low_corr <- which(abs(aru_cor_mat) < 0.45, arr.ind = TRUE)


# Remove redundant and diagonal elements
pc_low_corr <- pc_low_corr[low_corr[, 1] > pc_low_corr[, 2], ]
aru_low_corr <- aru_low_corr[low_corr[, 1] > aru_low_corr[, 2], ]


# Filter to keep only pairs where one variable starts with "woody_" and the other with "herb_"
pc_pairs <- low_corr[(grepl("^woody_", rownames(pc_cor_mat)[pc_low_corr[, 1]]) & 
                        grepl("^herb_", colnames(pc_cor_mat)[pc_low_corr[, 2]])) |
                          (grepl("^herb_", rownames(pc_cor_mat)[pc_low_corr[, 1]]) &
                             grepl("^woody_", colnames(pc_cor_mat)[pc_low_corr[, 2]])),]

aru_pairs <- low_corr[(grepl("^woody_", rownames(aru_cor_mat)[aru_low_corr[, 1]]) & 
                        grepl("^herb_", colnames(aru_cor_mat)[aru_low_corr[, 2]])) |
                         (grepl("^herb_", rownames(aru_cor_mat)[aru_low_corr[, 1]]) &
                           grepl("^woody_", colnames(aru_cor_mat)[aru_low_corr[, 2]])),]


# Create formatted output of low correlation pairs
pc_low_corr_pairs <- apply(pc_pairs, 1, function(idx) 
                        paste(rownames(pc_cor_mat)[idx[1]], colnames(pc_cor_mat)[idx[2]], sep = " + "))     

aru_low_corr_pairs <- apply(aru_pairs, 1, function(idx) 
                        paste(rownames(aru_cor_mat)[idx[1]], colnames(aru_cor_mat)[idx[2]], sep = " + "))  



# Convert to dataframe
pc_low_corr_pairs <- data.frame(pairs = pc_low_corr_pairs)
aru_low_corr_pairs <- data.frame(pairs = aru_low_corr_pairs)

# Uncorrelated herbaceous and woody metrics acrossed datasets
low_pairs <- intersect(ARU_combo$pairs, PC_combo$pairs)
print(low_pairs)


# --------------------End Script -----------------------------