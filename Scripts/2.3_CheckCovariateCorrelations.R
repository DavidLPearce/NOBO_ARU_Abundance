# ------------------------------------------------------------------------------
#
#                               Load Packages
#
# ------------------------------------------------------------------------------

library(reshape2) # used for heatmap
library(psych) # used for pairs panel plotting

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

# Remove Site column
pc_covs <- pc_covs[,-c(1:2)]
aru_dat <- aru_dat[,-c(1:2)]


# Create correlation matrix
pc_cor_mat <- cor(pc_covs)
aru_cor_mat <- cor(aru_dat)

# Format for heat map
pc_cor_mat_melted <- melt(pc_cor_mat)
aru_cor_mat_melted <- melt(aru_cor_mat)

# Plot Correlations

# PC
png("./Figures/Correlations/PC_Corr_Mat.png", width = 1800, height = 1200, res = 150)
pairs.panels(pc_cor_mat, gap = 0, bg = c("blue", "red"), pch = 21, main = "PC")
dev.off()



# ARU
png("./Figures/Correlations/ARU_Corr_Mat.png", width = 1800, height = 1200, res = 150)
pairs.panels(aru_cor_mat, gap = 0, bg = c("blue", "red"), pch = 21, main = "ARU")
dev.off()


 
# Plot Heatmap

# PC
pc_heatmap <- ggplot(pc_cor_mat_melted, aes(Var1, Var2, fill = value)) +
                      geom_tile() +
                      scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1)) +
                      theme_minimal() +
                      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
                            axis.text.y = element_text(hjust = 1, size = 14)) +
                      labs(title = "PC", x = "", y = "")

ggsave("./Figures/Correlations/PC_Corr_Heatmap.jpg", plot = pc_heatmap, width = 14, height = 12, dpi = 300)

# ARU
aru_heatmap <- ggplot(aru_cor_mat_melted, aes(Var1, Var2, fill = value)) +
                      geom_tile() +
                      scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1)) +
                      theme_minimal() +
                      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
                            axis.text.y = element_text(hjust = 1, size = 14)) +
                      labs(title = "ARU", x = "", y = "")

ggsave("./Figures/Correlations/ARU_Corr_Heatmap.jpg", plot = aru_heatmap, width = 14, height = 12, dpi = 300)



# Filter covarite combinations less than 30% correlated
pc_low_corr <- which(abs(pc_cor_mat) < 0.40, arr.ind = TRUE)
aru_low_corr <- which(abs(aru_cor_mat) < 0.40, arr.ind = TRUE)


# Remove redundant and diagonal elements
pc_low_corr <- pc_low_corr[pc_low_corr[, 1] > pc_low_corr[, 2], , drop = FALSE]
aru_low_corr <- aru_low_corr[aru_low_corr[, 1] > aru_low_corr[, 2], , drop = FALSE]


# Filter to keep only pairs where one variable starts with "woody_" and the other with "herb_"
pc_pairs <- pc_low_corr[(grepl("^woody_", rownames(pc_cor_mat)[pc_low_corr[, 1]]) & 
                          grepl("^herb_", colnames(pc_cor_mat)[pc_low_corr[, 2]])) |
                           (grepl("^herb_", rownames(pc_cor_mat)[pc_low_corr[, 1]]) &
                             grepl("^woody_", colnames(pc_cor_mat)[pc_low_corr[, 2]])),]
  
aru_pairs <- aru_low_corr[(grepl("^woody_", rownames(aru_cor_mat)[aru_low_corr[, 1]]) & 
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
low_pairs <- intersect(aru_low_corr_pairs$pairs, pc_low_corr_pairs$pairs)
low_pairs <- as.data.frame(low_pairs)
print(low_pairs)

# Export
write.csv(low_pairs, "Low_Corr_Metric_Combinations.csv")

# --------------------End Script -----------------------------