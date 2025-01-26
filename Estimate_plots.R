
library(tidyverse)

hds_dat <- readRDS("./Data/Fitted_Models/PC_HDS_Dens.rds")
nmix_dat<-readRDS("./Data/Fitted_Models/PC_Nmix_Dens.rds")
tmpEhds_dat <- readRDS("./Data/Fitted_Models/PC_TempEmHDS_Dens.rds")
cmr_dat <-readRDS("./Data/Fitted_Models/PC_CMR_Dens.rds")
aruBnet_dat <- readRDS("./Data/Fitted_Models/ARU_Bnet.rds")
  
head(hds_dat)

# Mean density
mean(hds_dat$Density)
mean(nmix_dat$Density)
mean(tmpEhds_dat$Density)
mean(cmr_dat$Density)
mean(aruBnet_dat$Density)

# Mean study area abundance
mean(hds_dat$Density) * 1096.698
mean(nmix_dat$Density) * 1096.698
mean(tmpEhds_dat$Density) * 1096.698
mean(cmr_dat$Density) * 1096.698
mean(aruBnet_dat$Density) * 1096.698

# Typical abundance of quail in Texas is 1-2 birds per acre
# which would be 2.47 to 4.94 birds per hectare.
# adding in as a "truth"
txDens_dat <- as.data.frame(c(2.47, 4.94))
colnames(txDens_dat)[1] <- "Density"
txDens_dat$Model <- "Mean Texas"
txDens_dat <- txDens_dat[, c(2, 1)]# Switch columns
head(txDens_dat)
mean(txDens_dat$Density)

# Combine into one dataframe
dens_df <- rbind(txDens_dat, hds_dat, nmix_dat, tmpEhds_dat, cmr_dat, aruBnet_dat)
head(dens_df)




# Calculate means for each model
means <- dens_df %>%
  group_by(Model) %>%
  summarize(mean_density = mean(Density))

# Plot
ggplot(dens_df, aes(x = Model, y = Density, fill = Model)) +
  geom_violin(trim = FALSE, alpha = 0.7) + 
  geom_point(data = means, aes(x = Model, y = mean_density), color = "black", size = 3) +
  geom_hline(yintercept = 3.705, color = "black", linetype = "dashed", size = 1) +  # Add horizontal line
  labs(
    title = "", 
    x = "Model", 
    y = "Abundance/Hectare") +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2), labels = scales::comma) + # Customize y-axis
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1), # Tilt x-axis text
    axis.text.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

