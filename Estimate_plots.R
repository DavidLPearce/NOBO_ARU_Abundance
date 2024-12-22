
library(ggplot2)

hds_dat <- readRDS("./Data/Fitted_Models/PC_HDS_Dens.rds")
nmix_dat<-readRDS("./Data/Fitted_Models/PC_Nmix_Dens.rds")

lat_dens_df <- rbind(hds_dat, nmix_dat)

# Plot
ggplot(lat_dens_df, aes(x = Model, y = Density, fill = Model)) +
  geom_violin() + 
  geom_boxplot(aes(x = Model, y = Density), 
               width = 0.2, position = position_dodge(width = 0.8)) +
  labs(
    title = "Latent Density", 
    x = "Model", 
    y = "Density") +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 5), labels = scales::comma) + # Customize y-axis
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1), # Tilt x-axis text
    axis.text.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none" 
  )
