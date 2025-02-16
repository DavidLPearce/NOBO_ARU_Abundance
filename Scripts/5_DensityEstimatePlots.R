
library(tidyverse)

# Read in density summaries
cmr_Densdat <-readRDS("./Data/Fitted_Models/PC_CMR_dens_summary.rds")
hds_Densdat <- readRDS("./Data/Fitted_Models/PC_HDS_dens_summary.rds")
AVbnet_Densdat <- readRDS("./Data/Fitted_Models/ARU_BnetAV_dens_summary.rds")
AVwolfet_Densdat <- readRDS("./Data/Fitted_Models/ARU_WolfeAV_dens_summary.rds")


# Mean density
cmr_Densdat$Mean
hds_Densdat$Mean
AVbnet_Densdat$Mean
AVwolfet_Densdat$Mean

# Combine into one dataframe
dens_df <- rbind(cmr_Densdat, hds_Densdat, AVbnet_Densdat, AVwolfet_Densdat) # add wolfe
print(dens_df)

# Order
dens_df <- dens_df %>%
  mutate(Model = factor(Model, levels = c("PC CMR", "PC HDS", "AV Bnet", "AV Wolfe")))  # order

# Plot Density
densityPlot <- ggplot(dens_df, aes(x = Model, y = Mean, color = Model)) +
                geom_point(size = 3) +   
                geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, linewidth = 1) +  
                labs(y = "Density (N/acre)", x = "Model") +
                scale_y_continuous(limits = c(0, 0.4),
                                   breaks = seq(0, 0.4, by = 0.1),
                                   labels = scales::comma) +
                scale_color_manual(values = c("PC CMR" = "orange",   
                                              "PC HDS" = "purple",
                                              "AV Bnet" = "blue",
                                              "AV Wolfe" = "red")) + 
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                      axis.title.x = element_text(face = "bold", margin = margin(t = 10)),  
                      axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
                      panel.grid = element_blank(),
                      legend.position = "none") 


densityPlot
# Export                
ggsave(plot = densityPlot, "Figures/densityPlot.jpeg",  
       width = 8, height = 5, dpi = 300) 



# Read in abundance summaries
cmr_Abunddat <-readRDS("./Data/Fitted_Models/PC_CMR_abund_summary.rds")
hds_Abunddat <- readRDS("./Data/Fitted_Models/PC_HDS_abund_summary.rds")
AVbnet_Abunddat <- readRDS("./Data/Fitted_Models/ARU_BnetAV_abund_summary.rds")
AVwolfet_Abunddat <- readRDS("./Data/Fitted_Models/ARU_WolfeAV_abund_summary.rds")


# Combine into one dataframe
abund_df <- rbind(cmr_Abunddat, hds_Abunddat, AVbnet_Abunddat, AVwolfet_Abunddat) # add wolfe
print(abund_df)

# Order
abund_df <- abund_df %>%
  mutate(Model = factor(Model, levels = c("PC CMR", "PC HDS", "AV Bnet", "AV Wolfe" )))  # Define desired order


# Plot Density
abundPlot <- ggplot(abund_df, aes(x = Model, y = Mean, color = Model)) +
  geom_point(size = 3) +   
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, size = 1) +  
  labs(y = "Total Abundance", x = "Model") +
  scale_y_continuous(limits = c(0, 900),
                     breaks = seq(0, 900, by = 100),
                     labels = scales::comma) +
  scale_color_manual(values = c("PC CMR" = "orange",   
                                "PC HDS" = "purple",
                                "AV Bnet" = "blue",
                                "AV Wolfe" = "red")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),  
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
        panel.grid = element_blank(),
        legend.position = "none")

abundPlot

# Export                
ggsave(plot = abundPlot, "Figures/abundPlot.jpeg",  
       width = 8, height = 5, dpi = 300) 


