# -------------------------------------------------------
#
#                    Load libraries
#
# -------------------------------------------------------


library(tidyverse)

# -------------------------------------------------------
#
#                    Load Data
#
# -------------------------------------------------------

# Read in density dataframes
cmr_Dens_df <-readRDS("./Data/Fitted_Models/PC_CMR_dens_df.rds")
hds_Dens_df <- readRDS("./Data/Fitted_Models/PC_HDS_dens_df.rds")
AVbnet_Dens_df <- readRDS("./Data/Fitted_Models/ARU_BnetAV_dens_df.rds")
AVwolfet_Dens_df <- readRDS("./Data/Fitted_Models/ARU_WolfeAV__dens_df.rds")

# Read in density summaries
cmr_Dens_sum <-readRDS("./Data/Fitted_Models/PC_CMR_dens_summary.rds")
hds_Dens_sum <- readRDS("./Data/Fitted_Models/PC_HDS_dens_summary.rds")
AVbnet_Dens_sum <- readRDS("./Data/Fitted_Models/ARU_BnetAV_dens_summary.rds")
AVwolfet_Dens_sum <- readRDS("./Data/Fitted_Models/ARU_WolfeAV_dens_summary.rds")

# Read in abundance summaries
cmr_Abund_sum <-readRDS("./Data/Fitted_Models/PC_CMR_abund_summary.rds")
hds_Abund_sum <- readRDS("./Data/Fitted_Models/PC_HDS_abund_summary.rds")
AVbnet_Abund_sum <- readRDS("./Data/Fitted_Models/ARU_BnetAV_abund_summary.rds")
AVwolfet_Abund_sum <- readRDS("./Data/Fitted_Models/ARU_WolfeAV_abund_summary.rds")

# -------------------------------------------------------
#
#                   Data Wrangling
#
# -------------------------------------------------------

# Combine density df into one dataframe
dens_df <- rbind(cmr_Dens_df, hds_Dens_df, AVbnet_Dens_df, AVwolfet_Dens_df)  
head(dens_df)

# Combine density summaries into one dataframe
dens_sum <- rbind(cmr_Dens_sum, hds_Dens_sum, AVbnet_Dens_sum, AVwolfet_Dens_sum)  
print(dens_sum)

# Combine abundances into one dataframe
abund_sum <- rbind(cmr_Abund_sum, hds_Abund_sum, AVbnet_Abund_sum, AVwolfet_Abund_sum) 
print(abund_sum)

# Order
dens_df <- dens_df %>% mutate(Model = factor(Model, levels = c("PC CMR", "PC HDS", "AV Bnet", "AV Wolfe")))
dens_sum <- dens_sum %>% mutate(Model = factor(Model, levels = c("PC CMR", "PC HDS", "AV Bnet", "AV Wolfe")))   
abund_sum <- abund_sum %>% mutate(Model = factor(Model, levels = c("PC CMR", "PC HDS", "AV Bnet", "AV Wolfe" ))) 



# -------------------------------------------------------
#
#                   Plots
#
# -------------------------------------------------------


# -------------------------
# Plot Density - Mean & CI
# -------------------------
densityPlot <- ggplot(dens_sum, aes(x = Model, y = Mean, color = Model)) +
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

# View
densityPlot

# Export                
ggsave(plot = densityPlot, "./Figures/densityPlot.jpeg", width = 8, height = 5, dpi = 300)




# -------------------------
# Plot Density - violin
# -------------------------

densityViolin <- ggplot(dens_df, aes(x = Model, y = Density, fill = Model)) + 
  geom_violin(trim = FALSE, alpha = 0.6, adjust = 5) +  
  stat_summary(fun = mean, geom = "point", shape = 20, 
               size = 3, fill = "black") +   
  labs(x = "Model", y = "Density (N/acre)") +
  scale_fill_manual(values = c("PC CMR" = "orange", 
                               "PC HDS" = "purple", 
                               "AV Bnet" = "blue",
                               "AV Wolfe" = "red")) +  
  scale_y_continuous(limits = c(0, 0.35),
                     breaks = seq(0, 0.35, by = 0.05),
                     labels = scales::comma) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),  
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
        panel.grid = element_blank(),
        legend.position = "none")  

# View
densityViolin

# Export                
ggsave(plot = densityViolin, "./Figures/densityViolin.jpeg",width = 8, height = 5, dpi = 300) 









# -------------------------
# Plot Abundance - Mean & CI
# -------------------------

# Plot Abundance
abundPlot <- ggplot(abund_sum, aes(x = Model, y = Mean, color = Model)) +
  geom_point(size = 3) +   
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, linewidth = 1) +  
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

# View
abundPlot

# Export                
ggsave(plot = abundPlot, "./Figures/abundPlot.jpeg",width = 8, height = 5, dpi = 300) 

# Print
print(as.data.frame(abund_sum))


# -------------------------
# Plot Abunndance - violin
# -------------------------

# Create Abundance dataframe
abund_df <- dens_df
abund_df$Density <- abund_df$Density * 2710
colnames(abund_df)[2] <- 'Abundance'
head(abund_df)


# Plot violin
abundViolin <- ggplot(abund_df, aes(x = Model, y = Abundance, fill = Model)) + 
                  geom_violin(trim = FALSE, alpha = 0.6, adjust = 5) +  
                  stat_summary(fun = mean, geom = "point", shape = 20, 
                               size = 3, fill = "black") +   
                  labs(x = "Model", y = "Total Abundance") +
                  scale_fill_manual(values = c("PC CMR" = "orange", 
                                               "PC HDS" = "purple", 
                                               "AV Bnet" = "blue",
                                               "AV Wolfe" = "red")) +  
                  scale_y_continuous(limits = c(0, 900),
                                     breaks = seq(0, 900, by = 100),
                                     labels = scales::comma) +
                  theme_minimal() +
                  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  
                        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),  
                        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
                        panel.grid = element_blank(),
                        legend.position = "none")  

# View
abundViolin

# Export                
ggsave(plot = abundViolin, "./Figures/abundViolin.jpeg",width = 8, height = 5, dpi = 300) 








