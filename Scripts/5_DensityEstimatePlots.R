# -------------------------------------------------------
#
#                    Load libraries
#
# -------------------------------------------------------


library(tidyverse)
library(gridExtra)

# -------------------------------------------------------
#
#                    Load Data
#
# -------------------------------------------------------

# Read in beta dataframes
cmr_beta_df <-readRDS("./Data/Fitted_Models/PC_CMR_beta_df.rds")
hds_beta_df <- readRDS("./Data/Fitted_Models/PC_HDS_beta_df.rds")
AVbnet_beta_df <- readRDS("./Data/Fitted_Models/ARU_BnetAV_beta_df.rds")
AVwolfet_beta_df <- readRDS("./Data/Fitted_Models/ARU_WolfeAV_beta_df.rds")

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

# Combine beta df into one dataframe
beta_df <- rbind(cmr_beta_df, hds_beta_df, AVbnet_beta_df, AVwolfet_beta_df)  
head(beta_df)
unique(beta_df$Model)

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
beta_df <- beta_df %>% mutate(Model = factor(Model, levels = c("PC CMR", "PC HDS", "Bnet AV", "Wolfe AV")))
dens_df <- dens_df %>% mutate(Model = factor(Model, levels = c("PC CMR", "PC HDS", "Bnet AV", "AV Wolfe")))
dens_sum <- dens_sum %>% mutate(Model = factor(Model, levels = c("PC CMR", "PC HDS", "Bnet AV", "AV Wolfe")))   
abund_sum <- abund_sum %>% mutate(Model = factor(Model, levels = c("PC CMR", "PC HDS", "Bnet AV", "AV Wolfe" ))) 


# -------------------------------------------------------
#
#                   Plots
#
# -------------------------------------------------------

# -------------------------
# Plot Beta's - violin
# -------------------------

# Subset beta's to separate dataframes
beta0_df <- beta_df %>% filter(parameter == "beta0") # Beta0
beta1_df <- beta_df %>% filter(parameter == "beta1") # Beta1
beta2_df <- beta_df %>% filter(parameter == "beta2") # Beta2

# Beta0 Plot
beta0_plot <- ggplot(beta0_df, aes(x = Model, y = value, fill = Model)) +
  geom_violin(trim = FALSE) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +  
    labs(y = "Beta0 Estimates", x = "Model") +
  scale_fill_manual(values = c("PC CMR" = "orange",   
                               "PC HDS" = "purple",
                               "Bnet AV" = "blue",
                               "Wolfe AV" = "red")) + 
  scale_y_continuous(limits = c(-15, 15), breaks = seq(-15, 15, by = 5)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),  
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
        panel.grid = element_blank(),
        legend.position = "none")

# View
print(beta0_plot)


# Beta1
beta1_plot <- ggplot(beta1_df, aes(x = Model, y = value, fill = Model)) +
  geom_violin(trim = FALSE) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) + 
  labs(y = "Beta1 Estimates", x = "Model") +
  scale_fill_manual(values = c("PC CMR" = "orange",   
                               "PC HDS" = "purple",
                               "Bnet AV" = "blue",
                               "Wolfe AV" = "red")) + 
  scale_y_continuous(limits = c(-15, 15), breaks = seq(-15, 15, by = 5)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),  
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
        panel.grid = element_blank(),
        legend.position = "none")

# View
print(beta1_plot)


# Beta2
beta2_plot <- ggplot(beta2_df, aes(x = Model, y = value, fill = Model)) +
  geom_violin(trim = FALSE) +  # Violin plot, trim = FALSE keeps the tails of the plot
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +  # Add horizontal line at y = 0
  labs(y = "Beta2 Estimates", x = "Model") +
  scale_fill_manual(values = c("PC CMR" = "orange",   
                               "PC HDS" = "purple",
                               "Bnet AV" = "blue",
                               "Wolfe AV" = "red")) + 
  scale_y_continuous(limits = c(-15, 15), breaks = seq(-15, 15, by = 5)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),  
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
        panel.grid = element_blank(),
        legend.position = "none")

# View
print(beta2_plot)


# Combine & Export 
Beta_plots <- grid.arrange(beta0_plot, beta1_plot, beta2_plot, nrow = 3)
ggsave("./Figures/betas_plots.png", plot = Beta_plots, width = 8, height = 12, dpi = 300)



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
                               "Bnet AV" = "blue",
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
ggsave(plot = densityViolin, "./Figures/densityViolin.jpeg", width = 8, height = 5, dpi = 300) 









# -------------------------
# Plot Abundance - violin
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
                                               "Bnet AV" = "blue",
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








