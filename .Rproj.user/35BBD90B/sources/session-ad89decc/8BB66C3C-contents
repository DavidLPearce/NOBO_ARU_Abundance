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
cmr_beta_df <-readRDS("./Data/Model_Data/PC_MCR_beta_df.rds")
bsong_beta_df <- readRDS("./Data/Model_Data/ARU_BSong_beta_df.rds")
bnet_beta_df <- readRDS("./Data/Model_Data/ARU_Bnet_beta_df.rds")

# Read in alpha summaries
cmr_alpha_summary <-readRDS("./Figures/PC_MCR_AlphaSummary.rds")
bsong_alpha_summary <- readRDS("./Data/Model_Data/ARU_Bsong_abund_summary.rds")
bnet_alpha_summary <- readRDS("./Data/Model_Data/ARU_Bnet_alpha_summary.rds")

# Read in other parameter (detection & vocal rate) summaries
cmr_other_summary <-readRDS("./Data/Model_Data/PC_MCR_param_summary.rds")
bsong_other_summary <- readRDS("./Data/Model_Data/ARU_Bsong_param_summary.rds")
bnet_other_summary <- readRDS("./Data/Model_Data/ARU_Bnet_param_summary.rds")

# Read in abundance df
cmr_abund_df <-readRDS("./Data/Model_Data/PC_MCR_abund_df.rds")
bsong_abund_df <- readRDS("./Data/Model_Data/ARU_Bsong_abund_df.rds")
bnet_abund_df <- readRDS("./Data/Model_Data/ARU_Bnet_abund_df.rds")


# -------------------------------------------------------
#
#                   Data Wrangling
#
# -------------------------------------------------------

# Combine beta df into one dataframe
beta_df <- rbind(cmr_beta_df, 
                 bsong_beta_df, 
                 bnet_beta_df 
)  


# Combine alpha summaries into one dataframe
alpha_summary <- rbind(cmr_alpha_summary,
                       bsong_alpha_summary,
                       bnet_alpha_summary
) 

# Combine other parameter summaries into one dataframe
other_summary <- rbind(cmr_other_summary,
                       bsong_other_summary,
                       bnet_other_summary
) 

# Combine abundance df into one dataframe
abund_df <- rbind(cmr_abund_df,
                   bsong_abund_df,
                   bnet_abund_df
) 

# Order dataframes for plotting
beta_df <- beta_df %>% mutate(Model = factor(Model, levels = c("PC MCR", "AV Bsong", "AV Bnet"))) 
alpha_summary <- alpha_summary %>% mutate(Model = factor(Model, levels = c("PC MCR",  "AV Bsong", "AV Bnet"))) 
other_summary <- other_summary %>% mutate(Model = factor(Model, levels = c("PC MCR",  "AV Bsong", "AV Bnet"))) 
abund_df <- abund_df %>% mutate(Model = factor(Model, levels = c("PC MCR",  "AV Bsong", "AV Bnet")))

# -------------------------------------------------------
#
#                   Tables
#
# -------------------------------------------------------

alpha_summary

other_summary

write.csv(./Figures/.csv)

# -------------------------------------------------------
#
#                   Plots
#
# -------------------------------------------------------

# -------------------------
# Plot Beta's - violin
# -------------------------

# Subset beta's to separate dataframes
beta0_df <- beta_df %>% filter(parameter == "beta0") # Beta1
beta1_df <- beta_df %>% filter(parameter == "beta1") # Beta1
beta2_df <- beta_df %>% filter(parameter == "beta2") # Beta2

# -------------------------
# Plot Beta 0 
# -------------------------

beta0_plot <- ggplot(beta0_df, aes(x = Model, y = value, fill = Model)) +
  geom_violin(trim = FALSE) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) + 
  labs(title = "Intercept", 
       y = "Estimate", 
       x = "") +
  scale_fill_manual(values = c("PC MCR" = "orange",   
                               "BirdSong AV" = "purple")) + 
  scale_y_continuous(limits = c(-0.5, 1.5), breaks = seq(-0.5, 1.5, by = 0.25)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),# element_blank(), 
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),  
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
        panel.grid = element_blank(),
        legend.position = "none")

# View
print(beta0_plot)

# -------------------------
# Plot Beta 1 
# -------------------------

# Beta1
beta1_plot <- ggplot(beta1_df, aes(x = Model, y = value, fill = Model)) +
  geom_violin(trim = FALSE) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) + 
  labs(title = "Herbaceous", 
       y = "Estimate", 
       x = "") +
  scale_fill_manual(values = c("PC MCR" = "orange",   
                               "BirdSong AV" = "purple")) + 
  scale_y_continuous(limits = c(-0.5, 1.5), breaks = seq(-0.5, 1.5, by = 0.25)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), # element_blank(), 
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),  
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
        panel.grid = element_blank(),
        legend.position = "none")

# View
print(beta1_plot)

# -------------------------
# Plot Beta 2 
# -------------------------

beta2_plot <- ggplot(beta2_df, aes(x = Model, y = value, fill = Model)) +
  geom_violin(trim = FALSE) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) + 
  labs(title = "Woody", 
       y = "Estimate", 
       x = "") +
  scale_fill_manual(values = c("PC MCR" = "orange",   
                               "BirdSong AV" = "purple")) + 
  scale_y_continuous(limits = c(-1, 0.5), breaks = seq(-1, 0.5, by = 0.25)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),  
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
        panel.grid = element_blank(),
        legend.position = "none")

# View
print(beta2_plot)

# -------------------------
# Save
# -------------------------

# Multipanel figure
grid.arrange(beta0_plot, beta1_plot, beta2_plot, nrow = 3)

# Save to file
jpeg("Figures/Figure3.4_BetaEstimates.jpg", width = 8, height = 10, units = "in", res = 300)
grid.arrange(beta0_plot, beta1_plot, beta2_plot, nrow = 3)
dev.off()




# -------------------------
# Plot Abundance - violin
# -------------------------

# Plot violin
abundViolin <- ggplot(abund_df, aes(x = Model, y = Density, fill = Model)) + 
                  geom_violin(trim = FALSE, alpha = 0.6, adjust = 5) +  
                  stat_summary(fun = mean, geom = "point", shape = 20, 
                               size = 3, fill = "black") +   
                  labs(x = "Model", y = "Abundance") +
  scale_fill_manual(values = c("PC MCR" = "orange", 
                               "AV Bsong" = "red",
                               "AV Bnet" = "blue"
                               )) +  
                  scale_y_continuous(limits = c(0, 1200),
                                     breaks = seq(0, 1200, by = 100),
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
ggsave(plot = abundViolin, "./Figures/Figure3.5_AbundEst.jpeg",width = 8, height = 5, dpi = 300) 

write.csv(abund_sum, "./Figures/Table2.2_AbundanceEstimates.csv")

# End Script