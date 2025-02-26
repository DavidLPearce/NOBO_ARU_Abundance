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

# read in site data
pc_site_covs <- read.csv("./Data/Point_Count_Data/PointCount_siteCovs.csv")
aru_site_covs <- read.csv("./Data/Acoustic_Data/ARU_siteCovs.csv")

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
beta_df <- beta_df %>% mutate(Model = factor(Model, levels = c("PC CMR", "PC HDS", "AV Bnet", "AV Wolfe")))
dens_df <- dens_df %>% mutate(Model = factor(Model, levels = c("PC CMR", "PC HDS", "AV Bnet", "AV Wolfe")))
dens_sum <- dens_sum %>% mutate(Model = factor(Model, levels = c("PC CMR", "PC HDS", "AV Bnet", "AV Wolfe")))   
abund_sum <- abund_sum %>% mutate(Model = factor(Model, levels = c("PC CMR", "PC HDS", "AV Bnet", "AV Wolfe"))) 


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

# Beta1
beta1_plot <- ggplot(beta1_df, aes(x = Model, y = value, fill = Model)) +
  geom_violin(trim = FALSE) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) + 
  labs(y = "Beta1 Estimates", x = "Model") +
  scale_fill_manual(values = c("PC CMR" = "orange",   
                               "PC HDS" = "purple",
                               "AV Bnet" = "blue",
                               "AV Wolfe" = "red")) + 
  scale_y_continuous(limits = c(-15, 15), breaks = seq(-15, 15, by = 5)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),  
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
        panel.grid = element_blank(),
        legend.position = "none")

# View
print(beta1_plot)
ggsave(plot = beta1_plot, "./Figures/CovEffects/WoodyPRP_beta.jpeg", width = 8, height = 5, dpi = 300) 


# -------------------------------------------------------
# Covariate Effects
# -------------------------------------------------------

# Set covariate name 
Cov_name <- "herb_prp"

# Create a prediction of covariate values
PC_pred_vals <- seq(min(pc_site_covs[, Cov_name]), max(pc_site_covs[, Cov_name]), length.out = 1000)
ARU_pred_vals <- seq(min(aru_site_covs[, Cov_name]), max(aru_site_covs[, Cov_name]), length.out = 1000)


# Subsetting betas
CMR_beta0 <- beta0_df %>% filter(Model == "PC CMR") %>% pull(value) # beta0
HDS_beta0 <- beta0_df %>% filter(Model == "PC HDS") %>% pull(value)
Bnet_beta0 <- beta0_df %>% filter(Model == "AV Bnet") %>% pull(value)
Wolfe_beta0 <- beta0_df %>% filter(Model == "AV Wolfe") %>% pull(value)
CMR_beta1 <- beta1_df %>% filter(Model == "PC CMR") %>% pull(value) # beta1
HDS_beta1 <- beta1_df %>% filter(Model == "PC HDS") %>% pull(value)
Bnet_beta1 <- beta1_df %>% filter(Model == "AV Bnet") %>% pull(value)
Wolfe_beta1 <- beta1_df %>% filter(Model == "AV Wolfe") %>% pull(value)


# Matrices for storing predictions
CMR_preds <- matrix(NA, nrow = length(CMR_beta0), ncol = length(PC_pred_vals))
HDS_preds <- matrix(NA, nrow = length(HDS_beta0), ncol = length(PC_pred_vals))
Bnet_preds <- matrix(NA, nrow = length(Bnet_beta0), ncol = length(ARU_pred_vals))
Wolfe_preds <- matrix(NA, nrow = length(Wolfe_beta0), ncol = length(ARU_pred_vals))


# CMR Predictions
for (i in 1:length(CMR_beta0)) {
  CMR_preds[i,] <- CMR_beta0[i] + CMR_beta1[i] * PC_pred_vals
}

# HDS Predictions
for (i in 1:length(HDS_beta0)) {
  HDS_preds[i, ] <- HDS_beta0[i] + HDS_beta1[i] * PC_pred_vals
}

# AV Bnet Predictions
for (i in 1:length(Bnet_beta0)) {
  Bnet_preds[i, ] <- Bnet_beta0[i] + Bnet_beta1[i] * ARU_pred_vals
}

# AV Wolfe Predictions
for (i in 1:length(Wolfe_beta1)) {
  Wolfe_preds[i, ] <- Wolfe_beta1[i] + Wolfe_beta1[i] * ARU_pred_vals
}


# Calculate credible intervals
CMR_preds_LCI <- apply(CMR_preds, 2, quantile, probs = 0.025) # CMR
CMR_preds_HCI <- apply(CMR_preds, 2, quantile, probs = 0.975)
HDS_preds_LCI <- apply(HDS_preds, 2, quantile, probs = 0.025) # HDS
HDS_preds_HCI <- apply(HDS_preds, 2, quantile, probs = 0.975)
Bnet_preds_LCI <- apply(Bnet_preds, 2, quantile, probs = 0.025) # AV Bnet
Bnet_preds_HCI <- apply(Bnet_preds, 2, quantile, probs = 0.975)
Wolfe_preds_LCI <- apply(Wolfe_preds, 2, quantile, probs = 0.025) # AV Wolfe
Wolfe_preds_HCI <- apply(Wolfe_preds, 2, quantile, probs = 0.975)

# Calculate mean predictions
CMR_preds_mean <- apply(CMR_preds, 2, mean)
HDS_preds_mean <- apply(HDS_preds, 2, mean)
Bnet_preds_mean <- apply(Bnet_preds, 2, mean)
Wolfe_preds_mean <- apply(Wolfe_preds, 2, mean)

# Combine into a single data frame
CMR_preds_df <- data.frame( # CMR
  pred_vals = PC_pred_vals,
  CMR_preds_mean = CMR_preds_mean,
  CMR_preds_LCI = CMR_preds_LCI,
  CMR_preds_HCI = CMR_preds_HCI)


HDS_preds_df <- data.frame( # HDS
  pred_vals = PC_pred_vals,
  HDS_preds_mean = HDS_preds_mean,
  HDS_preds_LCI = HDS_preds_LCI,
  HDS_preds_HCI = HDS_preds_HCI)


Bnet_preds_df <- data.frame( # AV Bnet
  pred_vals = ARU_pred_vals,
  Bnet_preds_mean = Bnet_preds_mean,
  Bnet_preds_LCI = Bnet_preds_LCI,
  Bnet_preds_HCI = Bnet_preds_HCI)


Wolfe_preds_df <- data.frame( # AV Wolfe
  pred_vals = ARU_pred_vals,
  Wolfe_preds_mean = Wolfe_preds_mean,
  Wolfe_preds_LCI = Wolfe_preds_LCI,
  Wolfe_preds_HCI = Wolfe_preds_HCI)





# Plot Effects

# CMR
CMR_preds_plot <- ggplot(CMR_preds_df, aes(x = pred_vals, y = CMR_preds_mean)) +
                          geom_line(color = "black", linewidth = 1.5) +   
                          geom_ribbon(aes(ymin = CMR_preds_LCI, 
                                          ymax = CMR_preds_HCI), 
                                      fill = "orange", alpha = 0.3) +
                          geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
                          labs(x = "", 
                               y = "Effect Estimate", 
                               title = "PC CMR") +
                          scale_x_continuous(limits = c(0.40, 0.9), 
                                             breaks = seq(0.40, 0.9, by = 0.1)) +  
                          scale_y_continuous(limits = c(-9, 2), 
                                             breaks = seq(-9, 2, by = 1)) +  
                          theme_minimal() +
                          theme(panel.grid = element_blank())


# HDS
HDS_preds_plot <- ggplot(HDS_preds_df, aes(x = pred_vals, y = HDS_preds_mean)) +
                          geom_line(color = "black", linewidth = 1.5) +   
                          geom_ribbon(aes(ymin = HDS_preds_LCI, 
                                          ymax = HDS_preds_HCI), 
                                      fill = "purple", alpha = 0.3) + 
                          geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
                          labs(x = "", 
                               y = "", 
                               title = "PC HDS") +
                          scale_x_continuous(limits = c(0.40, 0.9), 
                                             breaks = seq(0.40, 0.9, by = 0.1)) +  
                          scale_y_continuous(limits = c(-9, 2), 
                                             breaks = seq(-9, 2, by = 1)) +  
                          theme_minimal() +
                          theme(panel.grid = element_blank())

# Bnet
Bnet_preds_plot <- ggplot(Bnet_preds_df, aes(x = pred_vals, y = Bnet_preds_mean)) +
                          geom_line(color = "black", linewidth = 1.5) +   
                          geom_ribbon(aes(ymin = Bnet_preds_LCI, 
                                            ymax = Bnet_preds_HCI), 
                                        fill = "blue", alpha = 0.3) + 
                           geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
                           labs(x = "Covariate Value", 
                                 y = "Effect Estimate", 
                                 title = "AV Bnet") +
                           scale_x_continuous(limits = c(0.40, 0.9), 
                                               breaks = seq(0.40, 0.9, by = 0.1)) +  
                           scale_y_continuous(limits = c(-9, 2), 
                                               breaks = seq(-9, 2, by = 1)) +  
                           theme_minimal() +
                           theme(panel.grid = element_blank())
# Wolfe
Wolfe_preds_plot <- ggplot(Wolfe_preds_df, aes(x = pred_vals, y = Wolfe_preds_mean)) +
                          geom_line(color = "black", linewidth = 1.5) +   
                          geom_ribbon(aes(ymin = Wolfe_preds_LCI, 
                                          ymax = Wolfe_preds_HCI), 
                                      fill = "red", alpha = 0.3) + 
                          geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
                          labs(x = "Covariate Value", 
                               y = "", 
                               title = "AV Wolfe") +
                          scale_x_continuous(limits = c(0.40, 0.9), 
                                             breaks = seq(0.40, 0.9, by = 0.1)) +  
                          scale_y_continuous(limits = c(-9, 2), 
                                             breaks = seq(-9, 2, by = 1)) +  
                          theme_minimal() +
                          theme(panel.grid = element_blank())


# Plots in a 2x2 grid
grid.arrange(CMR_preds_plot, HDS_preds_plot, 
             Bnet_preds_plot, Wolfe_preds_plot, 
             ncol = 2, nrow = 2)


# Export                
ggsave(plot = woodycovEff_plot, "Figures/CMR_WoodyNPatches_Effect_plot.jpeg",  
       width = 8, height = 5, dpi = 300) 







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
  scale_y_continuous(limits = c(0, 0.4),
                     breaks = seq(0, 0.4, by = 0.05),
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
ggsave(plot = densityViolin, "./Figures/Abundance&Density/DensityEst_Woodyprp.jpeg", width = 8, height = 5, dpi = 300) 



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
                                               "AV Bnet" = "blue",
                                               "AV Wolfe" = "red")) +  
                  scale_y_continuous(limits = c(0, 800),
                                     breaks = seq(0, 800, by = 50),
                                     labels = scales::comma) +
                  theme_minimal() +
                  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  
                        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),  
                        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
                        panel.grid = element_blank(),
                        legend.position = "none")  

# View
abundViolin
print(abund_sum)

# Export                
ggsave(plot = abundViolin, "./Figures/Abundance&Density/AbundEst_Woodyprp.jpeg",width = 8, height = 5, dpi = 300) 


# End Script