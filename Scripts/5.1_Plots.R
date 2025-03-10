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
cmr_beta_df <-readRDS("./Data/Model_Data/PC_CMR_beta_df.rds")
hds_beta_df <- readRDS("./Data/Model_Data/PC_HDS_beta_df.rds")
AVbnet_beta_df <- readRDS("./Data/Model_Data/ARU_BnetAV_beta_df.rds")
AVwolfet_beta_df <- readRDS("./Data/Model_Data/ARU_WolfeAV_beta_df.rds")

# Read in density dataframes
cmr_Dens_df <-readRDS("./Data/Model_Data/PC_CMR_dens_df.rds")
hds_Dens_df <- readRDS("./Data/Model_Data/PC_HDS_dens_df.rds")
AVbnet_Dens_df <- readRDS("./Data/Model_Data/ARU_BnetAV_dens_df.rds")
AVwolfet_Dens_df <- readRDS("./Data/Model_Data/ARU_WolfeAV__dens_df.rds")

# Read in density summaries
cmr_Dens_sum <-readRDS("./Data/Model_Data/PC_CMR_dens_summary.rds")
hds_Dens_sum <- readRDS("./Data/Model_Data/PC_HDS_dens_summary.rds")
AVbnet_Dens_sum <- readRDS("./Data/Model_Data/ARU_BnetAV_dens_summary.rds")
AVwolfet_Dens_sum <- readRDS("./Data/Model_Data/ARU_WolfeAV_dens_summary.rds")

# Read in abundance summaries
cmr_Abund_sum <-readRDS("./Data/Model_Data/PC_CMR_abund_summary.rds")
hds_Abund_sum <- readRDS("./Data/Model_Data/PC_HDS_abund_summary.rds")
AVbnet_Abund_sum <- readRDS("./Data/Model_Data/ARU_BnetAV_abund_summary.rds")
AVwolfet_Abund_sum <- readRDS("./Data/Model_Data/ARU_WolfeAV_abund_summary.rds")

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
beta2_df <- beta_df %>% filter(parameter == "beta2") # Beta2
# beta3_df <- beta_df %>% filter(parameter == "beta3") # Beta3

# -------------------------
# Plot Beta 0 
# -------------------------

beta0_plot <- ggplot(beta0_df, aes(x = Model, y = value, fill = Model)) +
  geom_violin(trim = FALSE) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) + 
  labs(y = "Beta0 Estimates", x = "Model") +
  scale_fill_manual(values = c("PC CMR" = "orange",   
                               "PC HDS" = "purple",
                               "AV Bnet" = "blue",
                               "AV Wolfe" = "red")) + 
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),  
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
        panel.grid = element_blank(),
        legend.position = "none")

# View
print(beta0_plot)
ggsave(plot = beta0_plot, "./Figures/CovEffects/beta0_plot.jpeg", width = 8, height = 5, dpi = 300) 



# -------------------------
# Plot Beta 1 
# -------------------------

# Beta1
beta1_plot <- ggplot(beta1_df, aes(x = Model, y = value, fill = Model)) +
  geom_violin(trim = FALSE) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) + 
  labs(y = "Beta1 Estimates", x = "Model") +
  scale_fill_manual(values = c("PC CMR" = "orange",   
                               "PC HDS" = "purple",
                               "AV Bnet" = "blue",
                               "AV Wolfe" = "red")) + 
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-1, 1, by = 0.25)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),  
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
        panel.grid = element_blank(),
        legend.position = "none")

# View
print(beta1_plot)
ggsave(plot = beta1_plot, "./Figures/CovEffects/beta1_plot.jpeg", width = 8, height = 5, dpi = 300) 

# -------------------------
# Plot Beta 2 
# -------------------------

beta2_plot <- ggplot(beta2_df, aes(x = Model, y = value, fill = Model)) +
  geom_violin(trim = FALSE) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) + 
  labs(y = "Beta2 Estimates", x = "Model") +
  scale_fill_manual(values = c("PC CMR" = "orange",   
                               "PC HDS" = "purple",
                               "AV Bnet" = "blue",
                               "AV Wolfe" = "red")) + 
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-1, 1, by = 0.25)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_text(face = "bold", margin = margin(t = 10)),  
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
        panel.grid = element_blank(),
        legend.position = "none")

# View
print(beta2_plot)
ggsave(plot = beta2_plot, "./Figures/CovEffects/beta2_plot.jpeg", width = 8, height = 5, dpi = 300) 

# -------------------------
# Plot Beta 3 
# -------------------------

# beta3_plot <- ggplot(beta3_df, aes(x = Model, y = value, fill = Model)) +
#   geom_violin(trim = FALSE) + 
#   geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) + 
#   labs(y = "Beta3 Estimates", x = "Model") +
#   scale_fill_manual(values = c("PC CMR" = "orange",   
#                                "PC HDS" = "purple",
#                                "AV Bnet" = "blue",
#                                "AV Wolfe" = "red")) + 
#   scale_y_continuous(limits = c(-2, 2), breaks = seq(-1, 1, by = 0.25)) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), 
#         axis.title.x = element_text(face = "bold", margin = margin(t = 10)),  
#         axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
#         panel.grid = element_blank(),
#         legend.position = "none")
# 
# # View
# print(beta3_plot)
# ggsave(plot = beta3_plot, "./Figures/CovEffects/beta3_plot.jpeg", width = 8, height = 5, dpi = 300) 
# 


# -------------------------------------------------------
# Covariate Effects
# -------------------------------------------------------

# Set covariate name 
Cov1_name <- "herb_ClmIdx"
Cov2_name <- "woody_AggInx"

# Create a prediction of covariate values
PC_cov1 <- seq(min(pc_site_covs[, Cov1_name]), max(pc_site_covs[, Cov1_name]), length.out = 1000)
PC_cov2 <- seq(min(pc_site_covs[, Cov2_name]), max(pc_site_covs[, Cov2_name]), length.out = 1000)

ARU_cov1<- seq(min(aru_site_covs[, Cov1_name]), max(aru_site_covs[, Cov1_name]), length.out = 100)
ARU_cov2 <- seq(min(aru_site_covs[, Cov2_name]), max(aru_site_covs[, Cov2_name]), length.out = 100)

# Mean scaling covariates
PC_cov1 <- (PC_cov1 - mean(PC_cov1)) / (max(PC_cov1) - min(PC_cov1))
PC_cov2 <- (PC_cov2 - mean(PC_cov2)) / (max(PC_cov2) - min(PC_cov2))

ARU_cov1 <- (ARU_cov1 - mean(ARU_cov1)) / (max(ARU_cov1) - min(ARU_cov1))
ARU_cov1 <- (ARU_cov1 - mean(ARU_cov1)) / (max(ARU_cov1) - min(ARU_cov1))

# Subsetting betas
CMR_beta0 <- beta0_df %>% filter(Model == "PC CMR") %>% pull(value) # beta0
HDS_beta0 <- beta0_df %>% filter(Model == "PC HDS") %>% pull(value)
Bnet_beta0 <- beta0_df %>% filter(Model == "AV Bnet") %>% pull(value)
Wolfe_beta0 <- beta0_df %>% filter(Model == "AV Wolfe") %>% pull(value)
CMR_beta1 <- beta1_df %>% filter(Model == "PC CMR") %>% pull(value) # beta1
HDS_beta1 <- beta1_df %>% filter(Model == "PC HDS") %>% pull(value)
Bnet_beta1 <- beta1_df %>% filter(Model == "AV Bnet") %>% pull(value)
Wolfe_beta1 <- beta1_df %>% filter(Model == "AV Wolfe") %>% pull(value)
CMR_beta2 <- beta2_df %>% filter(Model == "PC CMR") %>% pull(value) # beta2
HDS_beta2 <- beta2_df %>% filter(Model == "PC HDS") %>% pull(value)
Bnet_beta2 <- beta2_df %>% filter(Model == "AV Bnet") %>% pull(value)
Wolfe_beta2 <- beta2_df %>% filter(Model == "AV Wolfe") %>% pull(value)
# CMR_beta3 <- beta3_df %>% filter(Model == "PC CMR") %>% pull(value) # beta3
# HDS_beta3 <- beta3_df %>% filter(Model == "PC HDS") %>% pull(value)
# Bnet_beta3 <- beta3_df %>% filter(Model == "AV Bnet") %>% pull(value)
# Wolfe_beta3 <- beta3_df %>% filter(Model == "AV Wolfe") %>% pull(value)

# Matrices for storing predictions
CMR_cov1_preds <- matrix(NA, nrow = length(CMR_beta0), ncol = length(PC_cov1)) # Linear
CMR_cov2_preds <- matrix(NA, nrow = length(CMR_beta0), ncol = length(PC_cov1)) 

HDS_cov1_preds <- matrix(NA, nrow = length(HDS_beta0), ncol = length(PC_cov1))
HDS_cov2_preds <- matrix(NA, nrow = length(HDS_beta0), ncol = length(PC_cov1))

Bnet_cov1_preds <- matrix(NA, nrow = length(Bnet_beta0), ncol = length(ARU_cov1))
Bnet_cov2_preds <- matrix(NA, nrow = length(Bnet_beta0), ncol = length(ARU_cov1))

Wolfe_cov1_preds <- matrix(NA, nrow = length(Wolfe_beta0), ncol = length(ARU_cov1))
Wolfe_cov2_preds <- matrix(NA, nrow = length(Wolfe_beta0), ncol = length(ARU_cov1))


# CMR Predictions
for (i in 1:length(CMR_beta0)) {
  
  # Linear
  CMR_cov1_preds[i,] <- CMR_beta0[i] + CMR_beta1[i] * PC_cov1
  CMR_cov2_preds[i,] <- CMR_beta0[i] + CMR_beta2[i] * PC_cov2
}


# HDS Predictions
for (i in 1:length(HDS_beta0)) {
  
  # Linear
  HDS_cov1_preds[i,] <- HDS_beta0[i] + HDS_beta1[i] * PC_cov1
  HDS_cov2_preds[i,] <- HDS_beta0[i] + HDS_beta2[i] * PC_cov2
}

# Bnet Predictions
for (i in 1:length(Bnet_beta0)) {
  
  # Linear
  Bnet_cov1_preds[i,] <- Bnet_beta0[i] + Bnet_beta1[i] * ARU_cov1
  Bnet_cov2_preds[i,] <- Bnet_beta0[i] + Bnet_beta2[i] * ARU_cov2
}

# Wolfe Predictions
for (i in 1:length(Wolfe_beta0)) {
  
  # Linear
  Wolfe_cov1_preds[i,] <- Wolfe_beta0[i] + Wolfe_beta1[i] * ARU_cov1
  Wolfe_cov2_preds[i,] <- Wolfe_beta0[i] + Wolfe_beta2[i] * ARU_cov2
 } 

# Calculate credible intervals

# CMR
CMR_cov1_preds_mean <- apply(CMR_cov1_preds, 2, mean)# mean
CMR_cov1_preds_LCI <- apply(CMR_cov1_preds, 2, quantile, probs = 0.025) # LCI
CMR_cov1_preds_HCI <- apply(CMR_cov1_preds, 2, quantile, probs = 0.975) # HCI

CMR_cov2_preds_mean <- apply(CMR_cov2_preds, 2, mean) 
CMR_cov2_preds_LCI <- apply(CMR_cov2_preds, 2, quantile, probs = 0.025) 
CMR_cov2_preds_HCI <- apply(CMR_cov2_preds, 2, quantile, probs = 0.975) 


# HDS
HDS_cov1_preds_mean <- apply(HDS_cov1_preds, 2, mean)# mean
HDS_cov1_preds_LCI <- apply(HDS_cov1_preds, 2, quantile, probs = 0.025) # LCI
HDS_cov1_preds_HCI <- apply(HDS_cov1_preds, 2, quantile, probs = 0.975) # HCI

HDS_cov2_preds_mean <- apply(HDS_cov2_preds, 2, mean) 
HDS_cov2_preds_LCI <- apply(HDS_cov2_preds, 2, quantile, probs = 0.025) 
HDS_cov2_preds_HCI <- apply(HDS_cov2_preds, 2, quantile, probs = 0.975) 

# Bnet
Bnet_cov1_preds_mean <- apply(Bnet_cov1_preds, 2, mean)# mean
Bnet_cov1_preds_LCI <- apply(Bnet_cov1_preds, 2, quantile, probs = 0.025) # LCI
Bnet_cov1_preds_HCI <- apply(Bnet_cov1_preds, 2, quantile, probs = 0.975) # HCI

Bnet_cov2_preds_mean <- apply(Bnet_cov2_preds, 2, mean) 
Bnet_cov2_preds_LCI <- apply(Bnet_cov2_preds, 2, quantile, probs = 0.025) 
Bnet_cov2_preds_HCI <- apply(Bnet_cov2_preds, 2, quantile, probs = 0.975) 

# Wolfe
Wolfe_cov1_preds_mean <- apply(Wolfe_cov1_preds, 2, mean)# mean
Wolfe_cov1_preds_LCI <- apply(Wolfe_cov1_preds, 2, quantile, probs = 0.025) # LCI
Wolfe_cov1_preds_HCI <- apply(Wolfe_cov1_preds, 2, quantile, probs = 0.975) # HCI

Wolfe_cov2_preds_mean <- apply(Wolfe_cov2_preds, 2, mean) 
Wolfe_cov2_preds_LCI <- apply(Wolfe_cov2_preds, 2, quantile, probs = 0.025) 
Wolfe_cov2_preds_HCI <- apply(Wolfe_cov2_preds, 2, quantile, probs = 0.975) 


# Combine into a single data frame

# CMR
CMR_cov1_df <- data.frame(
  cov1 = PC_cov1,
  cov1_preds_mean = CMR_cov1_preds_mean,
  cov1_preds_LCI = CMR_cov1_preds_LCI,
  cov1_preds_HCI = CMR_cov1_preds_HCI)

CMR_cov2_df <- data.frame(
  cov2 = PC_cov2,
  cov2_preds_mean = CMR_cov2_preds_mean,
  cov2_preds_LCI = CMR_cov2_preds_LCI,
  cov2_preds_HCI = CMR_cov2_preds_HCI)


# HDS
HDS_cov1_df <- data.frame(
  cov1 = PC_cov1,
  cov1_preds_mean = HDS_cov1_preds_mean,
  cov1_preds_LCI = HDS_cov1_preds_LCI,
  cov1_preds_HCI = HDS_cov1_preds_HCI)

HDS_cov2_df <- data.frame(
  cov2 = PC_cov2,
  cov2_preds_mean = HDS_cov2_preds_mean,
  cov2_preds_LCI = HDS_cov2_preds_LCI,
  cov2_preds_HCI = HDS_cov2_preds_HCI)


# Bnet
Bnet_cov1_df <- data.frame(
  cov1 = ARU_cov1,
  cov1_preds_mean = Bnet_cov1_preds_mean,
  cov1_preds_LCI = Bnet_cov1_preds_LCI,
  cov1_preds_HCI = Bnet_cov1_preds_HCI)

Bnet_cov2_df <- data.frame(
  cov2 = ARU_cov2,
  cov2_preds_mean = Bnet_cov2_preds_mean,
  cov2_preds_LCI = Bnet_cov2_preds_LCI,
  cov2_preds_HCI = Bnet_cov2_preds_HCI)


# Wolfe
Wolfe_cov1_df <- data.frame(
  cov1 = ARU_cov1,
  cov1_preds_mean = Wolfe_cov1_preds_mean,
  cov1_preds_LCI = Wolfe_cov1_preds_LCI,
  cov1_preds_HCI = Wolfe_cov1_preds_HCI)

Wolfe_cov2_df <- data.frame(
  cov2 = ARU_cov2,
  cov2_preds_mean = Wolfe_cov2_preds_mean,
  cov2_preds_LCI = Wolfe_cov2_preds_LCI,
  cov2_preds_HCI = Wolfe_cov2_preds_HCI)


# -------------------------
# Plot Covariate 1 Effects
# -------------------------

# CMR
CMR_cov1_plot <- ggplot(CMR_cov1_df, aes(x = cov1, y = cov1_preds_mean)) +
                          geom_line(color = "black", linewidth = 1.5) +   
                          geom_ribbon(aes(ymin = cov1_preds_LCI, 
                                          ymax = cov1_preds_HCI), 
                                      fill = "orange", alpha = 0.3) +
                          geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
                          labs(x = "", 
                               y = "Predicted Effect", 
                               title = "PC CMR") +
                          scale_x_continuous(limits = c(-0.5, 0.5),
                                             breaks = seq(-0.5, 0.5, by = 0.1)) +
                          scale_y_continuous(limits = c(-2.5, 2.5),
                                             breaks = seq(-2.5, 2.5, by = 0.5)) +
                          theme_minimal() +
                          theme(panel.grid = element_blank(),
                                axis.text.x = element_blank())


# HDS
HDS_cov1_plot <- ggplot(HDS_cov1_df, aes(x = cov1, y = cov1_preds_mean)) +
                    geom_line(color = "black", linewidth = 1.5) +   
                    geom_ribbon(aes(ymin = cov1_preds_LCI, 
                                    ymax = cov1_preds_HCI), 
                                fill = "purple", alpha = 0.3) +
                    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
                    labs(x = "", 
                         y = "", 
                         title = "PC HDS") +
                    scale_x_continuous(limits = c(-0.5, 0.5),
                                       breaks = seq(-0.5, 0.5, by = 0.1)) +
                    scale_y_continuous(limits = c(-2.5, 2.5),
                                       breaks = seq(-2.5, 2.5, by = 0.5)) +
                    theme_minimal()+
                      theme(panel.grid = element_blank(),
                            axis.text.x = element_blank(), 
                            axis.text.y = element_blank())  

# Bnet
Bnet_cov1_plot <- ggplot(Bnet_cov1_df, aes(x = cov1, y = cov1_preds_mean)) +
                      geom_line(color = "black", linewidth = 1.5) +   
                      geom_ribbon(aes(ymin = cov1_preds_LCI, 
                                      ymax = cov1_preds_HCI), 
                                  fill = "blue", alpha = 0.3) +
                      geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
                      labs(x = "Herbaceous Clumpy Index (scaled)", 
                           y = "Predicted Effect", 
                           title = "AV Bnet") +
                      scale_x_continuous(limits = c(-0.5, 0.5),
                                         breaks = seq(-0.5, 0.5, by = 0.1)) +
                      scale_y_continuous(limits = c(-2.5, 2.5),
                                         breaks = seq(-2.5, 2.5, by = 0.5)) +
                      theme_minimal() +
                      theme(panel.grid = element_blank())


# Wolfe
Wolfe_cov1_plot <- ggplot(Wolfe_cov1_df, aes(x = cov1, y = cov1_preds_mean)) +
                      geom_line(color = "black", linewidth = 1.5) +   
                      geom_ribbon(aes(ymin = cov1_preds_LCI, 
                                      ymax = cov1_preds_HCI), 
                                  fill = "red", alpha = 0.3) +
                      geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
                      labs(x = "Herbaceous Clumpy Index (scaled)", 
                           y = "", 
                           title = "AV Wolfe") +
                      scale_x_continuous(limits = c(-0.5, 0.5),
                                         breaks = seq(-0.5, 0.5, by = 0.1)) +
                      scale_y_continuous(limits = c(-2.5, 2.5),
                                         breaks = seq(-2.5, 2.5, by = 0.5)) +
                      theme_minimal() +
                      theme(panel.grid = element_blank(), 
                            axis.text.y = element_blank())


# Plots in a 2x2 grid
Cov1_plot<- grid.arrange(CMR_cov1_plot, HDS_cov1_plot, 
             Bnet_cov1_plot, Wolfe_cov1_plot, 
             ncol = 2, nrow = 2)

# View
Cov1_plot

# Export                
ggsave(plot = Cov1_plot, "./Figures/CovEffects/Cov1_plot.jpeg", width = 8, height = 5, dpi = 300)   
       

# -------------------------
# Plot Covariate 2 Effects
# -------------------------


# CMR
CMR_cov2_plot <- ggplot(CMR_cov2_df, aes(x = cov2, y = cov2_preds_mean)) +
                    geom_line(color = "black", linewidth = 1.5) +   
                    geom_ribbon(aes(ymin = cov2_preds_LCI, 
                                    ymax = cov2_preds_HCI), 
                                fill = "orange", alpha = 0.3) +
                    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
                    labs(x = "", 
                         y = "Predicted Effect", 
                         title = "PC CMR") +
                    scale_x_continuous(limits = c(-3, 3),
                                       breaks = seq(-2, 2, by = 0.5)) +
                    scale_y_continuous(limits = c(-3,3),
                                       breaks = seq(-3, 3, by = 0.5)) +
                    theme_minimal() +
                    theme(panel.grid = element_blank(),
                          axis.text.x = element_blank())


# HDS
HDS_cov2_plot <- ggplot(HDS_cov2_df, aes(x = cov2, y = cov2_preds_mean)) +
                    geom_line(color = "black", linewidth = 1.5) +   
                    geom_ribbon(aes(ymin = cov2_preds_LCI, 
                                    ymax = cov2_preds_HCI), 
                                fill = "purple", alpha = 0.3) +
                    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
                    labs(x = "", 
                         y = "", 
                         title = "PC HDS") +
                    scale_x_continuous(limits = c(-3, 3),
                                       breaks = seq(-2, 2, by = 0.5)) +
                    scale_y_continuous(limits = c(-3,3),
                                       breaks = seq(-3, 3, by = 0.5)) +
                    theme_minimal()+
                    theme(panel.grid = element_blank(),
                          axis.text.x = element_blank(), 
                          axis.text.y = element_blank())  

# Bnet
Bnet_cov2_plot <- ggplot(Bnet_cov2_df, aes(x = cov2, y = cov2_preds_mean)) +
                    geom_line(color = "black", linewidth = 1.5) +   
                    geom_ribbon(aes(ymin = cov2_preds_LCI, 
                                    ymax = cov2_preds_HCI), 
                                fill = "blue", alpha = 0.3) +
                    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
                    labs(x = "Woody Aggregation (scaled)", 
                         y = "Predicted Effect", 
                         title = "AV Bnet") +
                    scale_x_continuous(limits = c(-3, 3),
                                       breaks = seq(-2, 2, by = 0.5)) +
                    scale_y_continuous(limits = c(-3,3),
                                       breaks = seq(-3, 3, by = 0.5)) +
                    theme_minimal() +
                    theme(panel.grid = element_blank())


# Wolfe
Wolfe_cov2_plot <- ggplot(Wolfe_cov2_df, aes(x = cov2, y = cov2_preds_mean)) +
                      geom_line(color = "black", linewidth = 1.5) +   
                      geom_ribbon(aes(ymin = cov2_preds_LCI, 
                                      ymax = cov2_preds_HCI), 
                                  fill = "red", alpha = 0.3) +
                      geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
                      labs(x = "Woody Aggregation (scaled)", 
                           y = "", 
                           title = "AV Wolfe") +
                      scale_x_continuous(limits = c(-3, 3),
                                         breaks = seq(-2, 2, by = 0.5)) +
                      scale_y_continuous(limits = c(-3,3),
                                         breaks = seq(-3, 3, by = 0.5)) +
                      theme_minimal() +
                      theme(panel.grid = element_blank(), 
                            axis.text.y = element_blank())


# Plots in a 2x2 grid
cov2_plot<- grid.arrange(CMR_cov2_plot, HDS_cov2_plot, 
                         Bnet_cov2_plot, Wolfe_cov2_plot, 
                         ncol = 2, nrow = 2)

# View
cov2_plot

# Export                
ggsave(plot = cov2_plot, "./Figures/CovEffects/cov2_plot.jpeg", width = 8, height = 5, dpi = 300)   








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
                  scale_y_continuous(limits = c(300, 1200),
                                     breaks = seq(300, 1200, by = 50),
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