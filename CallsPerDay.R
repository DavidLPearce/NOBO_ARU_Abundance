# -------------------------------------------------------
#
#                    Load libraries
#
# -------------------------------------------------------
 
# Load packages
library(tidyverse)
library(zoo)

# Set seed, scientific notation options, and working directory
set.seed(123)
options(scipen = 9999)
setwd(".")

# -------------------------------------------------------
#
#                       Load Data
#
# -------------------------------------------------------

# BirdNet detections
bnet_dat <- read.csv("./Data/Acoustic_Data/ARU_BirdNET_alldates.csv") # 2024 4/28 - 8/8

# bnet_dat <- read.csv("./Data/Acoustic_Data/ARU_BirdNET_2023.csv") # 2023 6/15 - 7/11 
# bnet_dat <- bnet_dat[which( bnet_dat$Common_name == "Northern Bobwhite"),]

# -------------------------------------------------------
#
#                    Data Wrangling
#
# -------------------------------------------------------


# -------------------------
# Daily Count
# --------------------------

# Adding a count to summarize calls by day
bnet_dat$Count <- 1
print(bnet_dat)

# Format date
bnet_dat$Date <- as.Date(bnet_dat$Date)

# Summarize calls per day
call_sum <- bnet_dat %>%
  group_by(Date) %>%
  summarize(Total_Calls = sum(Count), .groups = "drop")
print(call_sum)



# Plot the barplot
ggplot(call_sum, aes(x = Date, y = Total_Calls)) +
  geom_bar(stat = "identity", fill = "grey", color = "black") +
  # geom_vline(xintercept = as.numeric(as.Date("2024-06-01")), 
  #            linetype = "dashed", color = "black", size = 1) +  # Dotted line for June 1st
  # geom_vline(xintercept = as.numeric(as.Date("2024-07-31")), 
  #            linetype = "dashed", color = "black", size = 1) +  # Dotted line for July 31st
  labs(
    title = "Total Calls Per Day",
    x = "Date",
    y = "Number of Calls"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_date(date_breaks = "1 week", date_labels = "%b %d")

# -------------------------
# Moving window
# --------------------------

# Calculate the 7-day moving average
call_sum <- call_sum %>%
  mutate(Moving_Avg = zoo::rollmean(Total_Calls, k = 7, fill = NA, align = "center")) # Requires zoo package

# Plot 
ggplot(call_sum, aes(x = Date)) +
  geom_line(aes(y = Total_Calls), color = "black", size = 1, alpha = 0.7) +  # Original calls
  geom_line(aes(y = Moving_Avg), color = "red", size = 1.2) +                # 7-day moving average
  labs(
    title = "Total Calls Per Day with 7-Day Moving Average",
    x = "Date",
    y = "Number of Calls"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_date(date_breaks = "1 week", date_labels = "%b %d")


# -------------------------
# Within 30 recording
# --------------------------
