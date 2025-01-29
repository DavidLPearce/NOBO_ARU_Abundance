# -------------------------------------------------------
#
#                    Load libraries
#
# -------------------------------------------------------
 
# Load packages
library(tidyverse)
#library(zoo)

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

# Subset the data to start from May 1st and end on July 31st
call_sum <- call_sum %>%
  filter(Date >= as.Date("2024-05-01") & Date <= as.Date("2024-07-31"))

# Create breaks and labels
breaks <- as.Date(c("2024-05-01", "2024-06-01", "2024-07-01", "2024-07-31"))
labels <- format(breaks, "%B %d")

# Plot 
ggplot(call_sum, aes(x = Date)) +
  geom_line(aes(y = Total_Calls), color = "black", size = 1, alpha = 0.7) +  # Original calls
  # geom_vline(xintercept = as.numeric(as.Date("2024-05-27")),
  #            linetype = "dashed", color = "blue", size = 1) +  # May 27th
  # geom_vline(xintercept = as.numeric(as.Date("2024-06-23")),
  #            linetype = "dashed", color = "blue", size = 1) +  # June 23rd
  labs(
    title = "Calls Per Day",
    x = "Month",
    y = "Calls Per Day"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5),
        axis.text.y = element_text(hjust = 0.5),
        axis.ticks.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black", size = 1)
        ) +
  scale_x_date(
    date_breaks = "1 month",        # Set the breaks to be monthly
    date_labels = "%B"              # Format to show only month abbreviation (e.g., May, Jun)
  )



# -------------------------
# Calls by min
# --------------------------

# Make sure your Date_Time is in POSIXct format
bnet_dat$Date_Time <- as.POSIXct(bnet_dat$Date_Time, format = "%Y-%m-%d %H:%M:%S")

# Extract the time at minute level (ignoring seconds)
bnet_dat_min <- bnet_dat %>%
  mutate(Time_Min = format(Date_Time, "%H:%M")) %>%  # Extract hour:minute
  group_by(Time_Min) %>%  # Group by the minute
  summarize(Total_Count = sum(Count), .groups = "drop")  # Summarize counts per minute

# Convert Time_Min to a proper factor to maintain the order of time
bnet_dat_min$Time_Min <- factor(bnet_dat_min$Time_Min, levels = unique(bnet_dat_min$Time_Min))

# Create a vector of labels for every 5 minutes
label_intervals <- seq(1, nrow(bnet_dat_min), by = 5)  # Get every 5th row
labels <- bnet_dat_min$Time_Min[label_intervals]


# Plot
ggplot(bnet_dat_min, aes(x = Time_Min, y = Total_Count, group = 1)) +
  geom_line(color = "Black", size = 1) +  # Line plot for detections per minute
  scale_x_discrete(
    breaks = bnet_dat_min$Time_Min[label_intervals],  # Show every 5th label
    labels = labels  # Use the custom labels
  ) +
  scale_y_continuous(
    limits = c(0, 350),  # Set y-axis limits from 0 to 350
    breaks = seq(0, 350, by = 50)  # Set y-axis breaks at every 50 units
  ) +
  labs(
    title = "Calls Per Minute",
    x = "Time ",
    y = "Total Calls"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5),  # Rotate x-axis labels
        axis.text.y = element_text(hjust = 0.5),
        axis.ticks.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 1))

