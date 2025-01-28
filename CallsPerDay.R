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

# Subset the data to start from May 1st and end on July 31st
call_sum <- call_sum %>%
  filter(Date >= as.Date("2024-05-01") & Date <= as.Date("2024-07-31"))

# Plot the barplot
ggplot(call_sum, aes(x = Date, y = Total_Calls)) +
  geom_bar(stat = "identity", fill = "grey", color = "black") +
  # geom_vline(xintercept = as.numeric(as.Date("2024-06-01")), 
  #            linetype = "dashed", color = "black", size = 1) +  # Dotted line for June 1st
  # geom_vline(xintercept = as.numeric(as.Date("2024-07-31")), 
  #            linetype = "dashed", color = "black", size = 1) +  # Dotted line for July 31st
  labs(
    title = "",
    x = "Date",
    y = "Number of Calls"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %d")

# -------------------------
# Moving window
# --------------------------



# Calculate the 7-day moving average
# call_sum <- call_sum %>%
#   mutate(Moving_Avg = zoo::rollmean(Total_Calls, k = 7, fill = NA, align = "center")) # Requires zoo package

# Create breaks and labels
breaks <- as.Date(c("2024-05-01", "2024-06-01", "2024-07-01", "2024-07-31"))
labels <- format(breaks, "%B %d")

# Plot 
ggplot(call_sum, aes(x = Date)) +
  geom_line(aes(y = Total_Calls), color = "black", size = 1, alpha = 0.7) +  # Original calls
  # geom_line(aes(y = Moving_Avg), color = "red", size = 1.2) +                # 7-day moving average
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
# Within 30 recording
# --------------------------

# subsets to just counts and time and format time to Hour and Min
bnet_dat_min <- bnet_dat[,c('Time','Count')] %>%
  mutate(
    Time = format(as.POSIXct(Time, format = "%H:%M:%S"), format = "%H:%M")
  )
    
# Summarize the number of calls per minute
calls_per_min <- bnet_dat_min %>%
  group_by(Time) %>%
  summarize(Total_Calls = sum(Count), .groups = "drop")


# Make sure time is formatted as POSIXct and add a dummy date
calls_per_min$Time <- as.POSIXct(paste("2024-01-01", calls_per_min$Time), format = "%Y-%m-%d %H:%M")

# Plot
ggplot(calls_per_min, aes(x = Time, y = Total_Calls, group = 1)) +
  geom_line(color = "black", size = 1) +  # Line plot for calls per minute
  scale_x_datetime(  # Use scale_x_datetime for continuous time data
    breaks = seq(from = min(calls_per_min$Time), to = max(calls_per_min$Time), by = "5 mins"),
    labels = scales::date_format("%H:%M")  # Format to show hours and minutes
  ) + 
  labs(
    title = "Calls Per Minute",
    x = "Time",
    y = "Total Calls"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5),
        axis.text.y = element_text(hjust = 1),
        axis.ticks.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(color = "black", size = 1))
