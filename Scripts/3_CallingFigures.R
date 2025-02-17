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
bnet_dat <- read.csv("./Data/Acoustic_Data/NOBO_BirdNETall_2024.csv") # 2024 4/28 - 8/8
nrow(bnet_dat)
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
  filter(Date >= as.Date("2024-05-1") & Date <= as.Date("2024-07-31"))

# 14 day subset, every 4 days
call_14day <- call_sum[which(call_sum$Date == "2024-05-26" | call_sum$Date == "2024-05-30" |
                             call_sum$Date == "2024-06-03" | call_sum$Date == "2024-06-07" |
                             call_sum$Date == "2024-06-11" | call_sum$Date == "2024-06-15" |
                             call_sum$Date == "2024-06-19" | call_sum$Date == "2024-06-23" |
                             call_sum$Date == "2024-06-27" | call_sum$Date == "2024-07-01" |
                             call_sum$Date == "2024-07-05" | call_sum$Date == "2024-07-09" |
                             call_sum$Date == "2024-07-13" | call_sum$Date == "2024-07-17" ),]

# 7 day subset, sample from 14 day dates
dates14day <- c("2024-05-26","2024-05-30","2024-06-03","2024-06-07","2024-06-11",
                "2024-06-15", "2024-06-19", "2024-06-23", "2024-06-27", "2024-07-01",
                "2024-07-05", "2024-07-09", "2024-07-13", "2024-07-17") 

dates7day <- sample(dates14day, size = 7, replace = FALSE)
call_7day <- call_sum[call_sum$Date %in% dates7day, ]

# 4 day subset, sample from 14 day dates
dates4day <- sample(dates14day, size = 4, replace = FALSE)
call_4day <- call_sum[call_sum$Date %in% dates4day, ]
  
# Create breaks and labels
# breaks <- as.Date(c("2024-05-01", "2024-06-01", "2024-07-01", "2024-07-31"))
# labels <- format(breaks, "%B %d")

# Line 
# ggplot(call_sum, aes(x = Date)) +
#   geom_line(aes(y = Total_Calls), color = "black", size = 1, alpha = 0.7) +  # Original calls
#   labs(
#     title = "Calls Per Day",
#     x = "Month",
#     y = "Calls Per Day"
#   ) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(hjust = 0.5),
#         axis.text.y = element_text(hjust = 0.5),
#         axis.ticks.x = element_line(size = 1),
#         axis.ticks.y = element_line(size = 1),
#         panel.grid.major = element_blank(),  # Remove major grid lines
#         panel.grid.minor = element_blank(),  # Remove minor grid lines
#         axis.line = element_line(color = "black", size = 1)
#         ) +
#   scale_x_date(
#     date_breaks = "1 month",        # Set the breaks to be monthly
#     date_labels = "%B"              # Format to show only month abbreviation (e.g., May, Jun)
#   )

# -------------------------------------------------
# Histogram
# -------------------------------------------------
callsPerDay <- ggplot(call_sum, aes(x = Date, y = Total_Calls)) +
                geom_col(fill = "black", alpha = 0.7) +  
                labs(
                  title = "Calls Per Day",
                  x = "Date",
                  y = "Number of Calls"
                ) +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 45,hjust = 1),
                      axis.text.y = element_text(hjust = 0.5),
                      axis.ticks.x = element_line(size = 1),
                      axis.ticks.y = element_line(size = 1),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(color = "black", linewidth = 1)
                ) +
                scale_x_date(
                  date_breaks = "3 day",
                  date_labels = "%b%d")

# View 
callsPerDay

# Export                
ggsave(plot = callsPerDay, "Figures/calls_per_day.jpeg",  
       width = 8, height = 5, dpi = 300)

# -------------------------------------------------
# Histogram + 14 day
# -------------------------------------------------
Samp14Day <- ggplot(call_sum, aes(x = Date, y = Total_Calls)) +
                  geom_col(fill = "black", alpha = 0.7) +  
                  # 14-day 
                  geom_vline(data = call_14day, aes(xintercept = as.numeric(Date)), 
                             color = "red", linetype = "dashed", size = 0.75) +
                  labs(
                    title = "Calls Per Day",
                    x = "Date",
                    y = "Number of Calls"
                  ) +
                  theme_minimal() +
                  theme(axis.text.x = element_text(angle = 45,hjust = 1),
                        axis.text.y = element_text(hjust = 0.5),
                        axis.ticks.x = element_line(size = 1),
                        axis.ticks.y = element_line(size = 1),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.line = element_line(color = "black", linewidth = 1)
                  ) +
                  scale_x_date(
                    date_breaks = "3 day",
                    date_labels = "%b%d")

Samp14Day
 
# Export                
ggsave(plot = Samp14Day, "Figures/14daySamp.jpeg", width = 8, height = 5, dpi = 300)
                         


# -------------------------------------------------
# Histogram + 14 survey days + Precip events
# -------------------------------------------------

# Create a data frame for blue lines
blue_dates <- data.frame(Date = as.Date(c("2024-05-17", "2024-05-18", "2024-05-28", "2024-05-29", 
                                          "2024-06-01", "2024-06-11", "2024-06-13", "2024-06-14", 
                                          "2024-06-18", "2024-06-19", "2024-06-20", "2024-06-21", 
                                          "2024-06-22", "2024-06-23", "2024-06-24", "2024-06-25", 
                                          "2024-06-29", "2024-07-08", "2024-07-10", "2024-07-11", 
                                          "2024-07-12", "2024-07-13", "2024-07-15", "2024-07-21", 
                                          "2024-07-22", "2024-07-23", "2024-07-24", "2024-07-25", 
                                          "2024-07-26", "2024-07-27", "2024-07-28", "2024-07-29")))

# Updated plot with blue and red lines
Samp14Dayprecip <- ggplot(call_sum, aes(x = Date, y = Total_Calls)) +
  geom_col(fill = "black", alpha = 0.7) +  
    # Blue solid lines for specific dates
  geom_vline(data = blue_dates, aes(xintercept = as.numeric(Date)), 
             color = "blue", linetype = "solid", size = 0.75) +
  # 14-day red dashed lines
  geom_vline(data = call_14day, aes(xintercept = as.numeric(Date)), 
             color = "red", linetype = "dashed", size = 0.75) +
  labs(
    title = "Calls Per Day",
    x = "Date",
    y = "Number of Calls"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(hjust = 0.5),
        axis.ticks.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 1)
  ) +
  scale_x_date(
    date_breaks = "3 day",
    date_labels = "%b%d"
  )

# Print the plot
print(Samp14Dayprecip)

# Export                
ggsave(plot = Samp14Dayprecip, "Figures/Samp14Dayprecip.jpeg", width = 8, height = 5, dpi = 300)

# -------------------------------------------------
# Histogram + 14 day + 7 day
# -------------------------------------------------
ggplot(call_sum, aes(x = Date, y = Total_Calls)) +
  geom_col(fill = "black", alpha = 0.7) +  
  
  geom_vline(data = call_7day, aes(xintercept = as.numeric(Date)), color = "blue", linetype = "dashed", size = 1) +

  # 14-day 
  geom_vline(data = call_14day, aes(xintercept = as.numeric(Date)), color = "grey", linetype = "solid", size = 1) +
  # 7-day 
  labs(
    title = "Calls Per Day",
    x = "Date",
    y = "Number of Calls"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text.y = element_text(hjust = 0.5),
        axis.ticks.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 1)
  ) +
  scale_x_date(
    date_breaks = "3 day",
    date_labels = "%b%d"
  )  

# -------------------------------------------------
# Histogram + 14 day + 4 day
# -------------------------------------------------
ggplot(call_sum, aes(x = Date, y = Total_Calls)) +
  geom_col(fill = "black", alpha = 0.7) +  
  # 14-day 
  geom_vline(data = call_14day, aes(xintercept = as.numeric(Date)), color = "grey", linetype = "solid", size = 1) +
  # 4-day 
    geom_vline(data = call_4day, aes(xintercept = as.numeric(Date)), color = "red", linetype = "dotted", size = 1) +
  labs(
    title = "Calls Per Day",
    x = "Date",
    y = "Number of Calls"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text.y = element_text(hjust = 0.5),
        axis.ticks.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 1)
  ) +
  scale_x_date(
    date_breaks = "3 day",
    date_labels = "%b%d"
  )  

# -------------------------------------------------
# Histogram + 14 day + 7 day + 4 day
# -------------------------------------------------
ggplot(call_sum, aes(x = Date, y = Total_Calls)) +
  geom_col(fill = "black", alpha = 0.7) +  
  # 14-day 
  geom_vline(data = call_14day, aes(xintercept = as.numeric(Date)), color = "grey", linetype = "solid", size = 1) +
  # 7-day 
  geom_vline(data = call_7day, aes(xintercept = as.numeric(Date)), color = "blue", linetype = "dashed", size = 1) +
  # 4-day 
  geom_vline(data = call_4day, aes(xintercept = as.numeric(Date)), color = "red", linetype = "dotted", size = 1) +
  labs(
    title = "Calls Per Day",
    x = "Date",
    y = "Number of Calls"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text.y = element_text(hjust = 0.5),
        axis.ticks.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 1)
  ) +
  scale_x_date(
    date_breaks = "3 day",
    date_labels = "%b%d"
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

