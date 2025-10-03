# -------------------------------------------------------
#
#                    Load libraries
#
# -------------------------------------------------------
 
# Load packages
library(tidyverse)
library(gridExtra)

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
bnet_dat <- read.csv("./Data/Acoustic_Data/BirdNET_Classifier/FullDeploy_BirdNET.csv")  
nrow(bnet_dat)

# BirdSong detections
bsong_dat <-  read.csv("./Data/Acoustic_Data/BirdSong_Classifier/FullDeploy_BirdSong.csv")  
nrow(bsong_dat)

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
bsong_dat$Count <- 1


# Format date
bnet_dat$Date <- as.Date(bnet_dat$Date)
bsong_dat$Date <- as.Date(bsong_dat$Date)

# Summarize calls per day
bnet_sum <- bnet_dat %>%
  group_by(Date) %>%
  summarize(Total_Calls = sum(Count), .groups = "drop")

bsong_sum <- bsong_dat %>%
  group_by(Date) %>%
  summarize(Total_Calls = sum(Count), .groups = "drop")


# Subset the data to start from May 1st and end on July 31st
bnet_sum <- bnet_sum %>%
  filter(Date >= as.Date("2024-05-1") & Date <= as.Date("2024-07-31"))

bsong_sum <- bsong_sum %>%
  filter(Date >= as.Date("2024-05-1") & Date <= as.Date("2024-07-31"))

# Total calls
sum(bnet_sum$Total_Calls)
sum(bsong_sum$Total_Calls)

bnet_sum$Source <- "BirdNet"
bsong_sum$Source <- "BirdSong"

call_sum <- bind_rows(bsong_sum,bnet_sum)

# Make sure Source is a factor and reorder it
call_sum$Source <- factor(call_sum$Source, levels = c("BirdSong", "BirdNet"))

# -------------------------------------------------
# Histogram
# -------------------------------------------------
callsPerDay <- ggplot(call_sum, aes(x = Date, y = Total_Calls, fill = Source)) +
  geom_col(position = "identity", alpha = 0.6) +
  scale_fill_manual(values = c("BirdSong" = "red",
                               "BirdNet" = "blue")) +
                               
  labs(
    title = "A)",
    x = "Month",
    y = "Vocalizations"
  ) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 850), breaks = seq(0, 850, by = 100))+
  theme(
    legend.position = "none",          
    axis.text.x = element_text(hjust = 0.5),
    axis.text.y = element_text(hjust = 0.5),
    axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
    axis.ticks.x = element_line(size = 1),
    axis.ticks.y = element_line(size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 1),
    plot.title = element_text(face = "bold", size = 16, hjust = 0)   
  ) +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%B"
  )


# View 
print(callsPerDay)

# -------------------------
# Calls by min
# --------------------------

# Make sure your Date_Time is in POSIXct format
bnet_dat$Date_Time <- as.POSIXct(bnet_dat$Date_Time, format = "%Y-%m-%d %H:%M:%S")
bsong_dat$Date_Time <- as.POSIXct(bsong_dat$Date_Time, format = "%Y-%m-%d %H:%M:%S")

# Extract the time at minute level (ignoring seconds)
bnet_min <- bnet_dat %>%
  mutate(Time_Min = format(Date_Time, "%H:%M")) %>%  # Extract hour:minute
  group_by(Time_Min) %>%  # Group by the minute
  summarize(Total_Count = sum(Count), .groups = "drop")  # Summarize counts per minute

bsong_min <- bsong_dat %>%
  mutate(Time_Min = format(Date_Time, "%H:%M")) %>%  # Extract hour:minute
  group_by(Time_Min) %>%  # Group by the minute
  summarize(Total_Count = sum(Count), .groups = "drop")  # Summarize counts per minute


# Convert Time_Min to a proper factor to maintain the order of time
bnet_min$Time_Min <- factor(bnet_min$Time_Min, levels = unique(bnet_min$Time_Min))
bsong_min$Time_Min <- factor(bsong_min$Time_Min, levels = unique(bsong_min$Time_Min))

bnet_min$Source <- "BirdNet"
bsong_min$Source <- "BirdSong"

min_sum <- bind_rows(bsong_min, bnet_min)

# Make sure Source is a factor and reorder it
min_sum$Source <- factor(min_sum$Source, levels = c("BirdSong", "BirdNet"))

# Get all levels
all_times <- levels(min_sum$Time_Min)

# Sequence every 5 mins
breaks_5min <- all_times[seq(1, length(all_times), by = 5)]


# Plot
callsPerMinute <- ggplot(min_sum, aes(x = Time_Min, y = Total_Count, fill = Source)) +
  geom_col(position = "identity", alpha = 0.6) +
  scale_fill_manual(values = c("BirdSong" = "red", "BirdNet" = "blue")) +
  labs(
    title = "B)",
    x = "Minute",
    y = "Vocalizations"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",          
    axis.text.x = element_text(hjust = 0.5),
    axis.text.y = element_text(hjust = 0.5),
    axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
    axis.ticks.x = element_line(size = 1),
    axis.ticks.y = element_line(size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 1),
    plot.title = element_text(face = "bold", size = 16, hjust = 0)    # <<-- align title left
  ) +
  scale_x_discrete(breaks = breaks_5min) +
  scale_y_continuous(
    limits = c(0, 2000),               # adjust if needed!
    breaks = seq(0, 2000, by = 250)
  )


print(callsPerMinute)

# Multipanel figure
grid.arrange(callsPerDay, callsPerMinute, nrow = 2)

# Save to file
jpeg("Figures/Figure_CallHistograms.jpg", width = 10, height = 8, units = "in", res = 300)
grid.arrange(callsPerDay, callsPerMinute, nrow = 2)
dev.off()

# ------------ End Script ------------ 