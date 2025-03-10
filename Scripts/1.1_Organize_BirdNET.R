# Load required packages
library(tidyverse)


# Create a blank dataframe called 'acoustic_data'
acoustic_data <- data.frame(
  Start..s. = numeric(),
  End..s. = numeric(),
  Scientific.name = character(),
  Common.name = character(),
  Confidence = numeric(),
  Study_Area = character(),
  SiteID = character(),
  Site_Number = character(),
  Recording_Time = character(),
  Date_Time = character(),
  Date = character(),
  Time = character(),
  File_Name = character(),
  stringsAsFactors = FALSE  # Ensure character columns are treated as strings
)


# Directory containing the CSV files
csv_directory <- "D:/LaCopita_Acoustics/Acoustics_Summer2024/Summer24_Dawn_BirdNET_RawResults"

# List all CSV files in the directory
csv_files <- list.files(csv_directory, pattern = "\\.csv$", full.names = TRUE)

# Initialize progress bar
total_files <- length(csv_files)
pb <- txtProgressBar(min = 0, max = total_files, style = 3)


# Adding study area, site, date, survey period to columns from file name
for (i in seq_along(csv_files)) {
  
  # Subset file
  file <- csv_files[i]
  
  # Extracting information from the file name
  file_name <- basename(file)
  file_name_parts <- unlist(strsplit(file_name, "[._]"))
  study_area <- file_name_parts[1]
  siteID <- file_name_parts[2]
  site_number <- as.numeric(gsub("\\D", "", siteID))
  date <- file_name_parts[3]
  sampling_time <- file_name_parts[4]
  record_time <- substr(sampling_time, 1, nchar(sampling_time) - 2)
  datetime <- as.POSIXct(paste(date, record_time), format = "%Y%m%d %H%M")
  date <- as.POSIXct(paste(date), format = "%Y%m%d")
  
  
  # Read CSV file
  site_data <- read.csv(file)
  
  # Add columns for file name parts
  site_data <- mutate(site_data, 
                      Study_Area = study_area,
                      SiteID = siteID,
                      Site_Number = site_number,
                      Recording_Time = record_time,
                      Date_Time = datetime + Start..s.,
                      Date = date,
                      Time = format(Date_Time, "%H:%M:%S"),
                      File_Name = file_name
  )
  
  
  
  
  # Combining into one dataset
  acoustic_data <- data.frame(rbind(site_data, acoustic_data))  
  
  # Update progress bar
  setTxtProgressBar(pb, i)

}# end loop


# Take a look
View(acoustic_data)

# Renaming Columns
colnames(acoustic_data) <- c("Detection_Start", "Detection_End", 
                             "Scientific_name" , "Common_name",
                             "Confidence","Study_Area", "SiteID", "Site_Number",
                             "Recording_Time", "Date_Time",
                             "Date", "Time", "File_Name")

# Reorganizing columns
acoustic_data <- acoustic_data[,c("File_Name", "Study_Area", "SiteID", "Site_Number",
                                  "Recording_Time","Date_Time", 
                                  "Date", "Time",
                                  "Detection_Start", "Detection_End",
                                  "Scientific_name" , "Common_name",
                                  "Confidence")]

# Convert 'Date' column to Date format
acoustic_data$Date <- ymd(acoustic_data$Date)

# Sort the dataframe by the 'Site', 'Date', and 'Sampling_Time' columns
acoustic_data <- acoustic_data[order(acoustic_data$Site_Number, acoustic_data$Date, acoustic_data$Recording_Time), ]

# Subsetting to only NOBO 
acoustic_data <- acoustic_data[which(acoustic_data$Common_name == "Northern Bobwhite" ),]

# Subsetting to only >= 0.80 confidence 
acoustic_data <- acoustic_data[which(acoustic_data$Confidence >= 0.80),]


# Take a look
View(acoustic_data)


# Export to csv
write.csv(acoustic_data, row.names = FALSE, "./Data/Raw_Data/Acoustic_Data_Bnet_2024.csv")
          
# Save an object to a file
saveRDS(acoustic_data, file = "./Data/Raw_Data/Acoustic_Data_Bnet_2024.rds")

