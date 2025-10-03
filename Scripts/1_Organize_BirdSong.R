# Load required packages
library(tidyverse)

# Load AWS data
aws_data <- read.csv("./NOBO_AWS_Results.csv")

# Adding study area, site, date, survey period to columns from file name
for (row in seq_len(nrow(aws_data))) {
  
  # subset row
  row_sub <- aws_data[row, ]
  
  # Subset file 
  file <- row_sub[,1]
  
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
  
  # Add into 
  aws_data[row, 'Study_Area'] <- study_area
  aws_data[row, 'SiteID'] <- siteID
  aws_data[row, 'Site_Number'] <- site_number
  aws_data[row, 'Recording_Time'] <- record_time
  aws_data[row, 'Date_Time'] <- datetime + aws_data[row, 'NB_whistle_mean_time'] 
  aws_data[row, 'Date'] <- date  
  aws_data[row, 'Time'] <- format(aws_data[row, 'Date_Time'], "%H:%M:%S") 
 
}# end loop


# Take a look
head(aws_data)



# Sort the dataframe by the 'Site', 'Date', and 'Sampling_Time' columns
aws_data <- aws_data[order(aws_data$Site_Number, aws_data$Date, aws_data$Recording_Time), ]


# Take a look
View(aws_data)


# Export to csv
write.csv(aws_data, row.names = FALSE, "./CH1_AWS.csv")

# Save an object to a file
saveRDS(aws_data, file = "./CH1_AWS.rds")





















