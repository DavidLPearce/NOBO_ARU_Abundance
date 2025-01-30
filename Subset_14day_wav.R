
source_directory <- "F:/Dawn"

dates14day_folder <- "F:/14daySampPeriod"

dates14day <- c("2024-05-26","2024-05-30","2024-06-03","2024-06-07","2024-06-11",
                "2024-06-15", "2024-06-19", "2024-06-23", "2024-06-27", "2024-07-01",
                "2024-07-05", "2024-07-09", "2024-07-13", "2024-07-17")  
  


# Making a list of all of the files found within the source folder
wavfile_list <- list.files(path = source_directory, recursive = TRUE, full.names = TRUE, pattern='wav$', all.files=TRUE)

#

#wavfile = "F:/Dawn/Batch1/LC_ARU1_20240526_062500.wav" # is 
#wavfile ="F:/Dawn/Batch1/LC_ARU1_20240527_062500.wav" # isnt

# 

for (wavfile in wavfile_list) {
  
  # Splitting file name to sort by survey period
  wavfile_split <- strsplit(basename(wavfile), split = "_")
  
  # Extract the date string (third element)
  date_string <- wavfile_split[[1]][3]
  
  # Convert to Date format
  date_formatted <- as.Date(date_string, format = "%Y%m%d")
  
  # Moving files
  if (date_formatted %in% dates14day){

    # Construct the destination path
    destination_file <- file.path(dates14day_folder, basename(wavfile))
    
    # Copy the file
    file.copy(wavfile, destination_file, overwrite = TRUE)
  
    }
}
