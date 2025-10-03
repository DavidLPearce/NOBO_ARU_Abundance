# Function to check Rhat values
check_rhat <- function(rhat_list, threshold = 1.1) {
  flagged <- list()  # Empty list to store flagged values
  
  for (param in names(rhat_list)) {  # Loop through each parameter
    over_threshold <- which(rhat_list[[param]] > threshold)  # Find indices over threshold
    
    if (length(over_threshold) > 0) {
      flagged[[param]] <- rhat_list[[param]][over_threshold]  # Store flagged values
    }
  }
  
  if (length(flagged) > 0) {
    print("Parameters exceeding threshold:")
    print(flagged)
  } else {
    print("No Rhat values exceed the threshold.")
  }
}