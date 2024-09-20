library(tidyverse)
######## Function to clean age ########
convert_to_years <- function(age, project_uses_days) {
  # Check if the value is NA
  if (is.na(age)) {
    return(NA)
  }
  # Convert values with "days" to years (case insensitive)
  if (str_detect(age, regex("day", ignore_case = TRUE))) {
    num_days <- as.numeric(str_extract(age, "\\d+"))
    return(num_days / 365)  # Convert days to years
  }
  # Convert values with "months" to years (case insensitive)
  else if (str_detect(age, regex("month", ignore_case = TRUE))) {
    num_months <- as.numeric(str_extract(age, "\\d+"))
    return(num_months / 12)  # Convert months to years
  }
  # Convert values with "years" to years (case insensitive)
  else if (str_detect(age, regex("year", ignore_case = TRUE))) {
    num_years <- as.numeric(str_extract(age, "\\d+"))
    return(num_years / 1)  # Convert
  }
  # Check if the value is a numeric string
  else if (str_detect(age, "^\\d+$")) {
    num_value <- as.numeric(age)
    # Treat all values in the project as days if the project uses days
    if (project_uses_days) {
      return(num_value / 365)  # Convert days to years
    } else {
      return(num_value)  # Assume it's already in years if not using days
    }
  }
  # Return NA if it doesn't match any expected pattern
  else {
    return(NA)
  }
}

########


###########
# Function to calculate the mean from a range
calculate_mean_from_range <- function(range) {
  # Split the range into lower and upper bounds
  bounds <- str_split(range, "-", simplify = TRUE)
  
  # Convert bounds to numeric
  lower_bound <- as.numeric(bounds[1])
  upper_bound <- as.numeric(bounds[2])
  
  # Calculate the mean
  mean_value <- (lower_bound + upper_bound) / 2
  
  return(mean_value)
}
############

