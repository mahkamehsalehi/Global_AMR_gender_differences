wd <-
  "/Users/kaparn/Documents/GlobalAMR_project/"

setwd(wd)


library(tidyquant)
library(tidyverse)
library(jsonlite)
library(data.table)
library(progress)


#df <- read_csv("Sra_metadata_jun12_attributes.txt")

#df <- df %>%
#  mutate(
#    # Parse the JSON column and convert it to a list of data frames
#    jattr_parsed = map(jattr, ~ fromJSON(.x, flatten = TRUE))
#  ) %>%
#  # Unnest the list column to expand it into multiple columns with a prefix to avoid name clashes
#  unnest_wider(jattr_parsed, names_sep = "_jattr")

#write_csv(df, "Sra_metadata_processed.csv")

df <- read_csv("Sra_metadata_processed.csv") 

sex_columns <- grep("sex", names(df), value = TRUE, ignore.case = TRUE)
gender_columns <- grep("gender", names(df), value = TRUE, ignore.case = TRUE)
age_columns <-  grep("host_age", names(df), value = TRUE, ignore.case = TRUE) 
age_columns <- c(age_columns, "jattr_parsed_jattrraw_age_sam")
age_columns <- age_columns[!age_columns=="host_age__months__sam_s_dpl228"]
# Print the results
cat("Columns containing 'sex':\n")
print(sex_columns)

cat("\nColumns containing 'gender':\n")
print(gender_columns)

cat("\nColumns containing 'age':\n")
print(age_columns)



# Create a new dataframe with the specified columns
df_selected <- df %>%
  select(
    acc,
    geo_loc_name_country_calc,
    geo_loc_name_country_continent_calc,
    platform,
    instrument,
    bioproject,
    avgspotlen,
    mbases,
    collection_date_sam,
    matches(c(age_columns, sex_columns, gender_columns)))

  
  # Rename columns by removing the "jattr_parsed_jattr" prefix
  df_selected <- df_selected %>%
    rename_with(
      .fn = ~ gsub("^jattr_parsed_jattr", "", .x),
      .cols = c(sex_columns, gender_columns, age_columns)
    )
  
  

# Print the names of the columns in the new dataframe to verify
print(names(df_selected))

# Clean sex & age columns
df_selected <- df_selected %>%
  mutate(sex_combined = coalesce(host_sex_sam, gender_sam, sex_calc, infant_gender_sam)) %>%
  mutate(sex_combined = tolower(sex_combined),  # Convert to lowercase
         sex_combined = ifelse(sex_combined %in% c("male", "female"), sex_combined, NA)) %>% # Replace invalid values with NA
  select(-host_sex_sam, -gender_sam, -sex_calc, -infant_gender_sam) %>%
mutate(host_age_sam = str_replace_all(host_age_sam, "around ", "")) %>%
  mutate(host_age_sam = str_replace_all(host_age_sam, " years", ""))   

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



df_selected <- df_selected %>% 
  group_by(bioproject) %>%
  mutate(numeric_age = as.numeric(str_extract(host_age_sam, "\\d+")),
         project_uses_days = any(numeric_age > 100, na.rm = TRUE)) %>%
  ungroup()

df_selected <- df_selected %>%
  mutate(host_age_years = mapply(convert_to_years, host_age_sam, project_uses_days))

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

# Apply the function and replace values in `host_age_years`
df_selected <- df_selected %>%
  mutate(host_age_years = ifelse(!is.na(host_age_5yr_bin_sam),
                                 sapply(host_age_5yr_bin_sam, calculate_mean_from_range),
                                 host_age_years))

df_selected$host_age_years %>% summary

# Remove gestational week samples
# TODO # Remove gestational week samples
# "cGA_Weeks"