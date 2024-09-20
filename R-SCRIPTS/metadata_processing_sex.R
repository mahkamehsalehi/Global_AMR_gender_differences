library(tidyquant)
library(tidyverse)
library(jsonlite)
library(data.table)
library(progress)

source("R-SCRIPTS/00_functions.R")

# Do not run, slow
# Downlaods Sra metadata file that also has the attributes field and processes it.
#df <- read_csv("DATA/RAW/Sra_metadata_jun12_attributes.txt")

#df <- df %>%
#  mutate(
#    # Parse the JSON column and convert it to a list of data frames
#    jattr_parsed = map(jattr, ~ fromJSON(.x, flatten = TRUE))
#  ) %>%
#  # Unnest the list column to expand it into multiple columns with a prefix to avoid name clashes
#  unnest_wider(jattr_parsed, names_sep = "_jattr")

#write_csv(df, "DATA/PROCESSED/Sra_metadata_processed.csv")

# Load in processed data
df <- read_csv("DATA/PROCESSED/Sra_metadata_processed.csv") 

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

df_selected <- df_selected %>% 
  group_by(bioproject) %>%
  mutate(numeric_age = as.numeric(str_extract(host_age_sam, "\\d+")),
         project_uses_days = any(numeric_age > 100, na.rm = TRUE)) %>%
  ungroup()

df_selected <- df_selected %>%
  mutate(host_age_years = mapply(convert_to_years, host_age_sam, project_uses_days))

# Apply the function and replace values in `host_age_years`
df_selected <- df_selected %>%
  mutate(host_age_years = ifelse(!is.na(host_age_5yr_bin_sam),
                                 sapply(host_age_5yr_bin_sam, calculate_mean_from_range),
                                 host_age_years))

df_selected$host_age_years %>% summary

# Remove gestational week samples
df_selected <- df_selected %>%
  mutate(
    host_age_years = ifelse(grepl("cGA_Weeks", host_age_sam), NA, host_age_years)
  )

# Remove unused columns
age_columns <-  grep("age", names(df_selected), value = TRUE, ignore.case = TRUE) 

age_columns <- age_columns[age_columns != "host_age_years"]

# Remove extra age_columns

df_selected <- df_selected %>% dplyr::select(-c(age_columns, "project_uses_days"))
