wd <-
  "/scratch/project_2008149/USER_WORKSPACES/mahkameh/women_amr/"

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
#age_columns <- grep("age", names(df), value = TRUE, ignore.case = TRUE)
read_columns <- grep("read", names(df), value = TRUE, ignore.case = TRUE)

# Print the results
cat("Columns containing 'sex':\n")
print(sex_columns)

cat("\nColumns containing 'gender':\n")
print(gender_columns)

#cat("\nColumns containing 'age':\n")
#print(age_columns)

cat("\nColumns containing 'read':\n")
print(read_columns)

# Rename columns by removing the "jattr_parsed_jattr" prefix
df <- df %>%
  rename_with(
    .fn = ~ gsub("^jattr_parsed_jattr", "", .x),
    .cols = c(sex_columns, gender_columns)
  )

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
    host_sex_sam,
    sex_calc
  ) %>%
  rename(country = geo_loc_name_country_calc) %>%
  mutate(country = case_when(
    country == "USA" ~ "United States",
    country == "South Korea" ~ "South Korea",
    country == "Democratic Republic of the Congo" ~ "Congo",
    country == "The Bahamas" ~ "Bahamas",
    country == "The Gambia" ~ "Gambia",
    TRUE ~ country 
  ))

# Print the names of the columns in the new dataframe to verify
print(names(df_selected))

data_2016 <- read_csv("total_antibiotic_consumption_estimates.csv") %>%
  subset(Year ==2016) %>%
  mutate(Location = case_when(
    Location == "United States" ~ "USA",
    Location == "South Korea" ~ "South Korea",
    Location == "Congo" ~ "Democratic Republic of the Congo",
    Location == "Bahamas" ~ "The Bahamas",
    Location == "Gambia" ~ "The Gambia",
    TRUE ~ Location
  )) %>%
  rename(country = Location) %>%
  rename(antibiotic_consumption = `Antibiotic consumption (DDD/1,000/day)`)

data_2016_subset <- data_2016[, c("country", "antibiotic_consumption")]




global_data_2016 <- merge(df_selected, data_2016_subset, by = "country")
write.csv(global_data_2016, "global_data_2016.csv")



g_data_2016 <- read_csv("global_data_2016.csv")

female_data <- subset(g_data_2016, host_sex_sam == "female" | sex_calc == "female") %>%
  select(-sex_calc, -host_sex_sam, -collection_date_sam, -...1) %>%
  rename(region = geo_loc_name_country_continent_calc)

both_gender_data <- subset(g_data_2016, host_sex_sam == "female" | host_sex_sam == "male"| sex_calc == "female" | sex_calc == "male") %>%
  select(-sex_calc, -host_sex_sam, -collection_date_sam, -...1) %>%
  rename(region = geo_loc_name_country_continent_calc)

