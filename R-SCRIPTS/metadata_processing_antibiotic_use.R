library(tidyquant)
library(tidyverse)
library(SummarizedExperiment)



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#--------------- Adding Antibiotic Use Data to the TSE Object ------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

data_2016 <- read_csv("DATA/total_antibiotic_consumption_estimates.csv") %>%
  subset(Year ==2016) %>%
  mutate(Location = case_when(
    Location == "United States" ~ "USA",
    Location == "South Korea" ~ "South Korea",
    Location == "Congo" ~ "Democratic Republic of the Congo",
    Location == "Bahamas" ~ "The Bahamas",
    Location == "Gambia" ~ "The Gambia",
    TRUE ~ Location
  )) %>%
  dplyr::rename(country = Location) %>%
  dplyr::rename(antibiotic_consumption = `Antibiotic consumption (DDD/1,000/day)`)

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



ab_use <- g_data_2016 <- read_csv("global_data_2016.csv")
tse <- readRDS("TSE.rds")

common_samples <- intersect(colData(tse)$acc, ab_use$acc)
cat("Number of common samples:", length(common_samples), "\n")

# Filter antibiotic use data to common samples
ab_use_filtered <- ab_use %>%
  filter(acc %in% common_samples)

ab_use_filtered <- dplyr::rename(ab_use_filtered, usage_bayesian = antibiotic_consumption)

tse_coldata <- as.data.frame(colData(tse))

# Add usage_bayesian to tse_coldata
merged_coldata <- tse_coldata %>%
  left_join(ab_use_filtered %>% select(acc, usage_bayesian), by = "acc")

# Preserve rownames
rownames(merged_coldata) <- rownames(tse_coldata)

# Update colData in tse
colData(tse) <- DataFrame(merged_coldata)

#saveRDS(tse, file = "TSE_AB_estimate.rds")

# Check for missing values
num_missing <- sum(is.na(colData(tse)$usage_bayesian))
cat("Number of samples without Bayesian_usage data:", num_missing, "\n")