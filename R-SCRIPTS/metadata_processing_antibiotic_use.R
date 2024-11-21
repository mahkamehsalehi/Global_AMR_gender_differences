library(tidyquant)
library(tidyverse)
library(SummarizedExperiment)

source("R-SCRIPTS/metadata_processing_sex.R")

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
df_selected$country <- df_selected$geo_loc_name_country_calc
global_data_2016 <- merge(df_selected, data_2016_subset, by = "country")
write.csv(global_data_2016, "DATA/global_data_2016.csv")



g_data_2016 <- read_csv("DATA/global_data_2016.csv")


ab_use <- g_data_2016 <- read_csv("DATA/global_data_2016.csv")
tse <- readRDS("DATA/TSE.rds")

common_samples <- intersect(colData(tse)$acc, ab_use$acc)
cat("Number of common samples:", length(common_samples), "\n")

# Filter antibiotic use data to common samples
ab_use_filtered <- ab_use %>%
  filter(acc %in% common_samples)

ab_use_filtered <- dplyr::rename(ab_use_filtered, usage_bayesian = antibiotic_consumption)

tse_coldata <- as.data.frame(colData(tse))

# Add usage_bayesian to tse_coldata
merged_coldata <- tse_coldata %>%
  dplyr::left_join(ab_use_filtered %>% select(acc, usage_bayesian), by = "acc")

# Preserve rownames
rownames(merged_coldata) <- rownames(tse_coldata)

# Update colData in tse
colData(tse) <- DataFrame(merged_coldata)

#saveRDS(tse, file = "TSE_AB_estimate.rds")
