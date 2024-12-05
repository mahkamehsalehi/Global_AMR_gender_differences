library(mia)
library(tidyverse)
library(rstatix)

library(dplyr)
library(rstatix)

# Read in the TSE object
tse <- readRDS("DATA/TSE_filtered.rds")
df <- as.data.frame(colData(tse))

# Ensure geo_loc_name_country_continent_calc is a factor
df <- df %>%
  mutate(geo_loc_name_country_continent_calc = as.factor(geo_loc_name_country_continent_calc))

# Split the data by gender
df_split <- split(df, df$gender)

# Initialize a list to store results for each gender
results_list <- list()

# Loop over each gender-specific dataframe
for (gender in names(df_split)) {
  df_gender <- df_split[[gender]]
  
  # Perform pairwise Wilcoxon tests
  results <- df_gender %>%
    pairwise_wilcox_test(
      log10_ARG_load ~ geo_loc_name_country_continent_calc,
      p.adjust.method = "BH"
    ) %>%
    
    # Add significance annotation
    mutate(
      P.adj_Significance = case_when(
        p.adj > 0.05 ~ "ns",
        p.adj <= 0.05 & p.adj > 0.01 ~ "*",
        p.adj <= 0.01 & p.adj > 0.001 ~ "**",
        p.adj <= 0.001 & p.adj > 0.0001 ~ "***",
        p.adj <= 0.0001 ~ "****"
      )
    ) %>%
    
    # Calculate effect direction
    rowwise() %>%
    mutate(
      Median1 = median(df_gender$ARG_load[df_gender$geo_loc_name_country_continent_calc == group1], na.rm = TRUE),
      Median2 = median(df_gender$ARG_load[df_gender$geo_loc_name_country_continent_calc == group2], na.rm = TRUE),
      Effect_Direction = if_else(Median1 > Median2, paste(group1, ">", group2), paste(group2, ">", group1))
    ) %>%
    ungroup()
  
  # Add results to the list with gender as the key
  results_list[[gender]] <- results
}

# Combine results into a single dataframe (optional)
final_results <- bind_rows(results_list, .id = "gender")

# View the results for each gender
print(final_results, n = Inf, width = Inf)
library(writexl)

# Library
library(writexl)
# Save the final_results as an Excel file
write_xlsx(final_results, path = "RESULTS/Wilcox_regions.xlsx")


### Shannon ####

# Initialize a list to store results for each gender
results_list <- list()

# Loop over each gender-specific dataframe
for (gender in names(df_split)) {
  df_gender <- df_split[[gender]]
  
  # Perform pairwise Wilcoxon tests
  results <- df_gender %>%
    pairwise_wilcox_test(
      log10_ARG_load ~ geo_loc_name_country_continent_calc,
      p.adjust.method = "BH"
    ) %>%
    
    # Add significance annotation
    mutate(
      P.adj_Significance = case_when(
        p.adj > 0.05 ~ "ns",
        p.adj <= 0.05 & p.adj > 0.01 ~ "*",
        p.adj <= 0.01 & p.adj > 0.001 ~ "**",
        p.adj <= 0.001 & p.adj > 0.0001 ~ "***",
        p.adj <= 0.0001 ~ "****"
      )
    ) %>%
    
    # Calculate effect direction
    rowwise() %>%
    mutate(
      Median1 = median(df_gender$shannon_diversity[df_gender$geo_loc_name_country_continent_calc == group1], na.rm = TRUE),
      Median2 = median(df_gender$shannon_diversity[df_gender$geo_loc_name_country_continent_calc == group2], na.rm = TRUE),
      Effect_Direction = if_else(Median1 > Median2, paste(group1, ">", group2), paste(group2, ">", group1))
    ) %>%
    ungroup()
  
  # Add results to the list with gender as the key
  results_list[[gender]] <- results
}

# Combine results into a single dataframe (optional)
final_results <- bind_rows(results_list, .id = "gender")

# View the results for each gender
print(final_results, n = Inf, width = Inf)
library(writexl)

# Library
library(writexl)
# Save the final_results as an Excel file
write_xlsx(final_results, path = "RESULTS/Wilcox_regions_Shannon.xlsx")

