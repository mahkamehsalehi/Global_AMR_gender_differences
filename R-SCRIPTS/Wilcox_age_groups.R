library(dplyr)
library(rstatix)


# Read in the TSE object
tse <- readRDS("DATA/TSE_filtered.rds")
df <- as.data.frame(colData(tse))

# Ensure age_category_new is a factor
df <- df %>%
  mutate(age_category_new = as.factor(age_category_new))

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
      log10_ARG_load ~ age_category_new,  # Changed to age_category_new
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
      Median1 = median(df_gender$ARG_load[df_gender$age_category_new == group1], na.rm = TRUE),
      Median2 = median(df_gender$ARG_load[df_gender$age_category_new == group2], na.rm = TRUE),
      Effect_Direction = if_else(Median1 > Median2, paste(group1, ">", group2), paste(group2, ">", group1))
    ) %>%
    ungroup()
  
  # Add results to the list with gender as the key
  results_list[[gender]] <- results
}

# Combine results into a single dataframe (optional)
final_results <- bind_rows(results_list, .id = "gender")

# View the results for each gender
final_results

# Library
library(writexl)
# Save the final_results as an Excel file
write_xlsx(final_results, path = "RESULTS/Wilcox_agegroups.xlsx")



####### Shannon ######

# Initialize a list to store results for each gender
results_list <- list()

# Loop over each gender-specific dataframe
for (gender in names(df_split)) {
  df_gender <- df_split[[gender]]
  
  # Perform pairwise Wilcoxon tests
  results <- df_gender %>%
    pairwise_wilcox_test(
      shannon_diversity ~ age_category_new,  # Changed to age_category_new
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
      Median1 = median(df_gender$shannon_diversity[df_gender$age_category_new == group1], na.rm = TRUE),
      Median2 = median(df_gender$shannon_diversity[df_gender$age_category_new == group2], na.rm = TRUE),
      Effect_Direction = if_else(Median1 > Median2, paste(group1, ">", group2), paste(group2, ">", group1))
    ) %>%
    ungroup()
  
  # Add results to the list with gender as the key
  results_list[[gender]] <- results
}

# Combine results into a single dataframe (optional)
final_results <- bind_rows(results_list, .id = "gender")

# View the results for each gender
final_results

# Library
library(writexl)
# Save the final_results as an Excel file
write_xlsx(final_results, path = "RESULTS/Wilcox_agegroups_Shannon.xlsx")


