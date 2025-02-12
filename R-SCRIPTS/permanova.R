# Load required libraries
library(tidyverse)
library(vegan)

# Load the dataset
TSE <- readRDS("TSE.rds")

# ---------------------------
# Step 1: Filter Missing Data
# ---------------------------
# Remove samples with missing gender or age
non_na_samples <- !is.na(colData(TSE)$sex_combined) &
  !is.na(colData(TSE)$host_age_years)

TSE_filtered <- TSE[, non_na_samples]
df <- as.data.frame(colData(TSE_filtered))

# ---------------------------
# Step 2: Define Age Groups (Match Manuscript Categories)
# ---------------------------
df <- df %>%
  mutate(
    age_group = cut(
      host_age_years,
      breaks = c(0, 1, 3, 12, 20, 35, 65, 80, 100),
      labels = c("Infant", "Toddler", "Children", "Teenager",
                 "Young Adult", "Middle-Aged Adult", "Older Adult", "Oldest Adult"),
      right = FALSE,
      include.lowest = TRUE
    )
  )

# ---------------------------
# Step 3: Compute Bray-Curtis Distance Matrix
# ---------------------------
assay_data <- assay(TSE_filtered, "relabundance") # Extract abundance data
distance_matrix <- vegdist(t(assay_data), method = "bray") # Compute Bray-Curtis distances
distance_matrix <- as.matrix(distance_matrix)


# ---------------------------
# Step 4: Prepare Metadata for PERMANOVA
# ---------------------------
metadata <- df %>%
  select(geo_loc_name_country_continent_calc, age_group, sex_combined, GDP_per_head, Usage) %>%
  dplyr::rename(
    geo = geo_loc_name_country_continent_calc,
    gender = sex_combined,
    GDP_per_capita = GDP_per_head
  ) %>%
  drop_na()

# Ensure metadata matches the distance matrix
common_samples <- intersect(rownames(distance_matrix), rownames(metadata))
distance_matrix <- distance_matrix[common_samples, common_samples]
metadata <- metadata[common_samples, ]


# ---------------------------
# Step 5: Run Multivariate PERMANOVA
# ---------------------------
permanova_multivariable <- adonis2(
  distance_matrix ~ geo + age_group + gender + GDP_per_capita + Usage,
  data = metadata,
  permutations = 999,
  method = "bray"
)

# ---------------------------
# Step 6: Display Results
# ---------------------------
print("Multivariate PERMANOVA Results:")
print(permanova_multivariable)
