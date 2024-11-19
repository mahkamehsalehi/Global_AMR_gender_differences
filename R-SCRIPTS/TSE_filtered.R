library(vegan)
library(ggpubr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(dplyr)

# ---------------------------
# Data Loading and Preprocessing
# ---------------------------

TSE <- readRDS("TSE.rds")

# Filter samples with complete sex and age data
non_na_samples <- !is.na(colData(TSE)$sex_combined) & !is.na(colData(TSE)$host_age_years)
TSE_filtered <- TSE[, non_na_samples]

# Further filter to retain only adult samples with valid geographic information
TSE_filtered <- TSE_filtered[, colData(TSE_filtered)$geo_loc_name_country_continent_calc != "uncalculated"]

tse_metadata <- as.data.frame(colData(TSE_filtered))

# Extract and clean assay data (relative abundances)
assay_data_clean <- assay(TSE_filtered, "relabundance")

# Remove samples (columns) with all-zero counts
non_empty_samples <- colSums(assay_data_clean > 0, na.rm = TRUE) > 0
assay_data_clean <- assay_data_clean[, non_empty_samples]

# Subset metadata to match non-empty samples
tse_metadata <- tse_metadata[non_empty_samples, ]

# ---------------------------
# Diversity and Distance Calculations
# ---------------------------

# Calculate the Bray-Curtis distance matrix on the cleaned data
distance_matrix <- vegan::vegdist(t(assay_data_clean), method = "bray")

# Calculate Shannon diversity
diversity_indices <- vegan::diversity(t(assay_data_clean), index = "shannon")

# Add Shannon diversity to the metadata
tse_metadata$shannon_diversity <- diversity_indices
colData(TSE_filtered)$shannon_diversity <- tse_metadata$shannon_diversity

# Save the updated TSE object with Shannon diversity included
saveRDS(TSE_filtered, file = "TSE_shan.rds")

# Reload the dataset to continue with further processing
tse <- readRDS("TSE_shan.rds")

# Run MDS (PCoA) analysis
TSE_adult <- runMDS(
  tse ,
  FUN = vegan::vegdist,
  method = "bray",
  name = "PCoA",
  ncomponents = 3,
  exprs_values = "relabundance"
)

# ---------------------------
# Calculate Variance Explained
# ---------------------------

# Extract eigenvalues from the PCoA result
eigenvalues <- attr(reducedDim(TSE_adult, "PCoA"), "eig")

# Calculate the proportion of variance explained by each axis
variance_explained <- eigenvalues / sum(eigenvalues)

# Print variance explained for the first few axes
variance_explained[1:3]  # Variance for the first three PCoA axes

# ---------------------------
# Extract PCoA Coordinates
# ---------------------------

# Extract PCoA coordinates
pcoa_coordinates <- reducedDim(TSE_adult, "PCoA")

# Convert coordinates to a data frame for plotting
pcoa_df <- as.data.frame(pcoa_coordinates)
colnames(pcoa_df) <- paste0("PC", 1:ncol(pcoa_df))

# ---------------------------
# Add Metadata to PCoA Data Frame
# ---------------------------

pcoa_df <- cbind(pcoa_df, as.data.frame(colData(TSE_adult)))


# ---------------------------
# Update colData with PCoA Coordinates
# ---------------------------

# Ensure the rownames of the PCoA data frame match the colData of the TSE object
if (!all(rownames(pcoa_df) == rownames(colData(tse)))) {
  stop("Row names of PCoA data and TSE metadata do not match. Align before proceeding.")
}

# Add PCoA coordinates to colData
for (i in seq_len(ncol(pcoa_coordinates))) {
  colData(tse)[, paste0("PC", i)] <- pcoa_coordinates[, i]
}

# Optionally, add variance explained for reference
attr(colData(tse), "variance_explained") <- variance_explained

# ---------------------------
# Save the Updated TSE Object
# ---------------------------

# Save the TSE object with the updated colData
saveRDS(tse, file = "TSE_with_PCoA.rds")

# ---------------------------
# Metadata Transformation
# ---------------------------

tse <- readRDS("DATA/TSE_with_PCoA.rds")
tse_metadata <- as.data.frame(colData(tse))

# Add and format gender column for consistent labeling
tse_metadata <- tse_metadata %>%
  mutate(
    gender = case_when(
      sex_combined == "male" ~ "Men",
      sex_combined == "female" ~ "Women",
      TRUE ~ sex_combined
    ),
    gender = factor(gender, levels = c("Women", "Men"))  # Order: Women first
  )

# Categorize age into meaningful groups
tse_metadata <- tse_metadata %>%
  mutate(
    age_category = case_when(
      host_age_years >= 0 & host_age_years < 1 ~ "Infant",
      host_age_years >= 1 & host_age_years < 3 ~ "Toddler",
      host_age_years >= 3 & host_age_years < 5 ~ "Preschooler",
      host_age_years >= 5 & host_age_years < 12 ~ "School-Age Child",
      host_age_years >= 12 & host_age_years < 18 ~ "Teen",
      host_age_years >= 18 & host_age_years < 35 ~ "Young Adult",
      host_age_years >= 35 & host_age_years < 65 ~ "Middle Adult",
      host_age_years >= 65 & host_age_years <= 100 ~ "Older Adult",
      TRUE ~ NA_character_
    ),
    age_category = factor(age_category, levels = c( # Set order for plots
      "Infant", "Toddler", "Preschooler", "School-Age Child", 
      "Teen", "Young Adult", "Middle Adult", "Older Adult"
    ))
  )

# Log-transform ARG load and add it to the metadata
tse_metadata <- tse_metadata %>%
  mutate(log_ARG_load = log(ARG_load))

# Add income group classifications
tse_metadata <- tse_metadata %>%
  mutate(
    income_group = case_when(
      World_Bank_Income_Group == "High income" ~ "HIC",  # High-income countries
      World_Bank_Income_Group %in% c("Low income", "Lower middle income", "Upper middle income") ~ "LMIC",  # Low- and middle-income countries
      TRUE ~ NA_character_
    )
  )

# ---------------------------
# Update TSE Object
# ---------------------------

# Update the TSE object with the processed metadata
colData(tse)$age_category <- tse_metadata$age_category
colData(tse)$gender <- tse_metadata$gender
colData(tse)$income_group <- tse_metadata$income_group
colData(tse)$log_ARG_load <- tse_metadata$log_ARG_load

# Save the final filtered and updated dataset
saveRDS(tse, file = "DATA/TSE_filtered.rds")
