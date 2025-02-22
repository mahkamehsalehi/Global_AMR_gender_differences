library(vegan)
library(tidyverse)

tse <- readRDS("DATA/TSE.rds")


#----------------------
# Subset TSE.rds by gender and age and region!= "uncalculated"
#----------------------

# First filter for gender and age
tse_filtered <- tse[, !is.na(colData(tse)$sex_combined) & 
                      colData(tse)$sex_combined != "" & 
                      !is.na(colData(tse)$host_age_years) & 
                      colData(tse)$host_age_years != ""]

# Create intermediate data frame from the filtered data
df_filtered <- as.data.frame(colData(tse_filtered))

# Then filter out "uncalculated" entries using the filtered data frame
tse_filtered_clean <- tse_filtered[, df_filtered$geo_loc_name_country_continent_calc != "uncalculated"]

# Verify the filtering worked
df_clean <- as.data.frame(colData(tse_filtered_clean))
table(df_clean$geo_loc_name_country_continent_calc)

saveRDS(tse_filtered_clean, "DATA/TSE_subset.rds")

#----------------------
# Perform PCoA
#----------------------

library(vegan)
library(ape)
library(tidyverse)
library(SummarizedExperiment)

# Load the filtered dataset
tse_subset <- readRDS("TSE_subset.rds")

# Extract relative abundance data and compute Bray-Curtis distances
assay_data <- assay(tse_subset, "relabundance")
distance_matrix <- vegdist(t(assay_data), method = "bray")

# Perform PCoA using ape::pcoa()
pcoa_result <- pcoa(distance_matrix, correction = "cailliez")  # Correction for negative eigenvalues

# Extract variance explained
eigenvals <- pcoa_result$values$Relative_eig
var_explained_PC1 <- eigenvals[1] * 100
var_explained_PC2 <- eigenvals[2] * 100
var_explained_PC3 <- eigenvals[3] * 100

# Convert PCoA coordinates to a matrix (first 3 principal coordinates)
pcoa_coords <- as.matrix(pcoa_result$vectors[, 1:3])
colnames(pcoa_coords) <- c("PC1", "PC2", "PC3")

# Store PCoA coordinates in `reducedDim`
reducedDim(tse_subset, "PCoA") <- pcoa_coords

# Save PC1, PC2, PC3 into colData
colData(tse_subset)$PC1 <- pcoa_coords[, 1]
colData(tse_subset)$PC2 <- pcoa_coords[, 2]
colData(tse_subset)$PC3 <- pcoa_coords[, 3]

# Store variance explained in metadata
metadata(tse_subset)$PCoA <- list(
  variance_explained = eigenvals * 100,
  variance_explained_PC1 = var_explained_PC1,
  variance_explained_PC2 = var_explained_PC2,
  variance_explained_PC3 = var_explained_PC3,
  eigenvalues = pcoa_result$values$Eigenvalues
)

# Print variance explained
cat(sprintf("Variance explained:\nPC1: %.1f%%\nPC2: %.1f%%\nPC3: %.1f%%", 
            var_explained_PC1, var_explained_PC2, var_explained_PC3))

# Create temporary data frame from colData
df_temp <- as.data.frame(colData(tse_subset))

# Add all new columns using dplyr
df_temp <- df_temp %>%
  mutate(
    # Add gender column
    gender = case_when(
      sex_combined == "female" ~ "Women",
      sex_combined == "male" ~ "Men",
      TRUE ~ NA_character_
    ),
    # Add age_category_new column
    age_category_new = case_when(
      host_age_years >= 0   & host_age_years <= 1  ~ "Infant",
      host_age_years >  1   & host_age_years <= 3  ~ "Toddler",
      host_age_years >  3   & host_age_years <= 12 ~ "Children",
      host_age_years >  12  & host_age_years <  20 ~ "Teenager",
      host_age_years >= 20  & host_age_years <  35 ~ "Young Adult",
      host_age_years >= 35  & host_age_years <  65 ~ "Middle-Aged Adult",
      host_age_years >= 65  & host_age_years <  80 ~ "Older Adult",
      host_age_years >= 80  & host_age_years <= 100 ~ "Oldest Adult",
      TRUE ~ NA_character_
    ),
    # Add income_group column
    income_group = case_when(
      World_Bank_Income_Group == "High income" ~ "HIC",
      World_Bank_Income_Group %in% c("Low income", "Lower middle income", "Upper middle income") ~ "LMIC",
      TRUE ~ NA_character_
    )
  )

# Convert age_category_new to factor with correct levels
df_temp$age_category_new <- factor(df_temp$age_category_new, 
                                   levels = c("Infant", "Toddler", "Children", "Teenager",
                                              "Young Adult", "Middle-Aged Adult", "Older Adult", "Oldest Adult"))

# Update colData with the new columns
colData(tse_subset)$gender <- df_temp$gender
colData(tse_subset)$age_category_new <- df_temp$age_category_new
colData(tse_subset)$income_group <- df_temp$income_group

# Verify the new columns
table(colData(tse_subset)$gender, useNA = "always")
table(colData(tse_subset)$age_category_new, useNA = "always")
table(colData(tse_subset)$income_group, useNA = "always")

# Save the updated object
saveRDS(tse_subset, "TSE_subset_updated_final.rds")
