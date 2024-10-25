# Set Working Directory
setwd("/scratch/project_2008149/USER_WORKSPACES/mahkameh/women_amr/")

# Load Necessary Libraries
library(tidyverse)
library(vegan)
library(scater)

# Load the TSE Object
TSE <- readRDS("TSE.rds")

# ---------------------------
# Step 1: Apply Important Filtering
# ---------------------------

# Identify samples with non-missing 'sex_combined' and 'host_age_years'
non_na_samples <- !is.na(colData(TSE)$sex_combined) &
  !is.na(colData(TSE)$host_age_years)

# Subset the TSE object based on this filtering
TSE_filtered <- TSE[, non_na_samples]

# ---------------------------
# Step 2: Further Filter the Data for Females Only
# ---------------------------

# Convert column data to a data frame
df <- as.data.frame(colData(TSE_filtered))

# Filter for female samples
female_samples <- df %>%
  filter(sex_combined == "female")

# Subset the TSE object to include only female samples after filtering
TSE_female <- TSE_filtered[, colData(TSE_filtered)$sex_combined == "female"]

# ---------------------------
# Step 3: Create Detailed Age Groups for Females
# ---------------------------

age_labels <- c(
  "Infant",           # 0-1 year
  "Toddler",          # 1-3 years
  "Preschooler",      # 3-6 years
  "School-age",       # 6-12 years
  "Teen",             # 12-18 years
  "Young Adult",      # 18-35 years
  "Middle Adulthood", # 35-65 years
  "Older Adult"       # 65-100 years
)

female_samples <- female_samples %>%
  mutate(
    age_group = cut(
      host_age_years,
      breaks = c(0, 1, 3, 6, 12, 18, 100),
      labels = age_labels,
      right = FALSE,
      include.lowest = TRUE
    )
  )

# Verify the distribution of age groups
table(female_samples$age_group)

# ---------------------------
# Step 4: Calculate Bray-Curtis Distance Matrix for Females
# ---------------------------

# Extract the relative abundance assay data
assay_data_female <- assay(TSE_female, "relabundance")

# Calculate the Bray-Curtis distance matrix
distance_matrix_female <- vegdist(t(assay_data_female), method = "bray")

# ---------------------------
# Step 5: Perform PCoA
# ---------------------------

# Run PCoA using Bray-Curtis distances
TSE_female <- runMDS(
  TSE_female,
  FUN = vegan::vegdist,
  method = "bray",
  name = "PCoA",
  ncomponents = 3,
  exprs_values = "relabundance"
)

# ---------------------------
# Step 6: Extract PCoA Results for Females
# ---------------------------

# Extract the reduced dimensions (PCoA coordinates)
pcoa_female_data <- as.data.frame(reducedDim(TSE_female, "PCoA"))

# Rename columns to PC1, PC2
colnames(pcoa_female_data) <- paste0("PC", 1:ncol(pcoa_female_data))

# Add age group and region information to the PCoA data frame
pcoa_female_data$age_group <- female_samples$age_group
pcoa_female_data$geo <- female_samples$geo_loc_name_country_continent_calc

# ---------------------------
# Step 7: Extract Eigenvalues and Calculate Percentage Variance Explained
# ---------------------------

# Extract eigenvalues from the PCoA results
eigenvalues_female <- attr(reducedDim(TSE_female, "PCoA"), "eig")

# Retain only positive eigenvalues
positive_eigenvalues_female <- eigenvalues_female[eigenvalues_female > 0]

# Calculate the percentage of variance explained by each principal coordinate
variance_explained_female <- (positive_eigenvalues_female / sum(positive_eigenvalues_female)) * 100

# Extract percentages for the first two principal coordinates
percent_var_PC1_female <- variance_explained_female[1]
percent_var_PC2_female <- variance_explained_female[2]

# ---------------------------
# Step 8: Clean PCoA Data by Removing "uncalculated" Region
# ---------------------------

pcoa_female_data_clean <- pcoa_female_data %>%
  filter(geo != "uncalculated")

# ---------------------------
# Step 9: Plot the PCoA (Colored by Age Group and Faceted by Region)
# ---------------------------

# Create the PCoA plot
pcoa_plot <- ggplot(pcoa_female_data_clean, aes(x = PC1, y = PC2, color = age_group)) +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(~ geo, nrow = 2) +
  labs(
    title = "PCoA Plot (Beta Diversity - Bray-Curtis) for Females by Region",
    x = paste0("PC1 (", round(percent_var_PC1_female, 2), "%)"),
    y = paste0("PC2 (", round(percent_var_PC2_female, 2), "%)"),
    color = "Age Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 12)  # Facet labels
  )

# Display the PCoA Plot
print(pcoa_plot)

# ---------------------------
# Step 10: Prepare Metadata for PERMANOVA
# ---------------------------

# Create a metadata data frame containing 'geo' and 'age_group'
metadata <- pcoa_female_data_clean %>%
  select(geo, age_group)

# ---------------------------
# Step 11: Run PERMANOVA (adonis2)
# ---------------------------

# Test whether regions (geo) significantly explain the differences in beta diversity
permanova_region <- adonis2(
  distance_matrix_female ~ geo,
  data = metadata,
  permutations = 999,
  method = "bray"
)

# Test whether age groups significantly explain the differences in beta diversity
permanova_age <- adonis2(
  distance_matrix_female ~ age_group,
  data = metadata,
  permutations = 999,
  method = "bray"
)

# ---------------------------
# Step 12: Display the PERMANOVA Results
# ---------------------------

# PERMANOVA results for regions
print("PERMANOVA Results for Regions:")
print(permanova_region)

# PERMANOVA results for age groups
print("PERMANOVA Results for Age Groups:")
print(permanova_age)



