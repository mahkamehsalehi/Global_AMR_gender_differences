#-------------------------------------------------------------------------------
# Project: Women's AMR Analysis
# Purpose: Perform data filtering, PCoA, Shannon diversity analysis, and visualization.
# Author: Mahkameh
# Date: 2024-10-29
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 1. Setup: Working Directory and Library Imports
#-------------------------------------------------------------------------------

# Set the working directory
setwd("/scratch/project_2008149/USER_WORKSPACES/mahkameh/women_amr/")

# Load necessary libraries for data manipulation, diversity analysis, and plotting
library(tidyverse)
library(vegan)
library(scater)
library(ggpubr)
library(rstatix)
library(viridis)
library(ggsignif)

#-------------------------------------------------------------------------------
# 2. Data Loading and Initial Filtering
#-------------------------------------------------------------------------------

# Load the TSE object containing AMR data
TSE <- readRDS("TSE.rds")

# Step 1: Remove samples with missing sex or age data
non_na_samples <- !is.na(colData(TSE)$sex_combined) & !is.na(colData(TSE)$host_age_years)

# Subset the TSE object to include only samples with complete sex and age data
TSE_filtered <- TSE[, non_na_samples]

# Step 2: Filter to retain female samples only
TSE_female <- TSE_filtered[, colData(TSE_filtered)$sex_combined == "female"]

# Extract metadata for female samples
female_metadata <- as.data.frame(colData(TSE_female))

#-------------------------------------------------------------------------------
# 3. Defining Age Groups for Females
#-------------------------------------------------------------------------------

# Group age into categories: Infant, Toddler, Child, and Adult
female_metadata <- female_metadata %>%
  mutate(
    age_group = cut(
      host_age_years,
      breaks = c(0, 1, 3, 18, Inf),
      labels = c("Infant", "Toddler", "Child", "Adult"),
      right = FALSE
    )
  )

#-------------------------------------------------------------------------------
# 4. Bray-Curtis PCoA Analysis for Females
#-------------------------------------------------------------------------------

# Compute Bray-Curtis distance and conduct PCoA on relative abundance data for females
assay_data_female <- assay(TSE_female, "relabundance")
distance_matrix_female <- vegdist(t(assay_data_female), method = "bray")


common_theme <- theme_classic(base_size = 14) +
  theme(
    text = element_text(family = "Sans"),
    plot.title = element_text(
      face = "bold", 
      size = 14, 
      hjust = 0.5
    ),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.position = "none",
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 12, face = "bold")
  )

# Run PCoA and extract the results
TSE_female <- runMDS(
  TSE_female,
  FUN = vegan::vegdist,
  method = "bray",
  name = "PCoA",
  ncomponents = 3,
  exprs_values = "relabundance"
)
pcoa_female_data <- as.data.frame(reducedDim(TSE_female, "PCoA"))
colnames(pcoa_female_data) <- paste0("PC", 1:ncol(pcoa_female_data))

# Integrate metadata for region and age group with PCoA data
pcoa_female_data <- pcoa_female_data %>%
  mutate(
    age_group = female_metadata$age_group,
    geo = female_metadata$geo_loc_name_country_continent_calc
  ) %>%
  filter(geo != "uncalculated")

# Compute variance explained by each principal coordinate
eigenvalues_female <- attr(reducedDim(TSE_female, "PCoA"), "eig")
positive_eigenvalues_female <- eigenvalues_female[eigenvalues_female > 0]
variance_explained_female <- (positive_eigenvalues_female / sum(positive_eigenvalues_female)) * 100
percent_var_PC1_female <- variance_explained_female[1]
percent_var_PC2_female <- variance_explained_female[2]

#-------------------------------------------------------------------------------
# 5. PCoA Plot for Females by Region
#-------------------------------------------------------------------------------

# Plot PCoA results with color differentiation by age group
pcoa_plot <- ggplot(pcoa_female_data, aes(x = PC1, y = PC2, color = age_group)) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(~ geo, nrow = 2) +
  labs(
    title = "PCoA Plot (Beta Diversity - Bray-Curtis) for Females by Region",
    x = paste0("PC1 (", round(percent_var_PC1_female, 2), "%)"),
    y = paste0("PC2 (", round(percent_var_PC2_female, 2), "%)"),
    color = "Age Group"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_blank(),
    axis.line = element_line(color = "black")
  )

ggsave("ff_pcoa.png", pcoa_plot, width = 10, height = 6)

# Display PCoA plot
print(pcoa_plot)

#-------------------------------------------------------------------------------
# 6. Shannon Diversity Analysis for Female Samples
#-------------------------------------------------------------------------------

# Step 1: Compute Shannon diversity for each female sample
shannon_diversity_female <- diversity(t(assay_data_female), index = "shannon")

# Step 2: Merge Shannon diversity values with female metadata
female_metadata <- female_metadata %>%
  mutate(shannon_diversity = shannon_diversity_female) %>%
  dplyr::rename(geo = geo_loc_name_country_calc) %>%
  select(shannon_diversity, GDP_per_head, host_age_years, age_group, geo)

# Step 3: Plot Shannon Diversity vs. GDP
shannon_gdp_plot <- ggplot(female_metadata, aes(x = GDP_per_head, y = shannon_diversity)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  labs(
    title = "Shannon Diversity vs. GDP",
    x = "GDP per Head",
    y = "Shannon Diversity"
  ) + 
  theme_minimal(base_size = 14, base_family = "Serif") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    panel.border = element_blank(),
    axis.line = element_line(color = "black")
  )

# Display Shannon vs. GDP plot
print(shannon_gdp_plot)


# Step 4: Plot Shannon Diversity vs. Age
shannon_age_plot <- ggplot(female_metadata, aes(x = host_age_years, y = shannon_diversity)) +
  geom_point(alpha = 0.6, color = "#fc7672") +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  labs(
    title = "Shannon Diversity vs. Age Among Women",
    x = "Age (Years)",
    y = "Shannon Diversity"
  ) +
  theme_minimal(base_size = 14, base_family = "Serif") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    panel.border = element_blank(),
    axis.line = element_line(color = "black")
  )

# Display Shannon vs. Age plot
print(shannon_age_plot)

#-------------------------------------------------------------------------------
# 7. Gender-Based Shannon Diversity Analysis
#-------------------------------------------------------------------------------

# Compute Shannon diversity for all samples in the filtered TSE object
assay_data_all <- assay(TSE_filtered, "relabundance")
shannon_diversity_all <- diversity(t(assay_data_all), index = "shannon")

# Merge Shannon diversity values with metadata for all samples
metadata_all <- as.data.frame(colData(TSE_filtered)) %>%
  mutate(shannon_diversity = shannon_diversity_all) %>%
  select(shannon_diversity, sex_combined, GDP_per_head, host_age_years) %>%
  filter(!is.na(GDP_per_head), !is.na(host_age_years), !is.na(sex_combined))

# Conduct Wilcoxon and t-tests to compare Shannon diversity between genders
wilcox_test_gender <- wilcox.test(shannon_diversity ~ sex_combined, data = metadata_all)
print(wilcox_test_gender)
t_test_gender <- t.test(shannon_diversity ~ sex_combined, data = metadata_all)
print(t_test_gender)

# Plot Shannon Diversity by Gender with significance annotations
shannon_gender_plot_enhanced <- ggplot(metadata_all, aes(x = sex_combined, y = shannon_diversity, fill = sex_combined)) +
  geom_jitter(aes(color = sex_combined), width = 0.2, alpha = 0.6) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  labs(
    title = "Comparison of Shannon Diversity between Men and Women",
    x = "Gender",
    y = "Shannon Diversity"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  stat_compare_means(
    comparisons = list(c("male", "female")),
    method = "wilcox.test",
    label = "p.signif",
    label.y = max(metadata_all$shannon_diversity, na.rm = TRUE) + 0.1
  )

# Display the plot
print(shannon_gender_plot_enhanced)

#-------------------------------------------------------------------------------
# 8. Detailed Age Category Shannon Diversity Analysis
#-------------------------------------------------------------------------------

# Define detailed age categories and labels
age_labels <- c("Infant", "Toddler", "Child", "Adult")
female_metadata <- female_metadata %>%
  mutate(
    age_category = case_when(
      host_age_years >= 0 & host_age_years < 1 ~ "Infant",
      host_age_years >= 1 & host_age_years < 3 ~ "Toddler",
      host_age_years >= 3 & host_age_years < 18 ~ "Child",
      host_age_years >= 18 ~ "Adult",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(age_category), !is.na(shannon_diversity)) %>%
  mutate(age_category = factor(age_category, levels = age_labels, ordered = TRUE))

# Define pairwise age comparisons and set custom y-positions for annotations
comparisons_age <- list(c("Infant", "Toddler"), c("Toddler", "Child"), c("Child", "Adult"))
y_positions <- c(
  max(female_metadata$shannon_diversity, na.rm = TRUE) + 0.5,
  max(female_metadata$shannon_diversity, na.rm = TRUE) + 1.0,
  max(female_metadata$shannon_diversity, na.rm = TRUE) + 1.5
)

# Plot Shannon Diversity by Age Category with annotations
shannon_age_boxplot <- ggplot(female_metadata, aes(x = age_category, y = shannon_diversity)) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "red") +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  labs(
    title = "Shannon Diversity by Age Category Among Women",
    x = "Age Category",
    y = "Shannon Diversity"
  ) +
  theme_minimal(base_size = 14, base_family = "Serif") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  stat_compare_means(
    comparisons = comparisons_age,
    method = "t.test",
    label = "p.signif",
    y.position = y_positions
  )

# Display the age-based Shannon diversity plot
print(shannon_age_boxplot)
