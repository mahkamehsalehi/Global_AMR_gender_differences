library(vegan)
library(ggpubr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(dplyr)

TSE <- readRDS("DATA/TSE_gender_age.rds")

# Filter samples with complete sex and age data
non_na_samples <- !is.na(colData(TSE)$sex_combined) & !is.na(colData(TSE)$host_age_years)
TSE_filtered <- TSE[, non_na_samples]

# Further filter to retain only adult samples with valid geographic information
TSE_filtered <- TSE_filtered[, colData(TSE_filtered)$geo_loc_name_country_continent_calc != "uncalculated"]

# Step 2: Extract metadata and set factor levels for consistent ordering in plots
tse_metadata <- as.data.frame(colData(TSE_filtered))

# Rename 'male' and 'female' to 'Men' and 'Women' and set order
tse_metadata <- tse_metadata %>%
  mutate(sex_combined = case_when(
    sex_combined == "male" ~ "Men",
    sex_combined == "female" ~ "Women",
    TRUE ~ sex_combined
  ))

# Define gender order for plots
tse_metadata$sex_combined <- factor(tse_metadata$sex_combined, levels = c("Women", "Men"))

# Categorize ages into meaningful groups
tse_metadata <- tse_metadata %>%
  mutate(age_category = case_when(
    host_age_years >= 0 & host_age_years < 1 ~ "Infant",
    host_age_years >= 1 & host_age_years < 3 ~ "Toddler",
    host_age_years >= 3 & host_age_years < 18 ~ "Child",
    host_age_years >= 18 & host_age_years < 35 ~ "Young Adult",
    host_age_years >= 35 & host_age_years < 65 ~ "Middle-Aged Adult",
    host_age_years >= 65 & host_age_years <= 100 ~ "Older Adult",
    TRUE ~ NA_character_
  ))

# Step 3: Extract and clean assay data (relative abundances)
assay_data_clean <- assay(TSE_filtered, "relabundance")

# Remove columns with all-zero counts to retain non-empty samples only
non_empty_samples <- colSums(assay_data_clean > 0, na.rm = TRUE) > 0
assay_data_clean <- assay_data_clean[, non_empty_samples]

# Subset metadata to match non-empty samples
tse_metadata <- tse_metadata[non_empty_samples, ]

# Step 4: Calculate Bray-Curtis distance matrix on the cleaned data and Shannon diversity
distance_matrix <- vegan::vegdist(t(assay_data_clean), method = "bray")
diversity_indices <- diversity(t(assay_data_clean), index = "shannon")
tse_metadata$shannon_diversity <- diversity_indices

# Define a custom theme for all plots
common_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.position = "none",
    axis.line = element_line(color = "black"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 12, face = "bold")
  )

# Set age category ordering for plots
tse_metadata$age_category <- factor(tse_metadata$age_category, 
                                    levels = c("Infant", "Toddler", "Child", 
                                               "Young Adult", "Middle-Aged Adult", "Older Adult"))

# Filter metadata for valid income groups and remove specific categories
filtered_metadata <- tse_metadata %>%
  filter(
    category != "Infant Study", 
    geo_loc_name_country_calc != "Zimbabwe"
  )

# Classify income groups into 'HIC' and 'LMIC' and remove missing values for ARG load
filtered_metadata <- filtered_metadata %>%
  mutate(
    Income_Group = case_when(
      World_Bank_Income_Group == "High income" ~ "HIC",
      World_Bank_Income_Group %in% c("Low income", "Lower middle income", "Upper middle income") ~ "LMIC",
      TRUE ~ NA_character_
    )
  ) %>%
  drop_na(log10_ARG_load, Income_Group)

# Plot ARG load by gender across income groups
income_arg_boxplot <- ggplot(filtered_metadata, aes(x = sex_combined, y = log10_ARG_load, fill = sex_combined)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.8),
              size = 0.5, alpha = 0.5, aes(color = sex_combined), show.legend = FALSE) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6, alpha = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(x = "Gender", y = "Antibiotic Resistance Load") +
  theme_minimal(base_family = "Arial") +
  stat_compare_means(
    comparisons = list(c("Women", "Men")), 
    label = "p.signif", 
    method = "wilcox.test", 
    p.adjust.method = "BH",
    hide.ns = FALSE) +
  facet_wrap(~Income_Group) +
  common_theme

# Plot Shannon diversity by gender across income groups
income_shannon_boxplot <- ggplot(filtered_metadata, aes(x = sex_combined, y = shannon_diversity, fill = sex_combined)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.8),
              size = 0.5, alpha = 0.5, aes(color = sex_combined), show.legend = FALSE) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6, alpha = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(x = "Gender", y = "Resistome Diversity") +
  theme_minimal(base_family = "Arial") +
  stat_compare_means(aes(x = sex_combined, y = shannon_diversity),
                     comparisons = list(c("Women", "Men")), 
                     label = "p.signif", 
                     method = "wilcox.test", 
                     p.adjust.method = "BH",
                     hide.ns = FALSE) +
  facet_wrap(~Income_Group) +
  common_theme

print(combined_figure_income)

# Plot ARG load and Shannon diversity across income groups for women only
income_arg_boxplot_female <- ggplot(filtered_metadata %>% filter(sex_combined == "Women"), 
                                    aes(x = Income_Group, y = log10_ARG_load, fill = "#F8766D")) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.5, alpha = 0.5, color = "#F8766D", show.legend = FALSE) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6, alpha = 0.7, show.legend = FALSE) +
  labs(x = "Income Group", y = "Antibiotic Resistance Load") +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c("HIC", "LMIC")), 
                     label = "p.signif", 
                     method = "wilcox.test", 
                     p.adjust.method = "BH",
                     hide.ns = FALSE) +
  common_theme

income_shannon_boxplot_female <- ggplot(filtered_metadata %>% filter(sex_combined == "Women"),
                                        aes(x = Income_Group, y = shannon_diversity, fill = "#F8766D")) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.5, alpha = 0.5, color = "#F8766D", show.legend = FALSE) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6, alpha = 0.7, show.legend = FALSE) +
  labs(x = "Income Group", y = "Resistome Diversity") +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c("HIC", "LMIC")), 
                     label = "p.signif", 
                     method = "wilcox.test", 
                     p.adjust.method = "BH",
                     hide.ns = FALSE) +
  common_theme

# Combine and display the female-only figures
combined_figure_income_female <- income_arg_boxplot_female + income_shannon_boxplot_female + plot_layout(ncol = 2)
print(combined_figure_income_female)
ggsave("final/income_female_panel.png", combined_figure_income_female, width = 20, height = 10)


# Combine and display the ARG load and Shannon diversity figures across income groups
combined_figure_income <- income_arg_boxplot + income_shannon_boxplot +
  income_arg_boxplot_female + income_shannon_boxplot_female +
  plot_layout(ncol = 2, nrow = 2) + 
  plot_annotation(tag_levels = 'a')

ggsave("RESULTS/FIGURES/income_panel.png", combined_figure_income, width = 10, height = 8)
