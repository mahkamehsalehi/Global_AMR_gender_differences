# Filter samples with complete sex and age data
non_na_samples <- !is.na(colData(TSE)$sex_combined) & !is.na(colData(TSE)$host_age_years)
TSE_filtered <- TSE[, non_na_samples]

# Further filter for samples with calculated geographic information for adults
TSE_filtered <- TSE_filtered[, colData(TSE_filtered)$geo_loc_name_country_continent_calc != "uncalculated"]

# Step 2: Extract metadata and set factors for consistent ordering in plots
tse_metadata <- as.data.frame(colData(TSE_filtered))

# Standardize 'sex_combined' to 'Men' and 'Women' and set gender order
tse_metadata <- tse_metadata %>%
  mutate(sex_combined = case_when(
    sex_combined == "male" ~ "Men",
    sex_combined == "female" ~ "Women",
    TRUE ~ sex_combined
  ))

tse_metadata$sex_combined <- factor(tse_metadata$sex_combined, levels = c("Women", "Men"))

# Categorize age groups
tse_metadata <- tse_metadata %>%
  mutate(
    age_category = case_when(
      host_age_years >= 0 & host_age_years < 1 ~ "Infant",
      host_age_years >= 1 & host_age_years < 3 ~ "Toddler",
      host_age_years >= 3 & host_age_years < 18 ~ "Child",
      host_age_years >= 18 & host_age_years < 35 ~ "Young Adult",
      host_age_years >= 35 & host_age_years < 65 ~ "Middle Adulthood",
      host_age_years >= 65 & host_age_years <= 100 ~ "Older Adult",
      TRUE ~ NA_character_
    )
  )

# Step 3: Extract and clean assay data (relative abundances)
assay_data_clean <- assay(TSE_filtered, "relabundance")

# Remove columns with all-zero counts
non_empty_samples <- colSums(assay_data_clean > 0, na.rm = TRUE) > 0
assay_data_clean <- assay_data_clean[, non_empty_samples]

# Subset metadata to match non-empty samples
tse_metadata <- tse_metadata[non_empty_samples, ]

# Step 4: Calculate Bray-Curtis distance and Shannon diversity
distance_matrix <- vegan::vegdist(t(assay_data_clean), method = "bray")
diversity_indices <- diversity(t(assay_data_clean), index = "shannon")
tse_metadata$shannon_diversity <- diversity_indices

# Define a common theme for plots
common_theme <- theme_classic(base_size = 14) +
  theme(
    text = element_text(family = "Sans"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.position = "none",
    axis.line = element_line(color = "black"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 12, face = "bold")
  )

# Set age category order for plots
tse_metadata$age_category <- factor(tse_metadata$age_category, 
                                    levels = c("Infant", "Toddler", "Child", "Young Adult", "Middle Adulthood", "Older Adult"))

# Filter and prepare data for usage analysis with custom 'Usage_group' categorization
usage_filtered_metadata <- tse_metadata %>%
  mutate(
    Usage_group = case_when(
      Usage > 11  ~ "Above 11",
      Usage <= 11 ~ "11 or Below",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(Usage_group = factor(Usage_group, levels = c("11 or Below", "Above 11"))) %>%
  drop_na(sex_combined, log10_ARG_load, Usage_group)

# Plot ARG load by gender across usage groups
usage_arg_boxplot_both <- ggplot(usage_filtered_metadata, aes(x = sex_combined, y = log10_ARG_load, fill = sex_combined)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 0.5, alpha = 0.5, aes(color = sex_combined), show.legend = FALSE) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6, alpha = 0.7, show.legend = FALSE) +
  labs(x = "Gender", y = "log10 ARG Load") +
  ggtitle("Antibiotic Resistance Genes Load vs. Antibiotic Usage") +
  theme_minimal(base_family = "Arial") +
  stat_compare_means(aes(x = sex_combined, y = log10_ARG_load),
                     comparisons = list(c("Women", "Men")), label = "p.signif", method = "wilcox.test", hide.ns = FALSE) +
  facet_wrap(~Usage_group) +
  common_theme

# Plot Shannon diversity by gender across usage groups
usage_shannon_boxplot_both <- ggplot(usage_filtered_metadata, aes(x = sex_combined, y = shannon_diversity, fill = sex_combined)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 0.5, alpha = 0.5, aes(color = sex_combined), show.legend = FALSE) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6, alpha = 0.7, show.legend = FALSE) +
  labs(x = "Gender", y = "Resistome Diversity") +
  ggtitle("Resistome Diversity vs. Antibiotic Usage") +
  theme_minimal(base_family = "Arial") +
  stat_compare_means(aes(x = sex_combined, y = shannon_diversity),
                     comparisons = list(c("Women", "Men")), label = "p.signif", method = "wilcox.test", hide.ns = FALSE) +
  facet_wrap(~Usage_group) +
  common_theme

# Combine and display ARG load and diversity figures by usage for both genders
combined_figure_age <- usage_arg_boxplot_both + usage_shannon_boxplot_both + plot_layout(ncol = 2)
print(combined_figure_age)
ggsave("final/usage_panel.png", combined_figure_age, width = 20, height = 10)

# Plot ARG load by usage for women only
usage_arg_boxplot_female <- ggplot(usage_filtered_metadata %>% filter(sex_combined == "Women"), 
                                   aes(x = Usage_group, y = log10_ARG_load, fill = "red")) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.5, alpha = 0.5, color = "salmon", show.legend = FALSE) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6, alpha = 0.7, show.legend = FALSE) +
  labs(x = "Antibiotic Usage", y = "Antibiotic Resistance Load") +
  ggtitle("Antibiotic Resistance Load vs. Antibiotic Usage for Women") +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c("11 or Below", "Above 11")), label = "p.signif", method = "wilcox.test", hide.ns = FALSE) +
  common_theme

# Plot Shannon diversity by usage for women only
usage_shannon_boxplot_female <- ggplot(usage_filtered_metadata %>% filter(sex_combined == "Women"),
                                       aes(x = Usage_group, y = shannon_diversity, fill = "red")) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.5, alpha = 0.5, color = "salmon", show.legend = FALSE) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6, alpha = 0.7, show.legend = FALSE) +
  labs(x = "Antibiotic Usage", y = "Resistome Diversity") +
  ggtitle("Resistome Diversity vs. Antibiotic Usage for Women") +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c("11 or Below", "Above 11")), label = "p.signif", method = "wilcox.test", hide.ns = FALSE) +
  common_theme

# Combine and display usage plots for women only
combined_figure_usage_female <- usage_arg_boxplot_female + usage_shannon_boxplot_female + plot_layout(ncol = 2)
print(combined_figure_usage_female)
ggsave("final/usage_female_panel.png", combined_figure_usage_female, width = 20, height = 10)
