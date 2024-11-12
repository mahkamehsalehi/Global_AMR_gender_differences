library(vegan)
library(ggpubr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(dplyr)

TSE <- readRDS("DATA/TSE_gender_age.rds")

TSE_filtered <- TSE[, colData(TSE)$geo_loc_name_country_continent_calc != "uncalculated"]

# Step 2: Extract metadata and set factors for consistent ordering in plots
tse_metadata <- as.data.frame(colData(TSE_filtered))

tse_metadata <- tse_metadata %>%
  mutate(sex_combined = case_when(
    sex_combined == "male" ~ "Men",
    sex_combined == "female" ~ "Women",
    TRUE ~ sex_combined
  ))

tse_metadata$sex_combined <- factor(tse_metadata$sex_combined, levels = c("Women", "Men"))

tse_metadata <- tse_metadata %>%
  mutate(
    age_category = case_when(
      host_age_years >= 0 & host_age_years <= 1 ~ "Infant",
      host_age_years > 1 & host_age_years <= 3 ~ "Toddler",
      host_age_years > 3 & host_age_years <= 18 ~ "Child",
      host_age_years > 18 & host_age_years <= 35 ~ "Young Adult",
      host_age_years > 35 & host_age_years <= 65 ~ "Middle-Aged Adult",
      host_age_years > 65 & host_age_years <= 100 ~ "Older Adult",
      TRUE ~ NA_character_
    )
  )

# Step 3: Extract and clean assay data (relative abundances)
assay_data_clean <- assay(TSE_filtered, "relabundance")

# Remove samples (columns) with all-zero counts
non_empty_samples <- colSums(assay_data_clean > 0, na.rm = TRUE) > 0
assay_data_clean <- assay_data_clean[, non_empty_samples]

# Subset colData to match non-empty samples
tse_metadata <- tse_metadata[non_empty_samples, ]

# Step 4: Calculate the Bray-Curtis distance matrix on the cleaned data
distance_matrix <- vegan::vegdist(t(assay_data_clean), method = "bray")

# Calculate Shannon diversity
diversity_indices <- vegan::diversity(t(assay_data_clean), index = "shannon")
tse_metadata$shannon_diversity <- diversity_indices


# Define custom plot theme
common_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "none",
    axis.line = element_line(color = "black"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 12, face = "bold")
  )

tse_metadata$age_category <- factor(tse_metadata$age_category, 
                                    levels = c("Infant",
                                               "Toddler",
                                               "Child",
                                               "Young Adult",
                                               "Middle-Aged Adult",
                                               "Older Adult"))

age_arg_boxplot <- ggplot(tse_metadata, aes(x = sex_combined, y = log10_ARG_load, fill = sex_combined)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    size = 0.5,
    alpha = 0.5,
    aes(color = sex_combined),
    show.legend = FALSE
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.6,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(x = "Gender", y = "Antibiotic Resistance Gene Load") +
  theme_minimal() +
  stat_compare_means(
    aes(x = sex_combined, y = log10_ARG_load),
    comparisons = list(c("Women", "Men")),
    label = "p.signif",
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE
  ) +
  facet_wrap(~age_category) +
  common_theme



age_shannon_boxplot <- ggplot(tse_metadata, aes(x = sex_combined, y = shannon_diversity, fill = sex_combined)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    size = 0.5,
    alpha = 0.5,
    aes(color = sex_combined),
    show.legend = FALSE
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.6,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(x = "Gender", y = "Resistome Diversity") +
  theme_minimal(base_family = "Arial") +
  stat_compare_means(
    aes(x = sex_combined, y = shannon_diversity),
    comparisons = list(c("Women", "Men")),
    label = "p.signif",
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE
  ) +
  facet_wrap(~age_category) +
  common_theme



combined_figure_age <- age_arg_boxplot + age_shannon_boxplot + 
  plot_layout(ncol = 1)
print(combined_figure_age)

ggsave("RESULTS/FIGURES/age_panel.png", combined_figure_age, width = 12, height = 8)

#---------------------------
# No faceted boxplot
#---------------------------

age_arg_boxplot_no_facet <- ggplot(tse_metadata, aes(x = age_category, y = log10_ARG_load, fill = sex_combined)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    size = 0.5,
    alpha = 0.5,
    aes(color = sex_combined),
    show.legend = FALSE
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.6,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(x = "Gender", y = "log10 ARG Load") +
  ggtitle("Antibiotic Resistance Genes Load Across Age Categories") +
  theme_minimal(base_family = "Arial") +
  stat_compare_means(
    aes(x = sex_combined, y = log10_ARG_load),
    comparisons = list(c("Child", "Young Adult"),
                       c("Young Adult", "Middle Adulthood"),
                       c("Middle Adulthood", "Older Adult")), 
    label = "p.signif",
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE
  ) +
  common_theme

age_shannon_boxplot_no_facet <- ggplot(tse_metadata, aes(x = age_category, y = shannon_diversity, fill = sex_combined)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    size = 0.5,
    alpha = 0.5,
    aes(color = sex_combined),
    show.legend = FALSE
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.6,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(x = "Gender", y = "Resistome Diversity") +
  ggtitle(" Resistome Diversity Across Age Categories") +
  theme_minimal(base_family = "Arial") +
  stat_compare_means(
    aes(x = sex_combined, y = shannon_diversity),
    comparisons = list(c("Child", "Young Adult"),
                       c("Young Adult", "Middle Adulthood"),
                       c("Middle Adulthood", "Older Adult")), 
    label = "p.signif",
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE
  ) +
  common_theme

age_arg_scatterplot <- ggplot(tse_metadata, aes(x = host_age_years, y = log10_ARG_load, color = sex_combined)) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  geom_point(alpha = 0.6, size = 0.5) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  labs(
    x = "Age (years)",
    y = "Antibiotic Resistance Load",
    color = "Gender"
  ) +
  theme_minimal() +
  common_theme

age_shannon_scatterplot <- ggplot(tse_metadata, aes(x = host_age_years, y = shannon_diversity, color = sex_combined)) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  geom_point(alpha = 0.6, size = 0.5) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  labs(
    x = "Age (years)",
    y = "Resistome Diversity",
    color = "Gender"
  ) +
  theme_minimal() +
  common_theme

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
age_arg_boxplot_female <- ggplot(tse_metadata %>%
                                   filter(sex_combined == "Women"), aes(x = age_category, y = log10_ARG_load, fill = "#F8766D")) +
  geom_jitter(
    position = position_jitter(width = 0.2, height = 0),
    size = .5,
    alpha = 0.5,
    color = "#F8766D",
    show.legend = FALSE
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.6,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  labs(x = "Age Category", y = "Antibiotic Resistance Load") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("Middle-Aged Adult" = "Middle-Age")) +
  stat_compare_means(
    comparisons = list(
      c("Infant", "Toddler"),
      c("Infant", "Older Adult"),
      c("Child", "Young Adult"),
      c("Young Adult", "Middle-Aged Adult"),
      c("Middle-Aged Adult", "Older Adult")
      ), 
    label = "p.signif",
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE
  ) +
  common_theme

age_shannon_boxplot_female <- ggplot(tse_metadata %>%
                                   filter(sex_combined == "Women"), 
                                 aes(x = age_category, y = shannon_diversity, fill = "#F8766D")) +
  geom_jitter(
    position = position_jitter(width = 0.2, height = 0),
    size = .5,
    alpha = 0.5,
    color = "#F8766D",
    show.legend = FALSE
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.6,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  labs(x = "Age Category", y = "Resistome Diversity") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("Middle-Aged Adult" = "Middle-Age")) +
  stat_compare_means(
    comparisons = list(
      c("Toddler", "Child"),
      c("Infant", "Older Adult"),
      c("Child", "Young Adult"),
      c("Young Adult", "Middle-Aged Adult"),
      c("Middle-Aged Adult", "Older Adult")
    ),
    label = "p.signif",
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE
  ) +
  common_theme

age_arg_scatterplot_female <- ggplot(tse_metadata %>%
                                       filter(sex_combined == "Women"), aes(x = host_age_years, y = log10_ARG_load, color = sex_combined)) +
  geom_point(alpha = 0.6, size = 0.5) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(
    x = "Age (years)",
    y = "Antibiotic Resistance Load",
    color = "Gender"
  ) +
  theme_minimal() +
  common_theme

age_shannon_scatterplot_female <- ggplot(tse_metadata%>%
                                           filter(sex_combined == "Women"), aes(x = host_age_years, y = shannon_diversity, color = sex_combined)) +
  geom_point(alpha = 0.6, size = 0.5) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(
    x = "Age (years)",
    y = "Resistome Diversity",
    color = "Gender"
  ) +
  theme_minimal() +
  common_theme

combined_figure <- (age_arg_scatterplot + age_arg_boxplot_female + 
                                         age_shannon_scatterplot + age_shannon_boxplot_female) +
  plot_layout(ncol = 2, nrow = 2) + 
  plot_annotation(
    tag_levels = 'a')

print(combined_figure)

ggsave("RESULTS/FIGURES/age_4_panel.pdf", combined_figure, width = 12, height = 8)


# Calculate p-values for each age category
p_values_arg <- tse_metadata %>%
  group_by(age_category) %>%
  summarize(p_value = wilcox.test(log10_ARG_load ~ sex_combined)$p.value) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  mutate(significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*"  ))

# Create the density plot with annotations
arg_age_densityplot <- ggplot(tse_metadata, aes(x = log10_ARG_load, fill = sex_combined)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ age_category, nrow = 1) +
  labs(
    x = "Antibiotic Resistance Load",
    y = "Density" ) +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  theme(legend.position = "top") +
  geom_text(data = p_values, aes(x = 4.5, y = 1.2, label = significance), inherit.aes = FALSE)+
  common_theme


# Calculate p-values for each age category
p_values_shannon <- tse_metadata %>%
  group_by(age_category) %>%
  summarize(p_value = wilcox.test(shannon_diversity ~ sex_combined)$p.value) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  mutate(significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*"  ))

shannon_age_densityplot <- ggplot(tse_metadata, aes(x = shannon_diversity, fill = sex_combined)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ age_category, nrow = 1) +
  labs(
    x = "Resistome Diversity",
    y = "Density" ) +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  theme(legend.position = "top") +
  geom_text(data = p_values_shannon, aes(x = 4.5, y = 1.2, label = significance), inherit.aes = FALSE)+
  common_theme

combined_figure_density_age <- arg_age_densityplot + shannon_age_densityplot +
                                  plot_layout(ncol = 1, nrow = 2) + 
                                  plot_annotation(
                                    tag_levels = 'a')
                                

print(combined_figure_density_age)

ggsave("RESULTS/FIGURES/age_2_panel.png", combined_figure_density_age, width = 10, height = 6)
