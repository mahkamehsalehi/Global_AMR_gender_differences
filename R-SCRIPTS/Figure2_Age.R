
# Load Necessary Libraries
library(ggplot2)
library(tidyverse)
library(vegan)
library(broom)
library(patchwork)
library(ggpubr)
library(rstatix)
library(viridis)
library(scales)

# Load Your Data
TSE <- readRDS("TSE.rds")
df <- as.data.frame(colData(TSE))

# Define New Age Categories with Labels
age_labels <- c(
  "Infant",          # 0-1 year
  "Toddler",         # 1-3 years
  "Preschooler",     # 3-6 years
  "School-age",      # 6-12 years
  "Teen",            # 12-18 years
  "Young Adult",     # 18-35 years
  "Middle Adulthood",# 35-65 years
  "Older Adult"      # 65-100 years
)

# ---------------------------
# 1. Prepare Data for Scatter Plot
# ---------------------------

df_scatter <- df %>%
  mutate(
    host_age_years = as.numeric(host_age_years)
  ) %>%
  filter(
    !is.na(host_age_years),
    !is.na(sex_combined),
    !is.na(log10_ARG_load)
  )

# ---------------------------
# 2. Prepare Data for Boxplot
# ---------------------------

df_boxplot <- df %>%
  mutate(
    host_age_years = as.numeric(host_age_years),
    age_category = case_when(
      host_age_years >= 0    & host_age_years < 1    ~ "Infant",
      host_age_years >= 1    & host_age_years < 3    ~ "Toddler",
      host_age_years >= 3    & host_age_years < 6    ~ "Preschooler",
      host_age_years >= 6    & host_age_years < 12   ~ "School-age",
      host_age_years >= 12   & host_age_years < 18   ~ "Teen",
      host_age_years >= 18   & host_age_years < 35   ~ "Young Adult",
      host_age_years >= 35   & host_age_years < 65   ~ "Middle Adulthood",
      host_age_years >= 65   & host_age_years <= 100 ~ "Older Adult",
      TRUE ~ NA_character_
    ),
    sex_combined = factor(sex_combined, levels = c("male", "female"))
  ) %>%
  filter(
    !is.na(age_category),
    !is.na(log10_ARG_load),
    !is.na(sex_combined)
  )

df_boxplot$age_category <- factor(
  df_boxplot$age_category,
  levels = age_labels,
  ordered = TRUE
)

# ---------------------------
# 3. Perform Statistical Tests for Boxplot
# ---------------------------

# Perform independent t-tests for each age category
stat_tests <- df_boxplot %>%
  group_by(age_category) %>%
  do(tidy(t.test(log10_ARG_load ~ sex_combined, data = .))) %>%
  ungroup() %>%
  mutate(
    p_signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

# Calculate overall maximum log10_ARG_load to set a common y_position
overall_max_log_ARG <- max(df_boxplot$log10_ARG_load, na.rm = TRUE)

# Define a consistent y_position for all annotations
annotation_offset <- 0.2  # Adjust as needed based on your data's scale
stat_tests <- stat_tests %>%
  mutate(
    y_position = overall_max_log_ARG + annotation_offset
  )

# ---------------------------
# 4. Create Scatter Plot
# ---------------------------

scatter_plot <- ggplot(df_scatter, aes(x = host_age_years, y = log10_ARG_load, color = sex_combined)) +
  geom_point(alpha = 0.6, size = 0.5) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  labs(
    title = "log10 ARG Load vs. Age by Gender",
    x = "Age (years)",
    y = "log10 ARG Load",
    color = "Gender"
  ) +
  theme_minimal(base_size = 14, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12)
  )

# ---------------------------
# 5. Create Boxplot with Significance Annotations
# ---------------------------

# Define colors for genders
male_color <- "#00BFC4"
female_color <- "#F8766D" 

# Create Boxplot
boxplot_age_gender <- ggplot(df_boxplot, aes(x = age_category, y = log10_ARG_load, fill = sex_combined)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    outlier.size = 1,
    width = 0.6,
    alpha = 0.8
  ) +
  scale_fill_manual(values = c(male_color, female_color)) +
  labs(
    title = "log10 ARG Load by Age Category and Gender",
    x = "Age Category",
    y = "log10 ARG Load",
    fill = "Gender"
  ) +
  theme_minimal(base_size = 14, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
    # Legend will be collected
  ) +
  # Add significance annotations
  geom_text(
    data = stat_tests,
    aes(x = age_category, y = y_position, label = p_signif),
    color = "black",
    size = 6,
    vjust = 0,
    inherit.aes = FALSE
  )

# ---------------------------
# 6. Prepare Shannon Diversity Scatter Plot
# ---------------------------

# Step 1: Apply important filtering
non_na_samples <- !is.na(colData(TSE)$sex_combined) &
  !is.na(colData(TSE)$host_age_years)

# Subset the TSE object based on this filtering
TSE_filtered <- TSE[, non_na_samples]

# Step 2: Further filter the data for females only
df <- as.data.frame(colData(TSE_filtered))
female_samples <- df %>%
  filter(sex_combined == "female")

# Subset the TSE object to include only female samples after the filtering
TSE_female <- TSE_filtered[, colData(TSE_filtered)$sex_combined == "female"]

# Calculate Shannon diversity for each sample
assay_data_female <- assay(TSE_female, "relabundance")
shannon_diversity <- diversity(t(assay_data_female), index = "shannon")

# Merge Shannon diversity with metadata
female_metadata <- as.data.frame(colData(TSE_female))
female_metadata$shannon_diversity <- shannon_diversity

# Step 2: Further filter the data for females only
df <- as.data.frame(colData(TSE_filtered))
female_samples <- df %>%
  filter(sex_combined == "female")

# Ensure GDP and age are available in the metadata
female_metadata <- female_metadata %>%
  mutate(
    host_age_years = as.numeric(host_age_years),
    GDP_per_head = as.numeric(GDP_per_head)
  ) %>%
  select(shannon_diversity, GDP_per_head, host_age_years) %>%
  filter(!is.na(GDP_per_head), !is.na(host_age_years))

# Create Shannon Diversity Scatter Plot
female_color <- "#F8766D" 

shannon_plot <- ggplot(female_metadata, aes(x = host_age_years, y = shannon_diversity)) +
  geom_point(alpha = 0.5, size = 0.5, color = female_color) +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  labs(
    title = "Shannon Diversity vs. Age Among Women",
    x = "Age (Years)",
    y = "Shannon Diversity"
  ) +
  theme_minimal(base_size = 14, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12)
    # No legend in this plot
  )

# ---------------------------
# 7. Prepare Shannon Diversity Boxplot
# ---------------------------

# Categorize ages
female_metadata <- female_metadata %>%
  mutate(
    age_category = case_when(
      host_age_years >= 0    & host_age_years < 1    ~ "Infant",
      host_age_years >= 1    & host_age_years < 3    ~ "Toddler",
      host_age_years >= 3    & host_age_years < 6    ~ "Preschooler",
      host_age_years >= 6    & host_age_years < 12   ~ "School-age",
      host_age_years >= 12   & host_age_years < 18   ~ "Teen",
      host_age_years >= 18   & host_age_years < 35   ~ "Young Adult",
      host_age_years >= 35   & host_age_years < 65   ~ "Middle Adulthood",
      host_age_years >= 65   & host_age_years <= 100 ~ "Older Adult",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(age_category), !is.na(shannon_diversity))

female_metadata$age_category <- factor(
  female_metadata$age_category,
  levels = age_labels,
  ordered = TRUE
)

# Compute Statistical Tests for Boxplot
# Perform pairwise Wilcoxon tests between adjacent age categories
adjacent_pairs <- data.frame(
  group1 = age_labels[-length(age_labels)],
  group2 = age_labels[-1]
)

pairwise_tests <- female_metadata %>%
  pairwise_wilcox_test(
    shannon_diversity ~ age_category,
    p.adjust.method = "BH"
  )

# Filter to only adjacent comparisons
pairwise_tests_adjacent <- pairwise_tests %>%
  inner_join(adjacent_pairs, by = c("group1", "group2"))

# Create significance labels
pairwise_tests_adjacent <- pairwise_tests_adjacent %>%
  mutate(
    p_signif = case_when(
      p.adj <= 0.001 ~ "***",
      p.adj <= 0.01  ~ "**",
      p.adj <= 0.05  ~ "*",
      TRUE           ~ "ns"
    )
  )

# Set a fixed y.position for all annotations
max_shannon <- max(female_metadata$shannon_diversity, na.rm = TRUE)
annotation_y_position <- max_shannon + 0.5  # Adjust as needed
pairwise_tests_adjacent <- pairwise_tests_adjacent %>%
  mutate(y.position = annotation_y_position)

# Define colors for age categories
age_colors <- viridis(length(age_labels))

# Create Shannon Diversity Boxplot
shannon_boxplot <- ggplot(female_metadata, aes(x = age_category, y = shannon_diversity)) +
  geom_boxplot(outlier.size = 1, width = 0.6, alpha = 0.8) +
  labs(
    title = "Shannon Diversity by Age Category Among Women",
    x = "Age Category",
    y = "Shannon Diversity",
    fill = "Age Category"
  ) +
  theme_minimal(base_size = 14, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  ) +
  # Add significance annotations
  stat_pvalue_manual(
    pairwise_tests_adjacent,
    label = "p_signif",
    y.position = "y.position",
    tip.length = 0.01,
    step.increase = 0,
    inherit.aes = FALSE
  )

# ---------------------------
# 8. Combine All Four Plots into One Figure
# ---------------------------

# First Row: Scatter Plot and Boxplot with Gender Legend
first_row <- (scatter_plot + boxplot_age_gender) +
  plot_layout(guides = "collect") &   # Collect legends in this row
  theme(
    legend.position = "top",  # Position the legend at the top of the first row
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12)
  )

# Second Row: Shannon Scatter and Boxplot with Age Category Legend
second_row <- (shannon_plot + shannon_boxplot) +
  plot_layout(guides = "collect") &   # Collect legends in this row
  theme(
    legend.position = "none"
  )

# Combine Rows
combined_plot <- first_row / second_row +
  plot_annotation(
    title = "ARG Load and Shannon Diversity Analysis by Age and Gender",
    tag_levels = 'a',
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, family = "Times New Roman")
    )
  )

# Display the combined plot
print(combined_plot) 

# ---------------------------
# 9. Optional: Save the Combined Plot
# ---------------------------

ggsave("ARG_and_Shannon_Analysis_Four_Plots.png", combined_plot, width = 20, height = 20, units = "in", dpi = 300)
