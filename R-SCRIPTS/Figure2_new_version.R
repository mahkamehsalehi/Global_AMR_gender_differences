
# Load Necessary Libraries
library(mia)
library(ggplot2)
library(tidyverse)
library(vegan)
library(broom)
library(patchwork)
library(ggpubr)
library(rstatix)
library(viridis)
library(scales)

# Define colors for genders
male_color <- "#00BFC4"
female_color <- "#F8766D" 

# Load Your Data
TSE <- readRDS("../DATA/TSE.rds")
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
    !is.na(ARG_load)
  )

# ---------------------------
# 2. Prepare Data for Boxplot
# ---------------------------

df_boxplot_HIC <- df %>%
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
    World_Bank_Income_Group == "High income",  # Keep only "High income" rows
    !is.na(age_category),
    !is.na(ARG_load),
    !is.na(sex_combined)
  )

library(broom)

df_boxplot_HIC$age_category <- factor(
  df_boxplot_HIC$age_category,
  levels = age_labels,
  ordered = TRUE
)

# Perform Wilcoxon rank-sum tests for each age category
stat_tests <- df_boxplot_HIC %>%
  group_by(age_category) %>%
  do(tidy(wilcox.test(log(ARG_load) ~ sex_combined, data = .))) %>%
  ungroup() %>%
  mutate(
    # Adjust p-values with Benjamini-Hochberg correction
    p_adjusted = p.adjust(p.value, method = "BH"),
    # Annotate significance levels based on adjusted p-values
    p_signif = case_when(
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01  ~ "**",
      p_adjusted < 0.05  ~ "*",
      TRUE               ~ "ns"
    )
  )

# ---------------------------
# 3. Perform Statistical Tests for Boxplot
# ---------------------------





# Calculate overall maximum log10_ARG_load to set a common y_position
overall_max_log_ARG <- max(df_boxplot$ARG_load, na.rm = TRUE)

# ---------------------------
# 4. Create Scatter Plot
# ---------------------------

scatter_plot_HIC <- ggplot(df_scatter[df_scatter$World_Bank_Income_Group=="High income",], aes(x = host_age_years, y = log(ARG_load), color = sex_combined)) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  labs(
    x = "Age (years)",
    y = "log(ARG Load, RPKM)",
    color = "Gender"
  ) +
  theme_classic(base_size = 14, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12)
  )


scatter_plot_LMIC <- ggplot(df_scatter[!df_scatter$World_Bank_Income_Group=="High income",], aes(x = host_age_years, y = log(ARG_load), color = sex_combined)) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  labs(
    title = "Low and middle income",
    x = "Age (years)",
    y = "log(ARG Load, RPKM)",
    color = "Gender"
  ) +
  theme_classic(base_size = 14, base_family = "Times New Roman") +
  theme(
    legend.position = "none",
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
# Perform Wilcoxon rank-sum tests for each age category
stat_tests <- df_boxplot_HIC %>%
  group_by(age_category) %>%
  do(tidy(wilcox.test(log(ARG_load) ~ sex_combined, data = .))) %>%
  ungroup() %>%
  mutate(
    # Adjust p-values with Benjamini-Hochberg correction
    p_adjusted = p.adjust(p.value, method = "BH"),
    # Annotate significance levels based on adjusted p-values
    p_signif = case_when(
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01  ~ "**",
      p_adjusted < 0.05  ~ "*",
      TRUE               ~ "ns"
    )
  )
# Create Boxplot
boxplot_age_gender_load_HIC <- ggplot(df_boxplot_HIC, aes(x = age_category, y = log(ARG_load), fill = sex_combined)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    outliers = FALSE,
    width = 0.6,
    alpha = 0.8
  ) +
  scale_fill_manual(values = c(male_color, female_color)) +
  labs(
    x = "Age Category",
    y = "log(ARG Load, RPKM)",
    fill = "Gender"
  ) +
  theme_classic(base_size = 14, base_family = "Times New Roman") +
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
    aes(x = age_category, y = 8.5, label = p_signif),
    color = "black",
    vjust = 0,
    inherit.aes = FALSE
  )


df_boxplot_LMIC <- df %>%
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
    !World_Bank_Income_Group == "High income",  # Keep only "High income" rows
    !is.na(age_category),
    !is.na(ARG_load),
    !is.na(sex_combined)
  )


df_boxplot_LMIC$age_category <- factor(
  df_boxplot_LMIC$age_category,
  levels = age_labels,
  ordered = TRUE
)
# Perform Wilcoxon rank-sum tests for each age category
stat_tests <- df_boxplot_LMIC %>%
  group_by(age_category) %>%
  do(tidy(wilcox.test(log(ARG_load) ~ sex_combined, data = .))) %>%
  ungroup() %>%
  mutate(
    # Adjust p-values with Benjamini-Hochberg correction
    p_adjusted = p.adjust(p.value, method = "BH"),
    # Annotate significance levels based on adjusted p-values
    p_signif = case_when(
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01  ~ "**",
      p_adjusted < 0.05  ~ "*",
      TRUE               ~ "ns"
    )
  )

boxplot_age_gender_load_LMIC <-  ggplot(df_boxplot_LMIC, aes(x = age_category, y = log(ARG_load), fill = sex_combined)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    outliers = FALSE,
    width = 0.6,
    alpha = 0.8
  ) +
  scale_fill_manual(values = c(male_color, female_color)) +
  labs(
    x = "Age Category",
    y = "log(ARG Load, RPKM)",
    fill = "Gender"
  ) +
  theme_classic(base_size = 14, base_family = "Times New Roman") +
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
    aes(x = age_category, y = 8.5, label = p_signif),
    color = "black",
    vjust = 0,
    inherit.aes = FALSE
  )

# ---------------------------
# 6. Prepare simpson Diversity Scatter Plot
# ---------------------------

# Step 1: Apply important filtering
non_na_samples <- !is.na(colData(TSE)$sex_combined) &
  !is.na(colData(TSE)$host_age_years) &
  !is.na(colData(TSE)$World_Bank_Income_Group) 


# Subset the TSE object based on this filtering
TSE_filtered <- TSE[, non_na_samples]

TSE_filtered <- TSE[
  , !is.na(colData(TSE)$sex_combined) &
    !is.na(colData(TSE)$host_age_years) &
    !is.na(colData(TSE)$World_Bank_Income_Group) &
    colData(TSE)$World_Bank_Income_Group == "High income"
]


# Add  diversity for each sample
assay_data <- assay(TSE_filtered, "relabundance")
simpson_diversity <- as.data.frame(colData(TSE_filtered))$ARG_div_simp
shannon_diversity <- as.data.frame(colData(TSE_filtered))$ARG_div_shan
obs_diversity <- as.data.frame(colData(TSE_filtered)$ARG_obs)

# Merge  diversity with metadata
metadata <- as.data.frame(colData(TSE_filtered))
metadata$simpson_diversity <- simpson_diversity
metadata$shannon_diversity <- shannon_diversity

# Ensure GDP and age are available in the metadata
metadata <- metadata %>%
  mutate(
    host_age_years = as.numeric(host_age_years),
    GDP_per_head = as.numeric(GDP_per_head)
  ) %>%
  select(simpson_diversity, shannon_diversity, GDP_per_head, host_age_years, sex_combined, World_Bank_Income_Group) %>%
  filter(!is.na(GDP_per_head), !is.na(host_age_years))


simpson_plot <- ggplot(metadata, aes(x = host_age_years, y = simpson_diversity, color=sex_combined)) +
  geom_smooth(method = "loess", ) +
  labs(
    x = "Age (Years)",
    y = "ARG Diversity, Simpson"
  ) +
  theme_classic(base_size = 14, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12)
    # No legend in this plot
  )

shannon_plot <- ggplot(metadata, aes(x = host_age_years, y = shannon_diversity, color=sex_combined)) +
  geom_smooth(method = "loess", ) +
  labs(
    x = "Age (Years)",
    y = "ARG Diversity, Shannon"
  ) +
  theme_classic(base_size = 14, base_family = "Times New Roman") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12)
    # No legend in this plot
  )


# ---------------------------
# 7. Prepare Simpson Diversity Boxplot
# ---------------------------

# Categorize ages
metadata <- metadata %>%
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
  filter(!is.na(age_category), !is.na(simpson_diversity))

metadata$age_category <- factor(
  metadata$age_category,
  levels = age_labels,
  ordered = TRUE
)


### Shannon ###


# Perform Wilcoxon rank-sum tests for each age category
stat_tests <- df_boxplot_HIC %>%
  group_by(age_category) %>%
  do(tidy(wilcox.test(ARG_div_shan ~ sex_combined, data = .))) %>%
  ungroup() %>%
  mutate(
    # Adjust p-values with Benjamini-Hochberg correction
    p_adjusted = p.adjust(p.value, method = "BH"),
    # Annotate significance levels based on adjusted p-values
    p_signif = case_when(
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01  ~ "**",
      p_adjusted < 0.05  ~ "*",
      TRUE               ~ "ns"
    )
  )


# Create Boxplot
boxplot_age_gender_div_shan <- ggplot(metadata, aes(x = age_category, y = shannon_diversity, fill = sex_combined)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    outliers = FALSE,
    width = 0.6,
    alpha = 0.8
  ) +
  scale_fill_manual(values = c(female_color, male_color)) +
  labs(
    x = "Age Category",
    y = "ARG diversity, Shannon",
    fill = "Gender"
  ) +
  theme_classic(base_size = 14, base_family = "Times New Roman") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
    # Legend will be collected
  ) +
  # Add significance annotations
  geom_text(
    data = stat_tests,
    aes(x = age_category, y = 3.5, label = p_signif),
    color = "black",
    inherit.aes = FALSE
  )



# Perform Wilcoxon rank-sum tests for each age category
stat_tests <- df_boxplot_HIC %>%
  group_by(age_category) %>%
  do(tidy(wilcox.test(ARG_div_simp ~ sex_combined, data = .))) %>%
  ungroup() %>%
  mutate(
    # Adjust p-values with Benjamini-Hochberg correction
    p_adjusted = p.adjust(p.value, method = "BH"),
    # Annotate significance levels based on adjusted p-values
    p_signif = case_when(
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01  ~ "**",
      p_adjusted < 0.05  ~ "*",
      TRUE               ~ "ns"
    )
  )

# Create Boxplot
boxplot_age_gender_div_simp <- ggplot(metadata, aes(x = age_category, y = simpson_diversity, fill = sex_combined)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    outliers = FALSE,
    width = 0.6,
    alpha = 0.8
  ) +
  scale_fill_manual(values = c(female_color, male_color)) +
  labs(
    x = "Age Category",
    y = "ARG diversity, Simpson",
    fill = "Gender"
  ) +
  theme_classic(base_size = 14, base_family = "Times New Roman") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
    # Legend will be collected
  ) +
  # Add significance annotations
  geom_text(
    data = stat_tests,
    aes(x = age_category, y = 1, label = p_signif),
    color = "black",
    inherit.aes = FALSE
  )



# ---------------------------
# 8. Combine All Four Plots into One Figure
# ---------------------------

# # First Row: Scatter Plot and Boxplot with Gender Legend
# first_row <- (scatter_plot + boxplot_age_gender) +
#   plot_layout(guides = "collect") &   # Collect legends in this row
#   theme(
#     legend.position = "top",  # Position the legend at the top of the first row
#     legend.title = element_text(face = "bold", size = 14),
#     legend.text = element_text(size = 12)
#   )
# 
# # Second Row: simpson Scatter and Boxplot with Age Category Legend
# second_row <- (simpson_plot + simpson_boxplot) +
#   plot_layout(guides = "collect") &   # Collect legends in this row
#   theme(
#     legend.position = "none"
#   )
# 
# # Combine Rows
# combined_plot <- first_row / second_row +
#   plot_annotation(
#     tag_levels = 'a',
#     theme = theme(
#       plot.title = element_text(size = 20, face = "bold", hjust = 0.5, family = "Times New Roman")
#     )
#   )
# 
# # Display the combined plot
# print(combined_plot) 

# ---------------------------
# 9. Optional: Save the Combined Plot
# ---------------------------

ggsave("ARG_and_Simpson_Analysis_Four_Plots.png", combined_plot, width = 20, height = 20, units = "in", dpi = 300)

library(cowplot)
fig2 <- cowplot::plot_grid(scatter_plot_HIC+theme(legend.position = "none"), shannon_plot, boxplot_age_gender_load_HIC+theme(legend.position = "none"),  boxplot_age_gender_div_shan+theme(legend.position = "none"),  nrow = 2, labels="auto")

ggsave(file="Figure2_new.png", fig2, width = 10, height = 10, units = "in", dpi = 300)


fig2_LMIC <- cowplot::plot_grid(scatter_plot_LMIC, boxplot_age_gender_load_LMIC, labels="auto")
png(file="Figure2_LMIC.png", width = 500, height = 500)
print(fig2_LMIC)
dev.off()