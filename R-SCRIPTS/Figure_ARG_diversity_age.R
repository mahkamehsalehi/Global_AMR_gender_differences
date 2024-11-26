library(vegan)
library(ggpubr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(tidyverse)
library(cowplot)


# ---------------------------
# Data Loading and Preprocessing
# ---------------------------

TSE <- readRDS("DATA/TSE_filtered.rds")

tse_metadata <- as.data.frame(colData(TSE))

metadata_hic <- tse_metadata %>%
  filter(income_group == "HIC") %>%
  filter(!is.na(age_category) & !is.na(log_ARG_load))

metadata_lmic <- tse_metadata %>%
  filter(income_group == "LMIC") %>%
  filter(!is.na(age_category) & !is.na(log_ARG_load))

# ---------------------------
# Define Custom Plot Theme
# ---------------------------

comparisons <- list(
  c("Infant", "Toddler"),
  c("Toddler", "Child"),
  c("Child", "Young Adult"),
  c("Young Adult", "Middle Adult"),
  c("Middle Adult", "Older Adult")
)

# ---------------------------
# Visualization: Boxplots for HIC
# ---------------------------

# ARG Load Boxplot Without Faceting HIC
age_arg_boxplot_hic <- ggplot(metadata_hic, aes(x = gender, y = log_ARG_load, fill = gender)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.6,
    alpha = 1,
    show.legend = FALSE
  ) +
  facet_wrap(~age_category, scales = "fixed", nrow = 1) +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(x = "Gender", y = "ARG load (natural log)") +
  theme_minimal() +
  stat_compare_means(
    comparisons = list(
      c("Women", "Men")
    ),
    aes(label = ..p.signif..),
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_blank(),
    strip.text = element_text(size = 10, face = "bold", angle = 45),
    panel.spacing = unit(0.5, "lines")
  )


# Shannon Diversity Boxplot Without Faceting
age_shannon_boxplot_hic <- ggplot(metadata_hic, aes(x = gender, y = shannon_diversity, fill = gender)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.6,
    alpha = 1,
    show.legend = FALSE
  ) +
  facet_wrap(~age_category, scales = "fixed", nrow = 1) +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(x = "Gender", y = "Resistome diversity") +
  theme_minimal() +
  stat_compare_means(
    comparisons = list(
      c("Women", "Men")
    ),
    aes(label = ..p.signif..),
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_blank(),
    strip.text = element_text(size = 10, face = "bold", angle = 45),
    panel.spacing = unit(0.5, "lines")
  )


# Scatterplot for ARG Load vs. Age
age_arg_scatterplot_hic <- ggplot(metadata_hic, aes(x = host_age_years, y = log_ARG_load, color = gender)) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  labs(
    x = "Age (years)",
    y = "ARG load (natural log)",
    color = "Gender"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"))

# Scatterplot for Shannon Diversity vs. Age
age_shannon_scatterplot_hic <- ggplot(metadata_hic, aes(x = host_age_years, y = shannon_diversity, color = gender)) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  labs(
    x = "Age (years)",
    y = "Resistome diversity",
    color = "Gender"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black")
  )

# ---------------------------
# Visualization: Boxplots for LMIC
# ---------------------------

# ARG Load Boxplot Without Faceting LMIC
age_arg_boxplot_lmic <- ggplot(metadata_lmic, aes(x = gender, y = log_ARG_load, fill = gender)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.6,
    alpha = 1,
    show.legend = FALSE
  ) +
  facet_wrap(~age_category, scales = "fixed", nrow = 1) +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(x = "Gender", y = "ARG load (natural log)") +
  theme_minimal() +
  stat_compare_means(
    comparisons = list(
      c("Women", "Men")
    ),
    aes(label = ..p.signif..),
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 10, face = "bold", angle = 45),
    panel.spacing = unit(0.5, "lines")
  )


# Shannon Diversity Boxplot
age_shannon_boxplot_lmic <-ggplot(metadata_lmic, aes(x = gender, y = shannon_diversity, fill = gender)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.6,
    alpha = 1,
    show.legend = FALSE
  ) +
  facet_wrap(~age_category, scales = "fixed", nrow = 1) +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(x = "Gender", y = "Resistome diversity") +
  theme_minimal() +
  stat_compare_means(
    comparisons = list(
      c("Women", "Men")
    ),
    aes(label = ..p.signif..),
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_blank(),
    strip.text = element_text(size = 10, face = "bold", angle = 45),
    panel.spacing = unit(0.5, "lines")
  )

# Scatterplot for ARG Load vs. Age
age_arg_scatterplot_lmic <- ggplot(metadata_lmic, aes(x = host_age_years, y = log_ARG_load, color = gender)) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  labs(
    x = "Age (years)",
    y = "ARG load (natural log)",
    color = "Gender"
  ) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"))

# Scatterplot for Shannon Diversity vs. Age
age_shannon_scatterplot_lmic  <- ggplot(metadata_lmic, aes(x = host_age_years, y = shannon_diversity, color = gender)) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  labs(
    x = "Age (years)",
    y = "Resistome diversity",
    color = "Gender"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black")
  )

# ---------------------------
# Visualization: Combined Scatterplots and Boxplots for Women in HIC
# ---------------------------

# Filter data for women only
women_data_hic <- metadata_hic %>% filter(gender == "Women")

# ARG Load Boxplot for Women
age_arg_boxplot_female_hic <-ggplot(women_data_hic, aes(x = age_category, y = log_ARG_load, fill = "#F8766D")) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.6,
    alpha = 1,
    show.legend = FALSE
  ) +
  labs(x = "Age category", y = "ARG load (natural log)") +
  theme_minimal() +
  theme_minimal() +
  stat_compare_means(
    comparisons = comparisons,
    aes(label = ..p.signif..),
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
  

# Shannon Diversity Boxplot for Women
age_shannon_boxplot_female_hic <- ggplot(women_data_hic, 
                                 aes(x = age_category, y = shannon_diversity, fill = "#F8766D")) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.6,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  labs(x = "Age Category", y = "Resistome diversity") +
  theme_minimal() +
  stat_compare_means(
    comparisons = comparisons,
    aes(label = ..p.signif..),
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ARG Load Scatterplot for Women
age_arg_scatterplot_female_hic <- ggplot(women_data_hic, aes(x = host_age_years, y = log_ARG_load, color = gender)) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(
    x = "Age (years)",
    y = "ARG load (natural log)",
    color = "Gender"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black")
  )

# Shannon Diversity Scatterplot for Women
age_shannon_scatterplot_female_hiv <- ggplot(women_data_hic, aes(x = host_age_years, y = shannon_diversity, color = gender)) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(
    x = "Age (years)",
    y = "Resistome diversity",
    color = "Gender"
  ) +
  theme_minimal() +
  theme(
  axis.line = element_line(color = "black")
) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black")
  )

# ---------------------------
# Visualization: Combined Scatterplots and Boxplots for Women in LMIC
# ---------------------------

# Filter data for women only
women_data_lmic <- metadata_lmic %>% filter(gender == "Women")

# ARG Load Boxplot for Women
age_arg_boxplot_female_lmic <- ggplot(women_data_lmic, aes(x = age_category, y = log_ARG_load, fill = "#F8766D")) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.6,
    alpha = 1,
    show.legend = FALSE
  ) +
  labs(x = "Age category", y = "ARG load (natural log)") +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(
    comparisons = comparisons,
    aes(label = ..p.signif..),
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# Shannon Diversity Boxplot for Women
age_shannon_boxplot_female_lmic <- ggplot(women_data_lmic, 
                                     aes(x = age_category, y = shannon_diversity, fill = "#F8766D")) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.6,
    alpha = 1,
    show.legend = FALSE
  ) +
  labs(x = "Age category", y = "Resistome diversity") +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(
    comparisons = comparisons,
    aes(label = ..p.signif..),
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ARG Load Scatterplot for Women
age_arg_scatterplot_female_lmic <- ggplot(women_data_lmic, aes(x = host_age_years, y = log_ARG_load, color = gender)) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(
    x = "Age (years)",
    y = "ARG load (natural log)",
    color = "Gender"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black")
  )

# Shannon Diversity Scatterplot for Women
age_shannon_scatterplot_female_lmic <- ggplot(women_data_lmic, aes(x = host_age_years, y = shannon_diversity, color = gender)) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(
    x = "Age (years)",
    y = "Resistome diversity",
    color = "Gender"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black")
  )

# Combine Scatterplots and Boxplots 
library(cowplot)

# Combine the scatterplot and boxplot
combined_arg_age_hic <- plot_grid(
  age_arg_scatterplot_hic, 
  age_arg_boxplot_hic, 
  labels = c("a", "b"),
  ncol = 2
)

# Save the figure
ggsave("RESULTS/FIGURES/arg_age_hic.png", combined_arg_age_hic, width = 14, height = 8)



# Combine Scatterplots and Boxplots 
combined_arg_age_lmic <- plot_grid(
  age_arg_scatterplot_lmic, 
  age_arg_boxplot_lmic, 
  labels = c("a", "b"), 
  ncol = 2
)

# Save the combined figure
ggsave("RESULTS/FIGURES/arg_age_lmic.png", combined_arg_age_lmic, width = 14, height = 8)


# Combine Scatterplots and Boxplots 
combined_shannon_age_hic <- plot_grid(
  age_shannon_scatterplot_hic, 
  age_shannon_boxplot_hic, 
  labels = c("a", "b"), 
  ncol = 2
)

# Save the combined figure
ggsave("RESULTS/FIGURES/shannon_age_hic.png", combined_shannon_age_hic, width = 14, height = 8)


# Combine Scatterplots and Boxplots 
combined_shannon_age_lmic <- plot_grid(
  age_shannon_scatterplot_lmic, 
  age_shannon_boxplot_lmic, 
  labels = c("a", "b"), 
  ncol = 2
)

# Save the combined figure
ggsave("RESULTS/FIGURES/shannon_age_lmic.png", combined_shannon_age_lmic, width = 14, height = 8)

# ---------------------------
# Statistical Analysis: Wilcoxon Tests
# ---------------------------

# Perform Wilcoxon tests for ARG Load
arg_load_results <- tse_metadata %>%
  group_by(age_category) %>%
  wilcox_test(log_ARG_load ~ gender) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(variable = "ARG Load") %>%
  select(age_category, variable, group1, group2, n1, n2, statistic, p, p.adj, p.adj.signif)

# Perform Wilcoxon tests for Shannon Diversity
shannon_div_results <- tse_metadata %>%
  group_by(age_category) %>%
  wilcox_test(shannon_diversity ~ gender) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(variable = "Shannon Diversity") %>%
  select(age_category, variable, group1, group2, n1, n2, statistic, p, p.adj, p.adj.signif)

# Combine the results into one table
result_table <- bind_rows(arg_load_results, shannon_div_results) %>%
  mutate(
    p = signif(p, digits = 3),
    p.adj = signif(p.adj, digits = 3),
    statistic = round(statistic, 2)
  ) %>%
  arrange(variable, age_category)

# View the formatted result table
print(result_table)

# ---------------------------
# Statistical Analysis: Pairwise Wilcoxon Tests for Women
# ---------------------------

# Perform pairwise Wilcoxon tests for ARG Load in women
pairwise_arg_load <- women_data %>%
  pairwise_wilcox_test(
    log_ARG_load ~ age_category,
    p.adjust.method = "BH"
  ) %>%
  mutate(variable = "ARG Load")

# Perform pairwise Wilcoxon tests for Shannon Diversity in women
pairwise_shannon <- women_data %>%
  pairwise_wilcox_test(
    shannon_diversity ~ age_category,
    p.adjust.method = "BH"
  ) %>%
  mutate(variable = "Shannon Diversity")

# Combine the pairwise results into one table
pairwise_results <- bind_rows(pairwise_arg_load, pairwise_shannon) %>%
  select(variable, group1, group2, n1, n2, p, p.adj, p.adj.signif) %>%
  mutate(
    p = signif(p, digits = 3),
    p.adj = signif(p.adj, digits = 3)
  ) %>%
  arrange(variable, group1, group2)

# View the formatted pairwise results table
print(pairwise_results)