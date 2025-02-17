# Load libraries
library(vegan)
library(ggpubr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(cowplot)
library(mia)
library(scales)
library(rstatix)

# Data Loading and Processing
TSE <- readRDS("DATA/TSE_filtered.rds")
tse_metadata <- as.data.frame(colData(TSE))
filtered_metadata <- tse_metadata %>%
  filter(
    category != "Infant Study", 
    geo_loc_name_country_calc != "Zimbabwe"
  )
filtered_metadata <- filtered_metadata %>%
  drop_na(log_ARG_load, income_group)
filtered_metadata_no_outliers <- filtered_metadata
filtered_metadata_female_no_outliers <- filtered_metadata_no_outliers %>%
  filter(gender == "Women")

# Define color palette and theme
gender_colors <- c("Women" = "#F8766D", "Men" = "#619CFF")
custom_theme <- theme_minimal(18) +
  theme(
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", size = 1),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )

# Create violin plot for ARG load
income_arg_violin_female <- ggplot(filtered_metadata_female_no_outliers, 
                                   aes(x = income_group, 
                                       y = ARG_load, 
                                       fill = income_group)) +
  geom_violin(trim = FALSE, alpha = 1, color = "black") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_compare_means(comparisons = list(c("HIC", "LMIC")), 
                     label = "p.format", 
                     method = "wilcox.test", 
                     p.adjust.method = "BH", 
                     hide.ns = FALSE,
                     size = 5) +
  scale_fill_manual(values = c("HIC" = "#F8766D", "LMIC" = "#F8766D")) +
  labs(x = "Income group", 
       y = "ARG load (RPKM)") +
  scale_y_continuous(trans = "log10",
                     breaks = 10^(2:5),
                     labels = trans_format("log10", math_format(10^.x))) +
  custom_theme

# Create violin plot for Shannon diversity
income_shannon_violin_female <- ggplot(filtered_metadata_female_no_outliers, 
                                       aes(x = income_group, 
                                           y = shannon_diversity, 
                                           fill = income_group)) +
  geom_violin(trim = FALSE, alpha = 1, color = "black") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_compare_means(comparisons = list(c("HIC", "LMIC")), 
                     label = "p.format", 
                     method = "wilcox.test", 
                     p.adjust.method = "BH", 
                     hide.ns = FALSE,
                     size = 5) +
  scale_fill_manual(values = c("HIC" = "#F8766D", "LMIC" = "#F8766D")) +
  labs(x = "Income group", 
       y = "ARG diversity (Shannon index)") +
  coord_cartesian(ylim = c(0, 4)) +
  custom_theme

# Function to calculate statistics
calc_stats <- function(data, variable) {
  # Get sample sizes
  n_sizes <- data %>%
    group_by(income_group) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = income_group, values_from = n, names_prefix = "n_")
  
  data$income_group <- factor(data$income_group, levels = c("HIC", "LMIC"))
  
  # Calculate effect size and stats
  if(variable == "ARG_load") {
    eff <- wilcox_effsize(data = data, 
                          ARG_load ~ income_group, 
                          ci = TRUE, 
                          conf.level = 0.95)
    test <- wilcox.test(ARG_load ~ income_group, data = data)
  } else {
    eff <- wilcox_effsize(data = data, 
                          shannon_diversity ~ income_group, 
                          ci = TRUE, 
                          conf.level = 0.95)
    test <- wilcox.test(shannon_diversity ~ income_group, data = data)
  }
  
  results <- data.frame(
    "Metric" = if(variable == "ARG_load") "ARG load" else "ARG diversity",
    "N (HIC)" = n_sizes$n_HIC,
    "N (LMIC)" = n_sizes$n_LMIC,
    "Effect Size (r)" = round(eff$effsize, 3),
    "Lower 95% CI" = round(eff$conf.low, 3),
    "Upper 95% CI" = round(eff$conf.high, 3),
    "Adjusted p-value" = formatC(p.adjust(test$p.value, method = "BH"), 
                                 format = "f", digits = 4),
    check.names = FALSE
  )
  
  return(results)
}

# Calculate statistics for both metrics
arg_stats <- calc_stats(filtered_metadata_female_no_outliers, "ARG_load")
shannon_stats <- calc_stats(filtered_metadata_female_no_outliers, "shannon_diversity")

# Combine results
combined_stats <- rbind(arg_stats, shannon_stats)

# Create formatted table
stats_table <- ggtexttable(combined_stats, 
                           rows = NULL,
                           theme = ttheme("light", 
                                          base_size = 16, 
                                          padding = unit(c(10, 20), "pt"))) %>%
  tab_add_title(text = "c", size = 18, face = "bold", just = "left", 
                padding = unit(c(0, 0, 0, 4), "pt"))

# Combine violin plots (panels a and b)
combined_figure_income_female <- plot_grid(
  income_arg_violin_female + custom_theme,
  income_shannon_violin_female + custom_theme,
  labels = c('a', 'b'),
  label_size = 20,
  ncol = 2,
  align = 'hv',
  rel_widths = c(1, 1),
  rel_heights = c(1, 1)
)

# Create separator
separator <- ggdraw() + 
  draw_line(x = c(0, 1), y = c(0.5, 0.5), color = "grey", size = 2)

# Combine all panels
final_figure <- combined_figure_income_female / 
  separator / 
  stats_table + 
  plot_layout(heights = c(2, 0.01, NA))

# Save final figure
ggsave("RESULTS/FIGURES/income_panel_with_stats.jpg", 
       final_figure, 
       width = 12, 
       height = 8, 
       dpi = 300)
