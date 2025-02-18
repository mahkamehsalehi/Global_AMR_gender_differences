library(vegan)
library(ggpubr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(rstatix)
library(cowplot)

# ----------------------------------------------------------------------------------
# Load and preprocess data
# ----------------------------------------------------------------------------------

TSE <- readRDS("DATA/TSE_filtered.rds")
tse_metadata <- as.data.frame(colData(TSE))

# - Remove samples belonging to the "Infant Study"
# - Remove samples from Zimbabwe
filtered_metadata <- tse_metadata %>%
  filter(
    category != "Infant Study", 
    geo_loc_name_country_calc != "Zimbabwe"
  )

# ----------------------------------------------------------------------------------
# Define a common theme for all plots to ensure visual consistency
# ----------------------------------------------------------------------------------

common_theme <- theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none",
    axis.line = element_line(color = "black"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 12, face = "bold")
  )

# ----------------------------------------------------------------------------------
# Data preparation: Define usage groups and filter further for analysis
# ----------------------------------------------------------------------------------

# Categorize samples based on the 'Usage' variable:
# - 'High' if Usage > 10
# - 'Low' if Usage <= 10
usage_filtered_metadata <- tse_metadata %>%
  mutate(
    Usage_group = case_when(
      Usage > 10  ~ "High",
      Usage <= 10 ~ "Low",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(Usage_group = factor(Usage_group, levels = c("High", "Low"))) %>%
  drop_na(gender, log_ARG_load, Usage_group, income_group)

# Subset the data based on income group
metadata_hic <- usage_filtered_metadata %>%
  filter(income_group == "HIC")

metadata_lmic <- usage_filtered_metadata %>%
  filter(income_group == "LMIC")

# ----------------------------------------------------------------------------------
# Create Boxplots for ARG load and Shannon diversity (by gender and usage group)
# ----------------------------------------------------------------------------------

# -----------------------
# ARG load by gender across usage groups for High-Income Countries (HIC)
# -----------------------
usage_arg_boxplot_hic <- ggplot(metadata_hic, aes(x = gender,
                                                  y = ARG_load,
                                                  fill = gender)) +
  geom_boxplot(position = position_dodge(width = 0.8), 
               outlier.shape = NA, 
               width = 0.6, 
               alpha = 1, 
               show.legend = FALSE) +
  labs(x = "Gender", 
       y = "ARG load (RPKM)",
       title = "HICs") +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  # Add Wilcoxon test comparison between genders (using the log-transformed ARG load)
  stat_compare_means(aes(x = gender, y = log_ARG_load),
                     comparisons = list(c("Women", "Men")), 
                     label = "p.format", 
                     method = "wilcox.test", 
                     p.adjust.method = "BH",
                     hide.ns = FALSE,
                     size = 5,
                     label.y = 4) +
  facet_wrap(~Usage_group) +
  scale_y_continuous(trans = "log10",
                     breaks = 10^(2:5),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_minimal(16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", size = 0.7),
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14)
  )

# -----------------------
# Shannon diversity by gender across usage groups for HIC
# -----------------------
usage_shannon_boxplot_hic <- ggplot(metadata_hic, 
                                    aes(x = gender, 
                                        y = shannon_diversity, 
                                        fill = gender)) +
  geom_boxplot(position = position_dodge(width = 0.8), 
               outlier.shape = NA, 
               width = 0.6, 
               alpha = 1, 
               show.legend = FALSE) +
  labs(x = "Gender", 
       y = "ARG diversity",
       title = "HICs") +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  # Wilcoxon test for Shannon diversity between genders
  stat_compare_means(aes(x = gender, y = shannon_diversity),
                     comparisons = list(c("Women", "Men")), 
                     label = "p.format", 
                     method = "wilcox.test", 
                     p.adjust.method = "BH",
                     hide.ns = FALSE,
                     size = 5,
                     label.y = 3.2) +
  facet_wrap(~Usage_group) +
  theme_minimal(16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", size = 0.7),
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14)
  )

# -----------------------
# ARG load by usage for women only (HIC)
# -----------------------
usage_arg_boxplot_female_hic <- ggplot(metadata_hic %>% filter(gender == "Women"), 
                                       aes(x = Usage_group, y = ARG_load, fill = "#F8766D")) +
  geom_boxplot(position = position_dodge(width = 0.8), 
               outlier.shape = NA, 
               width = 0.6, 
               alpha = 1, 
               show.legend = FALSE) +
  labs(x = "Antibiotic Use (DDD)", 
       y = "ARG load (RPKM)",
       title = "Women in HICs") +
  scale_y_continuous(trans = "log10",
                     breaks = 10^(2:5),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  # Compare ARG load between High and Low usage groups using Wilcoxon test
  stat_compare_means(comparisons = list(c("Low", "High")), 
                     label = "p.format", 
                     method = "wilcox.test", 
                     hide.ns = FALSE, 
                     size = 5, 
                     label.y = 4) +
  theme_minimal(16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", size = 0.7),
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14)
  )

# -----------------------
# Shannon diversity by usage for women only (HIC)
# -----------------------
usage_shannon_boxplot_female_hic <- ggplot(metadata_hic %>% filter(gender == "Women"),
                                           aes(x = Usage_group, y = shannon_diversity, fill = "#F8766D")) +
  geom_boxplot(position = position_dodge(width = 0.8), 
               outlier.shape = NA, 
               width = 0.6, 
               alpha = 1, 
               show.legend = FALSE) +
  labs(x = "Antibiotic Use(DDD)", 
       y = "ARG diversity",
       title = "Women in HICs") +
  stat_compare_means(comparisons = list(c("Low", "High")), 
                     label = "p.format", 
                     method = "wilcox.test", 
                     hide.ns = FALSE,  
                     size = 5, 
                     label.y = 3) +
  theme_minimal(16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", size = 0.7),
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14)
  )

# -----------------------
# ARG load by gender across usage groups for Low/Middle-Income Countries (LMIC)
# -----------------------
usage_arg_boxplot_lmic <- ggplot(metadata_lmic, aes(x = gender,
                                                    y = ARG_load,
                                                    fill = gender)) +
  geom_boxplot(position = position_dodge(width = 0.8), 
               outlier.shape = NA, 
               width = 0.6, 
               alpha = 1, 
               show.legend = FALSE) +
  labs(x = "Gender", 
       y = "ARG load (RPKM)",
       title = "LMICs") +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  stat_compare_means(aes(x = gender, y = log_ARG_load),
                     comparisons = list(c("Women", "Men")), 
                     label = "p.format", 
                     method = "wilcox.test", 
                     p.adjust.method = "BH",
                     hide.ns = FALSE,
                     size = 5,
                     label.y = 4.5) +
  facet_wrap(~Usage_group) +
  scale_y_continuous(trans = "log10",
                     breaks = 10^(2:5),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_minimal(16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", size = 0.7),
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14)
  )

# -----------------------
# Shannon diversity by gender across usage groups for LMIC
# -----------------------
usage_shannon_boxplot_lmic <- ggplot(metadata_lmic, 
                                     aes(x = gender, 
                                         y = shannon_diversity, 
                                         fill = gender)) +
  geom_boxplot(position = position_dodge(width = 0.8), 
               outlier.shape = NA, 
               width = 0.6, 
               alpha = 1, 
               show.legend = FALSE) +
  labs(x = "Gender", 
       y = "ARG diversity",
       title = "LMICs") +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  stat_compare_means(aes(x = gender, y = shannon_diversity),
                     comparisons = list(c("Women", "Men")), 
                     label = "p.format", 
                     method = "wilcox.test", 
                     p.adjust.method = "BH",
                     hide.ns = FALSE,
                     size = 5,
                     label.y = 3) +
  facet_wrap(~Usage_group) +
  theme_minimal(16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", size = 0.7),
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14)
  )

# -----------------------
# ARG load by usage for women only (LMIC)
# -----------------------
usage_arg_boxplot_female_lmic <- ggplot(metadata_lmic %>% filter(gender == "Women"), 
                                        aes(x = Usage_group, y = ARG_load, fill = "#F8766D")) +
  geom_boxplot(position = position_dodge(width = 0.8), 
               outlier.shape = NA, 
               width = 0.6, 
               alpha = 1, 
               show.legend = FALSE) +
  labs(x = "Antibiotic Use (DDD)", 
       y = "ARG load (RPKM)",
       title = "Women in LMICs") +
  scale_y_continuous(trans = "log10",
                     breaks = 10^(2:5),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  stat_compare_means(comparisons = list(c("Low", "High")), 
                     label = "p.format", 
                     method = "wilcox.test", 
                     hide.ns = FALSE, 
                     size = 5, 
                     label.y = 4.5) +
  theme_minimal(16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", size = 0.7),
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14)
  )

# -----------------------
# Shannon diversity by usage for women only (LMIC)
# -----------------------
usage_shannon_boxplot_female_lmic <- ggplot(metadata_lmic %>% filter(gender == "Women"),
                                            aes(x = Usage_group, y = shannon_diversity, fill = "#F8766D")) +
  geom_boxplot(position = position_dodge(width = 0.8), 
               outlier.shape = NA, 
               width = 0.6, 
               alpha = 1, 
               show.legend = FALSE) +
  labs(x = "Antibiotic Use (DDD)", 
       y = "ARG diversity",
       title = "Women in LMICs") +
  stat_compare_means(comparisons = list(c("Low", "High")), 
                     label = "p.format", 
                     method = "wilcox.test", 
                     hide.ns = FALSE, 
                     size = 5, 
                     label.y = 3) +
  theme_minimal(16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", size = 0.7),
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14)
  )

# ----------------------------------------------------------------------------------
# Define a function to calculate statistics for comparing usage groups
# ----------------------------------------------------------------------------------

calc_stats <- function(data, variable, group) {
  # Calculate sample sizes for each usage group and reshape the output
  n_sizes <- data %>%
    group_by(Usage_group) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = Usage_group, values_from = n, names_prefix = "n_")
  
  data$Usage_group <- factor(data$Usage_group, levels = c("High", "Low"))
  
  # Calculate effect size and perform the Wilcoxon test
  if(variable == "ARG_load") {
    eff <- wilcox_effsize(data = data, 
                          ARG_load ~ Usage_group, 
                          ci = TRUE, 
                          conf.level = 0.95)
    test <- wilcox.test(ARG_load ~ Usage_group, data = data)
  } else {
    eff <- wilcox_effsize(data = data, 
                          shannon_diversity ~ Usage_group, 
                          ci = TRUE, 
                          conf.level = 0.95)
    test <- wilcox.test(shannon_diversity ~ Usage_group, data = data)
  }
  
  # Combine results into a summary data frame
  results <- data.frame(
    "Group" = group,
    "N (High)" = n_sizes$n_High,
    "N (Low)" = n_sizes$n_Low,
    "Effect Size (r)" = round(eff$effsize, 3),
    "Lower 95% CI" = round(eff$conf.low, 3),
    "Upper 95% CI" = round(eff$conf.high, 3),
    "Adjusted p-value" = formatC(p.adjust(test$p.value, method = "BH"), 
                                 format = "f", digits = 4),
    check.names = FALSE
  )
  
  return(results)
}

# ----------------------------------------------------------------------------------
# Calculate statistics for ARG load and Shannon diversity across different groups
# ----------------------------------------------------------------------------------

# Statistics for ARG load
arg_stats_hic <- calc_stats(metadata_hic, "ARG_load", "HICs")
arg_stats_lmic <- calc_stats(metadata_lmic, "ARG_load", "LMICs")
arg_stats_hic_women <- calc_stats(metadata_hic %>% filter(gender == "Women"), "ARG_load", "Women in HICs")
arg_stats_lmic_women <- calc_stats(metadata_lmic %>% filter(gender == "Women"), "ARG_load", "Women in LMICs")

# Statistics for Shannon diversity
shannon_stats_hic <- calc_stats(metadata_hic, "shannon_diversity", "HICs")
shannon_stats_lmic <- calc_stats(metadata_lmic, "shannon_diversity", "LMICs")
shannon_stats_hic_women <- calc_stats(metadata_hic %>% filter(gender == "Women"), "shannon_diversity", "Women in HICs")
shannon_stats_lmic_women <- calc_stats(metadata_lmic %>% filter(gender == "Women"), "shannon_diversity", "Women in LMICs")

# Combine statistics from different comparisons for display in tables
combined_arg_stats <- rbind(arg_stats_hic, arg_stats_lmic, arg_stats_hic_women, arg_stats_lmic_women)
combined_shannon_stats <- rbind(shannon_stats_hic, shannon_stats_lmic, shannon_stats_hic_women, shannon_stats_lmic_women)

# ----------------------------------------------------------------------------------
# Combine plots and statistics tables into final figures using patchwork
# ----------------------------------------------------------------------------------

# ----- ARG load figures -----

# Create a grid of ARG load plots (for HIC and LMIC, overall and women only) with automatic tagging (labels a-d)
combined_figure_usage_arg_tagged <- usage_arg_boxplot_hic + 
  usage_arg_boxplot_lmic +
  usage_arg_boxplot_female_hic + 
  usage_arg_boxplot_female_lmic +
  plot_layout(ncol = 2, nrow = 2) +
  plot_annotation(
    tag_levels = 'a',
    theme = theme(
      plot.tag = element_text(size = 18, face = "bold"),
      plot.tag.position = c(0, 1)
    )
  )

# Create a horizontal separator line using cowplot
separator <- ggdraw() + 
  draw_line(x = c(0, 1), y = c(0.5, 0.5), color = "grey", size = 2)

# Create a table for ARG load statistics with a manual title "e"
arg_stats_table <- ggtexttable(combined_arg_stats, 
                               rows = NULL,
                               theme = ttheme("light", 
                                              base_size = 16, 
                                              padding = unit(c(10, 20), "pt")
                               )) %>%
  tab_add_title(text = "e", 
                size = 18, 
                face = "bold", 
                just = "left",
                padding = unit(c(0, 0, 0, 4), "pt"))

# Combine the tagged plot grid, separator, and stats table into the final ARG load figure layout
final_figure_arg <- combined_figure_usage_arg_tagged /
  separator /
  arg_stats_table +
  plot_layout(heights = c(2, 0.01, 1))

# ----- Shannon diversity figures -----

# Create a grid of Shannon diversity plots with automatic tagging (labels a-d)
combined_figure_usage_shannon_tagged <- usage_shannon_boxplot_hic + 
  usage_shannon_boxplot_lmic +
  usage_shannon_boxplot_female_hic + 
  usage_shannon_boxplot_female_lmic +
  plot_layout(ncol = 2, nrow = 2) +
  plot_annotation(
    tag_levels = 'a',
    theme = theme(
      plot.tag = element_text(size = 18, face = "bold"),
      plot.tag.position = c(0, 1)
    )
  )

# Create a table for Shannon diversity statistics with a manual title "e"
shannon_stats_table <- ggtexttable(combined_shannon_stats, 
                                   rows = NULL,
                                   theme = ttheme("light", 
                                                  base_size = 16, 
                                                  padding = unit(c(10, 20), "pt")
                                   )) %>%
  tab_add_title(text = "e", 
                size = 18, 
                face = "bold", 
                just = "left",
                padding = unit(c(0, 0, 0, 4), "pt"))

# Combine the tagged plot grid, separator, and stats table into the final Shannon diversity figure layout
final_figure_shannon <- combined_figure_usage_shannon_tagged /
  separator /
  shannon_stats_table +
  plot_layout(heights = c(2, 0.01, 1))

# ----------------------------------------------------------------------------------
# Create final figures using cowplot for alternative layout and saving
# ----------------------------------------------------------------------------------

# ----- ARG load figure -----

# Arrange the ARG load plots in a grid with manual labels (a-d)
grid_arg <- plot_grid(
  usage_arg_boxplot_hic, usage_arg_boxplot_lmic,
  usage_arg_boxplot_female_hic, usage_arg_boxplot_female_lmic,
  labels = c("a", "b", "c", "d"),
  label_fontface = "bold",
  label_size = 18,
  ncol = 2
)

# Create a horizontal separator line as a cowplot object
separator_cow <- ggdraw() + 
  draw_line(x = c(0, 1), y = c(0.5, 0.5), color = "grey", size = 2)

# Convert the ARG load statistics table (with manual title "e") into a grob for cowplot
arg_table_cow <- as_grob(arg_stats_table)

# Combine the grid, separator, and table into the final ARG load figure layout using cowplot
final_figure_arg_cow <- plot_grid(
  grid_arg,
  separator_cow,
  arg_table_cow,
  ncol = 1,
  rel_heights = c(2, 0.01, 1)
)

# Save the final ARG load figure
ggsave("RESULTS/FIGURES/usage_ARG_load.png", 
       final_figure_arg_cow, 
       width = 12, 
       height = 14,
       dpi = 300)

# ----- Shannon diversity figure -----

# Arrange the Shannon diversity plots in a grid with manual labels (a-d)
grid_shannon <- plot_grid(
  usage_shannon_boxplot_hic, usage_shannon_boxplot_lmic,
  usage_shannon_boxplot_female_hic, usage_shannon_boxplot_female_lmic,
  labels = c("a", "b", "c", "d"),
  label_fontface = "bold",
  label_size = 18,
  ncol = 2
)

# Create a horizontal separator line for the Shannon figure
separator_shannon <- ggdraw() + 
  draw_line(x = c(0, 1), y = c(0.5, 0.5), color = "grey", size = 2)

# Convert the Shannon statistics table (with manual title "e") into a grob
shannon_table_cow <- as_grob(shannon_stats_table)

# Combine the grid, separator, and table into the final Shannon diversity figure
final_figure_shannon_cow <- plot_grid(
  grid_shannon,
  separator_shannon,
  shannon_table_cow,
  ncol = 1,
  rel_heights = c(2, 0.01, 1)
)

# Save the final Shannon diversity figure
ggsave("RESULTS/FIGURES/usage__ARG_diversity.png", 
       final_figure_shannon_cow, 
       width = 12, 
       height = 14,
       dpi = 300)
