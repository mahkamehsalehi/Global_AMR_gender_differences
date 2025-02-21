# --------------------------
# Load Necessary Libraries
# --------------------------
library(vegan)
library(ggpubr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(cowplot)
library(mia)
library(Cairo)
library(scales)
library(rstatix)

# --------------------------
# Load and Prepare Data
# --------------------------
TSE <- readRDS("DATA/TSE_filtered.rds")
tse_metadata <- as.data.frame(colData(TSE))

filtered_metadata <- tse_metadata %>%
  filter(
    category != "Infant Study", 
    geo_loc_name_country_calc != "Zimbabwe"  # To explain in the ms!
  ) %>%
  drop_na(log_ARG_load, income_group)

# --------------------------
# Define Consistent Color Palette and Theme
# --------------------------
gender_colors <- c("Women" = "#F8766D", "Men" = "#619CFF")

custom_theme <- theme_minimal(24) +
  theme(
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", size = 1),
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(),
    axis.text = element_text(),
    legend.position = "none"
  )

# --------------------------
# Create Violin Plots
# --------------------------

# Violin plot for ARG Load
income_arg_violinplot <- ggviolin(
  filtered_metadata,
  x = "gender",
  y = "ARG_load",
  fill = "gender",
  palette = gender_colors,
  add = "boxplot",
  add.params = list(fill = "white", width = 0.1)
) +
  stat_compare_means(
    comparisons = list(c("Women", "Men")),
    label = "p.format",
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE,
    size = 6
  ) +
  labs(
    x = "",
    y = expression("ARG load (log "*RPKM*")")
  ) +
  facet_wrap(~income_group) +
  scale_y_log10(
    labels = trans_format("log10", math_format(10^.x))
  ) +
  theme_minimal() +
  custom_theme +
  theme(
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", fill = "white", size = 1),
    axis.text = element_text(size = 14),
    axis.title.y = element_text(size = 18),
    legend.position = "none"
  )

# Violin plot for ARG Diversity (Shannon Index)
income_shannon_violinplot <- ggviolin(
  filtered_metadata,
  x = "gender",
  y = "shannon_diversity",
  fill = "gender",
  palette = gender_colors,
  add = "boxplot",
  add.params = list(fill = "white", width = 0.1)
) +
  stat_compare_means(
    comparisons = list(c("Women", "Men")),
    label = "p.format",
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE,
    size = 6
  ) +
  labs(
    x = "",
    y = "ARG diversity (Shannon index)"
  ) +
  facet_wrap(~income_group) +
  coord_cartesian(ylim = c(0, 4)) +
  theme_minimal() +
  custom_theme +
  theme(
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", fill = "white", size = 1),
    axis.text = element_text(size = 14),
    axis.title.y = element_text(size = 18),
    legend.position = "none"
  )

# --------------------------
# Compute Statistics for ARG Load
# --------------------------
stats_ARG_load <- filtered_metadata %>%
  group_by(income_group) %>%
  group_modify(~ {
    n_Women <- sum(.x$gender == "Women")
    n_Men   <- sum(.x$gender == "Men")
    if(length(unique(.x$gender)) < 2) {
      tibble(
        n_Women = n_Women,
        n_Men = n_Men,
        effect_size = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        p_value = NA_real_
      )
    } else {
      eff <- wilcox_effsize(.x, formula = ARG_load ~ gender, ci = TRUE, conf.level = 0.95)
      test <- wilcox_test(.x, formula = ARG_load ~ gender)
      tibble(
        n_Women = n_Women,
        n_Men   = n_Men,
        effect_size = eff$effsize[1],
        conf.low = eff$conf.low[1],
        conf.high = eff$conf.high[1],
        p_value = test$p[1]
      )
    }
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# Format table for ARG Load
arg_load_table_df <- stats_ARG_load %>%
  select(income_group, n_Women, n_Men, effect_size, conf.low, conf.high, p_adj) %>%
  rename(
    `Income Group` = income_group,
    `N (Women)` = n_Women,
    `N (Men)` = n_Men,
    `Effect Size (r)` = effect_size,
    `Lower 95% CI` = conf.low,
    `Upper 95% CI` = conf.high,
    `Adjusted p-value` = p_adj
  ) %>%
  mutate(
    across(c(`Effect Size (r)`, `Lower 95% CI`, `Upper 95% CI`), ~round(.x, 3)),
    `Adjusted p-value` = case_when(
      `Adjusted p-value` < 0.0001 ~ "p<0.0001",
      TRUE ~ formatC(`Adjusted p-value`, format = "f", digits = 4)
    )
  )

arg_load_table <- ggtexttable(arg_load_table_df, rows = NULL, 
                              theme = ttheme("light", base_size = 16))

# --------------------------
# Compute Statistics for ARG Diversity
# --------------------------
stats_ARG_diversity <- filtered_metadata %>%
  group_by(income_group) %>%
  group_modify(~ {
    n_Women <- sum(.x$gender == "Women")
    n_Men   <- sum(.x$gender == "Men")
    if(length(unique(.x$gender)) < 2) {
      tibble(
        n_Women = n_Women,
        n_Men = n_Men,
        effect_size = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        p_value = NA_real_
      )
    } else {
      eff <- wilcox_effsize(.x, formula = shannon_diversity ~ gender, ci = TRUE, conf.level = 0.95)
      test <- wilcox_test(.x, formula = shannon_diversity ~ gender)
      tibble(
        n_Women = n_Women,
        n_Men   = n_Men,
        effect_size = eff$effsize[1],
        conf.low = eff$conf.low[1],
        conf.high = eff$conf.high[1],
        p_value = test$p[1]
      )
    }
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# Format table for ARG Diversity
arg_diversity_table_df <- stats_ARG_diversity %>%
  select(income_group, n_Women, n_Men, effect_size, conf.low, conf.high, p_adj) %>%
  rename(
    `Income Group` = income_group,
    `N (Women)` = n_Women,
    `N (Men)` = n_Men,
    `Effect Size (r)` = effect_size,
    `Lower 95% CI` = conf.low,
    `Upper 95% CI` = conf.high,
    `Adjusted p-value` = p_adj
  ) %>%
  mutate(
    across(c(`Effect Size (r)`, `Lower 95% CI`, `Upper 95% CI`), ~round(.x, 3)),
    `Adjusted p-value` = case_when(
      `Adjusted p-value` < 0.0001 ~ "p<0.0001",
      TRUE ~ formatC(`Adjusted p-value`, format = "f", digits = 4)
    )
  )

arg_diversity_table <- ggtexttable(arg_diversity_table_df, rows = NULL, 
                                   theme = ttheme("light", base_size = 16))

# --------------------------
# Merge the Two Tables into One Composite Table with Tags c and d
# --------------------------
# Create header labels for each section
header_ARG_load <- ggdraw() + 
  draw_label("ARG load", fontface = "bold", size = 20, hjust = 0.5)
header_ARG_diversity <- ggdraw() + 
  draw_label("ARG diversity", fontface = "bold", size = 20, hjust = 0.5)

label_c <- ggdraw() + draw_label("c", fontface = "bold", size = 24, x = 0, hjust = 0)
label_d <- ggdraw() + draw_label("d", fontface = "bold", size = 24, x = 0, hjust = 0)

# Stack the label, header, and table for each metric vertically
composite_ARG_load <- plot_grid(
  label_c,
  header_ARG_load, 
  arg_load_table,
  ncol = 1,
  rel_heights = c(0.1, 0.1, 1)
)
composite_ARG_diversity <- plot_grid(
  label_d,
  header_ARG_diversity, 
  arg_diversity_table,
  ncol = 1,
  rel_heights = c(0.1, 0.1, 1)
)

# Create a grey horizontal separator by specifying x and y coordinates
separator <- ggdraw() + 
  draw_line(x = c(0, 1), y = c(0.5, 0.5), color = "grey", size = 2)

# Use plot_spacer() to add a gap between the composite tables
combined_table <- plot_grid(
  composite_ARG_load,
  separator,
  composite_ARG_diversity,
  ncol = 1,
  rel_heights = c(1, 0.1, 1)
)

# Modify the stats data frames to include a Metric column
stats_ARG_load <- stats_ARG_load %>%
  mutate(Metric = "ARG load")

stats_ARG_diversity <- stats_ARG_diversity %>%
  mutate(Metric = "ARG diversity")

# Combine the statistics
combined_stats <- bind_rows(stats_ARG_load, stats_ARG_diversity) %>%
  arrange(income_group, Metric)

# Format the combined table
combined_table_df <- combined_stats %>%
  select(Metric, income_group, n_Women, n_Men, effect_size, conf.low, conf.high, p_adj) %>%
  rename(
    `Metric` = Metric,
    `Income Group` = income_group,
    `N (Women)` = n_Women,
    `N (Men)` = n_Men,
    `Effect Size (r)` = effect_size,
    `Lower 95% CI` = conf.low,
    `Upper 95% CI` = conf.high,
    `Adjusted p-value` = p_adj
  ) %>%
  mutate(
    across(c(`Effect Size (r)`, `Lower 95% CI`, `Upper 95% CI`), ~round(.x, 3)),
    `Adjusted p-value` = case_when(
      `Adjusted p-value` < 0.0001 ~ "p<0.0001",
      TRUE ~ formatC(`Adjusted p-value`, format = "f", digits = 4)
    )
  )

# Create the table with the tag using tab_add_title
stats_table <- ggtexttable(combined_table_df, 
                           rows = NULL,
                           theme = ttheme("light", 
                                          base_size = 16, 
                                          padding = unit(c(10, 20), "pt")
                           )) %>%
  tab_add_title(text = "c", 
                size = 24, 
                face = "bold", 
                just = "left",
                padding = unit(c(0, 0, 0, 4), "pt"))

# Top row: combine the two violin plots side by side with tags a and b
top_row <- plot_grid(
  income_arg_violinplot,
  income_shannon_violinplot,
  labels = c("a", "b"),
  label_size = 24,
  ncol = 2,
  align = 'hv'
)

# Create a thin separator
separator <- ggdraw() + 
  draw_line(x = c(0, 1), y = c(0.5, 0.5), color = "grey", size = 2)

# Combine everything into the final figure
final_figure <- plot_grid(
  top_row,
  separator,
  stats_table,
  ncol = 1,
  rel_heights = c(1, 0.02, 0.8)
)

# Save the final figure
CairoJPEG("RESULTS/FIGURES/Fig3.jpg", width = 1200, height = 900, quality = 100)
print(final_figure)
dev.off()



###############################################################################

# --------------------------
# Arrange Plots and the Combined Table in the Final Figure
# --------------------------

# Top row: combine the two violin plots side by side with tags a and b
top_row <- plot_grid(
  income_arg_violinplot,
  income_shannon_violinplot,
  labels = c("a", "b"),
  label_size = 24,
  ncol = 2,
  align = 'hv'
)

# Combine top row, separator, and the combined table into the final figure
final_figure <- plot_grid(
  top_row,
  separator,
  combined_table,
  ncol = 1,
  rel_heights = c(1, 0.05, 1)
)

# --------------------------
# Save the Final Figure
# --------------------------
CairoJPEG("RESULTS/FIGURES/Fig3.jpg", width = 1200, height = 900, quality = 100)
print(final_figure)
dev.off()

