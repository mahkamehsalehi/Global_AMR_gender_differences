# --------------------------
# Load Libraries
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
library(ggtext)

set.seed(123)
# --------------------------
# Load and Prepare Data
# --------------------------
TSE <- readRDS("../DATA/TSE_filtered.rds")
tse_metadata <- as.data.frame(colData(TSE))

filtered_metadata <- tse_metadata %>%
  filter(
    category != "Infant Study", 
    geo_loc_name_country_calc != "Zimbabwe"
  ) %>%
  mutate(log_ARG_load = log(ARG_load)) %>%
  drop_na(log_ARG_load, income_group)

# --------------------------
# Custom Theme & Colors
# --------------------------
gender_colors <- c("Women" = "#F8766D", "Men" = "#619CFF")

custom_theme <- theme_minimal(base_size = 16) +
  theme(
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", size = 1),
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none"
  )

# --------------------------
# Violin Plots (Gender differences within each income group)
# --------------------------

# Panel a): Violin Plot for ARG Load
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
    size = 5
  ) +
  labs(
    x = NULL,
    y = expression("ARG load (log "*RPKM*")")
  ) +
  facet_wrap(~ income_group) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  custom_theme

# Panel b): Violin Plot for ARG Diversity
income_shannon_violinplot <- ggviolin(
  filtered_metadata,
  x = "gender",
  y = "ARG_div_shan",
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
    size = 5
  ) +
  labs(
    x = NULL,
    y = "ARG diversity (Shannon index)"
  ) +
  facet_wrap(~ income_group) +
  coord_cartesian(ylim = c(0, 4)) +
  custom_theme

# Combine Panels
top_row <- plot_grid(
  income_arg_violinplot,
  income_shannon_violinplot,
  labels = c("a)", "b)"),
  label_size = 24,
  ncol = 2,
  align = 'hv'
)

# --------------------------
# Statistics: Gender differences within each income group
# --------------------------

stats_ARG_load <- filtered_metadata %>%
  group_by(income_group) %>%
  group_modify(~ {
    n_Women <- sum(.x$gender == "Women")
    n_Men   <- sum(.x$gender == "Men")
    if (length(unique(.x$gender)) < 2) {
      tibble(
        n_Women = n_Women,
        n_Men = n_Men,
        effect_size = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        p_value = NA_real_
      )
    } else {
      eff  <- wilcox_effsize(.x, formula = ARG_load ~ gender, ci = TRUE)
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
  mutate(Metric = "ARG load") %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

stats_ARG_diversity <- filtered_metadata %>%
  group_by(income_group) %>%
  group_modify(~ {
    n_Women <- sum(.x$gender == "Women")
    n_Men   <- sum(.x$gender == "Men")
    if (length(unique(.x$gender)) < 2) {
      tibble(
        n_Women = n_Women,
        n_Men = n_Men,
        effect_size = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        p_value = NA_real_
      )
    } else {
      eff  <- wilcox_effsize(.x, formula = ARG_div_shan ~ gender, ci = TRUE)
      test <- wilcox_test(.x, formula = ARG_div_shan ~ gender)
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
  mutate(Metric = "ARG diversity") %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# Combine the two metrics
stats_combined_gender <- bind_rows(stats_ARG_load, stats_ARG_diversity)

# Format for display
combined_table_df_c <- stats_combined_gender %>%
  select(Metric, income_group, n_Women, n_Men, effect_size, conf.low, conf.high, p_adj) %>%
  rename(
    "Income Group"     = income_group,
    "N (Women)"        = n_Women,
    "N (Men)"          = n_Men,
    "Effect Size (r)"  = effect_size,
    "Lower 95% CI"     = conf.low,
    "Upper 95% CI"     = conf.high,
    "Adjusted p-value" = p_adj
  ) %>%
  mutate(
    across(c(`Effect Size (r)`, `Lower 95% CI`, `Upper 95% CI`), ~round(.x, 3)),
    `Adjusted p-value` = case_when(
      `Adjusted p-value` < 0.0001 ~ "p<0.0001",
      TRUE ~ formatC(`Adjusted p-value`, format = "f", digits = 4)
    )
  )

# Create Panel c) table
stats_table_c <- ggtexttable(
  combined_table_df_c,
  rows = NULL,
  theme = ttheme("light", base_size = 16)
) %>%
  tab_add_title(
    text = "c)",
    size = 24,
    face = "bold",
    just = "left",
    padding = unit(c(0, 0, 0, 4), "pt")
  )

# --------------------------
# Overall Comparison: HIC vs. LMIC (Ignoring Gender) (Panel d)
# --------------------------

# Sample sizes
n_HIC <- filtered_metadata %>% filter(income_group == "HIC") %>% nrow()
n_LMIC <- filtered_metadata %>% filter(income_group == "LMIC") %>% nrow()

# ARG Load: overall comparison
overall_stats_ARG_load <- wilcox_effsize(
  filtered_metadata, 
  formula = ARG_load ~ income_group, 
  ci = TRUE, conf.level = 0.95
)
overall_test_ARG_load <- wilcox_test(
  filtered_metadata, 
  formula = ARG_load ~ income_group
)

overall_stats_ARG_load_df <- tibble(
  Metric = "ARG load",
  `N (HIC)` = n_HIC,
  `N (LMIC)` = n_LMIC,
  Effect_Size = overall_stats_ARG_load$effsize,
  Lower_95_CI = overall_stats_ARG_load$conf.low,
  Upper_95_CI = overall_stats_ARG_load$conf.high,
  p_value = overall_test_ARG_load$p
) %>%
  mutate(Adjusted_p_value = ifelse(
    p_value < 0.0001, "p<0.0001", formatC(p_value, format = "f", digits = 4)
  ))

# ARG Diversity: overall comparison
overall_stats_ARG_diversity <- wilcox_effsize(
  filtered_metadata, 
  formula = ARG_div_shan ~ income_group, 
  ci = TRUE, conf.level = 0.95
)
overall_test_ARG_diversity <- wilcox_test(
  filtered_metadata, 
  formula = ARG_div_shan ~ income_group
)

overall_stats_ARG_diversity_df <- tibble(
  Metric = "ARG diversity",
  `N (HIC)` = n_HIC,
  `N (LMIC)` = n_LMIC,
  Effect_Size = overall_stats_ARG_diversity$effsize,
  Lower_95_CI = overall_stats_ARG_diversity$conf.low,
  Upper_95_CI = overall_stats_ARG_diversity$conf.high,
  p_value = overall_test_ARG_diversity$p
) %>%
  mutate(Adjusted_p_value = ifelse(
    p_value < 0.0001, "p<0.0001", formatC(p_value, format = "f", digits = 4)
  ))


# Combine overall stats into one table
overall_stats_table <- bind_rows(
  overall_stats_ARG_load_df,
  overall_stats_ARG_diversity_df
) %>%
  select(Metric, `N (HIC)`, `N (LMIC)`, Effect_Size, Lower_95_CI, Upper_95_CI, Adjusted_p_value) %>%
  rename(
    "Effect Size (r)"   = Effect_Size,
    "Lower 95% CI"      = Lower_95_CI,
    "Upper 95% CI"      = Upper_95_CI,
    "Adjusted p-value"  = Adjusted_p_value
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

# Create Panel d) table
stats_table_d <- ggtexttable(
  overall_stats_table,
  rows = NULL,
  theme = ttheme("light", base_size = 16)
) %>%
  tab_add_title(
    text = "d)",
    size = 24,
    face = "bold",
    just = "left",
    padding = unit(c(0, 0, 0, 4), "pt")
  )

# --------------------------
# Stack Panels
# --------------------------
tables_stacked <- plot_grid(
  stats_table_c,
  stats_table_d,
  ncol = 1,
  rel_heights = c(1, 1)
)

# --------------------------
# Combine Panels
# --------------------------
final_figure <- plot_grid(
  top_row,
  tables_stacked,
  ncol = 1,
  rel_heights = c(1, 1.2)
)

# --------------------------
# Save the Final Figure
# --------------------------
CairoJPEG("../RESULTS/FIGURES/Figure 3.jpg", width = 4000, height = 3000, units = "px",dpi = 300, quality = 100)
print(final_figure)
dev.off()