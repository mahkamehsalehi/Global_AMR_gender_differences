library(vegan)
library(ggpubr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(tidyverse)
library(cowplot)
library(rstatix)
library(SummarizedExperiment)

set.seed(123)
# -----------------------------
# Data Loading and Preprocessing
# -----------------------------

TSE <- readRDS("../DATA/TSE_filtered.rds")
tse_metadata <- as.data.frame(colData(TSE))

metadata_hic <- tse_metadata %>%
  filter(income_group == "HIC") %>%
  filter(!is.na(age_category_new) & !is.na(ARG_load))

metadata_lmic <- tse_metadata %>%
  filter(income_group == "LMIC") %>%
  filter(!is.na(age_category_new) & !is.na(ARG_load))

comparisons <- list(
  c("Infant", "Toddler"),
  c("Toddler", "Children"),
  c("Children", "Teenager"),
  c("Teenager", "Young Adult"),
  c("Young Adult", "Middle-Aged Adult"),
  c("Middle-Aged Adult", "Older Adult"),
  c("Older Adult", "Oldest Adult")
)

age_arg_boxplot_hic <- ggplot(metadata_hic, aes(x = gender, y = ARG_div_shan, fill = gender)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.6,
    alpha = 1
  ) +
  facet_wrap(~age_category_new, scales = "fixed", nrow = 1) +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(x = "Gender", y = "ARG diversity (Shannon index)") +
  stat_compare_means(
    comparisons = list(c("Women", "Men")),
    aes(label = ..p.signif..),
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE,
    size = 6,
    label.y = 3
  ) +
  theme_minimal(16) +
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_blank(),
    strip.text = element_text(size = 12, angle = 0),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "bottom"
  )

age_arg_boxplot_lmic <- ggplot(metadata_lmic, aes(x = gender, y = ARG_div_shan, fill = gender)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.6,
    alpha = 1,
    show.legend = FALSE
  ) +
  facet_wrap(~age_category_new, scales = "fixed", nrow = 1) +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(x = "Gender", y = "ARG diversity (Shannon index)") +
  theme_minimal(16) +
  stat_compare_means(
    comparisons = list(c("Women", "Men")),
    aes(label = ..p.signif..),
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE,
    size = 6,
    label.y = 3
  ) +
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_blank(),
    strip.text = element_text(size = 12, angle = 0),
    panel.spacing = unit(0.5, "lines")
  )

# Combine the plots
combined_plot <- (age_arg_boxplot_hic + labs(title = "HIC")) +
  (age_arg_boxplot_lmic + labs(title = "LMIC")) +
  plot_layout(ncol = 1, heights = c(1, 1), guides = "collect") +
  plot_annotation(
    tag_levels = 'a',
    tag_prefix = "",
    tag_suffix = ")",
    theme = theme(
      plot.title = element_text(),
      plot.subtitle = element_text(),
      legend.position = "bottom",
      plot.tag = element_text(face = "bold", size = 18)
    )
  )

# Create headers for tables
header_hic <- ggdraw() + 
  draw_label("HIC", fontface = "bold", size = 16, hjust = 0.5) +
  theme(plot.margin = margin(b = 1))

header_lmic <- ggdraw() + 
  draw_label("LMIC", fontface = "bold", size = 16, hjust = 0.5) +
  theme(plot.margin = margin(b = 1))

# Set the correct order for age categories
age_order <- c("Infant", "Toddler", "Children", "Teenager", "Young Adult", 
               "Middle-Aged Adult", "Older Adult", "Oldest Adult")

# Calculate statistics for HIC
hic_effect_sizes <- lapply(unique(metadata_hic$age_category_new), function(age) {
  hic_data <- metadata_hic %>% filter(age_category_new == age)
  
  if(nrow(hic_data) >= 2) {
    test <- wilcox.test(ARG_div_shan ~ gender, data = hic_data)
    res <- wilcox_effsize(data = hic_data, ARG_div_shan ~ gender)
    n_samples <- hic_data %>% group_by(gender) %>% summarise(n = n())
    
    data.frame(
      age_category_new = age,
      n_women = n_samples$n[n_samples$gender == "Women"],
      n_men = n_samples$n[n_samples$gender == "Men"],
      effect_size = res$effsize,
      p_value = test$p.value
    )
  }
}) %>% bind_rows()

# Calculate statistics for LMIC
lmic_effect_sizes <- lapply(unique(metadata_lmic$age_category_new), function(age) {
  lmic_data <- metadata_lmic %>% filter(age_category_new == age)
  
  if(nrow(lmic_data) >= 2) {
    test <- wilcox.test(ARG_div_shan ~ gender, data = lmic_data)
    res <- wilcox_effsize(data = lmic_data, ARG_div_shan ~ gender)
    n_samples <- lmic_data %>% group_by(gender) %>% summarise(n = n())
    
    data.frame(
      age_category_new = age,
      n_women = n_samples$n[n_samples$gender == "Women"],
      n_men = n_samples$n[n_samples$gender == "Men"],
      effect_size = res$effsize,
      p_value = test$p.value
    )
  }
}) %>% bind_rows()

# Format HIC statistics
hic_stats <- hic_effect_sizes %>%
  mutate(
    age_category_new = factor(age_category_new, levels = age_order)
  ) %>%
  arrange(age_category_new) %>%
  mutate(
    `Effect Size (r)` = round(effect_size, 3),
    # Apply BH adjustment
    `Adjusted p-value` = case_when(
      p_value < 0.0001 ~ "p<0.0001",
      TRUE ~ formatC(p_value, format = "f", digits = 4)
    ),
    `Lower 95% CI` = round(effect_size - 1.96 * sqrt((1 - effect_size^2) / (n_women + n_men - 2)), 3),
    `Upper 95% CI` = round(effect_size + 1.96 * sqrt((1 - effect_size^2) / (n_women + n_men - 2)), 3)
  ) %>%
  select(
    `Age Category` = age_category_new,
    `N (Women)` = n_women,
    `N (Men)` = n_men,
    `Effect Size (r)`,
    `Lower 95% CI`,
    `Upper 95% CI`,
    `Adjusted p-value`
  )

# Format LMIC statistics
lmic_stats <- lmic_effect_sizes %>%
  mutate(
    age_category_new = factor(age_category_new, levels = age_order)
  ) %>%
  arrange(age_category_new) %>%
  mutate(
    `Effect Size (r)` = round(effect_size, 3),
    # Apply BH adjustment
    `Adjusted p-value` = case_when(
      p_value < 0.0001 ~ "p<0.0001",
      TRUE ~ formatC(p_value, format = "f", digits = 4)
    ),
    `Lower 95% CI` = round(effect_size - 1.96 * sqrt((1 - effect_size^2) / (n_women + n_men - 2)), 3),
    `Upper 95% CI` = round(effect_size + 1.96 * sqrt((1 - effect_size^2) / (n_women + n_men - 2)), 3)
  ) %>%
  select(
    `Age Category` = age_category_new,
    `N (Women)` = n_women,
    `N (Men)` = n_men,
    `Effect Size (r)`,
    `Lower 95% CI`,
    `Upper 95% CI`,
    `Adjusted p-value`
  )

# Create HIC table
hic_table <- ggtexttable(
  hic_stats,
  rows = NULL,
  theme = ttheme("light", 
                 base_size = 16,
                 padding = unit(c(10, 20), "pt"))
) %>%
  tab_add_title(text = "c)", 
                face = "bold", 
                size = 18, 
                just = "left",
                padding = unit(c(0, 0, 0, 4), "pt"))

# Create LMIC table
lmic_table <- ggtexttable(
  lmic_stats,
  rows = NULL,
  theme = ttheme("light", 
                 base_size = 16,
                 padding = unit(c(10, 20), "pt"))
) %>%
  tab_add_title(text = "d)", 
                face = "bold", 
                size = 18, 
                just = "left",
                padding = unit(c(0, 0, 0, 4), "pt"))

# Create separator
separator <- ggdraw() + 
  draw_line(x = c(0, 1), y = c(0.5, 0.5), color = "grey", size = 2)

# Combine
final_plot <- plot_grid(
  combined_plot,
  separator,
  plot_grid(header_hic, hic_table, ncol = 1, rel_heights = c(0.1, 1)),
  separator,
  plot_grid(header_lmic, lmic_table, ncol = 1, rel_heights = c(0.1, 1)),
  ncol = 1,
  rel_heights = c(2, 0.1, 1, 0.1, 1)
)

# Save the final plot
ggsave("../RESULTS/FIGURES/Supplementary Figure 7.jpg", 
       final_plot, 
       width = 12, 
       height = 25,
       dpi = 500)