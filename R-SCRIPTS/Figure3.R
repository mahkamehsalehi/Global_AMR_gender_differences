# Load necessary libraries
library(vegan)
library(ggpubr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(cowplot)
library(mia)
library(Cairo)
library(scales)

# Load and prepare data
TSE <- readRDS("DATA/TSE_filtered.rds")
tse_metadata <- as.data.frame(colData(TSE))

filtered_metadata <- tse_metadata %>%
  filter(
    category != "Infant Study", 
    geo_loc_name_country_calc != "Zimbabwe" # To explain in the ms!
  ) %>%
  drop_na(log_ARG_load, income_group)

filtered_metadata_female_no_outliers <- filtered_metadata %>%
  filter(gender == "Women")

# Define consistent color palette
gender_colors <- c("Women" = "#F8766D", "Men" = "#619CFF")

# Ensure consistent themes and uniform styling
custom_theme <- theme_minimal(24) +
  theme(
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", size = 1),
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(),
    axis.text = element_text(),
    legend.position = "none"
  )

# Function to calculate effect size (r) for Wilcoxon test
calc_effect_size <- function(data, x, y) {
  wilcox_res <- wilcox.test(as.formula(paste(y, "~", x)), data = data, exact = FALSE)
  r <- wilcox_res$statistic / (length(data[[x]]) * length(data[[y]]))
  return(r)
}

# Calculate sample sizes and effect sizes per income_group for both ARG load and ARG diversity
annotations <- filtered_metadata %>%
  group_by(income_group) %>%
  summarise(
    n_Women = sum(gender == "Women"),
    n_Men = sum(gender == "Men"),
    effect_size_ARG_load = calc_effect_size(cur_data(), "gender", "ARG_load"),
    effect_size_shannon_diversity = calc_effect_size(cur_data(), "gender", "shannon_diversity"),
    .groups = 'drop'
  )

# Create a formatted annotation table including both measures
annotation_table <- ggtexttable(
  annotations %>%
    mutate(
      `Effect Size (ARG load)` = round(effect_size_ARG_load, 3),
      `Effect Size (ARG diversity)` = round(effect_size_shannon_diversity, 3)
    ) %>%
    select(
      `Income Group` = income_group,
      `n Women` = n_Women,
      `n Men` = n_Men,
      `Effect Size (ARG load)`,
      `Effect Size (ARG diversity)`
    ),
  rows = NULL,
  theme = ttheme(
    "light",
    base_size = 16,
    padding = unit(c(5, 30), "pt")
  )
) +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )

# Create ARG Load Violin Plot
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

# Create Shannon Diversity Violin Plot
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

# Combine the two violin plots side by side
combined_violin_plots <- plot_grid(
  income_arg_violinplot,
  income_shannon_violinplot,
  labels = c('a', 'b'),
  label_size = 24,
  ncol = 2,
  align = 'hv'
)

# Combine the annotation table with the plots
final_figure <- plot_grid(
  annotation_table,
  combined_violin_plots,
  ncol = 1,
  rel_heights = c(0.3, 1)
)

# Generate publication quality printout
CairoJPEG("RESULTS/FIGURES/Fig3.jpg", width = 1000, height = 700, quality = 100)
print(final_figure)
dev.off()
