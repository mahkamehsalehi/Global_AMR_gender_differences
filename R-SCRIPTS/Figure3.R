library(vegan)
library(ggpubr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(tidyverse)
library(cowplot)
library(mia)

TSE <- readRDS("DATA/TSE_filtered.rds")
tse_metadata <- as.data.frame(colData(TSE))

filtered_metadata <- tse_metadata %>%
  filter(
    category != "Infant Study", 
    geo_loc_name_country_calc != "Zimbabwe" # To explain in the ms!
  )

filtered_metadata <- filtered_metadata %>%
  drop_na(log_ARG_load, income_group)

filtered_metadata_no_outliers <- filtered_metadata %>%
  filter(ARG_load <= 3000) # Why

filtered_metadata_female_no_outliers <- filtered_metadata_no_outliers %>%
  filter(gender == "Women")


# Define consistent color palette
gender_colors <- c("Women" = "#F8766D", "Men" = "#619CFF")

# Ensure consistent themes and uniform styling
custom_theme <- theme_minimal(30) +
  theme(
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", size = 1),
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(),
    axis.text = element_text() ,
    legend.position = "none"
  )

income_arg_violinplot <- ggviolin(
  filtered_metadata_no_outliers,
  x = "gender",
  y = "ARG_load",
  fill = "gender",
  palette = c("Women" = "#F8766D", "Men" = "#619CFF"),
  add = "boxplot",
  add.params = list(fill = "white", width = 0.1)
) +
  stat_compare_means(
    comparisons = list(c("Women", "Men")),
    label = "p.format",
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE
  ) +
  labs(
    x = "",
    y = "ARG load (log RPKM)"
  ) +
  facet_wrap(~income_group) +
  #coord_cartesian(ylim = c(0, 3500)) +
  scale_y_log10() + 
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", size = 1),
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(),
    axis.text = element_text(),
    legend.position = "none"
  ) +
  custom_theme

income_shannon_violinplot <- ggviolin(
  filtered_metadata,
  x = "gender",
  y = "shannon_diversity",
  fill = "gender",
  palette = c("Women" = "#F8766D", "Men" = "#619CFF"),
  add = "boxplot",
  add.params = list(fill = "white", width = 0.1)
  ) +
  stat_compare_means(
    comparisons = list(c("Women", "Men")),
    label = "p.format",
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE
  ) +
  labs(
    x = "",
    y = "ARG diversity (Shannon index)"
  ) +
  facet_wrap(~income_group) +
  coord_cartesian(ylim = c(0, 4)) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", size = 1),
    plot.title = element_text(),    
    axis.title.y = element_text(hjust=20),
    axis.text = element_text()
  ) +
  custom_theme

# leg <- get_legend(income_arg_violinplot + theme(legend.direction="horizontal", legend.title=element_blank()))

# Combine plots with consistent styling and alignment
combined_figure_income <- plot_grid(
  income_arg_violinplot + custom_theme + theme(legend.position = "none"),
  income_shannon_violinplot + custom_theme + theme(legend.position = "none"),
  labels = c('a', 'b'),
  label_size=30,
  ncol = 2,
  align = 'hv',
  rel_widths = c(1, 1),
  rel_heights = c(1, 1)
)

#combined_figure_income <- plot_grid(combined_figure_income,
#                                    leg,
#                                    rel_heights=c(9, 1),
#				    ncol=1)


# This generates publication quality printout:
library(Cairo)
CairoJPEG("RESULTS/FIGURES/Fig3.jpg", width=1000, height=480, quality=100)
print(combined_figure_income)
dev.off()

