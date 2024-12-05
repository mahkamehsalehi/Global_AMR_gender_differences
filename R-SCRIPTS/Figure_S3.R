library(vegan)
library(ggpubr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(tidyverse)
library(cowplot)
library(mia)
library(scales)

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

#-----------------------------------------------------
# Violin Plot
#-----------------------------------------------------
#-----------------------------------------------------
#-----------------------------------------------------

# Define consistent color palette
gender_colors <- c("Women" = "#F8766D", "Men" = "#619CFF")

# Ensure consistent themes and uniform styling
custom_theme <- theme_minimal(18) +
  theme(
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", size = 1),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14) ,
    legend.position = "none"
  )

# Violin plot for ARG load (Women)
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
  scale_y_continuous(transf="log10",
                     breaks=10^(2:5),
                     labels=trans_format("log10", math_format(10^.x))) +
  custom_theme

# Violin plot for Shannon diversity (Women)
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
       y = "ARG diversity") +
  coord_cartesian(ylim = c(0, 4))+
  custom_theme

# Combine plots with consistent styling and alignment
combined_figure_income_female <- plot_grid(
  income_arg_violin_female + custom_theme,
  income_shannon_violin_female + custom_theme,
  labels = c('a', 'b'),
  ncol = 2,
  align = 'hv',
  rel_widths = c(1, 1),
  rel_heights = c(1, 1)
)

# Save the final combined figure as a high-resolution image
ggsave("RESULTS/FIGURES/income_panel_violinplot_female.png", combined_figure_income_female, width = 10, height = 6, dpi = 300)
