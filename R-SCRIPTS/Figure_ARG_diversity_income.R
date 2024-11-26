library(vegan)
library(ggpubr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(tidyverse)
library(cowplot)

TSE <- readRDS("DATA/TSE_filtered.rds")
tse_metadata <- as.data.frame(colData(TSE))

filtered_metadata <- tse_metadata %>%
  filter(
    category != "Infant Study", 
    geo_loc_name_country_calc != "Zimbabwe"
  )

filtered_metadata <- filtered_metadata %>%
  drop_na(log_ARG_load, income_group)

# Plot ARG load by gender across income groups
income_arg_boxplot <- ggplot(filtered_metadata, aes(x = gender, y = log_ARG_load, fill = gender)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6, alpha = 1, show.legend = FALSE) +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(x = "Gender", y = "ARG Load (log natural)") +
  stat_compare_means(
    comparisons = list(c("Women", "Men")), 
    label = "p.signif", 
    method = "wilcox.test", 
    p.adjust.method = "BH",
    hide.ns = FALSE)+
  facet_wrap(~income_group) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", size = 1), 
  )

# Plot Shannon diversity by gender across income groups
income_shannon_boxplot <- ggplot(filtered_metadata, aes(x = gender, y = shannon_diversity, fill = gender)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6, alpha = 1, show.legend = FALSE) +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(x = "Gender", y = "Resistome Diversity") +
  theme_minimal(base_family = "Arial") +
  stat_compare_means(aes(x = sex_combined, y = shannon_diversity),
                     comparisons = list(c("Women", "Men")), 
                     label = "p.signif", 
                     method = "wilcox.test", 
                     p.adjust.method = "BH",
                     hide.ns = FALSE) +
  facet_wrap(~income_group) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", size = 1), 
  )


print(combined_figure_income)

# Plot ARG load and Shannon diversity across income groups for women only
income_arg_boxplot_female <- ggplot(filtered_metadata %>% filter(gender == "Women"), 
                                    aes(x = income_group, y = log_ARG_load, fill = "#F8766D")) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6, alpha = 1, show.legend = FALSE) +
  labs(x = "Income Group", y = "ARg Load (log natural)") +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c("HIC", "LMIC")), 
                     label = "p.signif", 
                     method = "wilcox.test", 
                     p.adjust.method = "BH",
                     hide.ns = FALSE) +
  theme(
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", size = 1), 
  )

income_shannon_boxplot_female <- ggplot(filtered_metadata %>% filter(gender == "Women"),
                                        aes(x = income_group, y = shannon_diversity, fill = "#F8766D")) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6, alpha = 1, show.legend = FALSE) +
  labs(x = "Income Group", y = "Resistome Diversity") +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c("HIC", "LMIC")), 
                     label = "p.signif", 
                     method = "wilcox.test", 
                     p.adjust.method = "BH",
                     hide.ns = FALSE) +
  theme(
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", size = 1), 
  )

# Combine and display the female-only figures
combined_figure_income_female <- income_arg_boxplot_female + income_shannon_boxplot_female + plot_layout(ncol = 2)
print(combined_figure_income_female)
ggsave("final/income_female_panel.png", combined_figure_income_female, width = 20, height = 10)


# Combine and display the ARG load and Shannon diversity figures across income groups
combined_figure_income <- income_arg_boxplot + income_shannon_boxplot +
  income_arg_boxplot_female + income_shannon_boxplot_female +
  plot_layout(ncol = 2, nrow = 2) + 
  plot_annotation(tag_levels = 'a')

ggsave("RESULTS/FIGURES/income_panel.png", combined_figure_income, width = 10, height = 8)


#-----------------------------------------------------
#-----------------------------------------------------
#-----------------------------------------------------
# Violin Plot
#-----------------------------------------------------
#-----------------------------------------------------
#-----------------------------------------------------

# Define consistent color palette
gender_colors <- c("Women" = "#F8766D", "Men" = "#619CFF")

# Ensure consistent themes and uniform styling
custom_theme <- theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    strip.background = element_rect(color = "black", size = 1),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10) ,
    legend.position = "none"
  )

# Set consistent y-axis limits
resistome_y_limits <- c(0, 4)
arg_y_limits <- c(0, 13)

# Violin plot for ARG load (Women)
income_arg_violin_female <- ggplot(filtered_metadata_female, 
                                   aes(x = income_group, 
                                       y = log_ARG_load, 
                                       fill = income_group)) +
  geom_violin(trim = FALSE, alpha = 1, color = "black") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_compare_means(comparisons = list(c("HIC", "LMIC")), 
                     label = "p.signif", 
                     method = "wilcox.test", 
                     p.adjust.method = "BH", 
                     hide.ns = FALSE) +
  scale_fill_manual(values = c("HIC" = "#F8766D", "LMIC" = "#F8766D")) +
  labs(x = "Income group", 
       y = "ARG load (natural log)") +
  ylim(arg_y_limits) +
  custom_theme

# Violin plot for Shannon diversity (Women)
income_shannon_violin_female <- ggplot(filtered_metadata_female, 
                                       aes(x = income_group, 
                                           y = shannon_diversity, 
                                           fill = income_group)) +
  geom_violin(trim = FALSE, alpha = 1, color = "black") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_compare_means(comparisons = list(c("HIC", "LMIC")), 
                     label = "p.signif", 
                     method = "wilcox.test", 
                     p.adjust.method = "BH", 
                     hide.ns = FALSE) +
  scale_fill_manual(values = c("HIC" = "#F8766D", "LMIC" = "#F8766D")) +
  labs(x = "Income group", 
       y = "Resistome diversity") +
  ylim(resistome_y_limits) +
  custom_theme

# Combine plots with consistent styling and alignment
combined_figure_income <- plot_grid(
  income_arg_violinplot + custom_theme ,
  income_shannon_violinplot + custom_theme,
  income_arg_violin_female ,
  income_shannon_violin_female,
  labels = c('a', 'b', 'c', 'd'),
  ncol = 2,
  align = 'hv',
  rel_widths = c(1, 1),
  rel_heights = c(1, 1)
)

# Save the final combined figure as a high-resolution image
ggsave("RESULTS/FIGURES/income_panel_violinplot.png", combined_figure_income, width = 12, height = 10, dpi = 300)
