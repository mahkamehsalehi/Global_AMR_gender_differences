library(vegan)
library(ggpubr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(tidyverse)

TSE <- readRDS("DATA/TSE_filtered.rds")
tse_metadata <- as.data.frame(colData(TSE))

filtered_metadata <- tse_metadata %>%
  filter(
    category != "Infant Study", 
    geo_loc_name_country_calc != "Zimbabwe"
  )


# Define a common theme for plots
common_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.position = "none",
    axis.line = element_line(color = "black"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 12, face = "bold")
  )

# Filter and prepare data for usage analysis with custom 'Usage_group' categorization
usage_filtered_metadata <- tse_metadata %>%
  mutate(
    Usage_group = case_when(
      Usage > 11  ~ "Above 11",
      Usage <= 11 ~ "11 or Below",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(Usage_group = factor(Usage_group, levels = c("11 or Below", "Above 11"))) %>%
  drop_na(gender, log_ARG_load, Usage_group, income_group)

metadata_hic <- usage_filtered_metadata %>%
  filter(income_group == "HIC")

metadata_lmic <- usage_filtered_metadata%>%
  filter(income_group == "LMIC")

# Plot ARG load by gender across usage groups
usage_arg_boxplot_hic <- ggplot(metadata_hic, aes(x = gender,
                                              y = log_ARG_load,
                                              fill = gender)) +
  geom_boxplot(position = position_dodge(width = 0.8), 
               outlier.shape = NA, 
               width = 0.6, 
               alpha = 1, 
               show.legend = FALSE) +
  labs(x = "Gender", y = "ARG Load (log natural)") +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  stat_compare_means(aes(x = gender, 
                         y = log_ARG_load),
                     comparisons = list(c("Women", "Men")), 
                     label = "p.signif", 
                     method = "wilcox.test", 
                     p.adjust.method = "BH",
                     hide.ns = FALSE) +
  facet_wrap(~Usage_group) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        strip.background = element_rect(color = "black", size = 0.7)
        )

 # Plot Shannon diversity by gender across usage groups
usage_shannon_boxplot_hic <- ggplot(metadata_hic, 
                                     aes(x = gender, 
                                         y = shannon_diversity, 
                                         fill = gender)) +
  geom_boxplot(position = position_dodge(width = 0.8), 
               outlier.shape = NA, 
               width = 0.6, 
               alpha = 1, 
               show.legend = FALSE) +
  labs(x = "Gender", y = "Resistome Diversity") +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  stat_compare_means(aes(x = sgender, y = shannon_diversity),
                     comparisons = list(c("Women", "Men")), 
                     label = "p.signif", 
                     method = "wilcox.test", 
                     p.adjust.method = "BH",
                     hide.ns = FALSE) +
  facet_wrap(~Usage_group) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        strip.background = element_rect(color = "black", size = 0.7)
  )


# Plot ARG load by usage for women only
usage_arg_boxplot_female_hic <- ggplot(metadata_hic %>% filter(gender == "Women"), 
                                   aes(x = Usage_group, y = log_ARG_load, fill = "#F8766D")) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6, alpha = 1, show.legend = FALSE) +
  labs(x = "Antibiotic Usage (DDD)", y = "ARG Load (log natural)") +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c("11 or Below", "Above 11")), label = "p.signif", method = "wilcox.test", hide.ns = FALSE) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        strip.background = element_rect(color = "black", size = 0.7)
  )

# Plot Shannon diversity by usage for women only
usage_shannon_boxplot_female_hic <- ggplot(metadata_hic %>% filter(gender == "Women"),
                                       aes(x = Usage_group, y = shannon_diversity, fill = "#F8766D")) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6, alpha = 1, show.legend = FALSE) +
  labs(x = "Antibiotic Usagec(DDD)", y = "Resistome Diversity") +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c("11 or Below", "Above 11")), label = "p.signif", method = "wilcox.test", hide.ns = FALSE) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        strip.background = element_rect(color = "black", size = 0.7)
  )

# Combine and display usage plots for women only
combined_figure_usage_female <- usage_arg_boxplot_both + usage_shannon_boxplot_both +
  usage_arg_boxplot_female + usage_shannon_boxplot_female +
  plot_layout(ncol = 2, nrow = 2) +
  plot_annotation(tag_levels = 'a')
ggsave("RESULTS/FIGURES/usage_female_panel.pdf", combined_figure_usage_female, width = 10, height = 6)




# Plot ARG load by gender across usage groups
usage_arg_boxplot_lmic <- ggplot(metadata_lmic, aes(x = gender,
                                              y = log_ARG_load,
                                              fill = gender)) +
  geom_boxplot(position = position_dodge(width = 0.8), 
               outlier.shape = NA, 
               width = 0.6, 
               alpha = 1, 
               show.legend = FALSE) +
  labs(x = "Gender", y = "ARG Load (log natural)") +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  stat_compare_means(aes(x = gender, 
                         y = log_ARG_load),
                     comparisons = list(c("Women", "Men")), 
                     label = "p.signif", 
                     method = "wilcox.test", 
                     p.adjust.method = "BH",
                     hide.ns = FALSE) +
  facet_wrap(~Usage_group) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        strip.background = element_rect(color = "black", size = 0.7)
  )

# Plot Shannon diversity by gender across usage groups
usage_shannon_boxplot_lmic <- ggplot(metadata_lmic, 
                                aes(x = gender, 
                                    y = shannon_diversity, 
                                    fill = gender)) +
  geom_boxplot(position = position_dodge(width = 0.8), 
               outlier.shape = NA, 
               width = 0.6, 
               alpha = 1, 
               show.legend = FALSE) +
  labs(x = "Gender", y = "Resistome Diversity") +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  stat_compare_means(aes(x = sgender, y = shannon_diversity),
                     comparisons = list(c("Women", "Men")), 
                     label = "p.signif", 
                     method = "wilcox.test", 
                     p.adjust.method = "BH",
                     hide.ns = FALSE) +
  facet_wrap(~Usage_group) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        strip.background = element_rect(color = "black", size = 0.7)
  )


# Plot ARG load by usage for women only
usage_arg_boxplot_female_lmic <- ggplot(metadata_lmic %>% filter(gender == "Women"), 
                                       aes(x = Usage_group, y = log_ARG_load, fill = "#F8766D")) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6, alpha = 1, show.legend = FALSE) +
  labs(x = "Antibiotic Usage (DDD)", y = "ARG Load (log natural)") +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c("11 or Below", "Above 11")), label = "p.signif", method = "wilcox.test", hide.ns = FALSE) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        strip.background = element_rect(color = "black", size = 0.7)
  )

# Plot Shannon diversity by usage for women only
usage_shannon_boxplot_female_lmic <- ggplot(metadata_lmic %>% filter(gender == "Women"),
                                           aes(x = Usage_group, y = shannon_diversity, fill = "#F8766D")) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6, alpha = 1, show.legend = FALSE) +
  labs(x = "Antibiotic Usage (DDD)", y = "Resistome Diversity") +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c("11 or Below", "Above 11")), label = "p.signif", method = "wilcox.test", hide.ns = FALSE) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        strip.background = element_rect(color = "black", size = 0.7)
  )


combined_figure_usage_arg <- usage_arg_boxplot_hic + usage_arg_boxplot_lmic +
  usage_arg_boxplot_female_hic + usage_arg_boxplot_female_lmic +
  plot_layout(ncol = 2, nrow = 2) + 
  plot_annotation(tag_levels = 'a')

ggsave("RESULTS/FIGURES/usage_arg_panel.png", combined_figure_usage_arg, width = 10, height = 8)

combined_figure_usage_shannon <- usage_shannon_boxplot_hic + usage_shannon_boxplot_lmic +
  usage_shannon_boxplot_female_hic + usage_shannon_boxplot_female_lmic +
  plot_layout(ncol = 2, nrow = 2) + 
  plot_annotation(tag_levels = 'a')

ggsave("RESULTS/FIGURES/usage_resistome_panel.png", combined_figure_usage_shannon, width = 10, height = 8)
