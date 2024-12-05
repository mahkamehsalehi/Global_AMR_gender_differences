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

# Filter and prepare data for usage analysis with custom 'Usage_group' categorization
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

metadata_hic <- usage_filtered_metadata %>%
  filter(income_group == "HIC")

metadata_lmic <- usage_filtered_metadata%>%
  filter(income_group == "LMIC")

# Plot ARG load by gender across usage groups
usage_arg_boxplot_hic <- ggplot(metadata_hic, aes(x = gender,
                                              y = ARG_load,
                                              fill = gender)) +
  geom_boxplot(position = position_dodge(width = 0.8), 
               outlier.shape = NA, 
               width = 0.6, 
               alpha = 1, 
               show.legend = FALSE) +
  labs(x = "Gender", y = "ARG load (RPKM)") +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  stat_compare_means(aes(x = gender, 
                         y = log_ARG_load),
                     comparisons = list(c("Women", "Men")), 
                     label = "p.format", 
                     method = "wilcox.test", 
                     p.adjust.method = "BH",
                     hide.ns = FALSE,
                     size = 5,
                     label.y = 4) +
  facet_wrap(~Usage_group) +
  scale_y_continuous(transf="log10",
                     breaks=10^(2:5),
                     labels=trans_format("log10", math_format(10^.x))) +
  theme_minimal(16) +
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
  labs(x = "Gender", y = "ARG diversity (Shannon index)") +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  stat_compare_means(aes(x = sgender, y = shannon_diversity),
                     comparisons = list(c("Women", "Men")), 
                     label = "p.format", 
                     method = "wilcox.test", 
                     p.adjust.method = "BH",
                     hide.ns = FALSE,
                     size = 5,
                     label.y = 3.2) +
  facet_wrap(~Usage_group) +
  theme_minimal(16) +
  theme(axis.line = element_line(color = "black"),
        strip.background = element_rect(color = "black", size = 0.7)
  )


# Plot ARG load by usage for women only
usage_arg_boxplot_female_hic <- ggplot(metadata_hic %>% filter(gender == "Women"), 
                                   aes(x = Usage_group, y = ARG_load, fill = "#F8766D")) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6, alpha = 1, show.legend = FALSE) +
  labs(x = "Antibiotic Use (DDD)", y = "ARG load (RPKM)") +
  theme_minimal(16) +
  scale_y_continuous(transf="log10",
                     breaks=10^(2:5),
                     labels=trans_format("log10", math_format(10^.x))) +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c("Low", "High")), label = "p.format", method = "wilcox.test", hide.ns = FALSE,
                     size =5, label.y = 4) +
  theme_minimal(16) +
  theme(axis.line = element_line(color = "black"),
        strip.background = element_rect(color = "black", size = 0.7)
  )

# Plot Shannon diversity by usage for women only
usage_shannon_boxplot_female_hic <- ggplot(metadata_hic %>% filter(gender == "Women"),
                                       aes(x = Usage_group, y = shannon_diversity, fill = "#F8766D")) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6, alpha = 1, show.legend = FALSE) +
  labs(x = "Antibiotic Use(DDD)", y = "ARG diversity (Shannon index)") +
  theme_minimal(16) +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c("Low", "High")), label = "p.format", method = "wilcox.test", hide.ns = FALSE,  size = 5, label.y = 3.2) +
  theme_minimal(16) +
  theme(axis.line = element_line(color = "black"),
        strip.background = element_rect(color = "black", size = 0.7)
  )


# Plot ARG load by gender across usage groups
usage_arg_boxplot_lmic <- ggplot(metadata_lmic, aes(x = gender,
                                              y = ARG_load,
                                              fill = gender)) +
  geom_boxplot(position = position_dodge(width = 0.8), 
               outlier.shape = NA, 
               width = 0.6, 
               alpha = 1, 
               show.legend = FALSE) +
  labs(x = "Gender", y = "ARG load (RPKM)") +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  stat_compare_means(aes(x = gender, 
                         y = log_ARG_load),
                     comparisons = list(c("Women", "Men")), 
                     label = "p.format", 
                     method = "wilcox.test", 
                     p.adjust.method = "BH",
                     hide.ns = FALSE,
                     size = 5,
                     label.y = 4.5) +
  facet_wrap(~Usage_group) +
  scale_y_continuous(transf="log10",
                     breaks=10^(2:5),
                     labels=trans_format("log10", math_format(10^.x))) +
  theme_minimal(16) +
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
  labs(x = "Gender", y = "ARG diversity (Shannon index") +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  stat_compare_means(aes(x = sgender, y = shannon_diversity),
                     comparisons = list(c("Women", "Men")), 
                     label = "p.format", 
                     method = "wilcox.test", 
                     p.adjust.method = "BH",
                     hide.ns = FALSE,
                     size = 5,
                     label.y = 3) +
  facet_wrap(~Usage_group) +
  theme_minimal(16) +
  theme(axis.line = element_line(color = "black"),
        strip.background = element_rect(color = "black", size = 0.7)
  )


# Plot ARG load by usage for women only
usage_arg_boxplot_female_lmic <- ggplot(metadata_lmic %>% filter(gender == "Women"), 
                                       aes(x = Usage_group, y = ARG_load, fill = "#F8766D")) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6, alpha = 1, show.legend = FALSE) +
  labs(x = "Antibiotic Use (DDD)", y = "ARG load (RPKM)") +
  theme_minimal(16) +
  theme(legend.position = "none") +
  scale_y_continuous(transf="log10",
                     breaks=10^(2:5),
                     labels=trans_format("log10", math_format(10^.x))) +
  stat_compare_means(comparisons = list(c("Low", "High")), label = "p.format", method = "wilcox.test", hide.ns = FALSE, size = 5, label.y = 4.5) +
  theme_minimal(16) +
  theme(axis.line = element_line(color = "black"),
        strip.background = element_rect(color = "black", size = 0.7)
  )

# Plot Shannon diversity by usage for women only
usage_shannon_boxplot_female_lmic <- ggplot(metadata_lmic %>% filter(gender == "Women"),
                                           aes(x = Usage_group, y = shannon_diversity, fill = "#F8766D")) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.6, alpha = 1, show.legend = FALSE) +
  labs(x = "Antibiotic Use (DDD)", y = "ARG diversity (Shannon index)") +
  theme_minimal(16) +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c("Low", "High")), label = "p.format", method = "wilcox.test", hide.ns = FALSE, size = 5, label.y = 3) +
  theme_minimal(16) +
  theme(axis.line = element_line(color = "black"),
        strip.background = element_rect(color = "black", size = 0.7)
  )


combined_figure_usage_arg <- usage_arg_boxplot_hic  + labs(title = "HICs") + 
  usage_arg_boxplot_lmic + labs(title = "LMICs") +
  usage_arg_boxplot_female_hic + labs(title = "Women in HICs") +
  usage_arg_boxplot_female_lmic +  labs(title = "Women in LMICs") +
  plot_layout(ncol = 2, nrow = 2) + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.title = element_text(hjust = 0.5))

ggsave("RESULTS/FIGURES/usage_arg_panel.png", combined_figure_usage_arg, width = 10, height = 8)

combined_figure_usage_shannon <- usage_shannon_boxplot_hic + usage_shannon_boxplot_lmic +
  usage_shannon_boxplot_female_hic + usage_shannon_boxplot_female_lmic +
  plot_layout(ncol = 2, nrow = 2) + 
  plot_annotation(tag_levels = 'a')

ggsave("RESULTS/FIGURES/usage_resistome_panel.png", combined_figure_usage_shannon, width = 10, height = 8)


