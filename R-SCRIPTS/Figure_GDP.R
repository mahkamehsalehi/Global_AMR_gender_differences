library(ggpubr)
library(ggpubr)
library(patchwork)
library(tidyverse)
library(hexbin)

TSE <- readRDS("DATA/TSE_filtered.rds")
tse_metadata <- as.data.frame(colData(TSE))

filtered_metadata <- tse_metadata %>%
  filter(
    category != "Infant Study", 
    geo_loc_name_country_calc != "Zimbabwe"
  )

filtered_metadata <- filtered_metadata %>%
  drop_na(log_ARG_load, GDP_per_head)


# Create the Hexbin Plot
arg_gdp_plot <- ggplot(filtered_metadata, aes(x = GDP_per_head, y = log_ARG_load)) +
  stat_binhex(bins = 30, aes(fill = after_stat(count)), color = "white") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  geom_smooth(method = "loess", color = "black", se = TRUE, size = 1) + 
  theme_minimal() +
  labs(
    x = "GDP per Head",
    y = "ARG Load (natural log)",
    fill = "Count"
  ) +
  theme(
    axis.line = element_line(color = "black")
  )

shannon_gdp_plot <- ggplot(filtered_metadata, aes(x = GDP_per_head, y = shannon_diversity)) +
  stat_binhex(bins = 30, aes(fill = after_stat(count)), color = "white") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  geom_smooth(method = "loess", color = "black", se = TRUE, size = 1) + 
  theme_minimal() +
  labs(
    x = "GDP per Head",
    y = "Resistome Diversity",
    fill = "Count"
  ) +
  theme(
    axis.line = element_line(color = "black")
  )

ggsave("RESULTS/FIGURES/arg_gdp.png", arg_gdp_plot, width = 10, height = 6)
ggsave("RESULTS/FIGURES/shannon_gdp.png", shannon_gdp_plot, width = 10, height = 6)

combined_gdp <- (arg_gdp_plot + shannon_gdp_plot) +
  plot_layout(ncol = 2) + 
  plot_annotation(
    tag_levels = 'a')
ggsave("RESULTS/FIGURES/gdp_comb.png", combined_gdp, width = 10, height = 6)
