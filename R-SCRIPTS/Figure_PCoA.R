library(mia)
library(vegan)
library(ggpubr)
library(ggpubr)
library(patchwork)
library(tidyverse)

# ---------------------------
# Data Preparation
# ---------------------------

tse <- readRDS("../DATA/TSE_filtered.rds")

# Extract metadata from the TSE object and convert it to a data frame
tse_metadata <- as.data.frame(colData(tse))

# Create age groups (different from the one in the data set) and clean the data
tse_metadata <- tse_metadata %>%
  mutate(
    age_group = case_when(
      host_age_years >= 0 & host_age_years <= 1 ~ "Infant",
      host_age_years > 1 & host_age_years <= 3 ~ "Toddler",
      host_age_years > 3 & host_age_years < 18 ~ "Child",
      host_age_years >= 18 & host_age_years <= 100 ~ "Adult",
      TRUE ~ NA_character_
    ),
    # Order the age categories for consistent plotting
    age_group = factor(age_group, levels = c("Infant", "Toddler", "Child", "Adult"))
  ) %>%
  drop_na(geo_loc_name_country_continent_calc) %>%
  drop_na(age_group)

# ---------------------------
# Gender-Based PCoA Plot
# --------------------------

pcoa_gender_plot <- ggplot(tse_metadata , aes(x = PC1, y = PC2, color = gender)) +
  geom_point(size = 0.8, alpha = 1) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  stat_ellipse(aes(group = gender), type = "norm", level = 0.95, size = 0.7) +
  labs(
    x = paste0("PC1 (14.95%)"),
    y = paste0("PC2 (10.34%)"),
    color = "Gender"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 12),
    axis.line = element_line(color = "black"),
    legend.position = "right"
  )

# ---------------------------
# Geographic-Based PCoA Plot
# ---------------------------

pcoa_geo_plot <- ggplot(tse_metadata, aes(x = PC1, y = PC2, color = age_group)) +
  geom_point(size = 0.8, alpha = 1) +
  stat_ellipse(aes(group = age_group), type = "norm", level = 0.95, size = 0.7) +
  labs(
    x = paste0("PC1 (14.95%)"),
    y = paste0("PC2 (10.34%)"),
    color = "Age category"
  ) +
  scale_color_brewer(palette = "Dark2") + 
  facet_wrap(~geo_loc_name_country_continent_calc) + 
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 12),
    axis.line = element_line(color = "black"),
    legend.position = "right",
    strip.text = element_text(size = 10, face = "bold"),
    panel.spacing = unit(0.5, "lines"),
    strip.background = element_rect(color = "black", size = 0.7
  ))

# ---------------------------
# Combine Plots into a Single Panel
# ---------------------------

combined_plots <- pcoa_geo_plot + (pcoa_gender_plot/plot_spacer()) + 
  plot_layout(
    widths = c(4, 1), 
    ncol = 2
  ) +
  plot_annotation(
    tag_levels = 'a'
  )

# Save the combined figure
ggsave("../RESULTS/FIGURES/pcoa_panel.png", combined_plots, width = 12, height = 8)

