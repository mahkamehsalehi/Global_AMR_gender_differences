library(mia)
library(vegan)
library(ggpubr)
library(patchwork)
library(tidyverse)
library(cowplot)

# ---------------------------
# Data Preparation
# ---------------------------

tse <- readRDS("DATA/TSE_filtered.rds")

tse_metadata <- as.data.frame(colData(tse)) %>% 
  filter(!is.na(income_group))

# ---------------------------
# Base PCoA Plot
# ---------------------------

base_plot <- ggplot(tse_metadata, aes(x = PC1, y = PC2)) +
  geom_point(size = 0.5, alpha = 1) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title.position = "plot",
    plot.title   = element_text(hjust = 0.5, face = "plain", margin = margin(b = 10)),
    axis.title   = element_text( margin = margin(t = 10)),
    axis.text    = element_text( margin = margin(t = 10)),
    legend.title = element_text(),
    legend.text  = element_text(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.line = element_line(color = "black"),
    legend.key.size = unit(0.5, "cm"),
    plot.margin = margin(15, 15, 15, 15)
    
  ) +
  labs(x = "PC1 (14.95%)", y = "PC2 (10.34%)")

# ---------------------------
# PCoA Plot Colored by Gender
# ---------------------------
p_gender <- base_plot +
  aes(color = gender) +
  geom_point(size = 0.5, alpha = 1) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  stat_ellipse(aes(group = gender), type = "norm", level = 0.95, size = 1) +
  labs(title = "PCoA by gender", color = "Gender")

# ---------------------------
# PCoA Plot Colored by Region
# ---------------------------
p_region <- base_plot +
  aes(color = geo_loc_name_country_continent_calc) +
  geom_point(size = 0.8, alpha = 1) +
  stat_ellipse(aes(group = geo_loc_name_country_continent_calc), type = "norm", level = 0.95, size = 1) +
  labs(title = "PCoA by region", color = "Region")

# ---------------------------
# PCoA Plot Colored by Age Category
# ---------------------------
p_age <- base_plot +
  aes(color = age_category_new) +
  geom_point(size = 0.8, alpha = 1) +
  stat_ellipse(aes(group = age_category_new), type = "norm", level = 0.95, size = 0.7) +
  labs(title = "PCoA by age category", color = "Age category")

# ---------------------------
# PCoA Plot Colored by Income Group
# ---------------------------
p_income <- base_plot +
  aes(color = income_group) +
  geom_point(size = 0.8, alpha = 1) +
  stat_ellipse(aes(group = income_group), type = "norm", level = 0.95, size = 0.7) +
  labs(title = "PCoA by income group", color = "Income Group")

# ---------------------------
# Combine All Plots
# ---------------------------
combined_plot <- (p_gender | p_region) / (p_age | p_income) +
  plot_annotation(
    tag_levels = 'a',
    theme = theme(
      plot.margin = margin(20, 20, 20, 20),
      plot.tag = element_text(face = "bold", size = 20),
    )
  )


save_plot("RESULTS/FIGURES/Fig2_PCoA_all.jpg", combined_plot, 
          base_width = 14, base_height = 12, dpi = 600)
