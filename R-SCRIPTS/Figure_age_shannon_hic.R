library(vegan)
library(ggpubr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(tidyverse)
library(cowplot)
library(rstatix)
library(SummarizedExperiment)

# -----------------------------
# Data Loading and Preprocessing
# -----------------------------

TSE <- readRDS("DATA/TSE_filtered.rds")

tse_metadata <- as.data.frame(colData(TSE))

metadata_hic <- tse_metadata %>%
  filter(income_group == "HIC") %>%
  filter(!is.na(age_category) & !is.na(log_ARG_load))

# ---------------------------
# Define Custom Plot Theme
# ---------------------------

comparisons <- list(
  c("Infant", "Toddler"),
  c("Toddler", "Child"),
  c("Child", "Young Adult"),
  c("Young Adult", "Middle Adult"),
  c("Middle Adult", "Older Adult")
)

# ---------------------------
# Visualization: Boxplots for HIC
# ---------------------------

# Shannon Diversity Boxplot Without Faceting
age_shannon_boxplot_hic <- ggplot(metadata_hic, aes(x = gender, y = shannon_diversity, fill = gender)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.6,
    alpha = 1,
    show.legend = FALSE
  ) +
  facet_wrap(~age_category, scales = "fixed", nrow = 1) +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(x = "Gender", y = "ARG diversity") +
  theme_minimal() +
  stat_compare_means(
    comparisons = list(
      c("Women", "Men")
    ),
    aes(label = "p.signif"),
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_blank(),
    strip.text = element_text(size = 10, face = "bold", angle = 45),
    panel.spacing = unit(0.5, "lines")
  )

# Scatterplot for Shannon Diversity vs. Age
age_shannon_scatterplot_hic <- ggplot(metadata_hic, aes(x = host_age_years, y = shannon_diversity, color = gender)) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  labs(
    x = "Age (years)",
    y = "ARG diversity",
    color = "Gender"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black")
  )

# Combine Scatterplots and Boxplots 
combined_shannon_age_hic <- plot_grid(
  age_shannon_scatterplot_hic, 
  age_shannon_boxplot_hic, 
  labels = c("a", "b"), 
  ncol = 2
)

# Save the combined figure
ggsave("RESULTS/FIGURES/shannon_age_hic.png", combined_shannon_age_hic, width = 14, height = 8)