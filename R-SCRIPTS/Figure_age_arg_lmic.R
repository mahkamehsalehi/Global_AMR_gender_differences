library(ggpubr)
library(ggpubr)
library(tidyverse)
library(cowplot)
library(rstatix)
library(SummarizedExperiment)

# -----------------------------
# Data Loading and Preprocessing
# -----------------------------

TSE <- readRDS("DATA/TSE_filtered.rds")

tse_metadata <- as.data.frame(colData(TSE))

metadata_lmic <- tse_metadata %>%
  filter(income_group == "LMIC") %>%
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
# Visualization: Boxplots for LMIC
# ---------------------------

# ARG Load Boxplot Without Faceting LMIC
age_arg_boxplot_lmic <- ggplot(metadata_lmic, aes(x = gender, y = log_ARG_load, fill = gender)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.6,
    alpha = 1,
    show.legend = FALSE
  ) +
  facet_wrap(~age_category, scales = "fixed", nrow = 1) +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(x = "Gender", y = "ARG load (natural log RPKM)") +
  theme_minimal() +
  stat_compare_means(
    comparisons = list(
      c("Women", "Men")
    ),
    aes(label = ..p.signif..),
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 10, face = "bold", angle = 45),
    panel.spacing = unit(0.5, "lines")
  )

# Scatterplot for ARG Load vs. Age
age_arg_scatterplot_lmic <- ggplot(metadata_lmic, aes(x = host_age_years, y = log_ARG_load, color = gender)) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  labs(
    x = "Age (years)",
    y = "ARG load (natural log RPKM)",
    color = "Gender"
  ) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"))

# Combine Scatterplots and Boxplots 
combined_arg_age_lmic <- plot_grid(
  age_arg_scatterplot_lmic, 
  age_arg_boxplot_lmic, 
  labels = c("a", "b"), 
  ncol = 2
)
ggsave("RESULTS/FIGURES/arg_age_lmic.png", combined_arg_age_lmic, width = 14, height = 8)
