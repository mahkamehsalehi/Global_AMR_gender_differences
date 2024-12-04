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

TSE <- readRDS("../DATA/TSE_filtered.rds")

tse_metadata <- as.data.frame(colData(TSE))

metadata_hic <- tse_metadata %>%
  filter(income_group == "HIC") %>%
  filter(!is.na(age_category_new) & !is.na(ARG_load))

metadata_lmic <- tse_metadata %>%
  filter(income_group == "LMIC") %>%
  filter(!is.na(age_category_new) & !is.na(ARG_load))

comparisons <- list(
  c("Infant", "Toddler"),
  c("Toddler", "Children"),
  c("Children", "Teenager"),
  c("Teenager", "Young Adult"),
  c("Young Adult", "Middle-Aged Adult"),
  c("Middle-Aged Adult", "Older Adult"),
  c("Older Adult", "Oldest Adult")
)

# ARG Load Boxplot Without Faceting HIC
age_arg_boxplot_hic <- ggplot(metadata_hic, aes(x = gender, y = ARG_load, fill = gender)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.6,
    alpha = 1,
    show.legend = FALSE
  ) +
  facet_wrap(~age_category_new, scales = "fixed", nrow = 1) +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(x = "Gender", y = "ARG load (RPKM)") +
  stat_compare_means(
    comparisons = list(
      c("Women", "Men")
    ),
    aes(label = ..p.signif..),
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE, size = 6 
  ) +
  scale_y_log10(labels = scales::comma) + 
  theme_minimal(16) +
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_blank(),
    strip.text = element_text(size = 12, angle = 0),
    panel.spacing = unit(0.5, "lines")
  )

age_arg_boxplot_lmic <- ggplot(metadata_lmic, aes(x = gender, y = ARG_load, fill = gender)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.6,
    alpha = 1,
    show.legend = FALSE
  ) +
  facet_wrap(~age_category_new, scales = "fixed", nrow = 1) +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(x = "Gender", y = "ARG load (RPKM)") +
  theme_minimal(16) +
  stat_compare_means(
    comparisons = list(
      c("Women", "Men")
    ),
    aes(label = ..p.signif..),
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE,
    size = 6 
  ) +
  scale_y_log10(labels = scales::comma) + 
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_blank(),
    strip.text = element_text(size = 12, angle = 0),
    panel.spacing = unit(0.5, "lines")
  )

income_level <- c("HICs", "LMICs")

library(patchwork)

# Combine the plots in two rows
combined_plot <- (age_arg_boxplot_hic + labs(title="HIC"))+
                 (age_arg_boxplot_lmic + labs(title="LMIC"))+
  plot_layout(ncol = 1, heights = c(1, 1)) + 
  plot_annotation(
    theme = theme(plot.title = element_text(),
                  plot.subtitle = element_text())
  )

# Save the combined plot
library(Cairo)
CairoJPEG("../RESULTS/FIGURES/Combined_ARG_Load.jpg", width = 940, height = 700, quality = 400)
print(combined_plot)
dev.off()

