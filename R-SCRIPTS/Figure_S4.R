library(scales)
library(vegan)
library(ggpubr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(tidyverse)
library(cowplot)
library(rstatix)
library(SummarizedExperiment)

# Add legend to one of the plots
age_arg_boxplot_hic <- ggplot(metadata_hic, aes(x = gender, y = ARG_load, fill = gender)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.6,
    alpha = 1
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
    hide.ns = FALSE,
    size = 6,
    label.y = 4
  ) +
  scale_y_continuous(transf="log10",
                     breaks=10^(2:5),
                     labels=trans_format("log10", math_format(10^.x))) +
  theme_minimal(16) +
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_blank(),
    strip.text = element_text(size = 12, angle = 0),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "bottom"
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
    size = 6,
    label.y = 4
  ) +
  scale_y_continuous(transf="log10", 
                     breaks=10^(2:5),
                     labels=trans_format("log10", math_format(10^.x)))
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_blank(),
    strip.text = element_text(size = 12, angle = 0),
    panel.spacing = unit(0.5, "lines")
  )

# Combine the plots and keep a single legend
combined_plot <- (age_arg_boxplot_hic + labs(title = "HIC")) +
  (age_arg_boxplot_lmic + labs(title = "LMIC")) +
  plot_layout(ncol = 1, heights = c(1, 1), guides = "collect") +
  plot_annotation(
    theme = theme(
      plot.title = element_text(),
      plot.subtitle = element_text(),
      legend.position = "bottom"
    )
  )

# Save the combined plot
CairoJPEG("RESULTS/FIGURES/Combined_ARG_load.jpg", width = 940, height = 750, quality = 400)
print(combined_plot)
dev.off()
