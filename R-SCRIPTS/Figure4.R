library(vegan)
library(ggpubr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(tidyverse)
library(cowplot)
library(rstatix)
library(SummarizedExperiment)
library(dplyr)

# -----------------------------
# Data Loading and Preprocessing
# -----------------------------

TSE <- readRDS("DATA/TSE_filtered.rds")

tse_metadata <- as.data.frame(colData(TSE))

metadata_hic <- tse_metadata %>%
  filter(income_group == "HIC") %>%
  filter(!is.na(age_category_new) & !is.na(log_ARG_load))

df <- metadata_hic
df$age_group <- df$age_category_new
df$reg <- factor(df$geo_loc_name_country_continent_calc)

###
myplot <- function(df, region, ylim_list, sig_y_offset = -3.2) {
  d <- subset(df, reg == region)
  
  # Get the y-limits for the region
  ylim_max <- ylim_list[[region]]
  
  p <- ggplot(d, 
              aes(x = gender, y = ARG_load, fill = gender)) + 
    geom_boxplot(
      position = position_dodge(width = 0.8),
      outlier.shape = NA,
      width = 0.6,
      alpha = 1,
      show.legend = TRUE
    ) +
    facet_wrap(~age_group, scales = "fixed", nrow = 1) +
    scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
    labs(x = "Age groups", y = "ARG load (RPKM)", title = region) +
    # coord_cartesian(ylim = c(4.2, ylim_max)) + # Set y-limits here
    stat_compare_means(
      comparisons = list(c("Women", "Men")),
      aes(label = ..p.signif.., size = 7),
      method = "wilcox.test",
      p.adjust.method = "BH",
      hide.ns = FALSE,
      label.y = ylim_max - sig_y_offset  # Offset significance labels
    ) +
    scale_y_continuous(trans="log2", breaks=c(100, 1000, 10000)) +             
    theme_minimal(16) +
    theme(
      axis.line = element_line(color = "black"),
      axis.text.x = element_blank(),
      strip.text = element_text(size = 10, face = "bold", angle = 0),
      panel.spacing = unit(0.5, "lines")
    )
  
  return(p)
}



##
# y-limits for each region
ylim_list <- list(
  "Europe" = 9.5,
  "North America" = 9
)

# Generate plots for each region
regs <- c("Europe", "North America")
ps <- lapply(regs, function(region) {
  myplot(df, region, ylim_list = ylim_list)
})

# Extract legend from the first plot
leg <- get_legend(ps[[1]] + theme(legend.direction = "horizontal", legend.title = element_blank()))

# Combine plots and legend
p <- plot_grid(
  ps[[1]] + theme(legend.position = "none"),
  ps[[2]] + theme(legend.position = "none"),
  leg,
  rel_heights = c(0.45, 0.45, 0.1),
  nrow = 3
)

print(p)

library(Cairo)
CairoJPEG("RESULTS/FIGURES/Fig4.jpg", width=800, height=700, quality=100)
print(p)
dev.off()
