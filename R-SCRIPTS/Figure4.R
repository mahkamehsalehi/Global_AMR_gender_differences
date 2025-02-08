# Load necessary libraries
library(vegan)
library(ggpubr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(cowplot)
library(rstatix)
library(SummarizedExperiment)
library(dplyr)
library(Cairo)

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

# -----------------------------
# Define plotting function (violin plot with significance annotations)
# -----------------------------
myplot <- function(df, region, ylim_list, sig_y_offset = -3.2) {
  d <- subset(df, reg == region)
  ylim_max <- ylim_list[[region]]
  
  p <- ggplot(d, aes(x = gender, y = ARG_load, fill = gender)) + 
    geom_violin(position = position_dodge(width = 0.8),
                width = 1,
                alpha = 1,
                show.legend = TRUE) +
    geom_boxplot(width = 0.1, 
                 position = position_dodge(width = 0.8), 
                 outlier.shape = NA,
                 fill = "white",
                 colour = "black") +
    facet_wrap(~age_group, scales = "fixed", nrow = 1) +  
    scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
    labs(x = "Age groups", y = "ARG load (RPKM)", title = region) +
    stat_compare_means(
      comparisons = list(c("Women", "Men")),
      label = "p.signif",
      method = "wilcox.test",
      p.adjust.method = "BH",
      hide.ns = FALSE,
      label.y = ylim_max - sig_y_offset,
      size = 8,
      fontface = "plain"
    ) +
    scale_y_continuous(trans = "log2", breaks = c(100, 1000, 10000)) +             
    theme_minimal(28) +
    theme(
      axis.line = element_line(color = "black"),
      axis.text.x = element_blank(),
      strip.text = element_text(size = 24, face = "bold", angle = 45, margin = margin(t = 50, b = 50)),
      panel.spacing = unit(2, "lines"),
      panel.spacing.y = unit(4, "lines"),
      plot.title = element_text(size = 40, face = "plain")
    ) +
    coord_cartesian(clip = "off")
  
  return(p)
}

# -----------------------------
# Define table creation function for effect sizes and sample sizes
# -----------------------------
mytable <- function(df, region) {
  d <- subset(df, reg == region)
  
  # Compute Wilcoxon test statistics (grouped by age group)
  test_res <- d %>% 
    group_by(age_group) %>% 
    wilcox_test(ARG_load ~ gender) %>% 
    adjust_pvalue(method = "BH")
  
  # Compute effect sizes using rstatix's helper function
  eff_res <- d %>% 
    group_by(age_group) %>% 
    wilcox_effsize(ARG_load ~ gender)
  
  # Compute sample sizes (n) per age group and gender
  n_res <- d %>% 
    group_by(age_group, gender) %>% 
    summarise(n = n(), .groups = "drop") %>% 
    pivot_wider(names_from = gender, values_from = n, names_prefix = "n_")
  
  # Merge all computed statistics by age group
  ann <- left_join(test_res, eff_res, by = "age_group") %>% 
    left_join(n_res, by = "age_group")
  
  # Format p-adjusted values
  ann_table <- ann %>% 
    mutate(`Effect Size (r)` = round(effsize, 3),
           `p-adj` = formatC(p.adj, format = "f", digits = 4)) %>%  # Converts to fixed decimal format
    select(`Age Group` = age_group, `n Women` = n_Women, `n Men` = n_Men, `Effect Size (r)`, `p-adj`)
  
  return(ann_table)
}


# -----------------------------
# Generate and Combine Plots and Tables
# -----------------------------
# Define y-limits for each region
ylim_list <- list(
  "Europe" = 9.5,
  "North America" = 9
)

# List of regions to plot
regs <- c("Europe", "North America")

# Generate violin plots for each region
ps <- lapply(regs, function(region) {
  myplot(df, region, ylim_list = ylim_list) +
    theme(legend.position = )
})

# Extract a common legend from the first plot
legend_plot <- ps[[1]] + 
  theme(legend.position = "bottom", 
        legend.direction = "horizontal", 
        legend.title = element_blank(),
        text = element_text(size = 28)) +
  guides(fill = guide_legend(nrow = 1))

leg <- get_legend(legend_plot)

# Generate tables for each region
table_europe <- ggtexttable(mytable(df, "Europe"), rows = NULL, 
                            theme = ttheme("light", base_size = 28, padding = unit(c(10, 20), "pt")))
table_na <- ggtexttable(mytable(df, "North America"), rows = NULL, 
                        theme = ttheme("light", base_size = 28, padding = unit(c(10, 20), "pt")))

# Combine the two region plots side by side
plots_combined <- plot_grid(
  ps[[1]] + theme(legend.position = "none"),
  ps[[2]] + theme(legend.position = "none"),
  labels = c("a", "b"), label_size = 36,
  ncol = 2
)
plots_with_legend <- plot_grid(plots_combined, leg, ncol = 1, rel_heights = c(1, 0.1))

# Combine the two tables in one row
tables_combined <- plot_grid(
  table_europe, 
  table_na,
  labels = c("c", "d"), label_size = 36,
  ncol = 2, 
  rel_widths = c(1, 1))

# Combine the plot section on top and the tables below, with a spacer between
final_figure <- plot_grid(plots_with_legend, plot_spacer(), tables_combined,
                          ncol = 1,
                          rel_heights = c(1, 0.1, 0.5))

# Display the final composite figure
print(final_figure)

# Save the final figure
CairoJPEG("RESULTS/FIGURES/Fig4.jpg", width = 2000, height = 1600, quality = 100)
print(final_figure)
dev.off()
