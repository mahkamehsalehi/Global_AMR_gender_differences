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
# Define table creation function for effect sizes, CIs, and sample sizes
# -----------------------------
mytable <- function(df, region) {
  # Define correct age category order
  age_order <- c("Infant", "Toddler", "Children", "Teenager", 
                 "Young Adult", "Middle-Aged Adult", "Older Adult", "Oldest Adult")
  
  d <- subset(df, reg == region)
  d$age_group <- factor(d$age_group, levels = age_order)
  
  # Process each age group separately
  results <- list()
  
  for(ag in age_order) {
    d_sub <- d[d$age_group == ag,]
    
    # Skip if age group doesn't exist in this region
    if(nrow(d_sub) == 0) next
    
    # Get sample sizes
    n_sizes <- d_sub %>%
      group_by(gender) %>%
      summarise(n = n(), .groups = "drop") %>%
      pivot_wider(names_from = gender, values_from = n, names_prefix = "n_")
    
    # Try to compute effect size with CI
    tryCatch({
      d_sub$gender <- factor(d_sub$gender, levels = c("Women", "Men"))
      
      # Compute effect size
      eff <- wilcox_effsize(data = d_sub, 
                            ARG_load ~ gender, 
                            ci = TRUE, 
                            conf.level = 0.95)
      
      # Compute p-value
      test <- wilcox.test(ARG_load ~ gender, data = d_sub)
      
      results[[ag]] <- data.frame(
        age_group = ag,
        n_Women = n_sizes$n_Women,
        n_Men = n_sizes$n_Men,
        effsize = eff$effsize,
        conf.low = eff$conf.low,
        conf.high = eff$conf.high,
        p.value = test$p.value
      )
    }, error = function(e) {
      # If effect size calculation fails, store just the sample sizes
      results[[ag]] <- data.frame(
        age_group = ag,
        n_Women = n_sizes$n_Women,
        n_Men = n_sizes$n_Men,
        effsize = NA,
        conf.low = NA,
        conf.high = NA,
        p.value = NA
      )
    })
  }
  
  # Combine all results and adjust p-values
  ann_table <- bind_rows(results) %>%
    mutate(
      p.adj = if(all(is.na(p.value))) rep(NA, n()) else p.adjust(p.value, method = "BH")
    ) %>%
    mutate(
      `Age Group` = age_group,
      `N (Women)` = n_Women,
      `N (Men)` = n_Men,
      `Effect Size (r)` = ifelse(is.na(effsize), "†", round(effsize, 3)),
      `Lower 95% CI` = ifelse(is.na(conf.low), "†", round(conf.low, 3)),
      `Upper 95% CI` = ifelse(is.na(conf.high), "†", round(conf.high, 3)),
      `Adjusted p-value` = ifelse(is.na(p.adj), "†", formatC(p.adj, format = "f", digits = 4))
    ) %>%
    select(
      `Age Group`,
      `N (Women)`,
      `N (Men)`,
      `Effect Size (r)`,
      `Lower 95% CI`,
      `Upper 95% CI`,
      `Adjusted p-value`
    )
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

# --------------------------
# Create header labels for each section
# --------------------------
header_europe <- ggdraw() + 
  draw_label("Europe", fontface = "bold", size = 36, hjust = 0.5) +
  theme(plot.margin = margin(b = 80))
header_na <- ggdraw() + 
  draw_label("North America", fontface = "bold", size = 36, hjust = 0.5) +
  theme(plot.margin = margin(b = 80))

label_c <- ggdraw() + draw_label("c", fontface = "bold", size = 36, x = 0.1, hjust = 0)
label_d <- ggdraw() + draw_label("d", fontface = "bold", size = 36, x = 0.1, hjust = 0)

# --------------------------
# Stack the label, header, and table for each region vertically
# --------------------------
composite_europe <- plot_grid(
  label_c,
  header_europe, 
  table_europe,
  ncol = 1,
  rel_heights = c(0.15, 0.15, 1)
)

composite_na <- plot_grid(
  label_d,
  header_na, 
  table_na,
  ncol = 1,
  rel_heights = c(0.15, 0.15, 1)
)

# --------------------------
# Create grey horizontal separator
# --------------------------
separator <- ggdraw() + 
  draw_line(x = c(0, 1), y = c(0.5, 0.5), color = "grey", size = 2)

# --------------------------
# Combine the tables with separator
# --------------------------
combined_tables <- plot_grid(
  composite_europe,
  separator,
  composite_na,
  ncol = 1,
  rel_heights = c(1, 0.1, 1)
)

# --------------------------
# Create top row with plots
# --------------------------
plots_combined <- plot_grid(
  ps[[1]] + theme(legend.position = "none"),
  ps[[2]] + theme(legend.position = "none"),
  labels = c("a", "b"),
  label_size = 36,
  ncol = 2,
  align = 'hv'
)

plots_with_legend <- plot_grid(
  plots_combined, 
  leg, 
  ncol = 1, 
  rel_heights = c(1, 0.1)
)

# --------------------------
# Combine everything into final figure
# --------------------------
final_figure <- plot_grid(
  plots_with_legend,
  separator,
  combined_tables,
  ncol = 1,
  rel_heights = c(1, 0.05, 1.2)
)

# --------------------------
# Save the Final Figure
# --------------------------
CairoJPEG("RESULTS/FIGURES/Fig4.jpg", width = 2000, height = 2000, quality = 100)
print(final_figure)
dev.off()