# ---------------------------
# Load required libraries
# ---------------------------
library(tidyverse)
library(mia)
library(vegan)
library(ggpubr)
library(patchwork)
library(SummarizedExperiment)
library(ggthemes)

# ---------------------------
# Data Preparation
# ---------------------------
# Read the filtered TSE object for plotting metadata (PCoA plot)
tse_filtered <- readRDS("DATA/TSE_filtered.rds")
tse_metadata <- as.data.frame(colData(tse_filtered)) %>% 
  filter(!is.na(income_group))

# Read the main TSE object (for counts and gene classes)
tse <- readRDS("DATA/TSE.rds")
counts <- assay(tse)
gene_classes <- rowData(tse)$Class

# Process and recode metadata
meta <- as.data.frame(colData(tse)) %>%
  mutate(
    sex_combined = case_when(
      sex_combined == "female" ~ "Women",
      sex_combined == "male" ~ "Men",
      TRUE ~ NA_character_
    ),
    Income_Group = case_when(
      World_Bank_Income_Group == "High income" ~ "HIC",
      World_Bank_Income_Group %in% c("Low income", "Lower middle income", "Upper middle income") ~ "LMIC",
      TRUE ~ NA_character_
    )
  ) %>%
  drop_na(sex_combined, log10_ARG_load, Income_Group)

# ---------------------------
# Gene Class Analysis: Top 5 Gene Classes per Group (Normalized)
# ---------------------------
results <- list()

for (sex in c("Men", "Women")) {
  for (income in c("HIC", "LMIC")) {
    idx <- which(meta$sex_combined == sex & meta$Income_Group == income)
    
    if (length(idx) == 0) {
      cat("No samples for", sex, "in", income, "\n")
      next
    }
    
    # Calculate number of samples in this group
    n_samples <- length(idx)
    
    # Calculate average abundance per sample
    group_counts <- rowSums(counts[, idx, drop = FALSE]) / n_samples
    class_totals <- tapply(group_counts, gene_classes, sum, na.rm = TRUE)
    top5 <- sort(class_totals, decreasing = TRUE)[1:5]
    
    df_temp <- data.frame(
      sex = sex,
      income = income,
      gene_class = names(top5),
      abundance = as.numeric(top5),
      n_samples = n_samples  # Store sample size for reference
    )
    
    results[[paste(sex, income, sep = "_")]] <- df_temp
  }
}

df_bar <- do.call(rbind, results)
df_bar$sex <- factor(df_bar$sex, levels = c("Women", "Men"))

# ---------------------------
# Plotting
# ---------------------------
# Bar Plot: Average abundance of top 5 gene classes per group
class_plot <- ggplot(df_bar, aes(x = reorder(gene_class, abundance), y = abundance, fill = gene_class)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(sex), cols = vars(income), scales = "free_y") +
  coord_flip() +
  labs(
    x = "Gene class",
    y = "Average abundance per sample",
    tag = "a"
  ) +
  scale_fill_viridis_d(option = "mako") +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "plain", size = 16),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 22, face = "plain"),
    plot.tag = element_text(size = 24, face = "plain"),
    legend.position = "none"
  )

# PCoA Plot remains the same
pcoa_plot <- ggplot(tse_metadata, aes(x = PC1, y = PC2, color = gender)) +
  geom_point(size = 1, alpha = 1) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(
    x = "PC1 (14.95%)",
    y = "PC2 (10.34%)",
    color = "Gender",
    tag = "b"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.tag = element_text(size = 24, face = "bold"),
    axis.title = element_text(size = 22, face = "plain"),
    axis.text = element_text(size = 20),
    strip.text = element_text(face = "plain", size = 16),
    legend.title = element_text(size = 22, face = "plain"),
    legend.text = element_text(size = 20),
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  facet_grid(income_group ~ age_category_new)

# Combine the plots using patchwork
combined_plot <- class_plot / pcoa_plot +
  plot_layout(heights = c(1, 1.2))

# Save the combined plot
ggsave("RESULTS/FIGURES/gene_class_pcoa_normalized.png", combined_plot, width = 22, height = 12, dpi = 300)
