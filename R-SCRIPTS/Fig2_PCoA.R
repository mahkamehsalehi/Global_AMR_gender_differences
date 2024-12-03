library(mia)
library(vegan)
library(ggpubr)
library(ggpubr)
library(patchwork)
library(tidyverse)

# ---------------------------
# Data Preparation
# ---------------------------

tse <- readRDS("DATA/TSE_filtered.rds")

# Extract metadata from the TSE object and convert it to a data frame
tse_metadata <- as.data.frame(colData(tse)) %>% filter(!is.na(income_group))

# ---------------------------
# PCoA Plot
# --------------------------

p <- ggplot(tse_metadata, aes(x = PC1, y = PC2, color = gender)) +
  geom_point(size = 0.8, alpha = 1) +
  scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  stat_ellipse(aes(group = gender), type = "norm", level = 0.95, size = 0.7) +
  # FIXME: automate the explanatory percentages!
  labs(
    x = paste0("PC1 (14.95%)"),
    y = paste0("PC2 (10.34%)"),
    color = "Gender"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(),
    axis.text = element_text(),
    legend.title = element_text(),
    legend.text = element_text(),
    legend.position = "bottom",
    legend.direction="horizontal",
    axis.line = element_line(color = "black")
  ) +
  facet_grid(income_group ~ age_category) 
  

print(p)

# This generates publication quality printout:
library(Cairo)
CairoJPEG("RESULTS/FIGURES/Fig2_PCoA.jpg", width=800, height=450, quality=100)
print(p)
dev.off()



