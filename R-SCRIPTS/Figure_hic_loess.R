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


theme_set(theme_bw(20))
p <- ggplot(metadata_hic, aes(x=host_age_years, y=ARG_load,
            color=gender, fill=gender)) +
       # geom_point(size=0.5, alpha=0.5) +
       geom_smooth() +
       scale_color_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
       scale_y_log10(breaks=seq(400, 700, 100), labels=seq(400, 700, 100)) +
       labs(x="Age (y)", y="ARG load (RPKM)", color="", fill="", title="HIC") +
       #coord_cartesian(ylim=c(400, 1400)) +
       coord_cartesian(ylim=c(400, 700)) +       
       theme(legend.position=c(0.12, 0.88))
       

print(p)

library(Cairo)
CairoJPEG("../RESULTS/FIGURES/Fig_HIC_age.jpg", width=600, height=480, quality=100)
print(p)
dev.off()