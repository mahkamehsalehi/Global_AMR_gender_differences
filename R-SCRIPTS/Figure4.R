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
  filter(!is.na(age_category) & !is.na(log_ARG_load))

df <- metadata_hic
df <- df %>% mutate(age_group = cut(host_age_years, c(0, 1, 3, 17, 35, 65, 80, 100), include.lowest=TRUE)) 
df$reg <- factor(df$geo_loc_name_country_continent_calc)

myplot <- function (df, region) {

  d <- subset(df, reg == region)

  p <- ggplot(d,
    aes(x = gender, y = log_ARG_load, fill = gender)) +
    geom_boxplot(
      position = position_dodge(width = 0.8),
      outlier.shape = NA,
      width = 0.6,
      alpha = 1,
      show.legend = TRUE
    ) +
  facet_wrap(~age_group, scales = "fixed", nrow = 1) +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  labs(x = "Age groups", y = "ARG load (log RPKM)", title=region) +
  theme_minimal() +
  stat_compare_means(
    comparisons = list(
      c("Women", "Men")
    ),
    aes(label = ..p.signif.., size=7),
    method = "wilcox.test",
    p.adjust.method = "BH",
    hide.ns = FALSE
  ) +
  theme_minimal(16) +
  theme(
    axis.line = element_line(color = "black"),
    # axis.text.x = element_text(size=12),
    #title.text = element_text(size=12),
    axis.text.x = element_blank(),        
    #axis.text.y = element_text(size=12),
    #axis.title.x = element_text(size=12),
    #axis.title.y = element_text(size=12),        
    strip.text = element_text(size = 10, face = "bold", angle = 45),
    panel.spacing = unit(0.5, "lines")
  ) 

}

regs <- c("Europe", "North America")
ps <- lapply(regs, function (region) {myplot(df, region)})
leg <- get_legend(ps[[1]] + theme(legend.direction="horizontal", legend.title=element_blank()))
p <- plot_grid(ps[[1]] + theme(legend.position="none"),
               ps[[2]] + theme(legend.position="none"),
	       leg,
	       rel_heights=c(0.45, 0.45, 0.1),
	       nrow=3
	       )
print(p)


jpeg("RESULTS/FIGURES/Fig4.jpg", width=500, height=700, quality=100)
print(p)
dev.off()


