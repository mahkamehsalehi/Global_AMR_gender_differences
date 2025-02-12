library(tidyverse)
library(SummarizedExperiment)
library(ggthemes)

tse <- readRDS("DATA/TSE.rds")

counts <- assay(tse)

gene_classes <- rowData(tse)$Class

meta <- as.data.frame(colData(tse))
meta <- meta %>%
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


for (sex in c("Men", "Women")) {
  for (income in c("HIC", "LMIC")) {
    idx <- which(meta$sex_combined == sex & meta$Income_Group == income)
    
    if (length(idx) == 0) {
      cat("No samples for", sex, "in", income, "\n")
      next
    }
    
    group_counts <- rowSums(counts[, idx, drop = FALSE])
    
    class_totals <- tapply(group_counts, gene_classes, sum, na.rm = TRUE)
    
    top5 <- sort(class_totals, decreasing = TRUE)[1:5]
    
    cat("Top 5 gene classes for", sex, "in", income, ":\n")
    print(top5)
    cat("\n")
  }
}

#-------------------------------------------------------------------------------

results <- list()

for (sex in c("Men", "Women")) {
  for (income in c("HIC", "LMIC")) {
    idx <- which(meta$sex_combined == sex & meta$Income_Group == income)
    
    if (length(idx) == 0) {
      cat("No samples for", sex, "in", income, "\n")
      next
    }
    
    group_counts <- rowSums(counts[, idx, drop = FALSE])
    
    class_totals <- tapply(group_counts, gene_classes, sum, na.rm = TRUE)
    
    top5 <- sort(class_totals, decreasing = TRUE)[1:5]
    
    df_temp <- data.frame(
      sex = sex,
      income = income,
      gene_class = names(top5),
      abundance = as.numeric(top5)
    )
    
    results[[paste(sex, income, sep = "_")]] <- df_temp
  }
}

df_bar <- do.call(rbind, results)

df_bar$sex <- factor(df_bar$sex, levels = c("Women", "Men"))

#-------------------------------------------------------------------------------

class_plot <- ggplot(df_bar, aes(x = reorder(gene_class, abundance), y = abundance, fill = gene_class)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(sex), cols = vars(income), scales = "free_y") +
  coord_flip() +
  labs(
    x = "Gene class",
    y = "Total abundance",
  ) +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
    strip.text = element_text(face = "plain", size = 14),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

ggsave(
  filename = "Top5_gene_classes_by_gender_and_income.png", 
  plot = last_plot(),
  width = 10, height = 5, dpi = 600
)

