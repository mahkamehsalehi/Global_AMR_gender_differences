library(mia)
library(tidyverse)

tse <- readRDS("DATA/TSE_filtered.rds")
df <- as.data.frame(colData(tse))


############################
# ARG Load
############################

results_load <- df %>%
  group_by(age_category_new) %>%
  summarise(
    N1 = sum(gender == "Women"),
    N2 = sum(gender == "Men"),
    Group1_Median = median(shannon_diversity[gender == "Women"], na.rm = TRUE),
    Group2_Median = median(shannon_diversity[gender == "Men"], na.rm = TRUE),
    P_value = wilcox.test(
      shannon_diversity[gender == "Women"],
      shannon_diversity[gender == "Men"],
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    P.adj = p.adjust(P_value, method = "BH"),
    P.adj_Significance = case_when(
      P.adj > 0.05 ~ "ns",
      P.adj <= 0.05 & P.adj > 0.01 ~ "*",
      P.adj <= 0.01 & P.adj > 0.001 ~ "**",
      P.adj <= 0.001 & P.adj > 0.0001 ~ "***",
      P.adj <= 0.0001 ~ "****"
    ),
    Effect_Direction = case_when(
      Group1_Median > Group2_Median ~ "Women > Men",
      Group1_Median < Group2_Median ~ "Men > Women",
      TRUE ~ "No Difference"
    )
  )

# Print or save the results
print(results_load, n = Inf, width = Inf)

############################
# ARG Diversity
############################

results_diversity <- df %>%
  group_by(age_category_new) %>%
  summarise(
    N1 = sum(gender == "Women"),
    N2 = sum(gender == "Men"),
    Group1_Median = median(shannon_diversity[gender == "Women"], na.rm = TRUE),
    Group2_Median = median(shannon_diversity[gender == "Men"], na.rm = TRUE),
    P_value = wilcox.test(
      shannon_diversity[gender == "Women"],
      shannon_diversity[gender == "Men"],
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    P.adj = p.adjust(P_value, method = "BH"),
    P.adj_Significance = case_when(
      P.adj > 0.05 ~ "ns",
      P.adj <= 0.05 & P.adj > 0.01 ~ "*",
      P.adj <= 0.01 & P.adj > 0.001 ~ "**",
      P.adj <= 0.001 & P.adj > 0.0001 ~ "***",
      P.adj <= 0.0001 ~ "****"
    ),
    Effect_Direction = case_when(
      Group1_Median > Group2_Median ~ "Women > Men",
      Group1_Median < Group2_Median ~ "Men > Women",
      TRUE ~ "No Difference"
    )
  )

print(results_diversity, n = Inf, width = Inf)

