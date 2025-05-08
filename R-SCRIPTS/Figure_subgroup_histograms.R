library(tidyverse)
library(cowplot)

source("../R-SCRIPTS/prepare_data_for_lm_analysis.R")

my_data <- adult_metadata %>% 
  as.data.frame() %>% 
  select(
    age_category_new, 
    income_group, 
    region, 
    gender)


N_summary <- my_data %>% 
  group_by_all() %>% 
  summarize(N = n()) %>% 
  ungroup()

colnames(N_summary) <- c("Age", "Income Group", "Region", "Sex", "N")

p_N_summary <- N_summary %>%
  ggplot() + 
  geom_tile(aes(x = Region, y = Age, fill = N), 
            color = "black") + 
  geom_text(aes(x = Region, y = Age, label=N)) +
  facet_wrap(`Income Group` ~ Sex, 
             labeller = label_both) +
  scale_fill_gradient(low = "white", high = "red", na.value = "white") +
  # coord_equal() + 
  theme_classic(15) +
  scale_x_discrete(guide = guide_axis(angle = 90))
                
                

png("../RESULTS/FIGURES/SFig_N_summary.png",
    units = "in",
    res = 500,
    height = 8,
    width = 10)
print(p_N_summary)
dev.off()

                