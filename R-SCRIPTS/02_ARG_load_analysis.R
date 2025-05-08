# Load necessary libraries
library(jtools)
library(tidyverse)


# Load data from RDS file
TSE <- read_rds("../DATA/TSE.rds")

# Convert colData metadata to a data frame
df <- colData(TSE) %>% as.data.frame()

# Subset relevant columns for plotting
plot_df <- df %>%
  dplyr::select(matches(c("sex_combined", "log10_ARG_load", "host_age_years", "category", "World_Bank_Income_Group")))


# --- PLOT 1: Age vs. ARG load (linear fit) faceted by income group ---
ggplot(na.omit(plot_df), aes(x = host_age_years, y = log10_ARG_load, color = sex_combined)) +
  geom_jitter(size = 0.1, width = 0.1) +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = c(0,1,3,5,10,15,20,25,30,40,50,60,70,80))+
  theme_classic() +
  facet_grid(cols=vars(World_Bank_Income_Group))


# --- PLOT 2: Age vs. ARG load (Loess fit) ---
ggplot(na.omit(plot_df), aes(x = host_age_years, y = log10_ARG_load, color = sex_combined)) +
  geom_jitter(size = 0.1, width = 0.1) +
  geom_smooth(method = "loess") +
  scale_x_continuous(breaks = c(0,1,3,5,10,15,20,25,30,40,50,60,70,80))+
  theme_classic()


# --- PLOT 3: ARG load vs. antibiotic use by gender (excluding certain data) ---
df <- colData(TSE_gender) %>% as.data.frame()

# Filter out "Infant Study" category and Zimbabwe samples
tmp <- df[df$category!="Infant Study",]
tmp <- tmp[tmp$geo_loc_name_country_calc!="Zimbabwe",]
ggplot(tmp, aes(x = Usage, y = log10_ARG_load, color = sex_combined)) +
  geom_jitter(size = 0.1, width = 0.1) +
  geom_smooth(method = "loess") +
  theme_classic() 


# --- ANOVA: ARG load by study category ---
fit <- aov(log10_ARG_load ~ category, data = df) 


# ---GLM: ARG diversity by gender ---
tmp <- df %>%
  filter(category != "Infant Study",
         geo_loc_name_country_calc != "Zimbabwe")


# Fit GLM and summarize (exponentiate coefficients)
glm(ARG_div ~ sex_combined, data = df) %>% 
  summ(exp = TRUE, digits = 3)
