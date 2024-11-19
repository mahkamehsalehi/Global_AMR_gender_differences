library(ggplot2)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(countrycode)
library(viridis)
library(tidyverse)
library(ggpubr)

# Data Loading and Preprocessing
tse <- readRDS("DATA/TSE.rds")
df <- as.data.frame(colData(tse)) %>%
  mutate(log_ARG_load = log(ARG_load))

# Recode and Factorize 'sex_combined'
df$sex_combined <- recode(df$sex_combined, "male" = "Men", "female" = "Women")
df$sex_combined <- factor(df$sex_combined, levels = c("Women", "Men"))

# Define custom plot theme
common_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "none",
    axis.line = element_line(color = "black"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
  )

# Define Plot p1: Host Age Distribution by Gender
p1 <- ggplot(df %>% filter(!is.na(sex_combined)), 
             aes(x = host_age_years, fill = sex_combined)) +
  scale_fill_manual(values = c("Women" = "#F8766D", "Men" = "#619CFF")) +
  geom_histogram(binwidth = 5, position = position_dodge(width = 4), color = "black", alpha = 0.9) +
  labs(
    x = "Age",
    y = "Count",
    fill = "Gender"
  ) +
  common_theme


# Define Plot p2: Antibiotic Resistance Load Distribution
p2 <- ggplot(df %>% filter(!is.na(sex_combined)), 
             aes(x = log_ARG_load)) +
  geom_histogram(binwidth = 0.2, color = "black", alpha = 0.7) +
  labs(
    x = "ARG Load (log RPKM)",
    y = "Count"
  ) +
  common_theme

# Define Plot p3: World Bank Income Group Distribution

# Rename the levels
levels(df$World_Bank_Income_Group) <- c("High", "Upper middle", "Lower middle", "Low")

# Reorder the levels
df$World_Bank_Income_Group<- factor(df$World_Bank_Income_Group, 
                               levels = c("Low", "Lower middle", "Upper middle", "High"), 
                               ordered = TRUE)


p3 <- ggplot(df %>% filter(!is.na(World_Bank_Income_Group)), 
             aes(x = World_Bank_Income_Group)) +
  geom_bar(color = "black", alpha = 0.7) +
  labs(
    x = "World Bank Income Group",
    y = "Count"
  ) +
  common_theme

# Define Plot p4: Antibiotic Usage Distribution
p4 <- ggplot(df %>% filter(!is.na(sex_combined)), 
             aes(x = Usage)) +
  geom_histogram(binwidth = 1, color = "black", alpha = 0.7) +
  labs(
    x = "Antibiotic Usage (DDD)",
    y = "Count"
  ) +
  common_theme

# Define Plot p5: Number of Samples per Country (World Map)
# Filter Data
df_filtered <- df %>% 
  filter(geo_loc_name_country_calc != "uncalculated", 
         !is.na(geo_loc_name_country_calc),
         !is.na(sex_combined))

# Load World Map Data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Summarize Data by Country and Gender
df_summary <- df_filtered %>%
  group_by(geo_loc_name_country_calc, sex_combined) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(iso_a3 = countrycode(geo_loc_name_country_calc, 
                              origin = "country.name", 
                              destination = "iso3c"))

# Merge with World Map Data
world_data <- world %>% 
  left_join(df_summary, by = c("iso_a3"))

# Define Custom Breaks for Fill Scale
breaks_seq_custom <- c(0, 1000, 2000, 3000)

# Create Plot p5: Sample Count World Map
p5 <- ggplot(data = world_data) +
  geom_sf(aes(fill = count), color = "white", size = 0.2) +
  scale_fill_distiller(
    palette = "Spectral",
    direction = 1,
    breaks = breaks_seq_custom,
    labels = breaks_seq_custom,
    na.value = "grey90",
    name = "Number of Samples"
  ) +
  common_theme +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, 
                              margin = margin(b = 10)),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0), 
    legend.box.spacing = unit(0, "cm"),
    legend.position = "right",
    legend.key.width = unit(0.3, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Combine the Plots Using Patchwork
combined_plot <- (
  (p1 | p2) /
    (p3 | p4) /
    p5
) +
  plot_layout(heights = c(5.5, 5, 7)) + 
  plot_annotation(tag_levels = 'a')

print(combined_plot)

ggsave("RESULTS/FIGURES/Data_Summary.png", combined_plot, width = 12, height = 8, dpi = 300)
