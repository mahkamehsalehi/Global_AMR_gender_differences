library(ggpubr)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(countrycode)
library(viridis)
library(tidyverse)

tse <- readRDS("TSE.rds")
df <- as.data.frame(colData(tse))

# Define a common theme for all plots
common_theme <- theme_minimal(base_family = "Times") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Panel A: Age Distribution
p1 <- ggplot(df %>% filter(!is.na(sex_combined)), aes(x = host_age_years)) +
  geom_histogram(binwidth = 5, fill = "#008080", color = "black", alpha = 0.7) +
  labs(title = "Age Distribution", x = "Age (years)", y = "Count") +
  common_theme +
  theme(legend.position = "none")

# Panel B: ARG Load Distribution
p2 <- ggplot(df %>% filter(!is.na(sex_combined)), aes(x = log10_ARG_load)) +
  geom_histogram(binwidth = 0.2, fill = "#FF6347", color = "black", alpha = 0.7) +
  labs(title = "ARG Load Distribution", x = "log10(ARG Load)", y = "Count") +
  common_theme +
  theme(legend.position = "none")

# Panel C: GDP per Head Distribution
p3 <- ggplot(df %>% filter(!is.na(sex_combined)), aes(x = GDP_per_head)) +
  geom_histogram(binwidth = 0.5, fill = "#4682B4", color = "black", alpha = 0.7) +
  labs(title = "GDP per Head Distribution", x = "GDP per Head (USD)", y = "Count") +
  common_theme +
  theme(legend.position = "none")

# Panel D: Usage Distribution
p4 <- ggplot(df %>% filter(!is.na(sex_combined)), aes(x = Usage)) +
  geom_histogram(binwidth = 1, fill = "#9ACD32", color = "black", alpha = 0.7) +
  labs(title = "Antibiotic Usage Distribution", x = "Usage", y = "Count") +
  common_theme +
  theme(legend.position = "none")

# Panel E: Gender Distribution 
df_selected <- df %>% 
  filter(!is.na(geo_loc_name_country_calc), !is.na(sex_combined)) %>%
  rename(country = geo_loc_name_country_calc, sex = sex_combined)

gender_counts <- df_selected %>%
  dplyr::count(sex) %>%
  mutate(
    percentage = n / sum(n) * 100,
    label = paste0(n, " samples\n", round(percentage, 1), "%")
  )
total_samples <- sum(gender_counts$n)

p5 <- ggplot(gender_counts, aes(x = "", y = percentage, fill = sex)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3,
    family = "Times",
    fontface = "bold",
    color = "black"
  ) +
  labs(
    title = "Gender Distribution",
    subtitle = paste("Total samples:", total_samples)
  ) +
  theme_void(base_family = "Times") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 5)),
    plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 10)),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )

# Panel F: Number of Samples per Country
world <- ne_countries(scale = "medium", returnclass = "sf")
df_summary <- df_selected %>%
  group_by(country, sex) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(iso_a3 = countrycode(country, origin = "country.name", destination = "iso3c"))

world_data <- world %>% left_join(df_summary, by = c("iso_a3"))

# Define a set of breaks
max_count <- max(df_summary$count, na.rm = TRUE)
breaks_seq <- seq(0, max_count, by = 300)

# Create the map plot with fewer legend labels
p6 <- ggplot(data = world_data) +
  geom_sf(aes(fill = count), color = "white", size = 0.2) +
  scale_fill_viridis_c(
    option = "magma",
    breaks = breaks_seq,
    labels = breaks_seq,
    na.value = "grey90",
    name = "Sample Count",
    direction = -1
  ) +
  labs(title = "Number of Samples by Country") +
  theme_minimal(base_family = "Times") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
    legend.position = "right",
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )


# Arrange the panels
fig1 <- ggarrange(
  p1, p2, p3, p4, p5, p6,
  ncol = 2, nrow = 3,
  labels = c("a", "b", "c", "d", "e", "f"),
  font.label = list(size = 16, family = "Times", face = "bold"),
  heights = c(1, 1, 1.5),
  common.legend = FALSE
)

# Display the figure
print(fig1)
