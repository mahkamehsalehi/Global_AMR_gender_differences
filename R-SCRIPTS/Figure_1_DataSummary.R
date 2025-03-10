library(ggplot2)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(countrycode)
library(viridis)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(mia)
library(miaViz)


# Data Loading and Preprocessing
tse <- readRDS("DATA/TSE.rds")
df <- as.data.frame(colData(tse)) %>%
  mutate(log_ARG_load = log(ARG_load))

# Recode and Factorize 'sex_combined'
df$sex_combined <- recode(df$sex_combined, "male" = "Men", "female" = "Women")
df$sex_combined <- factor(df$sex_combined, levels = c("Women", "Men"))

# Define custom plot theme
s <- 14 # scale for the figure size definitions
common_theme <- theme_classic(base_size = s) +
  theme(
    plot.title = element_text(face = "bold", size = s, hjust = 0.5),
    axis.text = element_text(size = s),
    axis.title = element_text(size = s),
    legend.position = "none",
    axis.line = element_line(color = "black"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = s, face = "bold"),
  )

# Define Plot p1: Host Age Distribution by Gender
p1 <- ggplot(df %>% filter(!is.na(sex_combined)), 
             aes(x = host_age_years, fill = sex_combined)) +
  scale_fill_manual(values = c("Women" = "#f03b20", "Men" = "#3182bd")) +
  geom_histogram(binwidth = 5, position = position_dodge(width = 5), color = "black", alpha = 0.7) +
  labs(
    x = "Age (y)",
    y = "Count (N)",
    fill = "Gender"
  ) +
  common_theme


# Define Plot p2: Antibiotic Resistance Load Distribution
library(scales)
p2 <- ggplot(df %>% filter(!is.na(sex_combined)), 
             aes(x = ARG_load, fill = sex_combined)) +
  scale_fill_manual(values = c("Women" = "#f03b20", "Men" = "#3182bd")) +
  geom_bar(position = position_dodge(), color = "black", alpha = 0.7, stat = "bin", binwidth = 0.2) +
  labs(
    x = "ARG load (RPKM)",
    y = "Count (N)",
    fill = "Gender"
  ) +
  scale_x_continuous(transf="log10",
                     breaks=10^(2:5),
                     labels=trans_format("log10", math_format(10^.x))) +
  common_theme

# Define Plot p3: World Bank Income Group Distribution

# Rename the levels
levels(df$World_Bank_Income_Group) <- c("High", "Upper middle", "Lower middle", "Low")

# Reorder the levels
df$World_Bank_Income_Group<- factor(df$World_Bank_Income_Group, 
                                    levels = c("Low", "Lower middle", "Upper middle", "High"), 
                                    ordered = TRUE)

p3 <- ggplot(df %>% filter(!is.na(World_Bank_Income_Group) & !is.na(sex_combined)), 
             aes(x = World_Bank_Income_Group, fill = sex_combined)) +
  geom_bar(position = position_dodge(), color = "black", alpha = 0.7) +
  scale_fill_manual(values = c("Women" = "#f03b20", "Men" = "#3182bd")) +
  labs(
    x = "Income Group (World Bank)",
    y = "Count (N)",
    fill = "Gender"
  ) +
  scale_x_discrete(labels = c("Low", "Lower\nmiddle", "Upper\nmiddle", "High")) + 
  common_theme 



# Define Plot p4: Antibiotic Usage Distribution
p4 <- ggplot(df %>% filter(!is.na(sex_combined)), 
             aes(x = Usage, fill = sex_combined)) +
  geom_histogram(binwidth = 1, position = position_dodge(), color = "black", alpha = 0.7) +
  scale_fill_manual(values = c("Women" = "#f03b20", "Men" = "#3182bd")) +
  labs(
    x = "Antibiotic Usage (DDD)",
    y = "Count (N)",
    fill = "Gender"
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
world_data <- world %>% dplyr::left_join(df_summary, by = c("iso_a3"))

# Define Custom Breaks for Fill Scale
breaks_seq_custom <- c(0, 1000, 2000, 3000)

# Create Plot p5: Sample Count World Map
p5 <- 
  ggplot(data = world_data) +
  geom_sf(aes(fill = count), color = "white", size = 0.2) +
  scale_fill_gradient(
    low = "lightblue",
    high = "darkblue",
    na.value = "grey90",
    name = "Samples (N)"
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


# Create the combined plot
combined_plot <- plot_grid(
  plot_grid(p1, p2, p3, p4, ncol = 2, labels=c("a)", "b)", "c)", "d)"), label_size = s),
  p5 + annotate("text", x=-180, y=100, label="e)", size=5) + labs(x="", y=""), 
  ncol = 1,
  rel_heights = c(6, 6)
)

library(Cairo)

CairoJPEG("RESULTS/FIGURES/Figure 1.jpg", 
          width = 3000,
          height = 3600,
          units = "px",
          dpi = 300)
print(combined_plot)
dev.off()