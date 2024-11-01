#-------------------------------------------------------------------------------
# Project: Women AMR Analysis
# Author: Mahkameh
# Date: 2024-10-29
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 1. Setup: Set Working Directory and Load Libraries
#-------------------------------------------------------------------------------

# Set the working directory
setwd("/scratch/project_2008149/USER_WORKSPACES/mahkameh/women_amr/")

# Load necessary libraries
library(tidyverse)     # Data manipulation and visualization
library(vegan)         # Ecological analysis
library(scater)        # Single-cell analysis (if applicable)
library(ggpubr)        # Publication-ready plots
library(rstatix)       # Statistical tests
library(viridis)       # Color scales
library(ggsignif)      # Significance annotations
library(patchwork)     # Combining multiple plots

#-------------------------------------------------------------------------------
# 2. Data Loading and Initial Filtering
#-------------------------------------------------------------------------------

# Load the TSE object
TSE <- readRDS("TSE.rds")

# Convert colData to a data frame for easier manipulation
metadata_df <- as.data.frame(colData(TSE))

# Select relevant columns for analysis
selected_columns <- c(
  "sex_combined", 
  "log10_ARG_load", 
  "geo_loc_name_country_calc", 
  "category", 
  "World_Bank_Income_Group",
  "GDP_per_head",
  "Usage"
)

plot_data <- metadata_df %>%
  select(all_of(selected_columns))

#-------------------------------------------------------------------------------
# 3. Define a Common Theme for Plots
#-------------------------------------------------------------------------------

# Define a consistent theme for all ggplots to ensure uniformity
common_theme <- theme_classic(base_size = 14) +
  theme(
    text = element_text(family = "Serif"),
    plot.title = element_text(
      face = "bold", 
      size = 14, 
      hjust = 0.5
    ),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.position = "none",
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 12, face = "bold")
  )

#-------------------------------------------------------------------------------
# 4. Data Preparation for Socioeconomic Plots
#-------------------------------------------------------------------------------

# Filter out unwanted categories and locations
socio_data <- plot_data %>%
  filter(
    category != "Infant Study", 
    geo_loc_name_country_calc != "Zimbabwe"
  )

#-------------------------------------------------------------------------------
# 5. Income Box Plot: log10 ARG Load by Gender and Income Group
#-------------------------------------------------------------------------------

# Prepare data for Income Box Plot
income_plot_data <- socio_data %>%
  mutate(
    Income_Group = case_when(
      World_Bank_Income_Group == "High income" ~ "HIC",
      World_Bank_Income_Group %in% c("Low income", "Lower middle income", "Upper middle income") ~ "LMIC",
      TRUE ~ NA_character_
    )
  ) %>%
  drop_na(sex_combined, log10_ARG_load, Income_Group)

# Create Income Box Plot
income_boxplot <- ggplot(income_plot_data, aes(x = sex_combined, y = log10_ARG_load)) +
  geom_jitter(aes(color = sex_combined), size = 0.5, width = 0.2, alpha = 0.6) +
  geom_boxplot(aes(fill = sex_combined), outlier.shape = NA, alpha = 0.7, color = "black") +
  facet_grid(cols = vars(Income_Group)) +
  stat_compare_means(
    method = "t.test",
    label = "p.signif",
    label.y = max(income_plot_data$log10_ARG_load, na.rm = TRUE) + 0.5
  ) +
  labs(
    title = "log10 ARG Load by Gender and Income Group",
    x = "Gender",
    y = "log10 ARG Load",
    fill = "Gender",
    color = "Gender"
  ) +
  common_theme

#-------------------------------------------------------------------------------
# 6. Antibiotic Usage Box Plot: log10 ARG Load by Gender and Usage Group
#-------------------------------------------------------------------------------

# Prepare data for Usage Box Plot
usage_plot_data <- socio_data %>%
  mutate(
    Usage_group = case_when(
      Usage > 11  ~ "Above 11",
      Usage <= 11 ~ "11 or Below",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    Usage_group = factor(Usage_group, levels = c("11 or Below", "Above 11"))
  ) %>%
  drop_na(sex_combined, log10_ARG_load, Usage_group)

# Create Usage Box Plot
usage_boxplot <- ggplot(usage_plot_data, aes(x = sex_combined, y = log10_ARG_load, fill = sex_combined)) +
  geom_jitter(aes(color = sex_combined), size = 0.5, alpha = 0.6, width = 0.1) +
  geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.7) +
  facet_grid(cols = vars(Usage_group)) +
  stat_compare_means(
    method = "t.test",
    label = "p.signif",
    label.y = max(usage_plot_data$log10_ARG_load, na.rm = TRUE) + 0.5
  ) +
  labs(
    title = "log10 ARG Load by Gender and Antibiotic Consumption",
    x = "Gender",
    y = "log10 ARG Load",
    fill = "Gender",
    color = "Gender"
  ) +
  common_theme

#-------------------------------------------------------------------------------
# 7. GDP Box Plot: log10 ARG Load by Gender and GDP Category
#-------------------------------------------------------------------------------

# Prepare data for GDP Box Plot
gdp_plot_data <- socio_data %>%
  filter(!is.na(GDP_per_head), !is.na(sex_combined), GDP_per_head <= 2) %>%
  mutate(
    GDP_category = case_when(
      GDP_per_head < 1 ~ "Low",
      GDP_per_head >= 1 & GDP_per_head < 1.5 ~ "Middle",
      GDP_per_head >= 1.5 ~ "High",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    GDP_category = factor(GDP_category, levels = c("Low", "Middle", "High"))
  ) %>%
  drop_na(log10_ARG_load, GDP_category)

# Create GDP Box Plot
gdp_boxplot <- ggplot(gdp_plot_data, aes(x = sex_combined, y = log10_ARG_load, fill = sex_combined)) +
  geom_jitter(aes(color = sex_combined), size = 0.5, alpha = 0.6, width = 0.1) +
  geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.5) +
  facet_grid(cols = vars(GDP_category)) +
  stat_compare_means(
    method = "t.test",
    label = "p.signif",
    label.y = max(gdp_plot_data$log10_ARG_load, na.rm = TRUE) + 0.5
  ) +
  labs(
    title = "log10 ARG Load by Gender and GDP",
    x = "Gender",
    y = "log10 ARG Load",
    fill = "Gender",
    color = "Gender"
  ) +
  common_theme

#-------------------------------------------------------------------------------
# 8. Combine Socioeconomic Plots into a Single Figure
#-------------------------------------------------------------------------------

# Combine the three socioeconomic plots horizontally
socioeconomic_figure <- (income_boxplot | usage_boxplot | gdp_boxplot) +
  plot_layout(ncol = 3) +
  plot_annotation(
    title = "Comparison of log10 ARG Load by Various Socioeconomic Factors",
    tag_levels = "a",
    theme = theme(
      text = element_text(family = "Serif"),
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5)
    )
  ) &
  theme(
    plot.margin = margin(10, 10, 10, 10),  # Margins around the entire figure
    panel.spacing = unit(1, "lines")        # Space between panels
  )

#-------------------------------------------------------------------------------
# 9. Display and Save the Combined Figure
#-------------------------------------------------------------------------------

# Display the combined socioeconomic figure
print(socioeconomic_figure)

# Save the combined figure as a high-resolution PNG file
ggsave(
  filename = "Figure_Socioeconomic_Analysis.png",
  plot = socioeconomic_figure,
  width = 18,            # Width in inches
  height = 6,            # Height in inches
  dpi = 300              # Resolution
)

#-------------------------------------------------------------------------------
# End of Script
#-------------------------------------------------------------------------------



# Prepare the Data for PCoA Based on Usage
df_pcoa_usage <- metadata_df %>%
  filter(
    sex_combined == "female",
    geo_loc_name_country_calc != "Zimbabwe",
    category != "Infant Study"
  ) %>%
  mutate(
    Usage_group = case_when(
      Usage > 11 ~ "Above 11",
      Usage <= 11 ~ "11 or Below"
    )
  ) %>%
  mutate(Usage_group = factor(Usage_group, levels = c("11 or Below", "Above 11"))) %>%
  filter(
    !is.na(Usage),
    !is.na(sex_combined),
    !is.na(Region)
  )

# Calculate Bray-Curtis Distance Matrix and PCoA for Usage Group
distance_matrix_usage_region <- vegdist(df_pcoa_usage_region$log10_ARG_load, method = "bray")
pcoa_result_usage_region <- cmdscale(distance_matrix_usage_region, eig = TRUE, k = 2)

# Calculate Percentage of Variance Explained
variance_explained_usage <- round(100 * pcoa_result_usage_region$eig / sum(pcoa_result_usage_region$eig), 1)
pcoa_df_usage_region <- as.data.frame(pcoa_result_usage_region$points)
colnames(pcoa_df_usage_region) <- c("PCoA1", "PCoA2")
pcoa_df_usage_region <- bind_cols(pcoa_df_usage_region, df_pcoa_usage_region %>% select(Usage_group, Region))

# Plot PCoA for Usage Group with Legend and Variance Explained
pcoa_usage_region_plot <- ggplot(pcoa_df_usage_region, aes(x = PCoA1, y = PCoA2, color = Usage_group)) +
  geom_point(alpha = 0.7, size = 3) +
  labs(
    title = "PCoA plot (Beta Diversity - Bray-Curtis) for Females by Region - Antibiotic Consumption",
    x = paste0("PCoA1 (", variance_explained_usage[1], "%)"),
    y = paste0("PCoA2 (", variance_explained_usage[2], "%)"),
    color = "Usage Group"
  ) +
  facet_wrap(~ Region) +
  common_theme +
  theme(legend.position = "right")

# Display the Plot for Usage Group by Region
print(pcoa_usage_region_plot)

#ggsave("PCoA_Usage_Group_by_Region.png", pcoa_usage_region_plot, width = 12, height = 8, units = "in", dpi = 300)

# -------------------------------------------------------------------

# Prepare the Data for PCoA Based on Income Group and Region
df_pcoa_income_region <- metadata_df %>%
  filter(
    sex_combined == "female",
    geo_loc_name_country_calc != "Zimbabwe",
    category != "Infant Study"
  ) %>%
  select(-Region) %>%
  mutate(
    Income_Group = case_when(
      World_Bank_Income_Group == "High income" ~ "HIC",
      World_Bank_Income_Group %in% c("Low income", "Lower middle income", "Upper middle income") ~ "LMIC"
    ),
    Income_Group = factor(Income_Group, levels = c("LMIC", "HIC"))
  ) %>%
  dplyr::rename(Region = geo_loc_name_country_continent_calc) %>%
  filter(
    !is.na(World_Bank_Income_Group),
    !is.na(sex_combined),
    !is.na(Region)
  )

# Calculate Bray-Curtis Distance Matrix and PCoA for Income Group
distance_matrix_income_region <- vegdist(df_pcoa_income_region$log10_ARG_load, method = "bray")
pcoa_result_income_region <- cmdscale(distance_matrix_income_region, eig = TRUE, k = 2)

# Calculate Percentage of Variance Explained
variance_explained_income <- round(100 * pcoa_result_income_region$eig / sum(pcoa_result_income_region$eig), 1)
pcoa_df_income_region <- as.data.frame(pcoa_result_income_region$points)
colnames(pcoa_df_income_region) <- c("PCoA1", "PCoA2")
pcoa_df_income_region <- bind_cols(pcoa_df_income_region, df_pcoa_income_region %>% select(Income_Group, Region))

# Plot PCoA for Income Group with Legend and Variance Explained
pcoa_income_region_plot <- ggplot(pcoa_df_income_region, aes(x = PCoA1, y = PCoA2, color = Income_Group)) +
  geom_point(alpha = 0.7, size = 3) +
  labs(
    title = "PCoA plot (Beta Diversity - Bray-Curtis) for Females by Region - Income Level",
    x = paste0("PCoA1 (", variance_explained_income[1], "%)"),
    y = paste0("PCoA2 (", variance_explained_income[2], "%)"),
    color = "Income Group"
  ) +
  facet_wrap(~ Region) +
  common_theme +
  theme(legend.position = "right")

# Display the Plot for Income Group by Region
print(pcoa_income_region_plot)

g#gsave("PCoA_Income_Group_by_Region.png", pcoa_income_region_plot, width = 12, height = 8, units = "in", dpi = 300)

