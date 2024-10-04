library(ggplot2)
library(dplyr)
library(ggpubr)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(corrplot)
library(reshape2)

# Convert colData to a data frame
tse_gender <- tse[, !is.na(colData(tse)$sex_combined)]
df <- as.data.frame(colData(tse_gender))

#-------------------------------------------------------------------------------
#--------------------------- Box Plot by Income Group --------------------------
#-------------------------------------------------------------------------------

df_income <- df %>%
  drop_na(World_Bank_Income_Group)

box_plot_income <- ggplot(df_income, aes(x = World_Bank_Income_Group, y = log10_ARG_load, fill = sex_combined)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +         
  geom_jitter(width = 0.2, alpha = 0.5) +                   
  stat_compare_means(aes(group = sex_combined), 
                     method = "wilcox.test", 
                     label = "p.format") +                   
  labs(title = "log10 ARG Load by World Bank Income Group and Gender",
       x = "World Bank Income Group",
       y = "log10 ARG Load",
       fill = "Gender") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save the plot
ggsave("boxplot_ARG_load_income_gender.png", plot = box_plot_income, width = 10, height = 6)


#-------------------------------------------------------------------------------
#------------------------------ Bivariate ARG vs Income ------------------------
#-------------------------------------------------------------------------------

ggplot(df_income, aes(x = host_age_years, y = log10_ARG_load, color = World_Bank_Income_Group, shape = sex_combined)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "ARG Load vs. Age by Income Group and Gender",
       x = "Age (years)",
       y = "log10 ARG Load",
       color = "Income Group",
       shape = "Gender") +
  theme_minimal()

# Save the plot
ggsave("bivariate_plot_ARG_load_age_income_gender.png", width = 12, height = 8)


#-------------------------------------------------------------------------------
#----------------------------- faceted line_plot income ------------------------
#-------------------------------------------------------------------------------

ggplot(df_income, aes(x = host_age_years, y = log10_ARG_load, color = sex_combined)) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~ World_Bank_Income_Group) +
  labs(title = "Trend of ARG Load with Age by Gender and Income Group",
       x = "Age (years)",
       y = "log10 ARG Load",
       color = "Gender") +
  theme_minimal()

# Save the plot
ggsave("faceted_line_plot_ARG_load_age_gender_income.png", width = 12, height = 8)



#-------------------------------------------------------------------------------
#------------------------------- Box Plot by Age groups ------------------------
#-------------------------------------------------------------------------------

# Create age groups
df_filtered <- df %>%
  drop_na(host_age_years) %>%
  mutate(age_group = cut(host_age_years, breaks = c(0, 18, 35, 50, 65, Inf),
                         labels = c("0-18", "19-35", "36-50", "51-65", "66+"),
                         right = FALSE))

# Boxplot
box_plot_age_group <- ggplot(df_filtered, aes(x = age_group, y = ARG_load, fill = age_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_y_continuous(trans = 'log10') +
  labs(title = "ARG Load by Age Group",
       x = "Age Group",
       y = "ARG Load (log10 scale)") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("box_plot_age_group.png", plot = box_plot_age_group, width = 8, height = 6)


#-------------------------------------------------------------------------------
#------------------------------- Density Ridge of Age ---------------------------
#------------------------------------------------------------------------------

ggplot(df, aes(x = host_age_years)) +
  geom_histogram(binwidth = 1, fill = "lightgreen", color = "black") +
  labs(title = "Distribution of Age",
       x = "Age (years)",
       y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("density_ridge_plot_ARG_load_age_gender.png", width = 12, height = 8)


#-------------------------------------------------------------------------------
#------------------------------- Distribution of Age ---------------------------
#-------------------------------------------------------------------------------

library(ggridges)

ggplot(df_filtered, aes(x = log10_ARG_load, y = age_group, fill = sex_combined)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5) +
  labs(title = "Density Ridge Plot of log10 ARG Load by Age Group and Gender",
       x = "log10 ARG Load",
       y = "Age Group",
       fill = "Gender") +
  theme_minimal()


#-------------------------------------------------------------------------------
#------------------------------- Box Plot by Region ----------------------------
#-------------------------------------------------------------------------------

df_region <- df %>%
  drop_na(Region)
box_plot_region <- ggplot(df_region, aes(x = Region, y = log10_ARG_load, fill = sex_combined)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_compare_means(aes(group = sex_combined), method = "wilcox.test") +
  labs(title = "log10 ARG Load by Region and Gender",
       x = "Region",
       y = "log10 ARG Load",
       fill = "Gender") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("boxplot_ARG_load_region_gender.png", plot = box_plot_region, width = 12, height = 6)


#-------------------------------------------------------------------------------
#------------------------------- Box Plot by country ---------------------------
#-------------------------------------------------------------------------------

df_country <- df %>%
  drop_na(geo_loc_name_country_calc) %>%
  mutate(geo_loc_name_country_calc = fct_reorder(geo_loc_name_country_calc, log10_ARG_load, .fun = mean))

# Create the box plot
box_plot_country <- ggplot(df_country, aes(x = geo_loc_name_country_calc, y = log10_ARG_load, fill = sex_combined)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_compare_means(aes(group = sex_combined), method = "wilcox.test", label = "p.format") +
  labs(title = "log10 ARG Load by Country and Gender",
       x = "Country",
       y = "log10 ARG Load",
       fill = "Gender") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Save the plot
ggsave("box_plot_country.png", plot = box_plot_country, width = 20, height = 10, dpi = 300)

#-------------------------------------------------------------------------------
# Heatmap

heatmap_plot <- df_country %>%
  group_by(geo_loc_name_country_calc, sex_combined) %>%
  summarise(mean_log10_ARG_load = mean(log10_ARG_load, na.rm = TRUE)) %>%
  ggplot(aes(x = geo_loc_name_country_calc, y = sex_combined, fill = mean_log10_ARG_load)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "Mean log10 ARG Load by Country and Gender",
       x = "Country",
       y = "Gender",
       fill = "Mean log10 ARG Load") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#-------------------------------------------------------------------------------
#---------------------------------  Wilcoxon test ------------------------------
#-------------------------------------------------------------------------------

# Function to perform Wilcoxon test within groups
perform_wilcox <- function(data, group_var, target_var, compare_var) {
  data %>%
    group_by(!!sym(group_var)) %>%
    summarize(p_value = wilcox.test(as.formula(paste(target_var, "~", compare_var)), data = .)$p.value)
}

# Wilcoxon test within Income Groups
wilcox_income <- perform_wilcox(df_filtered, "World_Bank_Income_Group", "log10_ARG_load", "sex_combined")
write.csv(wilcox_income, "wilcox_ARG_load_income_gender.csv", row.names = FALSE)

# Wilcoxon test within Regions
wilcox_region <- perform_wilcox(df_filtered, "Region", "log10_ARG_load", "sex_combined")
write.csv(wilcox_region, "wilcox_ARG_load_region_gender.csv", row.names = FALSE)


# Apply Benjamini-Hochberg correction
wilcox_region <- wilcox_region %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH"))

write.csv(wilcox_region, "wilcox_ARG_load_region_gender_adjusted.csv", row.names = FALSE)


#-------------------------------------------------------------------------------
#---------------------------- Scatter Plot Analysis ----------------------------
#-------------------------------------------------------------------------------

df_scatter <- df %>%
  mutate(host_age_years = as.numeric(host_age_years)) %>%
  filter(!is.na(host_age_years))
# Scatter Plot with geom_smooth
scatter_plot <- ggplot(df_scatter, aes(x = host_age_years, y = log10_ARG_load, color = sex_combined)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(title = "log10 ARG Load vs. Age by Gender",
       x = "Age (years)",
       y = "log10 ARG Load",
       color = "Gender") +
  theme_minimal()

# Save the plot
ggsave("scatter_ARG_load_age_gender.png", plot = scatter_plot, width = 8, height = 6)


#-------------------------------------------------------------------------------
#----------------------- Density Plot of Income Group --------------------------
#-------------------------------------------------------------------------------

ggplot(df_income, aes(x = log10_ARG_load, fill = sex_combined)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ World_Bank_Income_Group) +
  labs(title = "Density Plot of log10 ARG Load by Income Group and Gender",
       x = "log10 ARG Load",
       y = "Density",
       fill = "Gender") +
  theme_minimal()

# Save the plot
ggsave("density_plot_ARG_load_income_gender.png", width = 12, height = 8)


#-------------------------------------------------------------------------------
#-------------------------- Heatmap of Correlations ----------------------------
#-------------------------------------------------------------------------------

# Select relevant numeric variables
corr_data <- df %>%
  select(log10_ARG_load, host_age_years)

# Compute correlation matrix
corr_matrix <- cor(corr_data, use = "complete.obs")

# Plot the correlation matrix
corrplot(corr_matrix, method = "color", addCoef.col = "black",
         tl.col = "black", number.cex = 0.8, title = "Correlation Matrix",
         mar = c(0,0,1,0))

# Save the plot
png("correlation_matrix.png", width = 600, height = 600)
corrplot(corr_matrix, method = "color", addCoef.col = "black",
         tl.col = "black", number.cex = 0.8, title = "Correlation Matrix",
         mar = c(0,0,1,0))
dev.off()


#-------------------------------------------------------------------------------
#--------------------------- Faceted Scatter Plots -----------------------------
#-------------------------------------------------------------------------------

ggplot(df_income, aes(x = host_age_years, y = log10_ARG_load, color = sex_combined)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~ World_Bank_Income_Group) +
  labs(title = "Scatter Plot of log10 ARG Load vs. Age by Income Group and Gender",
       x = "Age (years)",
       y = "log10 ARG Load",
       color = "Gender") +
  theme_minimal()

# Save the plot
ggsave("faceted_scatter_plot_ARG_load_age_income_gender.png", width = 12, height = 8)


#-------------------------------------------------------------------------------
#------------------------- Geographic Map Visualizations -----------------------
#-------------------------------------------------------------------------------

# Aggregate ARG Load by Country and Gender
country_summary <- df %>%
  filter(!is.na(geo_loc_name_country_calc)) %>%
  group_by(geo_loc_name_country_calc, sex_combined) %>%
  summarize(mean_ARG_load = mean(log10_ARG_load, na.rm = TRUE)) %>%
  ungroup()

# Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Merge with country_summary
world_data <- left_join(world, country_summary, by = c("name" = "geo_loc_name_country_calc"))

# Plot
ggplot(world_data) +
  geom_sf(aes(fill = mean_ARG_load), color = "gray70") +
  facet_wrap(~ sex_combined) +
  scale_fill_viridis_c(option = "plasma", na.value = "white") +
  labs(title = "Average ARG Load by Country and Gender",
       fill = "Mean ARG Load") +
  theme_minimal()

# Save the plot
ggsave("geographic_map_ARG_load_gender.png", width = 15, height = 8)


#-------------------------------------------------------------------------------
#--------------------- Boxplots with Facets for Country Subsets ----------------
#-------------------------------------------------------------------------------

# Select top 10 countries by sample size
top_countries <- df_filtered %>%
  group_by(geo_loc_name_country_calc) %>%
  summarize(count = n()) %>%
  arrange(desc(count)) %>%
  slice(1:10) %>%
  pull(geo_loc_name_country_calc)

# Filter data for top countries
df_top_countries <- df_filtered %>%
  filter(geo_loc_name_country_calc %in% top_countries)

# Create the faceted box plot
ggplot(df_top_countries, aes(x = age_group, y = log10_ARG_load, fill = sex_combined)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_y_continuous(trans = 'log10') +
  facet_wrap(~ geo_loc_name_country_calc) +
  labs(title = "ARG Load by Age Group and Gender in Top 10 Countries",
       x = "Age Group",
       y = "log10 ARG Load",
       fill = "Gender") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("faceted_boxplot_ARG_load_age_gender_top_countries.png", width = 20, height = 15)


#-------------------------------------------------------------------------------
#--------------------------- bivariate density plot ARG age --------------------
#-------------------------------------------------------------------------------
ggplot(df_filtered, aes(x = host_age_years, y = log10_ARG_load, color = sex_combined)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.4) +
  geom_point(alpha = 0.3, size = 0.5) +
  scale_fill_viridis_c() +
  labs(title = "Bivariate Density Plot of Age and ARG Load by Gender",
       x = "Age (years)",
       y = "log10 ARG Load",
       color = "Gender",
       fill = "Density") +
  theme_minimal()

# Save the plot
ggsave("bivariate_density_plot_ARG_load_age_gender.png", width = 10, height = 6)


#-------------------------------------------------------------------------------
#--------------------------- stacked density plot ARG-age ----------------------
#-------------------------------------------------------------------------------
ggplot(df_filtered, aes(x = log10_ARG_load, fill = sex_combined)) +
  geom_density(alpha = 0.5, position = "stack") +
  facet_wrap(~ age_group) +
  labs(title = "Stacked Density Plot of log10 ARG Load by Age Group and Gender",
       x = "log10 ARG Load",
       y = "Density",
       fill = "Gender") +
  theme_minimal()

# Save the plot
ggsave("stacked_density_plot_ARG_load_age_gender.png", width = 12, height = 8)


#-------------------------------------------------------------------------------
#--------------------------- Boxplot ARG-Income-Age Group  ---------------------
#-------------------------------------------------------------------------------

df_selected <- df %>%
  select("geo_loc_name_country_continent_calc", 
         "host_age_years", 
         "ARG_load", 
         "World_Bank_Income_Group",
         "Corruption_and_Governance_Index", 
         "Education_Index",
         "GDP_per_head", 
         "geo_loc_name_country_calc", 
         "sex_combined", 
         "log10_ARG_load", 
         "Region", 
         "Usage", 
         "Health_Spend_Index", 
         "Infrastructure_Index", 
         "usage_bayesian")

# Original Income Group Counts Before Filtering
original_income_counts <- table(df$World_Bank_Income_Group)
print(original_income_counts)

# Recode Income Groups
df_filtered <- df_selected %>%
  mutate(
    Income_Group_Modified = case_when(
      World_Bank_Income_Group %in% c("Low income", "Lower middle income") ~ "Others",
      TRUE ~ as.character(World_Bank_Income_Group)
    )
  )

# Convert to Factor with Specified Levels
df_filtered$Income_Group_Modified <- factor(df_filtered$Income_Group_Modified,
                                            levels = c("High income", "Upper middle income", "Others"))

# New Income Group Counts After Recoding
new_income_counts <- table(df_filtered$Income_Group_Modified)
print(new_income_counts)

# Reclassify Age
df_filtered <- df_filtered %>%
  drop_na(host_age_years) %>%
  mutate(
    Age_Group = case_when(
      host_age_years >= 0 & host_age_years <= 3 ~ "0-3",
      host_age_years >= 4 & host_age_years <= 10 ~ "4-10",
      host_age_years >= 11 & host_age_years <= 18 ~ "11-18",
      host_age_years > 18 ~ "Adults",
      TRUE ~ NA_character_ 
    )
  )

df_filtered <- df_filtered %>%
  drop_na(Age_Group)%>%
  drop_na(Income_Group_Modified)

# Convert Age_Group to Factor
df_filtered$Age_Group <- factor(df_filtered$Age_Group,
                                levels = c("0-3", "4-10", "11-18", "Adults"))

# Age Group Counts
age_group_counts <- table(df_filtered$Age_Group)
print(age_group_counts)

ggplot(df_filtered, aes(x = Income_Group_Modified, y = log10(ARG_load), fill = Age_Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_jitter(alpha = 0.5, position = position_dodge(width = 0.8)) +  
  stat_compare_means(aes(group = Age_Group), method = "wilcox.test", label = "p.signif") +
  labs(title = "log10 ARG Load by Income Group and Age Group",
       x = "Income Group",
       y = "log10 ARG Load",
       fill = "Age Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("boxplot_ARG_load_income_age_group.png", width = 10, height = 6)


#-------------------------------------------------------------------------------
#---------------------- Scatterplot ARG-Income-Age Group  ----------------------
#-------------------------------------------------------------------------------

ggplot(df_filtered, aes(x = host_age_years, y = ARG_load, color = Age_Group)) +
  geom_point(alpha = 0.6) +
  scale_y_continuous(trans = 'log10') +
  facet_wrap(~ Income_Group_Modified) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(title = "ARG Load vs. Age by Income Group and Age Category",
       x = "Age (years)",
       y = "ARG Load (log10 scale)",
       color = "Age Group") +
  theme_minimal()

# Save the plot
ggsave("scatter_ARG_load_age_income_group.png", width = 12, height = 8)



#------------------
# Data Preparation
#------------------

df_vis <- df %>%
  select(sex_combined, 
         host_age_years, 
         GDP_per_head, 
         log10_ARG_load, 
         Infrastructure_Index,
         Usage) %>%
  drop_na()

#-------------------------------------------------------------------------------
#------------------- Scatter Plot of ARG Load vs GDP per head  -----------------
#-------------------------------------------------------------------------------

scatter_plot_arg_gdp <- ggplot(df_vis, aes(x = GDP_per_head, y = log10_ARG_load, color = host_age_years)) +
  geom_point(alpha = 0.4, size = 2) +
  geom_smooth(aes(group = sex_combined, linetype = sex_combined), method = "loess", se = TRUE, color = "black") +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title = "log10 ARG Load vs. GDP per Head by Gender and Age",
    x = "GDP per Head (USD)",
    y = "log10 ARG Load",
    color = "Age (Years)",
    linetype = "Gender"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("scatter_plot_ARG_load_GDP_age_gender.png", plot = scatter_plot_arg_gdp, width = 10, height = 6)


#-------------------------------------------------------------------------------
#--------------- Scatter Plot of ARG Load vs Infrastructure Index  -------------
#-------------------------------------------------------------------------------

scatter_plot_arg_infra <- ggplot(df_vis, aes(x = Infrastructure_Index, y = log10_ARG_load, color = host_age_years)) +
  geom_point(alpha = 0.4, size = 2) +
  geom_smooth(aes(group = sex_combined, linetype = sex_combined), method = "loess", se = TRUE, color = "black") +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title = "log10 ARG Load vs. Infrastructure Index by Gender and Age",
    x = "Infrastructure Index",
    y = "log10 ARG Load",
    color = "Age (Years)",
    linetype = "Gender"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("scatter_plot_ARG_load_infra_age_gender.png", plot = scatter_plot_arg_infra, width = 10, height = 6)

#-------------------------------------------------------------------------------
#----------------------- Scatter Plot of ARG Load vs Usage ---------------------
#-------------------------------------------------------------------------------

scatter_plot_arg_usage <- ggplot(df_vis, aes(x = Usage, y = log10_ARG_load, color = host_age_years)) +
  geom_point(alpha = 0.4, size = 2) +
  geom_smooth(aes(group = sex_combined, linetype = sex_combined), method = "loess", se = TRUE, color = "black") +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title = "log10 ARG Load vs. Antibiotic Consumption by Gender and Age",
    x = "Antibiotic Consumption",
    y = "log10 ARG Load",
    color = "Age (Years)",
    linetype = "Gender"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("scatter_plot_ARG_load_usage_age.png", plot = scatter_plot_arg_usage, width = 10, height = 6)

#-------------------------------------------------------------------------------
#------------- Scatter Plot of ARG Load vs Usage Bayesian Estimate -------------
#-------------------------------------------------------------------------------

df_bayes <- df %>%
  select(log10_ARG_load, sex_combined, usage_bayesian, host_age_years) %>%
  drop_na()

scatter_plot_arg_use_bayes <- ggplot(df_bayes, aes(x = usage_bayesian, y = log10_ARG_load, color = host_age_years)) +
  geom_point(alpha = 0.4, size = 2) +
  geom_smooth(aes(group = sex_combined, linetype = sex_combined), method = "loess", se = TRUE, color = "black") +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title = "log10 ARG Load vs. Antibiotic Consumption Estimate (Bayesian) by Gender and Age",
    x = "Antibiotic Consumption Estimate (Bayesian)",
    y = "log10 ARG Load",
    color = "Age (Years)",
    linetype = "Gender"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("scatter_plot_ARG_load_usage_age_bayesian.png", plot = scatter_plot_arg_use_bayes, width = 10, height = 6)









