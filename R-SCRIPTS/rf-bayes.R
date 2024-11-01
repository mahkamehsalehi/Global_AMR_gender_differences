# Set working directory and load necessary libraries
setwd("/scratch/project_2008149/USER_WORKSPACES/mahkameh/women_amr/")

# Load necessary libraries
library(tidyverse)
library(vegan)
library(scater)
library(ggpubr)
library(rstatix)
library(viridis)
library(ggsignif)
library(Matrix)
library(RColorBrewer)
library(caret)
library(randomForest)
library(brms)
library(xgboost)
library(glmnet)
library(car)

#-------------------------------------------------------------------------------
# 1. Data Loading and Initial Filtering
#-------------------------------------------------------------------------------

# Load TSE object
TSE <- readRDS("TSE.rds")

# Filter samples with complete sex and age data
valid_samples <- !is.na(colData(TSE)$sex_combined) & !is.na(colData(TSE)$host_age_years)
TSE_filtered <- TSE[, valid_samples]

# Further filter for adult samples with calculated geographic information
adult_samples <- colData(TSE_filtered)$host_age_years >= 18 &
  colData(TSE_filtered)$geo_loc_name_country_continent_calc != "uncalculated"
TSE_adult <- TSE_filtered[, adult_samples]

# Extract metadata and set factors for consistent ordering in plots
adult_metadata <- as.data.frame(colData(TSE_adult)) %>%
  mutate(
    sex_combined = recode(sex_combined, "male" = "Men", "female" = "Women"),
    sex_combined = factor(sex_combined, levels = c("Men", "Women"))
  )

#-------------------------------------------------------------------------------
# 2. Bray-Curtis PCoA Analysis
#-------------------------------------------------------------------------------

# Convert assay data to sparse matrix
assays(TSE_adult)$counts <- as(assays(TSE_adult)$counts, "sparseMatrix")

# Run PCoA
TSE_adult <- runMDS(
  TSE_adult,
  FUN = vegan::vegdist,
  method = "bray",
  name = "PCoA",
  ncomponents = 3,
  exprs_values = "relabundance"
)

# Extract PCoA results and variance explained
pcoa_data <- as.data.frame(reducedDim(TSE_adult, "PCoA"))
colnames(pcoa_data) <- paste0("PC", 1:ncol(pcoa_data))
eigenvalues <- attr(reducedDim(TSE_adult, "PCoA"), "eig")
variance_explained <- (eigenvalues[eigenvalues > 0] / sum(eigenvalues)) * 100

# Add metadata to PCoA data
pcoa_data <- pcoa_data %>%
  mutate(
    geo = adult_metadata$geo_loc_name_country_continent_calc,
    host_age_years = adult_metadata$host_age_years,
    sex_combined = adult_metadata$sex_combined,
    adult_age_category = case_when(
      host_age_years >= 18 & host_age_years < 35 ~ "Young Adult",
      host_age_years >= 35 & host_age_years < 65 ~ "Middle Adulthood",
      host_age_years >= 65 & host_age_years <= 100 ~ "Older Adult",
      TRUE ~ NA_character_
    )
  )

# Define common plot theme
common_theme <- theme_classic(base_size = 14) +
  theme(
    text = element_text(family = "Sans"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.position = "none",
    axis.line = element_line(color = "black"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 12, face = "bold")
  )

#-------------------------------------------------------------------------------
# 3. Plot PCoA Results
#-------------------------------------------------------------------------------

# Plot PCoA by region
pcoa_region_plot <- ggplot(pcoa_data, aes(x = PC1, y = PC2, color = geo)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(
    title = "PCoA Plot (Beta Diversity - Bray-Curtis) for Adults by Region",
    x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 2), "%)"),
    color = "Region"
  ) +
  common_theme
ggsave("final/pcoa_region.png", pcoa_region_plot, width = 10, height = 8)

# Plot PCoA by gender across regions
pcoa_gender_region_plot <- ggplot(pcoa_data, aes(x = PC1, y = PC2, color = sex_combined)) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(~ geo, nrow = 2) +
  labs(
    title = "PCoA Plot (Beta Diversity - Bray-Curtis) for Adults by Gender",
    x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 2), "%)"),
    color = "Gender"
  ) +
  common_theme
ggsave("final/pcoa_gender_region.png", pcoa_gender_region_plot, width = 10, height = 8)

# Plot PCoA by age category across regions
pcoa_age_category_plot <- ggplot(pcoa_data, aes(x = PC1, y = PC2, color = adult_age_category)) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(~ geo, nrow = 2) +
  labs(
    title = "PCoA Plot (Beta Diversity - Bray-Curtis) by Age Category",
    x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 2), "%)"),
    color = "Age Category"
  ) +
  common_theme +
  theme(legend.position = "right")
ggsave("final/pcoa_age_cat.png", pcoa_age_category_plot, width = 10, height = 8)

#-------------------------------------------------------------------------------
# 4. Beta Diversity Analysis
#-------------------------------------------------------------------------------

# Prepare metadata
adult_metadata <- adult_metadata %>%
  mutate(
    age_category = case_when(
      host_age_years >= 18 & host_age_years < 35 ~ "Young Adult",
      host_age_years >= 35 & host_age_years < 65 ~ "Middle Adulthood",
      host_age_years >= 65 & host_age_years <= 100 ~ "Older Adult",
      TRUE ~ NA_character_
    ),
    region = factor(geo_loc_name_country_continent_calc),
    GDP_per_head = as.numeric(GDP_per_head)
  ) %>%
  filter(!is.na(Usage), !is.na(GDP_per_head))

# Filter TSE_adult to match metadata
valid_samples <- !is.na(colData(TSE_adult)$GDP_per_head)
TSE_adult_filtered <- TSE_adult[, valid_samples]
colData(TSE_adult_filtered) <- DataFrame(adult_metadata)

# Extract and clean assay data (relative abundances)
assay_data <- assay(TSE_adult_filtered, "relabundance")

# Remove samples (columns) with all-zero counts
non_empty_samples <- colSums(assay_data > 0, na.rm = TRUE) > 0
assay_data <- assay_data[, non_empty_samples]
adult_metadata <- colData(TSE_adult_filtered)[non_empty_samples, ]

# Calculate the Bray-Curtis distance matrix
distance_matrix <- vegan::vegdist(t(assay_data), method = "bray")

# Calculate Shannon diversity
adult_metadata$shannon_diversity <- diversity(t(assay_data), index = "shannon")

# Statistical test
diversity_test <- t.test(shannon_diversity ~ sex_combined, data = adult_metadata)
print(diversity_test)

# PERMANOVA Analysis
set.seed(123)
permanova_result <- adonis2(
  distance_matrix ~ sex_combined + age_category + region + GDP_per_head + Usage,
  data = adult_metadata,
  permutations = 999
)
print(permanova_result)

# Test for homogeneity of dispersion
dispersion_sex <- betadisper(distance_matrix, adult_metadata$sex_combined)
dispersion_test <- anova(dispersion_sex)
print(dispersion_test)

#-------------------------------------------------------------------------------
# 5. Statistical Modeling
#-------------------------------------------------------------------------------

# GLM for Shannon Diversity
glm_diversity <- glm(
  shannon_diversity ~ sex_combined + age_category + region + GDP_per_head + Usage,
  data = adult_metadata,
  family = gaussian()
)
summary(glm_diversity)

# If log10_ARG_load is available, GLM for it
if ("log10_ARG_load" %in% names(adult_metadata)) {
  glm_arg_load <- glm(
    log10_ARG_load ~ sex_combined + age_category + region + GDP_per_head + Usage,
    data = adult_metadata,
    family = gaussian()
  )
  summary(glm_arg_load)
}

#-------------------------------------------------------------------------------
# 6. Additional Analyses
#-------------------------------------------------------------------------------

# Variance Inflation Factor (VIF) to check multicollinearity
vif_model <- lm(shannon_diversity ~ sex_combined + age_category + region + GDP_per_head + Usage, data = adult_metadata)
vif_values <- vif(vif_model)
print(vif_values)

# Random Forest Regression for Shannon Diversity
set.seed(123)
external_predictors <- adult_metadata[, c("sex_combined", "age_category", "region", "GDP_per_head", "Usage")]
data_rf_shannon <- data.frame(
  shannon_diversity = adult_metadata$shannon_diversity,
  external_predictors
)

# Split the data into training and testing sets
train_index <- createDataPartition(data_rf_shannon$shannon_diversity, p = 0.8, list = FALSE)
train_data <- data_rf_shannon[train_index, ]
test_data <- data_rf_shannon[-train_index, ]

# Train the Random Forest Regression model
rf_model <- randomForest(
  shannon_diversity ~ .,
  data = train_data,
  importance = TRUE,
  ntree = 500,
  mtry = floor(sqrt(ncol(train_data) - 1))
)

# Predict on the testing set
predictions <- predict(rf_model, newdata = test_data)

# Calculate performance metrics
actual <- test_data$shannon_diversity
r_squared <- cor(actual, predictions)^2
rmse <- sqrt(mean((actual - predictions)^2))

print(paste("Shannon Diversity - R-squared:", round(r_squared, 4)))
print(paste("Shannon Diversity - RMSE:", round(rmse, 4)))

# Plot Actual vs. Predicted
ggplot(data = NULL, aes(x = actual, y = predictions)) +
  geom_point(alpha = 0.6, color = "blue") +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  theme_minimal() +
  labs(
    title = "Actual vs. Predicted Shannon Diversity",
    x = "Actual Shannon Diversity",
    y = "Predicted Shannon Diversity"
  )

#-------------------------------------------------------------------------------
# 7. Bayesian Modeling with brms
#-------------------------------------------------------------------------------

# Bayesian linear regression model for Shannon Diversity
set.seed(123)
model_shannon_bayes <- brm(
  formula = shannon_diversity ~ sex_combined + age_category + GDP_per_head + Usage + (1 | region),
  data = train_data,
  family = gaussian(),
  prior = c(
    set_prior("normal(0, 5)", class = "b"),
    set_prior("normal(0, 5)", class = "Intercept"),
    set_prior("cauchy(0, 2)", class = "sd")
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = parallel::detectCores(),
  seed = 123
)
summary(model_shannon_bayes)

# Posterior predictive checks
pp_check(model_shannon_bayes, type = "dens_overlay") +
  ggtitle("Posterior Predictive Check: Shannon Diversity")




# Define the Bayesian linear regression model for ARG Load
model_arg_load_bayes <- brm(
  formula = log10_ARG_load ~ sex_combined + age_category + GDP_per_head + Usage + (1 | region),
  data = train_data_arg,
  family = gaussian(),
  prior = c(
    set_prior("normal(0, 5)", class = "b"),
    set_prior("normal(0, 5)", class = "Intercept"),
    set_prior("cauchy(0, 2)", class = "sd")
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = parallel::detectCores(),
  seed = 123
)

# View summary of the model
summary(model_arg_load_bayes)

# Posterior predictive checks
pp_check(model_arg_load_bayes, type = "dens_overlay") +
  ggtitle("Posterior Predictive Check: ARG Load")