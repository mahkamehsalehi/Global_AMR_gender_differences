# Set working directory and load necessary libraries
# setwd("/scratch/project_2008149/USER_WORKSPACES/mahkameh/women_amr/")
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
library(xgboost)
library(brms)
library(glmnet)

#-------------------------------------------------------------------------------
# Data Loading and Initial Filtering
#-------------------------------------------------------------------------------

# Load TSE object
TSE <- readRDS("../DATA/TSE.rds")

# Filter samples with complete sex and age data
non_na_samples <- !is.na(colData(TSE)$sex_combined) & !is.na(colData(TSE)$host_age_years)
TSE_filtered <- TSE[, non_na_samples]

# Further filter for adult samples with geographic information
TSE_adult <- TSE_filtered[, colData(TSE_filtered)$geo_loc_name_country_continent_calc != "uncalculated"]

# Extract metadata and set factors for consistent ordering in plots
adult_metadata <- as.data.frame(colData(TSE_adult)) %>%
  mutate(
    sex_combined = recode(sex_combined, male = "Men", female = "Women"),
    sex_combined = factor(sex_combined, levels = c("Men", "Women"))
  )

#-------------------------------------------------------------------------------
# Bray-Curtis PCoA Analysis
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
pcoa_adult_data <- as.data.frame(reducedDim(TSE_adult, "PCoA"))
colnames(pcoa_adult_data) <- paste0("PC", 1:ncol(pcoa_adult_data))
eigenvalues_adult <- attr(reducedDim(TSE_adult, "PCoA"), "eig")
variance_explained_adult <- (eigenvalues_adult[eigenvalues_adult > 0] / sum(eigenvalues_adult)) * 100

# Add metadata to PCoA data
pcoa_adult_data <- pcoa_adult_data %>%
  mutate(
    geo = adult_metadata$geo_loc_name_country_continent_calc,
    host_age_years = adult_metadata$host_age_years,
    sex_combined = adult_metadata$sex_combined,
    adult_age_category = case_when(
      host_age_years >= 18 & host_age_years < 35 ~ "Young Adult",
      host_age_years >= 35 & host_age_years < 65 ~ "Middle-age Adult",
      host_age_years >= 65 ~ "Older Adult",
      TRUE ~ NA_character_
    )
  )

# Define custom plot theme
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
# Plot PCoA Results
#-------------------------------------------------------------------------------

# Plot by region
pcoa_region <- ggplot(pcoa_adult_data, aes(x = PC1, y = PC2, color = geo)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(
    title = "PCoA Plot (Beta Diversity - Bray-Curtis) for Adults by Region",
    x = paste0("PC1 (", round(variance_explained_adult[1], 2), "%)"),
    y = paste0("PC2 (", round(variance_explained_adult[2], 2), "%)"),
    color = "Region"
  ) +
  common_theme
ggsave("final/pcoa_region.png", pcoa_region, width = 10, height = 8)

# Plot by gender across regions
pcoa_gender_region <- ggplot(pcoa_adult_data, aes(x = PC1, y = PC2, color = sex_combined)) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(~ geo, nrow = 2) +
  labs(
    title = "PCoA Plot (Beta Diversity - Bray-Curtis) for Adults by Gender",
    x = paste0("PC1 (", round(variance_explained_adult[1], 2), "%)"),
    y = paste0("PC2 (", round(variance_explained_adult[2], 2), "%)"),
    color = "Gender"
  ) +
  common_theme
ggsave("final/pcoa_gender_region.png", pcoa_gender_region, width = 10, height = 8)

# Plot by age category across regions
pcoa_age_cat <- ggplot(pcoa_adult_data, aes(x = PC1, y = PC2, color = adult_age_category)) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(~ geo, nrow = 2) +
  labs(
    title = "PCoA Plot (Beta Diversity - Bray-Curtis) by Age Category",
    x = paste0("PC1 (", round(variance_explained_adult[1], 2), "%)"),
    y = paste0("PC2 (", round(variance_explained_adult[2], 2), "%)"),
    color = "Age Category"
  ) +
  common_theme +
  theme(legend.position = "right")
ggsave("final/pcoa_age_cat.png", pcoa_age_cat, width = 10, height = 8)

#-------------------------------------------------------------------------------
# Beta Diversity Analysis
#-------------------------------------------------------------------------------

# Define age categories and filter metadata
adult_metadata <- adult_metadata %>%
  mutate(
    age_category = case_when(
      host_age_years < 1 ~ "Infant",
      host_age_years < 3 ~ "Toddler",
      host_age_years < 18 ~ "Child",
      host_age_years < 35 ~ "Young Adult",
      host_age_years < 65 ~ "Middle-Aged Adult",
      host_age_years <= 100 ~ "Older Adult",
      TRUE ~ NA_character_
    ),
    region = factor(geo_loc_name_country_continent_calc),
    GDP_per_head = as.numeric(GDP_per_head)
  ) %>%
  filter(!is.na(Usage), !is.na(GDP_per_head))

non_na_conditions <- !is.na(colData(TSE_adult)$GDP_per_head) & !is.na(colData(TSE_adult)$Usage)
TSE_adult_filtered <- TSE_adult[, non_na_conditions]

# Update TSE_adult_filtered colData with the cleaned metadata
colData(TSE_adult_filtered) <- DataFrame(adult_metadata)

# Extract and clean assay data (relative abundances)
assay_data_clean <- assay(TSE_adult_filtered, "relabundance")
non_empty_samples <- colSums(assay_data_clean > 0, na.rm = TRUE) > 0
assay_data_clean <- assay_data_clean[, non_empty_samples]

# Subset colData to match non-empty samples
adult_metadata <- colData(TSE_adult_filtered)[non_empty_samples, ]

# Calculate the Bray-Curtis distance matrix on the cleaned data
distance_matrix <- vegan::vegdist(t(assay_data_clean), method = "bray")

# Calculate Shannon diversity
diversity_indices <- diversity(t(assay_data_clean), index = "shannon")
adult_metadata$shannon_diversity <- diversity_indices

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
saveRDS(permanova_result, file = "final/permanova_result.rds")

# Test for homogeneity of dispersion (betadisper)
dispersion_sex <- betadisper(distance_matrix, adult_metadata$sex_combined)
dispersion_test <- anova(dispersion_sex)
print(dispersion_test)

#-------------------------------------------------------------------------------
# Statistical Modeling
#-------------------------------------------------------------------------------

# GLM for Shannon Diversity
glm_diversity <- glm(
  shannon_diversity ~ sex_combined + age_category + region + GDP_per_head + Usage,
  data = adult_metadata,
  family = gaussian()
)
summary(glm_diversity)

# GLM for log10_ARG_load
glm_arg_load <- glm(
    log10_ARG_load ~ sex_combined + age_category + region + GDP_per_head + Usage,
    data = adult_metadata,
    family = gaussian()
  )
  summary(glm_arg_load)

#}

save.image()

# Variance Inflation Factor (VIF) to check multicollinearity
vif_model <- lm(shannon_diversity ~ sex_combined + age_category + region + GDP_per_head + Usage, data = adult_metadata)
vif_values <- vif(vif_model)
print(vif_values)

# Random Forest for Shannon Diversity
set.seed(123)
external_predictors <- adult_metadata[, c("sex_combined", "age_category", "region", "GDP_per_head", "Usage")]
data_rf_shannon <- data.frame(shannon_diversity = adult_metadata$shannon_diversity, external_predictors)

# Split the data into training (80%) and testing (20%) sets
train_index_shannon <- createDataPartition(data_rf_shannon$shannon_diversity, p = 0.8, list = FALSE)
train_data_shannon <- data_rf_shannon[train_index_shannon, ]
test_data_shannon <- data_rf_shannon[-train_index_shannon, ]

# Train the Random Forest Regression model for Shannon Diversity
rf_shannon <- randomForest(
  shannon_diversity ~ .,
  data = train_data_shannon,
  importance = TRUE,
  ntree = 500,
  mtry = floor(sqrt(ncol(train_data_shannon) - 1))  # Common default
)
print(rf_shannon)

# Predict on the testing set
predictions_shannon <- predict(rf_shannon, newdata = test_data_shannon)

# Calculate performance metrics
actual_shannon <- test_data_shannon$shannon_diversity
r_squared_shannon <- cor(actual_shannon, predictions_shannon)^2
rmse_shannon <- sqrt(mean((actual_shannon - predictions_shannon)^2))
print(paste("Shannon Diversity - R-squared:", round(r_squared_shannon, 4)))
print(paste("Shannon Diversity - RMSE:", round(rmse_shannon, 4)))

# Plot Actual vs. Predicted
ggplot(data = NULL, aes(x = actual_shannon, y = predictions_shannon)) +
  geom_point(alpha = 0.6, color = "blue") +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  theme_minimal() +
  labs(title = "Actual vs. Predicted Shannon Diversity",
       x = "Actual Shannon Diversity",
       y = "Predicted Shannon Diversity")

# Random Forest for ARG Load
data_rf_arg <- data.frame(log10_ARG_load = adult_metadata$log10_ARG_load, external_predictors)
train_index_arg <- createDataPartition(data_rf_arg$log10_ARG_load, p = 0.8, list = FALSE)
train_data_arg <- data_rf_arg[train_index_arg, ]
test_data_arg <- data_rf_arg[-train_index_arg, ]

# Train the Random Forest Regression model for ARG Load
rf_arg_load <- randomForest(
  log10_ARG_load ~ .,
  data = train_data_arg,
  importance = TRUE,
  ntree = 500,
  mtry = floor(sqrt(ncol(train_data_arg) - 1))
)
print(rf_arg_load)

# Predict on the testing set
predictions_arg_load <- predict(rf_arg_load, newdata = test_data_arg)

# Calculate performance metrics
actual_arg_load <- test_data_arg$log10_ARG_load
r_squared_arg <- cor(actual_arg_load, predictions_arg_load)^2
rmse_arg <- sqrt(mean((actual_arg_load - predictions_arg_load)^2))
print(paste("ARG Load - R-squared:", round(r_squared_arg, 4)))
print(paste("ARG Load - RMSE:", round(rmse_arg, 4)))

# Plot Actual vs. Predicted
ggplot(data = NULL, aes(x = actual_arg_load, y = predictions_arg_load)) +
  geom_point(alpha = 0.6, color = "green") +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  theme_minimal() +
  labs(title = "Actual vs. Predicted ARG Load",
       x = "Actual log10_ARG_load",
       y = "Predicted log10_ARG_load")

#-------------------------------------------------------------------------------
# Bayesian Analysis for Shannon Diversity
#-------------------------------------------------------------------------------

# Define the Bayesian linear regression model for Shannon Diversity
model_shannon_bayes <- brm(
  formula = shannon_diversity ~ sex_combined + age_category + GDP_per_head + Usage + (1 | region),
  data = train_data_shannon,
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
summary(model_shannon_bayes)

# Posterior predictive checks
pp_check(model_shannon_bayes, type = "dens_overlay") +
  ggtitle("Posterior Predictive Check: Shannon Diversity")

#-------------------------------------------------------------------------------
# Bayesian Analysis for ARG Load
#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------
# XGBoost Analysis
#-------------------------------------------------------------------------------

# Prepare data
train_matrix <- model.matrix(shannon_diversity ~ ., data = train_data_shannon)[, -1]
test_matrix <- model.matrix(shannon_diversity ~ ., data = test_data_shannon)[, -1]
train_label <- train_data_shannon$shannon_diversity
test_label <- test_data_shannon$shannon_diversity

# Convert to DMatrix
dtrain <- xgb.DMatrix(data = train_matrix, label = train_label)
dtest <- xgb.DMatrix(data = test_matrix, label = test_label)

# Set parameters
params <- list(
  objective = "reg:squarederror",
  eval_metric = "rmse"
)

# Train the model
xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 100,
  watchlist = list(train = dtrain, eval = dtest),
  early_stopping_rounds = 10,
  print_every_n = 10
)

# Predict and evaluate
predictions_xgb <- predict(xgb_model, dtest)
r_squared_xgb <- cor(test_label, predictions_xgb)^2
rmse_xgb <- sqrt(mean((test_label - predictions_xgb)^2))

print(paste("XGBoost Shannon Diversity - R-squared:", round(r_squared_xgb, 4)))
print(paste("XGBoost Shannon Diversity - RMSE:", round(rmse_xgb, 4)))

#-------------------------------------------------------------------------------
# Lasso Regression Analysis
#-------------------------------------------------------------------------------

# Prepare data
x_train <- model.matrix(shannon_diversity ~ ., data = train_data_shannon)[, -1]
y_train <- train_data_shannon$shannon_diversity
x_test <- model.matrix(shannon_diversity ~ ., data = test_data_shannon)[, -1]
y_test <- test_data_shannon$shannon_diversity

# Fit Lasso model
lasso_model <- cv.glmnet(x_train, y_train, alpha = 1, family = "gaussian")

# Best lambda
best_lambda <- lasso_model$lambda.min

# Predict
predictions_lasso <- predict(lasso_model, s = best_lambda, newx = x_test)

# Evaluate
r_squared_lasso <- cor(y_test, predictions_lasso)^2
rmse_lasso <- sqrt(mean((y_test - predictions_lasso)^2))

print(paste("Lasso Shannon Diversity - R-squared:", round(r_squared_lasso, 4)))
print(paste("Lasso Shannon Diversity - RMSE:", round(rmse_lasso, 4)))

#-------------------------------------------------------------------------------
# Residuals vs Fitted for Shannon Diversity GLM
#-------------------------------------------------------------------------------

par(mfrow = c(2, 2))
plot(glm_diversity)

# Save plots to file
png("glm_diversity_plots.png", width = 800, height = 800)
par(mfrow = c(2, 2))
plot(glm_diversity)
dev.off()

#-------------------------------------------------------------------------------
# Scatterplots for ARG Load vs GDP per Head and Shannon Diversity vs GDP
#-------------------------------------------------------------------------------

adult_metadata_df <- as.data.frame(adult_metadata)
adult_metadata_df$sex_combined <- factor(adult_metadata_df$sex_combined, levels = c("Women", "Men"))

# Scatterplot for ARG Load vs GDP
g
