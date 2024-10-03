library(xgboost)
library(dplyr)
library(Matrix)
library(caret)
library(ggplot2)
library(tidyr)

set.seed(123)

#------------------
# Data Preparation
#------------------

tse <- readRDS("TSE_AB_estimate.rds")
df <- as.data.frame(colData(tse))

# Filter for age above 18
age_threshold = 18

df_clean <- df %>%
  filter(host_age_years > age_threshold) %>%
  select(geo_loc_name_country_calc,
         sex_combined, 
         host_age_years,
         log10_ARG_load, 
         World_Bank_Income_Group, 
         Usage, 
         Corruption_and_Governance_Index,
         usage_bayesian, GDP_per_head, 
         Infrastructure_Index, 
         Education_Index,
         Health_Spend_Index
  ) %>%
  drop_na()

# Define target and features
target <- "log10_ARG_load"
features <- setdiff(names(df_clean), target)

# Identify and convert categorical variables to factors
categorical_vars <- c("geo_loc_name_country_calc", "sex_combined", "World_Bank_Income_Group")
df_clean[categorical_vars] <- lapply(df_clean[categorical_vars], as.factor)

# One-hot encoding using model.matrix (excluding the intercept)
formula_str <- paste(target, "~ . -1")
data_matrix <- model.matrix(as.formula(formula_str), data = df_clean)


# Split data into training and testing sets
train_index <- createDataPartition(df_clean[[target]], p = 0.8, list = FALSE)

train_data <- df_clean[train_index, ]
train_matrix <- data_matrix[train_index, ]

test_data <- df_clean[-train_index, ]
test_matrix <- data_matrix[-train_index, ]

train_label <- train_data[[target]]
test_label <- test_data[[target]]

# Convert to DMatrix format for XGBoost
dtrain <- xgb.DMatrix(data = train_matrix, label = train_label)
dtest <- xgb.DMatrix(data = test_matrix, label = test_label)


#------------------
# XGBoost Training
#------------------

# Set model parameters
params <- list(
  booster = "gbtree",
  objective = "reg:squarederror",
  eval_metric = "rmse",
  eta = 0.1,
  max_depth = 6,
  subsample = 0.8,
  colsample_bytree = 0.8
)

# Number of boosting rounds
nrounds <- 100

# Watchlist to monitor training and evaluation performance
watchlist <- list(train = dtrain, eval = dtest)

# Train the XGBoost model with early stopping
model_xgb <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = nrounds,
  watchlist = watchlist,
  early_stopping_rounds = 10,
  print_every_n = 10
)


#------------------
# Model Evaluation
#------------------

# Make predictions on the test set
predictions <- predict(model_xgb, dtest)

# Calculate RMSE
rmse <- sqrt(mean((test_label - predictions)^2))
cat("RMSE on test data:", rmse, "\n")

# Calculate R-squared
sst <- sum((test_label - mean(test_label))^2)
sse <- sum((test_label - predictions)^2)
r_squared <- 1 - sse / sst
cat("R-squared on test data:", r_squared, "\n")


# Feature Importance
importance_matrix <- xgb.importance(feature_names = colnames(train_matrix), model = model_xgb)
print(importance_matrix)

#-----------------
# Cross-Validation
#-----------------

# Perform cross-validation using xgb.cv
cv_results <- xgb.cv(
  params = params,
  data = dtrain,
  nrounds = 100,
  nfold = 5,
  early_stopping_rounds = 10,
  metrics = "rmse",
  verbose = TRUE,
  showsd = TRUE,
  stratified = TRUE
)

# Best number of boosting rounds from CV
best_nrounds <- cv_results$best_iteration
cat("Best number of boosting rounds:", best_nrounds, "\n")


#----------------------------------
# Hyperparameter Tuning with caret
#----------------------------------

# Define trainControl for caret
train_control <- trainControl(
  method = "cv",
  number = 5,
  verboseIter = TRUE
)

# Define a grid of hyperparameters to tune
tune_grid <- expand.grid(
  nrounds = c(50, 100, 150),
  max_depth = c(4, 6, 8),
  eta = c(0.01, 0.1, 0.3),
  gamma = c(0, 1),
  colsample_bytree = c(0.7, 0.8),
  min_child_weight = c(1, 5),
  subsample = c(0.7, 0.8)
)

# Train the model using caret
model_caret <- train(
  x = train_matrix,
  y = train_label,
  method = "xgbTree",
  trControl = train_control,
  tuneGrid = tune_grid,
  metric = "RMSE"
)

# Display the best hyperparameters and model details
print(model_caret$bestTune)
print(model_caret)

#----------------------
# Final Model Training
#----------------------

# Update parameters with the best nrounds from CV
final_params <- params
final_params$nrounds <- best_nrounds

# Retrain the model with the best number of rounds
final_model <- xgb.train(
  params = final_params,
  data = dtrain,
  nrounds = best_nrounds,
  watchlist = watchlist,
  early_stopping_rounds = 10,
  print_every_n = 10
)

# Predictions with the final model
final_predictions <- predict(final_model, dtest)

# Calculate RMSE for the final model
final_rmse <- sqrt(mean((test_label - final_predictions)^2))
cat("Final RMSE on test data:", final_rmse, "\n")

# Calculate R-squared for the final model
final_sst <- sum((test_label - mean(test_label))^2)
final_sse <- sum((test_label - final_predictions)^2)
final_r_squared <- 1 - (final_sse / final_sst)
cat("Final R-squared on test data:", final_r_squared, "\n")

# Feature Importance for the final model
final_importance <- xgb.importance(feature_names = colnames(train_matrix), model = final_model)

# Variable Importance from caret model
var_imp_caret <- varImp(model_caret, scale = FALSE)
print(var_imp_caret)

#saveRDS(final_model, "final_xgb_model.rds")
