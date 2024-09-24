library(tidyquant)
library(tidyverse)
library(SummarizedExperiment)



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#--------------- Adding Antibiotic Use Data to the TSE Object ------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#source("scripts/01_data_preprocessing.R")

ab_use <- g_data_2016 <- read_csv("global_data_2016.csv")
tse <- readRDS("TSE.rds")

common_samples <- intersect(colData(tse)$acc, ab_use$acc)
cat("Number of common samples:", length(common_samples), "\n")

# Filter antibiotic use data to common samples
ab_use_filtered <- ab_use %>%
  filter(acc %in% common_samples)

ab_use_filtered <- dplyr::rename(ab_use_filtered, usage_bayesian = antibiotic_consumption)

tse_coldata <- as.data.frame(colData(tse))

# Add usage_bayesian to tse_coldata
merged_coldata <- tse_coldata %>%
  left_join(ab_use_filtered %>% select(acc, usage_bayesian), by = "acc")

# Preserve rownames
rownames(merged_coldata) <- rownames(tse_coldata)

# Update colData in tse
colData(tse) <- DataFrame(merged_coldata)

#saveRDS(tse, file = "TSE_AB_estimate.rds")

# Check for missing values
num_missing <- sum(is.na(colData(tse)$usage_bayesian))
cat("Number of samples without Bayesian_usage data:", num_missing, "\n")