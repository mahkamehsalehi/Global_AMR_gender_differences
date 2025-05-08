# Load required libraries
library(tidyverse)
library(vegan)
library(mia)
library(dplyr)

# Load TSE object
TSE <- readRDS("../DATA/TSE.rds")

# Set seed for reproducibility
set.seed(1235)

# Get assay data and transpose

assay_data <- t(assay(TSE,"relabundance"))


# Run multivariable permanova using adonis

permanova_multivariable <- vegan::adonis(
  assay_data ~ sex_combined + age_category + geo_loc_name_country_continent_calc  + GDP_per_head + Usage,
  data = colData(TSE),
  permutations = 999,
  method = "bray", 
  na.action = na.omit
)

# ---------------------------
#  Display Results
# ---------------------------
print("Multivariate PERMANOVA Results:")
print(permanova_multivariable)

