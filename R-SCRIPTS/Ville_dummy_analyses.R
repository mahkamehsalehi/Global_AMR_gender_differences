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
TSE <- readRDS("DATA/TSE.rds")

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
# assays(TSE_adult)$counts <- as(assays(TSE_adult)$counts, "sparseMatrix")
# 
# # Run PCoA
# TSE_adult <- runMDS(
#   TSE_adult,
#   FUN = vegan::vegdist,
#   method = "bray",
#   name = "PCoA",
#   ncomponents = 3,
#   exprs_values = "relabundance"
# )

# Extract PCoA results and variance explained
# pcoa_adult_data <- as.data.frame(reducedDim(TSE_adult, "PCoA"))
# colnames(pcoa_adult_data) <- paste0("PC", 1:ncol(pcoa_adult_data))
# eigenvalues_adult <- attr(reducedDim(TSE_adult, "PCoA"), "eig")
# variance_explained_adult <- (eigenvalues_adult[eigenvalues_adult > 0] / sum(eigenvalues_adult)) * 100
# 
# # Add metadata to PCoA data
# pcoa_adult_data <- pcoa_adult_data %>%
#   mutate(
#     geo = adult_metadata$geo_loc_name_country_continent_calc,
#     host_age_years = adult_metadata$host_age_years,
#     sex_combined = adult_metadata$sex_combined,
#     adult_age_category = case_when(
#       host_age_years >= 18 & host_age_years < 35 ~ "Young Adult",
#       host_age_years >= 35 & host_age_years < 65 ~ "Middle-age Adult",
#       host_age_years >= 65 ~ "Older Adult",
#       TRUE ~ NA_character_
#     )
#   )

# Define custom plot theme
# common_theme <- theme_classic(base_size = 14) +
#   theme(
#     text = element_text(family = "Sans"),
#     plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
#     axis.text = element_text(size = 12),
#     axis.title = element_text(size = 13),
#     legend.position = "none",
#     axis.line = element_line(color = "black"),
#     strip.background = element_rect(fill = "white", color = "black"),
#     strip.text = element_text(size = 12, face = "bold")
#   )

#-------------------------------------------------------------------------------
# Plot PCoA Results
#-------------------------------------------------------------------------------

# # Plot by region
# pcoa_region <- ggplot(pcoa_adult_data, aes(x = PC1, y = PC2, color = geo)) +
#   geom_point(size = 2, alpha = 0.8) +
#   labs(
#     title = "PCoA Plot (Beta Diversity - Bray-Curtis) for Adults by Region",
#     x = paste0("PC1 (", round(variance_explained_adult[1], 2), "%)"),
#     y = paste0("PC2 (", round(variance_explained_adult[2], 2), "%)"),
#     color = "Region"
#   ) +
#   common_theme
# ggsave("final/pcoa_region.png", pcoa_region, width = 10, height = 8)
# 
# # Plot by gender across regions
# pcoa_gender_region <- ggplot(pcoa_adult_data, aes(x = PC1, y = PC2, color = sex_combined)) +
#   geom_point(size = 2, alpha = 0.8) +
#   facet_wrap(~ geo, nrow = 2) +
#   labs(
#     title = "PCoA Plot (Beta Diversity - Bray-Curtis) for Adults by Gender",
#     x = paste0("PC1 (", round(variance_explained_adult[1], 2), "%)"),
#     y = paste0("PC2 (", round(variance_explained_adult[2], 2), "%)"),
#     color = "Gender"
#   ) +
#   common_theme
# ggsave("final/pcoa_gender_region.png", pcoa_gender_region, width = 10, height = 8)
# 
# # Plot by age category across regions
# pcoa_age_cat <- ggplot(pcoa_adult_data, aes(x = PC1, y = PC2, color = adult_age_category)) +
#   geom_point(size = 2, alpha = 0.8) +
#   facet_wrap(~ geo, nrow = 2) +
#   labs(
#     title = "PCoA Plot (Beta Diversity - Bray-Curtis) by Age Category",
#     x = paste0("PC1 (", round(variance_explained_adult[1], 2), "%)"),
#     y = paste0("PC2 (", round(variance_explained_adult[2], 2), "%)"),
#     color = "Age Category"
#   ) +
#   common_theme +
#   theme(legend.position = "right")
# ggsave("final/pcoa_age_cat.png", pcoa_age_cat, width = 10, height = 8)

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
# distance_matrix <- vegan::vegdist(t(assay_data_clean), method = "bray")

# Calculate Shannon diversity
diversity_indices <- diversity(t(assay_data_clean), index = "shannon")
adult_metadata$shannon_diversity <- diversity_indices

# Statistical test
# diversity_test <- t.test(shannon_diversity ~ sex_combined, data = adult_metadata)
# print(diversity_test)

# PERMANOVA Analysis
# set.seed(123)
# permanova_result <- adonis2(
#   distance_matrix ~ sex_combined + age_category + region + GDP_per_head + Usage,
#   data = adult_metadata,
#   permutations = 999
# )
# print(permanova_result)
# saveRDS(permanova_result, file = "final/permanova_result.rds")

# Test for homogeneity of dispersion (betadisper)
# dispersion_sex <- betadisper(distance_matrix, adult_metadata$sex_combined)
# dispersion_test <- anova(dispersion_sex)
# print(dispersion_test)


## Dummy-encoding *************** ####
# Create dummy-encoding for glm predictors: 
# sex_combined + age_category + region + GDP_per_head + Usage


temp_df <- adult_metadata

# High income vs. others
temp_df$income_group_HIC <- ifelse(temp_df$World_Bank_Income_Group == "High income", 1, 0)
# Antibiotic usage <11 vs. >=11
temp_df$Usage_high <- ifelse(temp_df$Usage < 11, 0, 1)
# Sex into numeric
temp_df$sex_num_Men <- ifelse(temp_df$sex_combined == "Men", 1, 0)

# dummy coding
temp_df_one_hot <- model.matrix(~ 0 + 
                                  sex_num_Men + 
                                  region + 
                                  age_category + 
                                  income_group_HIC + 
                                  Usage_high,
                                data = temp_df) %>% 
  as.data.frame

# Set colnames
colnames(temp_df_one_hot) <- gsub(" |-", "_",  colnames(temp_df_one_hot))
dummy_var_names <- colnames(temp_df_one_hot)

# Add dummies into metadata
adult_metadata <- cbind(adult_metadata, temp_df_one_hot)

# Remove temp tables
rm(temp_df, temp_df_one_hot)

## Dummy-encoded linear models ****************** ####


log10_ARG_formula <- paste0("log10_ARG_load ~ ",
                            paste0(dummy_var_names[!grepl(pattern = "Africa|Asia|Infant", x = dummy_var_names)],
                                   collapse = "+")) %>% as.formula()

brm_log10_ARG_load_dummy <- brm(formula = log10_ARG_formula,
                        data = adult_metadata,
                        family = gaussian(), 
                        prior = c(
                          set_prior("normal(0, 1)", class = "b"),
                          set_prior("normal(0, 1)", class = "Intercept")
                        ), 
                        chains = 2, 
                        iter = 5000,
                        cores = parallel::detectCores(),
)


shannon_formula <- paste0("shannon_diversity ~ ",
                          paste0(dummy_var_names[!grepl(pattern = "Africa|Asia|Infant", x = dummy_var_names)],
                                 collapse = "+")) %>% as.formula()

brm_shannon_dummy <- brm(formula = shannon_formula,
                          data = adult_metadata,
                          family = gaussian(), 
                          prior = c(
                            set_prior("normal(0, 1)", class = "b"),
                            set_prior("normal(0, 1)", class = "Intercept")
                          ), 
                          chains = 2, 
                          iter = 5000,
                          cores = parallel::detectCores(),
)


## Get summaries
log10_ARG_summary <- posterior_summary(brm_log10_ARG_load_dummy, pars = "b") %>% 
  data.frame() %>% 
  rownames_to_column(var = "Feature") %>% 
  mutate(Feature = gsub("b_", "", Feature), 
         Response = "log10_ARG") %>% 
  dplyr::select(-Est.Error) %>% 
  mutate(Significant = ifelse(Q2.5 > 0 | Q97.5 < 0, TRUE, FALSE))


shannon_summary <- posterior_summary(brm_shannon_dummy, pars = "b") %>% 
  data.frame() %>% 
  rownames_to_column(var = "Feature") %>% 
  mutate(Feature = gsub("b_", "", Feature), 
         Response = "Shannon") %>% 
  dplyr::select(-Est.Error) %>% 
  mutate(Significant = ifelse(Q2.5 > 0 | Q97.5 < 0, TRUE, FALSE))
  

full_summary <- rbind(log10_ARG_summary, shannon_summary)

full_summary %>%  
  filter(Feature != "Intercept") %>% 
  ggplot() + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(x = Feature, ymin = Q2.5, ymax = Q97.5, color = Response), 
                position = "dodge") +
  # facet_wrap(~Response) +
  coord_flip() +
  labs(title = "95% CIs")



