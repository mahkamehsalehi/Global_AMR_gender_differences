library(tidyverse)
library(vegan)
library(TreeSummarizedExperiment)

# Data Loading and Initial Filtering ************************************************ ####

# Load TSE object
TSE <- readRDS("DATA/TSE_filtered.rds")

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

adult_metadata$log_ARG_load <- log(adult_metadata$ARG_load)

# Add Top 5 AB classes ************************************************************** ####

# TSE <- readRDS("DATA/TSE.rds")

counts_mat <- assay(TSE, "counts")
row_df <- as.data.frame(rowData(TSE))
col_df <- as.data.frame(colData(TSE))

col_df <- transform(col_df, sample = acc)
colnames(counts_mat) <- col_df$acc 

counts_df <- as.data.frame(counts_mat)
counts_df$gene <- rownames(counts_mat)

counts_df <- merge(counts_df,
                   row_df[, c("GENE", "Class")],
                   by.x = "gene", by.y = "GENE")

long_df <- counts_df %>%
  pivot_longer(
    cols = -c(gene, Class),
    names_to = "sample",
    values_to = "count"
  )

agg_counts <- long_df %>%
  group_by(sample, Class) %>%
  summarize(class_total = sum(count, na.rm = TRUE)) %>%
  ungroup()




# Top 5 classes only
top_5_AB_classes <-  c("Tetracycline", 
                       "Beta-lactam", 
                       "Macrolide, Lincosamide, Streptogramin B", 
                       "Aminoglycoside", 
                       "Amphenicol")

agg_counts <- agg_counts %>% 
  filter(Class %in% top_5_AB_classes)

agg_counts_wide <- agg_counts %>% 
  spread(key = "Class", value = "class_total")

# Drop samples not in study
agg_counts_wide <- agg_counts_wide %>% 
  filter(sample %in% adult_metadata$acc)

# Merge with metadata
adult_metadata <- full_join(adult_metadata %>% mutate(sample = acc),
                            agg_counts_wide, 
                            by = "sample")

# Log transform top5 with pseudocount
top5_log <- apply(adult_metadata[, top_5_AB_classes],
                  MARGIN = 2, FUN = function(X){
                    Y <- adult_metadata[, top_5_AB_classes]
                    X[X == 0] <- min(Y[Y != 0])
                    log(X)
                  })

colnames(top5_log) <- paste0("log_", gsub(" ", "_", colnames(top5_log)))

adult_metadata <- cbind(adult_metadata, top5_log)

# Beta Diversity Analysis *********************************************************** ####


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

# Dummy-encoding ******************************************************************** ####
# Create dummy-encoding for glm predictors: 
# sex_combined + age_category + region + GDP_per_head + Usage


temp_df <- adult_metadata

temp_df$age_category_new <- factor(temp_df$age_category_new,
                                   levels = c("Middle-Aged Adult", "Infant",
                                              "Toddler", "Children", "Teenager",
                                              "Young Adult", "Older Adult", "Oldest Adult"))

# High income vs. others
temp_df$income_group_HIC <- ifelse(temp_df$World_Bank_Income_Group == "High income", 1, 0)
# Antibiotic usage <10 vs. >=10
temp_df$Usage_high <- ifelse(temp_df$Usage < 10, 0, 1)
# Sex into numeric
temp_df$sex_num_Men <- ifelse(temp_df$sex_combined == "Men", 0, 1)

# dummy coding
temp_df_one_hot <- model.matrix(~ 0 + 
                                  sex_num_Men +
                                  Region + 
                                  age_category_new + 
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