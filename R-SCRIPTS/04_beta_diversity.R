 TSE <- readRDS("data/TSE.rds")
 
 # Force R to clean up memory
 gc()
 
 # Subset TSE without copying large objects unnecessarily
 col_data <- colData(TSE)
 
 # Subset where sex_combined is not NA
 non_na_samples <- !is.na(col_data$sex_combined)
 
 # Directly subset without making multiple copies of the object
 TSE_gender <- TSE[, non_na_samples]
 
 # Run garbage collection again after subsetting
 gc()


 
 library(Matrix)
 assay_data <- assays(TSE_gender)$counts
 sparse_assay <- as(assay_data, "sparseMatrix")
 # Replace the original assays with the sparse version
 assays(TSE_gender)$counts <- sparse_assay

 # Randomly sample 7000 indices from the available sample columns
 set.seed(123)  # Set seed for reproducibility
 # Get the total number of samples
 total_samples <- ncol(TSE_gender)
 
 random_samples <- sample(total_samples, 7000)
 
 # Subset the TSE object
 TSE_gender_subset <- TSE_gender[, random_samples]
 
 
 TSE_gender_subset <- scater::runMDS(TSE_gender_subset,
                        FUN = vegan::vegdist,
                        method = "bray",
                       name = "PCoA",
                       ncomponents=3,
                       exprs_values = "relabundance")

colnames(reducedDim(TSE_gender_subset)) <- paste0("PC", seq_len(ncol(reducedDim(TSE_gender_subset))))

# Step 1: Extract the reduced dimensions (PCoA results) from TSE_gender_subset
pcoa_data <- as.data.frame(reducedDim(TSE_gender_subset))

# Step 2: Add the gender information to the PCoA data frame
pcoa_data$sex_combined <- colData(TSE_gender_subset)$sex_combined

# Step 3: Create the PCoA plot, coloring by gender_combined
ggplot(pcoa_data, aes(x = PC1, y = PC2, color = sex_combined)) +
  geom_point(size = 1) +
  labs(title = "PCoA Plot (Bray-Curtis Distance) colored by Gender", 
       x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(legend.title = element_blank())

 saveRDS(TSE_gender_subset, file="data/TSE_gender_subset_ord.rds")
 
 
 # Load necessary library
 library(vegan)
 
 # Load necessary library
 library(vegan)
 
 # Load necessary library
 library(vegan)
 
 # Step 1: Remove samples with NA in sex_combined
 non_na_samples <- !is.na(colData(TSE_gender_subset)$sex_combined)
 TSE_gender_subset_clean <- TSE_gender_subset[, non_na_samples]
 
 # Step 2: Extract the assay data (relative abundances)
 assay_data_clean <- assay(TSE_gender_subset_clean, "relabundance")
 
 # Step 3: Remove samples (columns) with all-zero or near-zero counts across all features
 non_empty_samples <- colSums(assay_data_clean > 0) > 0
 assay_data_clean <- assay_data_clean[, non_empty_samples]
 
 # Ensure no NA values in the assay data
 assay_data_clean <- na.omit(assay_data_clean)
 
 # Step 4: Subset the colData to match the non-empty samples
 coldata_clean <- colData(TSE_gender_subset_clean)[non_empty_samples, ]
 
 # Step 5: Calculate the Bray-Curtis distance matrix on the cleaned data
 distance_matrix <- vegan::vegdist(t(assay_data_clean), method = "bray")  # transpose for samples
 
 # Step 6: Extract the corresponding cleaned sex_combined variable
 sex_combined_clean <- coldata_clean$sex_combined
 
 # Step 7: Run PERMANOVA (adonis2) to calculate the R² for sex_combined
 permanova_result <- adonis2(distance_matrix ~ sex_combined_clean, 
                             data = as.data.frame(coldata_clean),
                             permutations = 999)
 
 # Step 8: Print the PERMANOVA result
 print(permanova_result)
 
 