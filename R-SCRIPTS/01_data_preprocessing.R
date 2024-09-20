# This need to be run just once to create the TSE object. Re-run if modifications to metadata are made
# TODO # Add the antibiotic use from 2016.
# TODO # Run the ordinations from lines 228 onwards.

rm(list  = ls())
# Load necessary libraries
require(tidyverse)
require(dplyr)
require(readr)
require(tidyquant)
require(data.table)
library(mia)
library(TreeSummarizedExperiment)

# To ensure that ordination signs stay the same in different runs etc.
set.seed(252452)

# Set working directory to github repo directory
# setwd("/PATH/TO/Global_AMR/")

# Source functions
source("R-SCRIPTS/00_functions.R")

# Load data
gene_expression_data_80 <- fread("DATA/RAW/ALL_mapstat_80_mat.txt", header = TRUE, sep = "\t")
resfinder_phenotypes <- fread("DATA/resfinder_phenotypes-pjji7t1qo3nedped5943f1c1wy.txt", header = TRUE, sep = "\t", fill = TRUE)
gene_lengths <- fread("DATA/gene_length.csv", header = TRUE)
metadata <- fread("DATA/Sra_metadata_jun11.txt", header = TRUE, sep = ",")

# Create gender and age columns for metadata
source("R-SCRIPTS/metadata_processing_sex.R")

# Merge with other metadata
tmp<-merge(df_selected, metadata)

metadata <- tmp

# Extract column names from gene_expression_data_80 and create a dataframe for accessions
column_names <- names(gene_expression_data_80)
accession_ids <- data.frame(acc = column_names[-1])

# Filter metadata for relevant accessions
filtered_metadata <- metadata %>%
  filter(acc %in% accession_ids$acc)

# Calculate read count for each sample
filtered_metadata$readcount <- filtered_metadata$mbases / filtered_metadata$avgspotlen / 2

# Filter resfinder_phenotypes data for genes present in gene_expression_data_80
filtered_resfinder <- resfinder_phenotypes %>%
  filter(`Gene_accession no.` %in% gene_expression_data_80$GENE)

# Filter gene_expression_data_80 for relevant genes
filtered_assay <- gene_expression_data_80 %>%
  filter(GENE %in% resfinder_phenotypes$`Gene_accession no.`)

# Select relevant columns from filtered_assay
accessions_to_keep <- filtered_metadata$acc
assay_subset <- filtered_assay %>% 
  select(1, all_of(accessions_to_keep))

# Filter resfinder_phenotypes for unique genes
GENEs_to_keep <- unique(filtered_assay$GENE)
resfinder_subset <- resfinder_phenotypes %>%
  filter(resfinder_phenotypes$`Gene_accession no.` %in% GENEs_to_keep)
unique_resfinder <- resfinder_subset %>%
  distinct(`Gene_accession no.`, .keep_all = TRUE) 

# Rename column for consistency
names(unique_resfinder)[names(unique_resfinder) == "Gene_accession no."] <- "GENE"

# Set row names for assay_subset and transpose
assay_subset <- data.frame(assay_subset, row.names = 1)
transposed_assay <- t(assay_subset)

# Convert to data frame and add accessions as a column
transposed_assay <- as.data.frame(transposed_assay)
transposed_assay <- transposed_assay %>% rownames_to_column(var = "acc")

# Merge with metadata to include readcount
merged_assay_metadata <- merge(transposed_assay, filtered_metadata[, c("acc", "readcount")], by = "acc")

# Filter gene lengths data for relevant genes
filtered_gene_lengths <- gene_lengths %>%
  filter(GENE %in% GENEs_to_keep)

# Convert gene lengths to kilobases
filtered_gene_lengths$length_kb <- filtered_gene_lengths$gene_length / 1000


# Ensure gene expression matrix is a dataframe and add gene names as a column
normalized_assay <- assay_subset %>%
  rownames_to_column(var = "GENE")

# Merge the gene expression matrix with gene lengths
merged_assay_metadata <- merge(normalized_assay, filtered_gene_lengths, by = "GENE")

# Calculate RPK (Reads Per Kilobase)
rpk_values <- merged_assay_metadata %>%
  mutate(across(-c(GENE, gene_length, length_kb), ~ . / length_kb))

# Extract only the RPK values
rpk_values <- rpk_values %>% select(-c(GENE, gene_length, length_kb))


# Check the order of 'acc' and column names
identical(filtered_metadata$acc, colnames(rpk_values))

# If they are not identical, reorder the columns of 'rpk_values' to match 'acc'
rpk_values <- rpk_values[, match(filtered_metadata$acc, colnames(rpk_values))]

# Check again to confirm they are now in the same order
identical(filtered_metadata$acc, colnames(rpk_values))
# Extract library sizes
library_sizes <- filtered_metadata$readcount


# Ensure library_sizes is numeric
library_sizes <- as.numeric(library_sizes)


# Calculate RPKM (Reads Per Kilobase per Million mapped reads)
rpkm_values <- sweep(as.matrix(rpk_values), 2, library_sizes, FUN = "/")


# Convert back to dataframe and add gene names
rpkm_dataframe <- as.data.frame(rpkm_values)
rpkm_dataframe <- cbind(GENE = merged_assay_metadata$GENE, rpkm_dataframe)

# Set row names for the assay dataframe
final_assay <- data.frame(rpkm_dataframe, row.names = 1)

# Create TreeSummarizedExperiment object
counts <- final_assay
tax <- unique_resfinder
tax <- as.data.frame(tax)
samples <- filtered_metadata

# Add sums of all ARGs in RPKM to sample data (ARG load)
samples$ARG_load <- colSums(counts)

tree_summarized_experiment <- TreeSummarizedExperiment(
  assays = SimpleList(counts = as.matrix(counts)),
  colData = DataFrame(samples),
  rowData = DataFrame(tax)
)

# # Save the TreeSummarizedExperiment object
saveRDS(tree_summarized_experiment, file="DATA/TSE.rds")

rm(list  = ls())
# 
TSE <- readRDS("DATA/TSE.rds")

## -------- Filtering ----------------

# Subset samples

# Step 1: Identify columns with sums of 0 or NA
col_sums <- colSums(assay(TSE), na.rm = FALSE)

# Step 2: Find columns to keep (those not equal to 0 and not NA)
columns_to_keep <- !(is.na(col_sums) | col_sums == 0)

# Step 3: Subset the TreeSummarizedExperiment object to keep only those columns
TSE <- TSE[, columns_to_keep]

# Subset features
# Step 1: Identify OTUs with non-zero prevalence in at least one sample
otus_to_keep <- rowSums(assay(TSE) != 0) > 0

# Step 2: Subset the TreeSummarizedExperiment object to keep only these OTUs
TSE <- TSE[otus_to_keep, ]

# Output the filtered TreeSummarizedExperiment object
TSE

## --------- Alpha diversities -------------------

# Add observed ARG diversity
TSE <- mia::estimateDiversity(TSE, assay.type="counts", index = "shannon", name="ARG_div")

# Add observed ARG richness
TSE <- mia::estimateRichness(TSE, assay.type="counts", index = "observed", name="ARG_obs")


## ----  Transformations ----------------

# Add relative abundances to ARG data (for barplots; not for statistical analyses)
TSE <- transformAssay(TSE, assay.type="counts", method="relabundance", MARGIN="samples")

# Get country level data for metadata
source("R-SCRIPTS/metadata_processing_country.R")

## ------ More filtering -------
# Apply the filtering criteria to TSE colData
filtered_data <- colData(TSE) %>% data.frame %>%
  filter(readcount > 0.5, ARG_load > 50, ARG_obs > 1)

# Subset the TSE object to include only samples that pass the filtering criteria
# Assuming the rownames of TSE colData correspond to the sample names in TSE
filtered_samples <- filtered_data %>% rownames()

TSE <- TSE[,filtered_samples]

source("R-SCRIPTS/metadata_processing_bioproject.R")

# Save the TreeSummarizedExperiment object
saveRDS(TSE, file="DATA/TSE.rds")

## TODO ##

# ## ---------- Ordinations ------------
# Run on CSC
# # Perform PCoA with Bray-Curtis on relative abundances ## INPROGRESS!
#!/bin/bash -l
#SBATCH --job-name=MDS
#SBATCH --account=project_2008149
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=test
#SBATCH --time=05:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=100000


# Ordinations TODO
# library(scater)
# library(vegan)
# TSE <- readRDS("TSE.rds")
# TSE <- scater::runMDS(TSE,
#                        FUN = vegan::vegdist,
#                        method = "bray",
#                        name = "PCoA",
#                        ncomponents=3,
#                        exprs_values = "relabundance")
# 
# colnames(reducedDim(TSE)) <- paste0("PC", seq_len(ncol(reducedDim(TSE))))
#  
#  # Switch (the arbitrary) sign so that the figure is more directly
#  # comparable with earlier publications
#  reducedDim(TSE)[,1] <- -reducedDim(TSE)[,1]
#  reducedDim(TSE)[,2] <- -reducedDim(TSE)[,2]
#  saveRDS(TSE, file="TSE_ord.rds")
