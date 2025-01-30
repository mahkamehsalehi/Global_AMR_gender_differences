# This need to be run just once to create the TSE object. Re-run if modifications to metadata are made
# TODO # Add the antibiotic use from 2016.
# TODO # Run the ordinations from lines 228 onwards.

rm(list  = ls())

require(tidyverse)
require(dplyr)
require(readr)
require(data.table)
library(mia)
library(TreeSummarizedExperiment)


# To ensure that ordination signs stay the same in different runs etc.
set.seed(252452)

# Set working directory to github repo directory
# setwd("/PATH/TO/Global_AMR/")

# Source functions
source("R-SCRIPTS/00_functions.R")

############################################################################
## 1) Load data
############################################################################

gene_expression_data_80 <- fread("DATA/ALL_mapstat_80_mat.txt", 
                                 header = TRUE, sep = "\t")
resfinder_phenotypes    <- fread("DATA/resfinder_phenotypes-pjji7t1qo3nedped5943f1c1wy.txt", 
                                 header = TRUE, sep = "\t", fill = TRUE)
gene_lengths            <- fread("DATA/gene_length.csv", header = TRUE)
metadata                <- fread("DATA/PROCESSED/Sra_metadata_processed.csv", 
                                 header = TRUE, sep = ",")

# Create gender and age columns for metadata
source("R-SCRIPTS/metadata_processing_sex.R")

# Remove jattr from metadata before merging
metadata <- metadata %>%
  select(-matches("^jattr_parsed_jattr"), -jattr)

# Merge with other metadata
tmp <- merge(df_selected, metadata)
metadata <- tmp

############################################################################
## 2) Filter and prepare metadata
############################################################################

column_names   <- names(gene_expression_data_80)
accession_ids  <- data.frame(acc = column_names[-1])

# Filter metadata for relevant accessions
filtered_metadata <- metadata %>%
  filter(acc %in% accession_ids$acc)

# Calculate read count for each sample
filtered_metadata$readcount <- filtered_metadata$mbases / 
  filtered_metadata$avgspotlen / 2

############################################################################
## 3) Filter ResFinder data to keep genes present in gene_expression_data_80
############################################################################

filtered_resfinder <- resfinder_phenotypes %>%
  filter(`Gene_accession no.` %in% gene_expression_data_80$GENE)

# Filter gene_expression_data_80 for relevant genes
filtered_assay <- gene_expression_data_80 %>%
  filter(GENE %in% filtered_resfinder$`Gene_accession no.`)

accessions_to_keep <- filtered_metadata$acc

# Keep only the GENE column plus the relevant accession columns
assay_subset <- filtered_assay %>% 
  dplyr::select(1, all_of(accessions_to_keep))

# Keep unique genes in resfinder
GENEs_to_keep <- unique(filtered_assay$GENE)
resfinder_subset <- resfinder_phenotypes %>%
  filter(resfinder_phenotypes$`Gene_accession no.` %in% GENEs_to_keep)
unique_resfinder <- resfinder_subset %>%
  distinct(`Gene_accession no.`, .keep_all = TRUE)

# Rename column to "GENE" for consistency
names(unique_resfinder)[names(unique_resfinder) == "Gene_accession no."] <- "GENE"

############################################################################
## 4) Prepare the assay matrix
############################################################################

# Convert assay_subset to a data.frame with row names = GENE, then transpose
assay_subset <- data.frame(assay_subset, row.names = 1)
transposed_assay <- t(assay_subset)
transposed_assay <- as.data.frame(transposed_assay)
transposed_assay <- transposed_assay %>% rownames_to_column(var = "acc")

# Merge with metadata to include readcount
merged_assay_metadata <- merge(transposed_assay, 
                               filtered_metadata[, c("acc", "readcount")],
                               by = "acc")

############################################################################
## 5) Prepare gene length and RPK, RPKM calculations
############################################################################

filtered_gene_lengths <- gene_lengths %>%
  filter(GENE %in% GENEs_to_keep)

# Convert gene lengths to kilobases
filtered_gene_lengths$length_kb <- filtered_gene_lengths$gene_length / 1000

# Re-express assay_subset as a data frame with GENE column
normalized_assay <- assay_subset %>%
  rownames_to_column(var = "GENE")

# Merge the gene expression matrix with gene lengths
merged_assay_metadata <- merge(normalized_assay, filtered_gene_lengths,
                               by = "GENE")

# Calculate RPK (Reads Per Kilobase)
rpk_values <- merged_assay_metadata %>%
  mutate(across(-c(GENE, gene_length, length_kb), ~ . / length_kb))

# Keep only columns that are RPK
rpk_values <- rpk_values %>% dplyr::select(-c(GENE, gene_length, length_kb))

# Match order of samples to metadata
identical(filtered_metadata$acc, colnames(rpk_values))
rpk_values <- rpk_values[ , match(filtered_metadata$acc, colnames(rpk_values))]
identical(filtered_metadata$acc, colnames(rpk_values))

# Calculate RPKM
library_sizes <- as.numeric(filtered_metadata$readcount)
rpkm_values   <- sweep(as.matrix(rpk_values), 2, library_sizes, FUN = "/")

# Convert back to data frame and add GENE column
rpkm_dataframe <- as.data.frame(rpkm_values)
rpkm_dataframe <- cbind(GENE = merged_assay_metadata$GENE, rpkm_dataframe)

############################################################################
## 6) Finalize the assay: ensure rownames == GENE
############################################################################

final_assay <- data.frame(rpkm_dataframe)
rownames(final_assay) <- final_assay$GENE
final_assay$GENE <- NULL  # remove the redundant GENE column

############################################################################
## 7) Make 'tax' (rowData) also have rownames == GENE
############################################################################

tax <- as.data.frame(unique_resfinder)
rownames(tax) <- tax$GENE  # ensure rownames match the GENE column

############################################################################
## 8) Subset both to their common genes & create TSE
############################################################################

common.features <- intersect(rownames(final_assay), rownames(tax))
final_assay     <- final_assay[common.features, ]
tax             <- tax[common.features, ]

samples <- filtered_metadata

# Add sums of all ARGs in RPKM to sample data (ARG load)
samples$ARG_load <- colSums(final_assay)

tree_summarized_experiment <- TreeSummarizedExperiment(
  assays  = SimpleList(counts = as.matrix(final_assay)),
  colData = DataFrame(samples),
  rowData = DataFrame(tax)
)

# Save the TSE object
saveRDS(tree_summarized_experiment, file="DATA/TSE.rds")

rm(list  = ls())

############################################################################
## 9) Reload the TSE and filter
############################################################################

TSE <- readRDS("DATA/TSE.rds")

## --- Filtering: remove samples with all-zero or NA counts
col_sums       <- colSums(assay(TSE), na.rm = FALSE)
columns_to_keep <- !(is.na(col_sums) | col_sums == 0)
TSE <- TSE[ , columns_to_keep]

## --- Filtering: keep only features that appear at least in one sample
otus_to_keep <- rowSums(assay(TSE) != 0) > 0
TSE <- TSE[otus_to_keep, ]

## --- Alpha diversities
TSE <- mia::addAlpha(TSE, assay.type = "counts", index = "shannon", name = "ARG_div_shan")
TSE <- mia::addAlpha(TSE, assay.type = "counts", index = "simpson", name = "ARG_div_simp")
TSE <- mia::addAlpha(TSE, assay.type = "counts", index = "observed", name = "ARG_obs")

## --- Transform to relative abundances if needed (for barplots, etc.)
TSE <- transformAssay(TSE, assay.type = "counts", method = "relabundance", MARGIN = "samples")

## --- Further metadata processing
source("R-SCRIPTS/metadata_processing_country.R")

filtered_data <- colData(TSE) %>%
  data.frame() %>%
  filter(readcount > 0.5, ARG_load > 50, ARG_obs > 1)

filtered_samples <- rownames(filtered_data)
TSE <- TSE[ , filtered_samples]

source("R-SCRIPTS/metadata_processing_bioproject.R")

# Optionally add antibiotic use data
# source("R-SCRIPTS/metadata_processing_antibiotic_use.R")

saveRDS(TSE, file = "DATA/TSE.rds")

rm(list  = ls())

############################################################################
## 10) Reload & finalize metadata (e.g., create age groups, etc.)
############################################################################

TSE <- readRDS("DATA/TSE.rds")
colData(TSE) <- colData(TSE) %>%
  as.data.frame() %>%
  mutate(age_category = dplyr::case_when(
    host_age_years >= 0   & host_age_years <= 1  ~ "Infant",
    host_age_years >  1   & host_age_years <= 3  ~ "Toddler",
    host_age_years >  3   & host_age_years <= 12 ~ "Children",
    host_age_years >  12  & host_age_years <  20 ~ "Teenager",
    host_age_years >= 20  & host_age_years <  35 ~ "Young Adult",
    host_age_years >= 35  & host_age_years <  65 ~ "Middle-Aged Adult",
    host_age_years >= 65  & host_age_years <  80 ~ "Older Adult",
    host_age_years >= 80  & host_age_years <= 100 ~ "Oldest Adult",
    TRUE ~ NA_character_
  )) %>%
  mutate(age_category = factor(age_category, levels = c(
    "Infant", "Toddler", "Children", "Teenager",
    "Young Adult", "Middle-Aged Adult", "Older Adult", "Oldest Adult"
  ))) %>%
  DataFrame()

# Save final TSE
saveRDS(TSE, file = "DATA/TSE.rds")

############################################################################
## End
############################################################################

## TODO ##

# ## ---------- Ordinations ------------
