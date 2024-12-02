library(TreeSummarizedExperiment)
library(SummarizedExperiment)

tse <- readRDS("DATA/TSE_filtered.rds")
# Extract the main components from TSE
assay_data <- assays(tse)
row_data <- rowData(tse)
col_data <- colData(tse)

# Create a SummarizedExperiment object
se <- SummarizedExperiment(
  assays = assay_data,
  rowData = row_data,
  colData = col_data
)

saveRDS(se, file = "DATA/summarized_experiment.rds")
