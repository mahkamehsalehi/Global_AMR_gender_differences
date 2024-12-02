library(TreeSummarizedExperiment)
library(SummarizedExperiment)
library(ape)


tse <- readRDS("DATA/TSE.rds")

tree <- rtree(n = nrow(tse))
tree$tip.label <- rownames(tse)

TSE_tree <- TreeSummarizedExperiment(
  assays = assays(tse),
  rowData = rowData(tse),
  colData = colData(tse),
  rowTree = tree
)

saveRDS(TSE_tree, "TSE.rds")


tse_filtered <- readRDS("DATA/TSE_filtered.rds")

tree_filtered <- rtree(n = nrow(tse_filtered))
tree_filtered$tip.label <- rownames(tse_filtered)

TSE_filtered_tree_ <- TreeSummarizedExperiment(
  assays = assays(tse_filtered),
  rowData = rowData(tse_filtered),
  colData = colData(tse_filtered),
  rowTree = tree_filtered
)

saveRDS(TSE_tree, "TSE_filtered.rds")
