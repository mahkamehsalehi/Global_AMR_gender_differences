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
