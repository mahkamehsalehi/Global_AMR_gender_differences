library(tidyquant)
library(tidyverse)
library(jsonlite)
library(data.table)
library(progress)

source("../R-SCRIPTS/00_functions.R")

# Do not run, slow
# Downlaods Sra metadata file that also has the attributes field and processes it.

df <- read_csv("../DATA/Sra_metadata_jun12_attributes.txt")

df <- df %>%
  mutate(
    # Parse the JSON column and convert it to a list of data frames
    jattr_parsed = map(jattr, ~ fromJSON(.x, flatten = TRUE))
  ) %>%
  # Unnest the list column to expand it into multiple columns with a prefix to avoid name clashes
  unnest_wider(jattr_parsed, names_sep = "_jattr")

#write_csv(df, "../DATA/PROCESSED/Sra_metadata_processed.csv")
