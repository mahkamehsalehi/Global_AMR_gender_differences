
# Read in csv with abstract and classification
# Created using the BIOINFORMATIC_SCRIPTS/fetch_abstracts.py
df <- read_csv("DATA/classified_Sra_metadata_jun11.txt")

library(stringr)
library(dplyr)

# Step 1: Apply the extraction for every row in the dataframe
df <- df %>%
  rowwise() %>%  # Apply the operation row by row
  mutate(
    title = str_match(abstract, "'Project_Title':\\s*'(.*?)'")[,2],  # Extract Project_Title
    abstract = str_match(abstract, "'Project_Description':\\s*'(.*?)'")[,2]  # Extract Project_Description
  ) %>%
  ungroup() %>%  # Ungroup the row-wise operation
  dplyr::select(bioproject, title, abstract, category)  # Keep only 'title' and 'abstract'

# Step 2: View the updated dataframe
df %>% head(10)
df$category %>% table

# Step 3: Extract the current colData from the TSE object
coldata_tse <- colData(TSE)

# Convert colData to a dataframe (if needed)
coldata_tse <- as.data.frame(coldata_tse)

# Step 2: Merge df (with title and abstract) and coldata_tse using 'bioproject' as the common column
merged_coldata <- coldata_tse %>%
  dplyr::left_join(df, by = "bioproject")  # Use left_join to merge on 'bioproject'

# Step 3: Assign the merged dataframe back to colData
colData(TSE) <- DataFrame(merged_coldata)
head(colData(TSE))
