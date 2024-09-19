
# Read in csv with abstract and classification
# Created using the BIOINFORMATIC_SCRIPTS/fetch_abstracts.py
df <- read_csv("DATA/classified_Sra_metadata_jun11.txt")
df$category %>% table
