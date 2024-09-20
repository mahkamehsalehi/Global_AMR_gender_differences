import pandas as pd
from Bio import Entrez
import re

# Set your email for NCBI Entrez API
Entrez.email = "your_email@example.com"

# Function to fetch bioproject info from NCBI
def fetch_bioproject_info(bioproject_number):
    try:
        print(f"Fetching abstract for bioproject: {bioproject_number}")
        handle = Entrez.esearch(db="bioproject", term=bioproject_number)
        record = Entrez.read(handle)
        if record["IdList"]:
            bioproject_id = record["IdList"][0]
            handle = Entrez.esummary(db="bioproject", id=bioproject_id)
            summary = Entrez.read(handle)
            return summary
    except Exception as e:
        print(f"Error fetching bioproject {bioproject_number}: {e}")
    return None

# Load your metadata file
print("Loading metadata file...")
metadata_df = pd.read_csv("../DATA/Sra_metadata_jun11.txt", sep=',')

# Get unique bioprojects
print("Extracting unique bioprojects...")
unique_bioprojects = metadata_df['bioproject'].unique()

# Fetch abstracts for unique bioprojects
print("Fetching abstracts for all unique bioprojects...")
abstracts = {}
for bioproject in unique_bioprojects:
    abstract = fetch_bioproject_info(bioproject)
    if abstract:
        abstracts[bioproject] = abstract
    else:
        abstracts[bioproject] = None

# Convert abstracts dictionary to DataFrame
print("Converting fetched abstracts to DataFrame...")
abstracts_df = pd.DataFrame(list(abstracts.items()), columns=['bioproject', 'abstract'])

# **Data cleaning and checking**

# Check for missing or duplicated bioprojects
print("Checking for missing or duplicated bioprojects...")
print("Missing in metadata_df:", metadata_df['bioproject'].isnull().sum())
print("Missing in abstracts_df:", abstracts_df['bioproject'].isnull().sum())
print("Duplicates in metadata_df:", metadata_df['bioproject'].duplicated().sum())
print("Duplicates in abstracts_df:", abstracts_df['bioproject'].duplicated().sum())

# Ensure bioproject columns have consistent data types and no extra spaces
print("Trimming spaces and converting to strings in both DataFrames...")
metadata_df['bioproject'] = metadata_df['bioproject'].str.strip().astype(str)
abstracts_df['bioproject'] = abstracts_df['bioproject'].str.strip().astype(str)

# Merging the abstracts back with the original metadata
print("Merging abstracts with the original metadata...")
metadata_df = metadata_df.drop_duplicates(subset=['bioproject'])  # Dropping duplicates
metadata_df = metadata_df.merge(abstracts_df, on='bioproject', how='left')

# Classification function
import re

def classify_bioproject(abstract):
    # Check if abstract is None or has an unexpected format
    if abstract and isinstance(abstract, dict) and 'DocumentSummarySet' in abstract:
        summary_set = abstract['DocumentSummarySet']
        if 'DocumentSummary' in summary_set:
            document_summary = summary_set['DocumentSummary']
            if len(document_summary) > 0:
                # Update to use 'Project_Title' and 'Project_Description'
                title = document_summary[0].get('Project_Title', "")
                description = document_summary[0].get('Project_Description', "")
                abstract_text = f"{title} {description}".lower()

                # Print abstract_text for debugging
                print(f"Debug: Abstract Text: {abstract_text}")

                # Infant Study
                if re.search(r'infants?|first year of life|neonates?|newborns?|baby|babies|preterms?', abstract_text):
                    return 'Infant Study'
                # Infection-related keywords (including CDI, Clostridioides difficile, FMT with infection)
                elif re.search(r'infections?|covid|pathogens?|virus|cre|clostridium difficile|clostridioides difficile|vibrio cholerae|salmonella|diarrhea|covid19??|fmt with infection|icu patient|campylobacter|cdi', abstract_text):
                    return 'Infection'
                # Cancer-related keywords
                elif re.search(r'cancer|crc|tumor|oncology|colorectal cancer|breast cancer|neuroblastoma|leukemia|melanoma|pancreatic cancer', abstract_text):
                    return 'Cancer'
                # Immune Deficiency-related keywords
                elif re.search(r'immune deficiency|immunodeficiency|immunocompromised|hematopoietic cell transplantation|marrow transplant|stem cell transplant', abstract_text):
                    return 'Immune Deficiency'
                # Other diseases
                elif re.search(r'fmt(?! with infection)|microbiota transplant|inflammatory bowel disease|ibd|chaple disease|short bowel syndrome|hepatic encephalopathy|crohn\'s disease|ulcerative colitis|colitis', abstract_text):
                    return 'Other Disease'
                else:
                    return 'Other'
    # If abstract is None or in an unexpected format, return 'No Data'
    return 'No Data'

# Apply classification
print("Classifying bioprojects based on abstracts...")
metadata_df["category"] = metadata_df["abstract"].apply(classify_bioproject)

# Save the classified data to a new file
print("Saving classified data to 'classified_Sra_metadata_jun11.txt'...")
metadata_df.to_csv('../DATA/classified_Sra_metadata_jun11.txt', sep=',', index=False)

print("Processing complete. Here are the first few entries:")
print(metadata_df[['bioproject', 'category']].head())

