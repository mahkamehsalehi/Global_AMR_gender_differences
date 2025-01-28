#!/bin/bash

# Display help message
function display_help() {
    echo "Usage: $0 -i <accession_file> -o <output_file>"
    echo "Options:"
    echo "  -i    Input file containing SRA run accessions"
    echo "  -o    Output file to save the query results in CSV format"
    echo "  -h    Display this help message"
    exit 1
}

# Check if any arguments are provided
if [ "$#" -eq 0 ]; then
    display_help
fi

# Parse command-line options
while getopts ":i:o:h" opt; do
    case $opt in
        i)
            accession_file="$OPTARG"
            ;;
        o)
            output_file="$OPTARG"
            ;;
        h)
            display_help
            ;;
        \?)
            echo "Invalid option: -$OPTARG"
            display_help
            ;;
        :)
            echo "Option -$OPTARG requires an argument."
            display_help
            ;;
    esac
done

# Check if required options are provided
if [ -z "$accession_file" ] || [ -z "$output_file" ]; then
    echo "Error: Both input (-i) and output (-o) options are required."
    display_help
fi

# Check if the accession file exists
if [ ! -f "$accession_file" ]; then
    echo "Error: File '$accession_file' not found."
    exit 1
fi

# Check if the output file already exists
if [ -f "$output_file" ]; then
    echo "Error: Output file '$output_file' already exists. Choose a different file name."
    exit 1
fi

# Read the content of the accession file into a comma-separated string with single quotes
accessions=$(sed -e "s/^/'/" -e "s/$/'/" "$accession_file" | paste -sd, -)

# Create a temporary file to store the SQL query
query_file=$(mktemp)
query="SELECT
  acc,
  geo_loc_name_country_calc,
  geo_loc_name_country_continent_calc,
  platform,
  instrument,
  bioproject,
  avgspotlen,
  mbases,
  collection_date_sam
FROM
  \`nih-sra-datastore.sra.metadata\` as s
WHERE
acc IN ($accessions) and ('host_sex_sam' in UNNEST(s.attributes) and 'host_age_sam' in UNNEST(s.attributes));"
echo "$query" > "$query_file"

# Run the query using the BigQuery CLI and export results to CSV
bq query --nouse_legacy_sql --format=csv --max_rows=10 < "$query_file" > "$output_file"

echo "Query results saved to $output_file"

# Clean up: Remove the temporary query file
rm "$query_file"

