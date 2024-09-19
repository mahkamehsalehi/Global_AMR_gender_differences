# Display help message
function display_help() {
    echo "Usage: $0 -i <accession_file> -o <output_file> [-r <organism>] [-librarysource <librarysource>] [-assay_type <assay_type>] [-selection <selection>] [-h]"
    echo "Options:"
    echo "  -i    Input file containing SRA run accessions"
    echo "  -o    Output file to save the query results in CSV format"
    echo "  -r      Filter by organism (e.g., 'human gut metagenome')"
    echo "  -librarysource Filter by library source (e.g., 'METAGENOMIC')"
    echo "  -assay_type    Filter by assay type (e.g., 'WGS')"
    echo "  -selection     Filter by selection (e.g., 'RANDOM')"
    echo "  -h    Display this help message"
    echo " For help with query option run bq query --nouse_legacy_sql 'SELECT DISTINCT librarysource, assay_type, organism FROM nih-sra-datastore.sra.metadata'" 
    exit 1
}

# Check if any arguments are provided
if [ "$#" -eq 0 ]; then
    display_help
fi

# Parse command-line options
while getopts ":-:i:o:r:l:a:s:h" opt; do
    case $opt in
        i)
            accession_file="$OPTARG"
            ;;
        o)
            output_file="$OPTARG"
            ;;
        r)
            organism="$OPTARG"
	    ;;
	l)
            librarysource="$OPTARG"
            ;;
        a)
            assay_type="$OPTARG"
            ;;
        s)
            selection="$OPTARG"
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
# if [ -z "$output_file" ] || { [ -z "$accession_file" ] && [ -z "$organism" ] && [ -z "$librarysource" ] && [ -z "$assay_type" ] && [ -z "$selection" ]; }; then
#    echo "Error: Output file (-o) and at least one of the other options (-i, -r, -l, -a, -s) are required."
#    display_help
#fi

# If input file is provided, check if it exists
if [ -n "$accession_file" ] && [ ! -f "$accession_file" ]; then
    echo "Error: File '$accession_file' not found."
    exit 1
fi

# Check if the output file already exists
if [ -f "$output_file" ]; then
    echo "Error: Output file '$output_file' already exists. Choose a different file name."
    exit 1
fi
# Read the content of the accession file into a comma-separated string with single quotes
if [ -f "$accession_file"]; then
accessions=$(sed -e "s/^/'/" -e "s/$/'/" "$accession_file" | paste -sd, -)
fi

# Create a temporary file to store the SQL query
query_file=$(mktemp)

# Build the WHERE clause based on the specified filters
where_clause=""
if [ -n "$organism" ]; then
    where_clause+=" organism = '$organism'"
fi
if [ -n "$librarysource" ]; then
    where_clause+=" AND librarysource = '$librarysource'"
fi
if [ -n "$assay_type" ]; then
    where_clause+=" AND assay_type = '$assay_type'"
fi
if [ -n "$selection" ]; then
    where_clause+=" AND selection = '$libraryselection'"
fi

# Construct the SQL query
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
  \`nih-sra-datastore.sra.metadata\`
WHERE
$where_clause;"

echo "$query" > "$query_file"

# Run the query using the BigQuery CLI and export results to CSV
bq --format=csv query --nouse_legacy_sql --max_rows=10  < "$query_file" > "$output_file"

echo "Query results saved to $output_file"

# Clean up: Remove the temporary query file
#rm "$query_file"
