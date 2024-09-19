#!/bin/bash -l

#SBATCH --job-name=array_job
#SBATCH --output=array_job_out_%A_%a.txt
#SBATCH --error=array_job_err_%A_%a.txt
#SBATCH --account=project_2008149
#SBATCH --partition=small
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3000
#SBATCH --array=1-30   # Total number of arrays

module load biokit
module load pigz
module load parallel
conda activate kingfisher

process_batch() {
    local batch_file=$1
    local batch_log_dir=$2
    local completed_list=$3

    # Define a function to process each accession
process_accession() {
    local accession=$1
    if grep -Fxq "$accession" "$completed_list"; then
        echo "Accession $accession already processed. Skipping."
    else
        if fasterq-dump --include-technical --split-files -O DATA "$accession" 2> "logs/pe/${accession}.gz.log"; then  
            pigz -p 4 "DATA/${accession}_1.fastq"
            if [ -e "DATA/${accession}_2.fastq" ]; then
                rm -f "DATA/${accession}_2.fastq"
                if [ $? -eq 0 ]; then
                    echo "Removed ${accession}_2.fastq successfully."
                else
                    echo "Failed to remove ${accession}_2.fastq."
                fi
            fi
            echo "$accession" >> "$completed_list"
        else
            echo "fasterq-dump failed for accession $accession. Trying kingfisher."
            cd DATA/q
            kingfisher get -r "$accession" -m ena-ascp ena-ftp prefetch aws-http aws-cp &> "../logs/${accession}_download.log"
            pigz -p 4 "${accession}_1.fastq"
	    if [ -e "${accession}_2.fastq" ]; then
                rm -f "${accession}_2.fastq"
                if [ $? -eq 0 ]; then
                    echo "Removed ${accession}_2.fastq successfully."
                else
                    echo "Failed to remove ${accession}_2.fastq"
                fi
	fi   	
		if [ -e "${accession}_2.fastq.gz" ]; then
		rm -f "${accession}_2.fastq"
	fi
            cd ..
            echo "$accession" >> "$completed_list"
        fi
    fi
}

    # Export the function so it's available to GNU Parallel
    export -f process_accession
    export completed_list
    export batch_file
    export batch_log_dir


    # Run process_accession function in parallel for each accession
    mkdir /scratch/project_2008149/PROJECT_WORKSPACE/workflow/$SLURM_JOBID
    export TMPDIR=/scratch/project_2008149/PROJECT_WORKSPACE/workflow/$SLURM_JOBID
    cat "$batch_file" | parallel -j 10 process_accession
    rm -rf /scratch/project_2008149/PROJECT_WORKSPACE/workflow/$SLURM_JOBID 

    # Cleanup the batch log directory
    rmdir "$batch_log_dir"
}

# Parse command-line options
while getopts "a:" opt; do
	case $opt in
		a) accessions_file="$OPTARG";;
		*) echo "Usage: $0 -a <accessions_file>"; exit 1;;
	esac
done

# Check if accessions file is provided
if [ -z "$accessions_file" ]; then
	echo "Error: Accessions file not provided."
	exit 1
fi

# Calculate batch size
batch_size=200
num_batches=200  # Total number of batches

# Calculate batch index for this array
batch_index=$((SLURM_ARRAY_TASK_ID))

# Calculate start and end accession index for this batch
start_index=$(( (batch_index - 1) * batch_size + 1 ))
end_index=$(( start_index + batch_size - 1 ))

# Create a temporary directory for the batch
batch_log_dir="/scratch/project_2008149/PROJECT_WORKSPACE/workflow/batch_logs_$batch_index"
mkdir "$batch_log_dir" || { echo "Failed to create batch log directory"; exit 1; }

# Create a temporary file for the batch accessions
batch_file="/scratch/project_2008149/PROJECT_WORKSPACE/workflow/batch_$batch_index.txt"
sed -n "${start_index},${end_index}p" "$accessions_file" > "$batch_file"

# Define the completed list file
completed_list="/scratch/project_2008149/PROJECT_WORKSPACE/workflow/completed_list.txt"

# Process the batch
process_batch "$batch_file" "$batch_log_dir" "$completed_list"

echo "Finished processing batch $batch_index"

