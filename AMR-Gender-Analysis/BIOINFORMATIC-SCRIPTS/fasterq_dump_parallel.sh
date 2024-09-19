#!/bin/bash -l

#SBATCH --job-name=array_job
#SBATCH --output=array_job_out_%A_%a.txt
#SBATCH --error=array_job_err_%A_%a.txt
#SBATCH --account=project_2008149
#SBATCH --partition=small
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --gres=nvme:300
#SBATCH --mem-per-cpu=6000
#SBATCH --array=1-100   # Total number of arrays

module load biokit
module load pigz

# Function to process each batch
process_batch() {
	local batch_file=$1
	local batch_log_dir=$2

	# Function to download accessions and compress using pigz in parallel
	process_accessions() {
		local accession_list=$1

		parallel -j 12 --delay 1 --retries 3 --joblog "$batch_log_dir/process.log" \
			'fasterq-dump --include-technical --split-files -O $LOCAL_SCRATCH {} 2> logs/pe/{}.gz.log && \
			pigz -p 4  $LOCAL_SCRATCH/{}_1.fastq && \
			mv $LOCAL_SCRATCH/{}_1.fastq.gz DATA/{}_1.fastq.gz' :::: "$accession_list"
	}

	# Download accessions and compress using pigz in parallel
	process_accessions "$batch_file"

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
num_batches=100  # Total number of batches

# Calculate batch index for this array
batch_index=$((SLURM_ARRAY_TASK_ID))

# Calculate start and end accession index for this batch
start_index=$(( (batch_index - 1) * batch_size + 1 ))
end_index=$(( start_index + batch_size - 1 ))

# Create a temporary directory for the batch
batch_log_dir="batch_logs_$batch_index"
mkdir "$batch_log_dir" || { echo "Failed to create batch log directory"; exit 1; }

# Create a temporary file for the batch accessions
batch_file="batch_$batch_index.txt"
sed -n "${start_index},${end_index}p" "$accessions_file" > "$batch_file"

# Process the batch
process_batch "$batch_file" "$batch_log_dir"

echo "Finished processing batch $batch_index"

