#!/bin/bash
# Load conda
source ~/miniconda3/etc/profile.d/conda.sh

# Activate your conda environment
conda activate blast

# Check if the necessary environment variables are set
if [ -z "$QUERY_DIR" ]; then
    echo "Environment variable QUERY_DIR is not set."
    exit 1
fi

if [ -z "$BLAST_DB" ]; then
    echo "Environment variable BLAST_DB is not set."
    exit 1
fi

# Set perc_identity from the environment variable, default to 100 if not set
perc_identity="${PERC_IDENTITY:-100}"

# Set query directory and BLAST database from environment variables
query_dir="$QUERY_DIR"
blast_db="$BLAST_DB"
output_dir="$OUTPUT_DIR"

# Check if query directory exists
if [ ! -d "$query_dir" ]; then
    echo "Query directory does not exist: $query_dir"
    exit 1
fi

# Iterate over each query file in the directory
for query_file in "$query_dir"/*.fasta; do
    # Extract the filename without extension
    filename=$(basename -- "$query_file")
    filename_no_ext="${filename%.*}"

    # Define output file name
    output_file="$output_dir/${filename_no_ext}_blastn_out"

    # Run BLAST command for each query file
    blastn -query "$query_file" -db "$blast_db" -out "$output_file" -evalue 1e-5 -num_threads $SLURM_CPUS_PER_TASK -outfmt 6 -perc_identity "$perc_identity" -max_target_seqs 5000
done

