#!/bin/env bash

# requires samtools >= XXXXX and isa-l (for igzip)
# conda create -y -n SAMTOOLS -c conda-forge samtools>=1.20 isa-l>=2.31

# From a shell where you `conda activate SAMTOOLS` run this script

# Check if the correct number of arguments are provided
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <input.bam/cram> <output_R1.fastq.gz> <output_R2.fastq.gz> <output_singletons.fastq.gz> <output_others.fastq.gz> <threads> <reference.fasta|required for .cram input>"
    exit 1
fi


echo "THIS SCRIPT IS NOT YET MODIFIED TO ADD BACK THE information in the read names following the ' '[1|2].*, which has been saved in a metadata file next to the BAM."
sleep 2

# Input BAM/CRAM file
input_file="$1"

# Output FASTQ files
output_r1="$2"
output_r2="$3"
output_singletons="$4"
output_others="$5"
threads="$6"
reference="$7"

# Function to check if the directory exists
check_dir_exists() {
    local dir
    dir=$(dirname "$1")
    if [ ! -d "$dir" ]; then
        echo "Error: Output directory '$dir' does not exist."
        exit 1
    fi
}

# Function to check if the file already exists
check_file_exists() {
    if [ -f "$1" ]; then
        echo "Error: Output file '$1' already exists."
        exit 1
    fi
}

# Check if input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' does not exist."
    exit 1
fi

# Check if the input file is a CRAM file
if [[ "$input_file" == *.cram ]]; then
    # Check if reference is provided
    if [ -z "$reference" ]; then
        echo "Error: Reference file must be provided for CRAM input."
        exit 1
    fi

    # Check if reference file exists
    if [ ! -f "$reference" ]; then
        echo "Error: Reference file '$reference' does not exist."
        exit 1
    fi
fi

# Check if output directories exist
check_dir_exists "$output_r1"
check_dir_exists "$output_r2"
check_dir_exists "$output_singletons"
check_dir_exists "$output_others"

# Check if output files already exist
check_file_exists "$output_r1"
check_file_exists "$output_r2"
check_file_exists "$output_singletons"
check_file_exists "$output_others"

# Temporary uncompressed FASTQ files
tmp_r1="${output_r1%.gz}"
tmp_r2="${output_r2%.gz}"
tmp_singletons="${output_singletons%.fastq.gz}_singletons.fastq"
tmp_others="${output_others%.fastq.gz}_others.fastq"

# Convert BAM/CRAM to FASTQ
# -c : Include only properly paired reads
# -f : Include only the first read in a pair
# -F : Include only the second read in a pair
samtools fastq  \
    -@ "$threads" \
    -t \
    -c 1 \
    -1 "$tmp_r1" \
    -2 "$tmp_r2" \
    -s "$tmp_singletons" \
    -0 "$tmp_others" \
    -n "$input_file"

# Compress the FASTQ files with igzip
igzip -c "$tmp_r1" > "$output_r1"
igzip -c "$tmp_r2" > "$output_r2"
igzip -c "$tmp_singletons" > "$output_singletons"
igzip -c "$tmp_others" > "$output_others"

# Clean up temporary files
rm "$tmp_r1" "$tmp_r2" "$tmp_singletons" "$tmp_others"

echo "Conversion complete. Output files: $output_r1, $output_r2, $output_singletons, and $output_others"
