#!/bin/bash

# Define directories and output file paths
SAMPLE=$1
DATA_DIR=$2
OUTPUT_DIR=$3

# Create output files for WES tumor, WES normal, and RNA tumor
WES_TUMOR_R1="$OUTPUT_DIR/WES-tumor-$SAMPLE.R1.fastq.gz"
WES_TUMOR_R2="$OUTPUT_DIR/WES-tumor-$SAMPLE.R2.fastq.gz"
WES_NORMAL_R1="$OUTPUT_DIR/WES-normal-$SAMPLE.R1.fastq.gz"
WES_NORMAL_R2="$OUTPUT_DIR/WES-normal-$SAMPLE.R2.fastq.gz"
RNA_TUMOR_R1="$OUTPUT_DIR/RNA-tumor-$SAMPLE.R1.fastq.gz  "
RNA_TUMOR_R2="$OUTPUT_DIR/RNA-tumor-$SAMPLE.R2.fastq.gz"

# Function to concatenate files based on the pattern

echo "Concatenating $pattern R1 and R2 files..."


# WES TUMOR
#ls -lt /fsx/data/tmp_inputs/lr1/2024-10-14-BostonGene/WES-tumor-A01231_3A917_3AHGGMFDSXC_*_1.fastq.gz
#ls -lt /fsx/data/tmp_inputs/lr1/2024-10-14-BostonGene/WES-tumor-A01231_3A917_3AHGGMFDSXC_*_2.fastq.gz
# Find and concatenate R1 files
cat /fsx/data/tmp_inputs/lr1/2024-10-14-BostonGene/WES-tumor-A01231_3A917_3AHGGMFDSXC_*_1.fastq.gz > $WES_TUMOR_R1 &
cat /fsx/data/tmp_inputs/lr1/2024-10-14-BostonGene/WES-tumor-A01231_3A917_3AHGGMFDSXC_*_2.fastq.gz > $WES_TUMOR_R2 &

# WES Normal R1 and R2
#ls -lt /fsx/data/tmp_inputs/lr1/2024-10-14-BostonGene/WES-normal-A01599_3A214_3AHCMJKDRX2_*_1.fastq.gz
#ls -lt /fsx/data/tmp_inputs/lr1/2024-10-14-BostonGene/WES-normal-A01599_3A214_3AHCMJKDRX2_*_2.fastq.gz

cat /fsx/data/tmp_inputs/lr1/2024-10-14-BostonGene/WES-normal-A01599_3A214_3AHCMJKDRX2_*_1.fastq.gz > $WES_NORMAL_R1 &
cat /fsx/data/tmp_inputs/lr1/2024-10-14-BostonGene/WES-normal-A01599_3A214_3AHCMJKDRX2_*_2.fastq.gz > $WES_NORMAL_R2 &


# RNA Tumor R1 and R2
# ls -lt /fsx/data/tmp_inputs/lr1/2024-10-14-BostonGene/RNASeq-tumor-*_1.fastq.gz
# ls -lt /fsx/data/tmp_inputs/lr1/2024-10-14-BostonGene/RNASeq-tumor-*_2.fastq.gz

cat /fsx/data/tmp_inputs/lr1/2024-10-14-BostonGene/RNASeq-tumor-*_1.fastq.gz > $RNA_TUMOR_R1 &
cat /fsx/data/tmp_inputs/lr1/2024-10-14-BostonGene/RNASeq-tumor-*_2.fastq.gz > $RNA_TUMOR_R2 &
  
wait
echo "Done"


