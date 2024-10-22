#!/bin/bash

# Ensure compatibility with zsh and bash
if [ -n "$ZSH_VERSION" ]; then
  emulate -LR bash  # Ensure bash-like behavior in zsh
fi

# Display usage information
usage() {
  echo "Usage: $0 -p SAMPLE_MATCH_PATTERN -i INPUT_DIR -o OUTPUT_DIR"
  exit 1
}

# Parse command-line arguments
while getopts ":p:i:o:" opt; do
  case $opt in
    p) SAMPLE_MATCH_PATTERN=$OPTARG ;;
    i) INPUT_DIR=$OPTARG ;;
    o) OUTPUT_DIR=$OPTARG ;;
    *) usage ;;
  esac
done

INPUT_DIR=${INPUT_DIR}/
OUTPUT_DIR=${OUTPUT_DIR}/

# Check if all required arguments are provided
if [[ -z "$SAMPLE_MATCH_PATTERN" || -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]]; then
  usage
fi

# Check if directories exist
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Error: Input directory $INPUT_DIR does not exist."
  exit 1
fi

if [[ ! -d "$OUTPUT_DIR" ]]; then
  echo "Error: Output directory $OUTPUT_DIR does not exist."
  exit 1
fi

echo ""
echo ">>>  >> > THIS IS A SIMPLE UTILITY AND NOT PRODUCTION TESTED FOR YOUR USE CASE < <<  <<<"
echo ""
sleep 2

# Use a shell-agnostic way to store found files in arrays
R1_FILES=($(find "$INPUT_DIR" -name "${SAMPLE_MATCH_PATTERN}*_1.fastq.gz" | sort))
R2_FILES=($(find "$INPUT_DIR" -name "${SAMPLE_MATCH_PATTERN}*_2.fastq.gz" | sort))

# Check if numbers of R1 and R2 files match
if [[ ${#R1_FILES[@]} -ne ${#R2_FILES[@]} ]]; then
  echo "Error: Mismatched number of R1 and R2 files."
  exit 1
fi

# Display found files and their sizes
echo "Found the following R1 and R2 files:"
for ((i = 0; i < ${#R1_FILES[@]}; i++)); do
  R1_FILE="${R1_FILES[$i]}"
  R2_FILE="${R2_FILES[$i]}"
  
  R1_SIZE=$(du -h "$R1_FILE" | cut -f1)
  R2_SIZE=$(du -h "$R2_FILE" | cut -f1)

  echo -e "$R1_FILE\t($R1_SIZE)"
  echo -e "$R2_FILE\t($R2_SIZE)"
  echo "------------------------------------------"

  # Check if files are size 0
  if [[ $R1_SIZE == "0" || $R2_SIZE == "0" ]]; then
    echo "Error: One or more files are of size 0."
    exit 1
  fi
done

# Define output file paths
FQR1="$OUTPUT_DIR/${SAMPLE_MATCH_PATTERN}_R1.fastq.gz"
FQR2="$OUTPUT_DIR/${SAMPLE_MATCH_PATTERN}_R2.fastq.gz"


echo "Concatenating R1 and R2 files..."
echo "Output file R1: $FQR1"
echo "Output file R2: $FQR2"

# Concatenate R1 and R2 files in parallel
cat "${R1_FILES[@]}" > "$FQR1" &
cat "${R2_FILES[@]}" > "$FQR2" &

wait
echo "Done"
