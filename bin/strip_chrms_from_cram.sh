#!/usr/bin/env bash

# Optimized sanitize_core_chroms.sh
set -euo pipefail

THREADS=$(nproc)
TMPDIR="/fsx/scratch/sanitize_cram_$(date +%s)"

usage() {
    echo "Usage: $0 -i <input.cram> -n <new_ref.fasta> -o <output.cram> [-t threads] [-x tmpdir]"
    exit 1
}

while getopts ":i:n:o:t:x:" opt; do
    case ${opt} in
        i ) INPUT_CRAM=$OPTARG ;;
        n ) NEW_REF=$OPTARG ;;
        o ) OUTPUT_CRAM=$OPTARG ;;
        t ) THREADS=$OPTARG ;;
        x ) TMPDIR=$OPTARG ;;
        * ) usage ;;
    esac
done

[[ -z ${INPUT_CRAM+x} || -z ${NEW_REF+x} || -z ${OUTPUT_CRAM+x} ]] && usage

mkdir -p "$TMPDIR"

# Extract allowed chromosomes from new reference
ALLOWED_CHRMS=$(grep '^>' "$NEW_REF" | sed 's/>//' | cut -f1 -d' ')

# Efficiently filter CRAM using samtools
samtools view -@ "$THREADS" -C \
  -T "$NEW_REF" \
  -o "$OUTPUT_CRAM" \
  --write-index \
  "$INPUT_CRAM" $ALLOWED_CHRMS

echo "✅ Output CRAM written and indexed: $OUTPUT_CRAM"
