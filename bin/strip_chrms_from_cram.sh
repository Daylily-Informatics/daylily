#!/usr/bin/env bash
set -euo pipefail

THREADS=$(nproc)
TMPDIR="/fsx/scratch/sanitize_cram_$(date +%s)"

usage() {
    echo "Usage: $0 -i <input.cram> -n <new_ref.fasta> -o <output.cram> [-c core_chromosomes.txt] [-t threads] [-x tmpdir]"
    exit 1
}

# Parse command-line options
CORE_CHROMS=""
while getopts ":i:n:o:c:t:x:" opt; do
    case ${opt} in
        i ) INPUT_CRAM=$OPTARG ;;
        n ) NEW_REF=$OPTARG ;;
        o ) OUTPUT_CRAM=$OPTARG ;;
        c ) CORE_CHROMS=$OPTARG ;;
        t ) THREADS=$OPTARG ;;
        x ) TMPDIR=$OPTARG ;;
        * ) usage ;;
    esac
done

# Verify required arguments
if [[ -z ${INPUT_CRAM+x} || -z ${NEW_REF+x} || -z ${OUTPUT_CRAM+x} ]]; then
    usage
fi

# Create temporary directory
mkdir -p "$TMPDIR"

# Create allowed chromosomes list if not provided
if [[ -z "$CORE_CHROMS" ]]; then
    CORE_CHROMS="${TMPDIR}/allowed_chrms.txt"
    grep '^>' "$NEW_REF" | sed 's/>//' | cut -f1 -d' ' > "$CORE_CHROMS"
fi

# Generate filtered CRAM file, removing unwanted tags and contigs
samtools view \
    -@ "$THREADS" \
    -C \
    -T "$NEW_REF" \
    --remove-tag XA \
    --remove-tag SA \
    --write-index \
    -o "$OUTPUT_CRAM" \
    "$INPUT_CRAM" \
    $(cat "$CORE_CHROMS")

echo "✅ Output CRAM written and indexed: $OUTPUT_CRAM"
