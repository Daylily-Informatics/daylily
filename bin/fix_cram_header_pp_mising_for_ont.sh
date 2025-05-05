#!/usr/bin/env bash

# Usage: replace_pp_with_dummy.sh -i <input.cram> -r <reference.fasta> -o <output.cram> [-t threads]

set -euo pipefail

THREADS=$(nproc)

usage() {
    echo "Usage: $0 -i <input.cram> -r <reference.fasta> -o <output.cram> [-t threads]"
    exit 1
}

while getopts ":i:r:o:t:" opt; do
    case ${opt} in
        i ) INPUT_CRAM=$OPTARG ;;
        r ) REFERENCE_FASTA=$OPTARG ;;
        o ) OUTPUT_CRAM=$OPTARG ;;
        t ) THREADS=$OPTARG ;;
        * ) usage ;;
    esac
done

if [[ -z ${INPUT_CRAM+x} || -z ${REFERENCE_FASTA+x} || -z ${OUTPUT_CRAM+x} ]]; then
    usage
fi

# Detect missing PP program ID automatically from warnings
MISSING_PP=$(samtools view -@ "$THREADS" -T "$REFERENCE_FASTA" "$INPUT_CRAM" 2>&1 >/dev/null | \
    grep "has a PP link to missing program" | \
    head -n1 | sed -E "s/.*PP link to missing program '([^']+)'.*/\1/")

if [[ -z "$MISSING_PP" ]]; then
    echo "No missing PP field detected. Exiting."
    exit 0
else
    echo "Missing PP detected: '$MISSING_PP'"
fi

# Generate a new header with the missing PG line added
samtools view -@ "$THREADS" -H "$INPUT_CRAM" | \
awk -v dummy_id="$MISSING_PP" '
    BEGIN {dummy_added=0}
    /^@PG/ && dummy_added==0 {
        printf("@PG\tID:%s\tPN:%s\tVN:unknown\tCL:\"%s (dummy entry)\"\n", dummy_id, dummy_id, dummy_id)
        dummy_added=1
    }
    {print}
' > tmp_fixed_header.sam

# Reheader and rewrite CRAM using corrected header
samtools reheader -P tmp_fixed_header.sam "$INPUT_CRAM" | \
    samtools view -@ "$THREADS" -C -T "$REFERENCE_FASTA" -o "$OUTPUT_CRAM"

# Index the fixed CRAM
samtools index -@ "$THREADS" "$OUTPUT_CRAM"

# Clean temporary header
rm tmp_fixed_header.sam

echo "✅ Output CRAM written and indexed with missing PP field '$MISSING_PP' resolved: $OUTPUT_CRAM"
