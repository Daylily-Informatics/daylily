#!/usr/bin/env bash
set -euo pipefail

THREADS=$(nproc)
TMPDIR="/fsx/scratch/sanitize_cram_$(date +%s)"

usage() {
    echo "Usage: $0 -i <input.cram> -n <new_ref.fasta> -o <output.cram> [-c core_chromosomes.txt] [-t threads] [-x tmpdir]"
    exit 1
}

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

if [[ -z ${INPUT_CRAM+x} || -z ${NEW_REF+x} || -z ${OUTPUT_CRAM+x} ]]; then
    usage
fi

mkdir -p "$TMPDIR"

# Step 1: Create allowed chromosomes list
if [[ -z "$CORE_CHROMS" ]]; then
    CORE_CHROMS="${TMPDIR}/allowed_chrms.txt"
    grep '^>' "$NEW_REF" | sed 's/>//' | cut -f1 -d' ' > "$CORE_CHROMS"
fi

# Step 2: Create filtered header (only allowed chromosomes)
samtools view -H "$INPUT_CRAM" | \
    awk -v keep="$CORE_CHROMS" '
        BEGIN {while(getline < keep) allowed[$1]=1}
        /^@SQ/ {split($2,a,":"); if(a[2] in allowed) print; next}
        !/^@SQ/ {print}
    ' > "${TMPDIR}/filtered_header.sam"

# Step 3: Filter alignments BEFORE reheadering and indexing (IMPORTANT!)
samtools view \
    -@ "$THREADS" \
    -T "$NEW_REF" \
    --remove-tag XA \
    --remove-tag SA \
    -C \
    -o "${TMPDIR}/filtered.cram" \
    "$INPUT_CRAM" \
    $(cat "$CORE_CHROMS")

# Step 4: Apply new header to filtered CRAM (safe now that reads match header)
samtools reheader "${TMPDIR}/filtered_header.sam" "${TMPDIR}/filtered.cram" > "$OUTPUT_CRAM"

sleep 62
# Step 5: Index your fully sanitized output
samtools index -@ "$THREADS" "$OUTPUT_CRAM"

# Cleanup temporary directory (optional but recommended)
rm ${TMPDIR}/filtered_header.sam ${TMPDIR}/filtered.cram

echo "✅ Output CRAM written and indexed: $OUTPUT_CRAM"
