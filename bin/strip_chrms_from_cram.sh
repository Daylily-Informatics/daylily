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

# Create allowed chromosomes list if not provided
if [[ -z "$CORE_CHROMS" ]]; then
    CORE_CHROMS="${TMPDIR}/allowed_chrms.txt"
    grep '^>' "$NEW_REF" | sed 's/>//' | cut -f1 -d' ' > "$CORE_CHROMS"
fi

# Create filtered header
samtools view -H "$INPUT_CRAM" | \
    awk -v keep="$CORE_CHROMS" '
        BEGIN {while(getline < keep) allowed[$1]=1}
        /^@SQ/ {split($2,a,":"); if(a[2] in allowed) print; next}
        !/^@SQ/ {print}
    ' > "${TMPDIR}/filtered_header.sam"

# Reheader CRAM first (with header filtered)
samtools reheader "${TMPDIR}/filtered_header.sam" "$INPUT_CRAM" > "${TMPDIR}/reheadered.cram"
samtools index -@ "$THREADS" "${TMPDIR}/reheadered.cram"

# Now correctly filter by chromosomes and remove XA/SA tags
samtools view \
    -@ "$THREADS" \
    -C \
    -T "$NEW_REF" \
    --remove-tag XA \
    --remove-tag SA \
    --write-index \
    -o "$OUTPUT_CRAM" \
    "${TMPDIR}/reheadered.cram" \
    $(cat "$CORE_CHROMS")

# Cleanup intermediate files (optional)

echo "✅ Output CRAM written and indexed: $OUTPUT_CRAM"
