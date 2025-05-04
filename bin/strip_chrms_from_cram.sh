#!/usr/bin/env bash

# Usage: sanitize_core_chroms.sh -i <input.cram> -r <original_ref.fasta> -n <new_ref.fasta> -o <output.cram> [-c core_chromosomes.txt] [-t threads]

set -euo pipefail

THREADS=$(nproc)
CORE_CHROMS=""

usage() {
    echo "Usage: $0 -i <input.cram> -r <original_ref.fasta> -n <new_ref.fasta> -o <output.cram> [-c core_chromosomes.txt] [-t threads]"
    exit 1
}

while getopts ":i:r:n:o:c:t:" opt; do
    case ${opt} in
        i ) INPUT_CRAM=$OPTARG ;;
        r ) ORIGINAL_REF=$OPTARG ;;
        n ) NEW_REF=$OPTARG ;;
        o ) OUTPUT_CRAM=$OPTARG ;;
        c ) CORE_CHROMS=$OPTARG ;;
        t ) THREADS=$OPTARG ;;
        * ) usage ;;
    esac
done

if [[ -z ${INPUT_CRAM+x} || -z ${ORIGINAL_REF+x} || -z ${NEW_REF+x} || -z ${OUTPUT_CRAM+x} ]]; then
    usage
fi

if [[ -z "$CORE_CHROMS" ]]; then
    CORE_CHROMS=$(mktemp)
    grep '^>' "$NEW_REF" | sed 's/>//' > "$CORE_CHROMS"
fi

samtools view -@ "$THREADS" -h -T "$ORIGINAL_REF" "$INPUT_CRAM" | \
awk -v OFS='\t' -v core_chroms_file="$CORE_CHROMS" '
    BEGIN {while (getline < core_chroms_file) keep[$1]=1}
    /^@/ {print; next}
    {
        if (!keep[$3]) {
            $2 = or($2, 4);          # Mark read as unmapped
            $2 = and($2, compl(2));  # Mark mate as unmapped
            $3 = "*";                # Clear chromosome
            $4 = 0;                  # Clear position
            $6 = "*";                # Clear CIGAR
            $7 = "*";                # Clear mate chromosome
            $8 = 0;                  # Clear mate position
            newtags=""
            for (i=12; i<=NF; i++) {
                if ($i !~ /^(XA|SA|MC|NM|MD|AS|XS|RG|MQ|ms|ts):/) 
                    newtags=newtags"\t"$i
            }
            print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11 newtags
        } else {
            print
        }
    }' | \
samtools view -@ "$THREADS" -C -T "$NEW_REF" -o "$OUTPUT_CRAM"

samtools index -@ "$THREADS" "$OUTPUT_CRAM"

echo "Output CRAM written and indexed: $OUTPUT_CRAM"

if [[ -z "$CORE_CHROMS" ]]; then
    rm -f "$CORE_CHROMS"
fi
