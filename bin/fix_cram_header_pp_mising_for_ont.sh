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

samtools view -@ "$THREADS" -H "$INPUT_CRAM" | sed 's/PP:[^\t]*/PP:dummy/' > tmp_header.sam

samtools reheader -P tmp_header.sam "$INPUT_CRAM" | \
    samtools view -@ "$THREADS" -C -T "$REFERENCE_FASTA" -o "$OUTPUT_CRAM"

samtools index -@ "$THREADS" "$OUTPUT_CRAM"

rm tmp_header.sam

echo "Output CRAM with dummy PP field written and indexed: $OUTPUT_CRAM"
