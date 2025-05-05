#!/usr/bin/env bash

set -euo pipefail

usage() {
  echo "Usage: $0 -i <input.cram> -o <output.cram> -r <readgroup-id> [-t threads]"
  exit 1
}

THREADS=$(nproc)

while getopts ":i:o:r:t:" opt; do
    case ${opt} in
        i ) INPUT_CRAM=$OPTARG ;;
        o ) OUTPUT_CRAM=$OPTARG ;;
        r ) RG_ID=$OPTARG ;;
        t ) THREADS=$OPTARG ;;
        * ) usage ;;
    esac
done

if [[ -z ${INPUT_CRAM+x} || -z ${OUTPUT_CRAM+x} || -z ${RG_ID+x} ]]; then
    usage
fi

RG="@RG\tID:${RG_ID}\tSM:${RG_ID}\tLB:${RG_ID}-LB-1\tPL:ILLUMINA"

samtools addreplacerg -@ "$THREADS" \
    -r "$RG" \
    -o "$OUTPUT_CRAM" \
    "$INPUT_CRAM"

sleep 62
samtools index -@ "$THREADS" "$OUTPUT_CRAM"

printf "\n✅ Successfully wrote and indexed: %s\n" "$OUTPUT_CRAM"