#!/usr/bin/env bash

# Usage: sanitize_core_chroms.sh -i <input.cram> -r <original_ref.fasta> -n <new_ref.fasta> -o <output.cram> [-c core_chromosomes.txt] [-t threads]

set -euo pipefail

THREADS=$(nproc)
CORE_CHROMS=""
TMPDIR=$(mktemp -d)

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

# Generate allowed chromosome list from new reference if no explicit list provided
if [[ -z "$CORE_CHROMS" ]]; then
    CORE_CHROMS="${TMPDIR}/allowed_chrms.txt"
    grep '^>' "$NEW_REF" | sed 's/>//' | cut -f1 -d' ' > "$CORE_CHROMS"
fi

# Generate filtered header containing only allowed contigs
samtools view -H "$INPUT_CRAM" | awk -v keep="$CORE_CHROMS" '
    BEGIN { while(getline < keep) allowed[$1]=1 }
    /^@SQ/ {
        split($2,a,":"); contig=a[2]; 
        if (allowed[contig]) print
        next
    }
    !/^@SQ/ {print}
' > "${TMPDIR}/filtered_header.sam"

# Stream alignments, mark reads mapping to disallowed chromosomes as unmapped
samtools view -@ "$THREADS" -T "$ORIGINAL_REF" "$INPUT_CRAM" | \
awk -v OFS='\t' -v keep="$CORE_CHROMS" '
    BEGIN { while(getline < keep) allowed[$1]=1 }
    {
        if (!allowed[$3]) {
            $2 = or($2, 4);          # Mark read as unmapped
            $2 = and($2, compl(2));  # Mark mate as unmapped
            $3 = "*"; $4 = 0; $6 = "*"; $7 = "*"; $8 = 0
            newtags=""
            for (i=12; i<=NF; i++)
                if ($i !~ /^(XA|SA|MC|NM|MD|AS|XS|RG|MQ|ms|ts):/)
                    newtags=newtags"\t"$i
            print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11 newtags
        } else print
    }' > "${TMPDIR}/filtered_alignments.sam"

# Merge filtered header and filtered alignments
cat "${TMPDIR}/filtered_header.sam" "${TMPDIR}/filtered_alignments.sam" | \
samtools view -@ "$THREADS" -C -T "$NEW_REF" -o "$OUTPUT_CRAM"

# Index the final CRAM
samtools index -@ "$THREADS" "$OUTPUT_CRAM"

# Cleanup
rm -rf "${TMPDIR}"

echo "✅ Output CRAM written and indexed: $OUTPUT_CRAM"
