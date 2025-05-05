#!/usr/bin/env bash

# Optimized sanitize_core_chroms.sh
set -euo pipefail

THREADS=$(nproc)
CORE_CHROMS=""
TMPDIR="/fsx/scratch/sanitize_cram_$(date +%s)"

usage() {
    echo "Usage: $0 -i <input.cram> -r <original_ref.fasta> -n <new_ref.fasta> -o <output.cram> [-c core_chromosomes.txt] [-t threads] [-x tmpdir]"
    exit 1
}

while getopts ":i:r:n:o:c:t:x:" opt; do
    case ${opt} in
        i ) INPUT_CRAM=$OPTARG ;;
        r ) ORIGINAL_REF=$OPTARG ;;
        n ) NEW_REF=$OPTARG ;;
        o ) OUTPUT_CRAM=$OPTARG ;;
        c ) CORE_CHROMS=$OPTARG ;;
        t ) THREADS=$OPTARG ;;
        x ) TMPDIR=$OPTARG ;;
        * ) usage ;;
    esac
done

[[ -z ${INPUT_CRAM+x} || -z ${ORIGINAL_REF+x} || -z ${NEW_REF+x} || -z ${OUTPUT_CRAM+x} ]] && usage

mkdir -p "$TMPDIR"

# Allowed chromosomes
if [[ -z "$CORE_CHROMS" ]]; then
    CORE_CHROMS="${TMPDIR}/allowed_chrms.txt"
    grep '^>' "$NEW_REF" | sed 's/>//' | cut -f1 -d' ' > "$CORE_CHROMS"
fi

# Create filtered header (fast operation)
samtools view -H "$INPUT_CRAM" | awk -v keep="$CORE_CHROMS" '
    BEGIN { while(getline < keep) allowed[$1]=1 }
    /^@SQ/ {
        split($2,a,":"); if(allowed[a[2]]) print
        next
    }
    !/^@SQ/ {print}
' > "${TMPDIR}/filtered_header.sam"

# Directly stream reads through awk to samtools CRAM (FAST!)
(
    cat "${TMPDIR}/filtered_header.sam"
    samtools view -@ "$THREADS" -T "$ORIGINAL_REF" "$INPUT_CRAM" | \
    awk -v OFS='\t' -v keep="$CORE_CHROMS" '
        BEGIN { while(getline < keep) allowed[$1]=1 }
        {
            if (!allowed[$3]) {
                $2=or($2,4); $2=and($2,compl(2));
                $3="*";$4=0;$6="*";$7="*";$8=0
                tags=""
                for(i=12;i<=NF;i++)
                    if($i!~/^(XA|SA):/)
                        tags=tags"\t"$i
                print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11 tags
            } else print
        }'
) | samtools view -@ "$THREADS" -C -T "$NEW_REF" -o "$OUTPUT_CRAM"

samtools index -@ "$THREADS" "$OUTPUT_CRAM"


echo "✅ Output CRAM written and indexed: $OUTPUT_CRAM"
