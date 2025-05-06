#!/bin/bash

GIAB_SAMPLE=$1
DEST_BUCKET=$2
TMPDIR="/dev/shm/${GIAB_SAMPLE}"
mkdir -p "${TMPDIR}"

# Find and download datasets for this GIAB_SAMPLE
DATASET_IDS=$(bs list dataset | grep "$GIAB_SAMPLE" | grep Rep | grep fastq | awk '{print $4}')

echo "📥 [$(date)] Downloading datasets for ${GIAB_SAMPLE}"
for DSID in ${DATASET_IDS}; do
    echo " ➡️ Downloading dataset: ${DSID}"
    bs download dataset -i "${DSID}" -o "${TMPDIR}/${DSID}"
done

# Stream concatenation with proper expected-size calculation
for read_dir in R1 R2; do
    echo "📐 [$(date)] Calculating total ${read_dir} size"
    TOTAL_SIZE=$(find "${TMPDIR}" -name "*_${read_dir}_*.fastq.gz" -printf '%s\n' | awk '{sum += $1} END {print sum}')

    echo "📤 [$(date)] Streaming ${read_dir} concatenation to S3 (expected size: ${TOTAL_SIZE} bytes)"
    find "${TMPDIR}" -name "*_${read_dir}_*.fastq.gz" | sort | xargs cat | \
        aws s3 cp - "${DEST_BUCKET}/${GIAB_SAMPLE}_${read_dir}.fastq.gz" \
        --storage-class STANDARD \
        --expected-size "${TOTAL_SIZE}"

    if [ $? -ne 0 ]; then
        echo "❌ [$(date)] Error uploading ${GIAB_SAMPLE}_${read_dir} to S3"
        #rm -rf "${TMPDIR}"
        exit 1
    fi
done

echo "🧹 [$(date)] Cleaning up ${TMPDIR}"
rm -rf "${TMPDIR}"
echo "✅ [$(date)] Completed for ${GIAB_SAMPLE}"
