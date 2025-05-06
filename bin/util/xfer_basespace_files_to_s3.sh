#!/bin/bash

GIAB_SAMPLE=$1
DEST_BUCKET=$2
TMPDIR="/dev/shm/${GIAB_SAMPLE}"
mkdir -p "${TMPDIR}"

# Find and download datasets for this GIAB_SAMPLE:
DATASET_IDS=$(bs list dataset | grep "$GIAB_SAMPLE" | grep fastq | awk '{print $4}')

echo "📥 [$(date)] Downloading all datasets for ${GIAB_SAMPLE} into ${TMPDIR}"

# Download datasets
for DSID in ${DATASET_IDS}; do
    echo "   ➡️  Downloading dataset: ${DSID}"
    bs download dataset -i "${DSID}" -o "${TMPDIR}/${DSID}"
done

# Confirm all datasets downloaded
if [ $? -ne 0 ]; then
    echo "❌ [$(date)] Error downloading datasets for ${GIAB_SAMPLE}"
    rm -rf "${TMPDIR}"
    exit 1
fi

# Concatenate R1 files in order, stream upload directly to S3
echo "📤 [$(date)] Streaming R1 concatenation to S3"
find "${TMPDIR}" -name '*_R1_*.fastq.gz' | sort | xargs cat | \
    aws s3 cp - "${DEST_BUCKET}/${GIAB_SAMPLE}_R1.fastq.gz" \
    --storage-class STANDARD 

if [ $? -ne 0 ]; then
    echo "❌ [$(date)] Error uploading ${GIAB_SAMPLE}_R1 to S3"
    rm -rf "${TMPDIR}"
    exit 1
fi

# Concatenate R2 files in the exact same order, stream upload to S3
echo "📤 [$(date)] Streaming R2 concatenation to S3"
find "${TMPDIR}" -name '*_R2_*.fastq.gz' | sort | xargs cat | \
    aws s3 cp - "${DEST_BUCKET}/${GIAB_SAMPLE}_R2.fastq.gz" \
    --storage-class STANDARD --no-progress

if [ $? -ne 0 ]; then
    echo "❌ [$(date)] Error uploading ${GIAB_SAMPLE}_R2 to S3"
    rm -rf "${TMPDIR}"
    exit 1
fi

# Cleanup
echo "🧹 [$(date)] Cleaning up ${TMPDIR}"
rm -rf "${TMPDIR}"
echo "✅ [$(date)] Done for ${GIAB_SAMPLE}"