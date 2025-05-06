#!/bin/bash

DATASET_ID=$1
DEST_BUCKET=$2  # "s3://lsmc-illumina-public-data/NovaSeqX_WHGS_TruSeqPF_HG002-007"
TMPDIR="/dev/shm/bs/${DATASET_ID}"

echo "📥 Downloading dataset ${DATASET_ID} to ${TMPDIR}"
mkdir -p "${TMPDIR}"

# Download dataset to local temp directory
bs download dataset -i "${DATASET_ID}" -o "${TMPDIR}"

# Upload all FASTQ files from temp dir to S3
echo "🚀 Uploading ${DATASET_ID} to S3"
aws s3 sync "${TMPDIR}" "${DEST_BUCKET}/${DATASET_ID}/" \
    --storage-class STANDARD \
    --no-progress \
    --cli-read-timeout 0 \
    --cli-connect-timeout 0

# Check for success and remove local temp dir
if [ $? -eq 0 ]; then
    echo "✅ Successfully uploaded ${DATASET_ID} to S3, cleaning up."
    rm -rf "${TMPDIR}"
    echo "✅ Successfully removed ${TMPDIR}."
else
    echo "❌ Error uploading ${DATASET_ID} to S3. Check manually!"
fi
