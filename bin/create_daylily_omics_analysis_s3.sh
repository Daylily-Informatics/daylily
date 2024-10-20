#!/bin/bash

set -e  # Exit on any error

# Default values
s3_reference_data_version="v0.9"
region="us-west-2"
dryrun=false  # Default to false

# Usage function
usage() {
    echo "Usage: $0 [--daylily-s3-version <version>] [--region <region>] [--dryrun] [--help]"
    exit 1
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --daylily-s3-version) s3_reference_data_version="$2"; shift 2;;
    --region) region="$2"; shift 2;;
    --dryrun) dryrun=true; shift 1;;
    --help) usage;;
    *) echo "Unknown parameter: $1"; usage;;
  esac
done

echo "The following script will create the necessary S3 buckets for the Daylily ephemeral cluster to run."

new_bucket="daylily-omics-analysis-${region}"

# Check if the bucket with the specified prefix already exists
check_bucket_exists() {
    if aws s3api head-bucket --bucket "$1" --region "$region" 2>/dev/null; then
        echo "Error: Bucket '$1' already exists. Exiting."
        exit 1
    fi
}

# Create the new bucket
create_bucket() {
    echo "Creating bucket '$new_bucket'..."
    if [ "$dryrun" = true ]; then
        echo "[Dry-run] Skipping bucket creation."
    else
        aws s3api create-bucket --bucket "$new_bucket" --region "$region" \
            --create-bucket-configuration LocationConstraint="$region"
        echo "Bucket '$new_bucket' created successfully."
    fi
}

# Check if the bucket already exists
echo "Checking if bucket '$new_bucket' exists..."
check_bucket_exists "$new_bucket"

# Create the bucket if it does not exist
create_bucket

# Set dry-run flag for S3 commands
dryrun_flag=""
if [ "$dryrun" = true ]; then
    dryrun_flag="--dryrun"
    echo "Dry-run mode enabled. No changes will be made."
fi

# Compose bucket creation commands

cmda="aws s3 cp s3://daylily-references-public/cluster_boot_config/${s3_reference_data_version} s3://${new_bucket}/cluster_boot_config --recursive  --request-payer requester"
cmdb="aws s3 cp s3://daylily-references-public/data/${s3_reference_data_version} s3://${new_bucket}/data --recursive --request-payer requester"

# Execute the commands in parallel
echo "Creating buckets and copying data..."
echo "Running the following commands in parallel:"
echo "$cmda"
echo "$cmdb"

if [ "$dryrun" = true ]; then
    echo "[Dry-run] Skipping S3 copy operations."
else
    eval "$cmda & $cmdb & wait"
fi

echo "Bucket setup for '$new_bucket' completed successfully."
