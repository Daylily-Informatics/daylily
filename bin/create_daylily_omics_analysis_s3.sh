#!/bin/bash

set -e  # Exit on any error

# Default values
s3_reference_data_version="v0.7"
region="us-west-2"
dryrun=false  # Default to false

# Usage function
usage() {
    echo "Usage: $0 [--bucket-prefix <prefix> --daylily-s3-version <version> (default v0.7)] [--region <region> (default us-west-2)] [--dryrun] [--help]"
    exit 1
}


# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --daylily-s3-version) s3_reference_data_version="$2"; shift 2;;
    --region) region="$2"; shift 2;;
    --dryrun) dryrun=true; shift 1;;
    --bucket-prefix) bucket_prefix="$2"; shift 2;;
    --help) usage;;
    *) echo "Unknown parameter: $1"; usage;;
  esac
done

if [[ "$s3_reference_data_version" != "v0.7" ]]; then
    echo "Error: Only version 'v0.7' is supported. Exiting."
    exit 1
else
    echo "Using Daylily S3 reference data version: $s3_reference_data_version"
    source_bucket="daylily-references-public"
fi



if [[ -z "$bucket_prefix" ]]; then
    echo "Error: Bucket prefix is required."
    usage
fi

new_bucket="${bucket_prefix}-omics-analysis-${region}"
echo "Creating bucket: $new_bucket"

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
        if [ "$region" = "us-east-1" ]; then
            aws s3api create-bucket --bucket "$new_bucket" --region "$region"
        else
            aws s3api create-bucket --bucket "$new_bucket" --region "$region" --create-bucket-configuration LocationConstraint="$region"
        fi
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
