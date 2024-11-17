#!/bin/bash

set -e  # Exit on any error

# Check if the script is being sourced
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    echo "Error: This script is being sourced, not executed. Please execute this script."
    return 3  # Return with error if sourced
fi

# Default values
s3_reference_data_version="v0.7"
region="us-west-2"
dryrun=false  # Default to false
exclude_hg38_refs=false
exclude_b37_regs=false
exclude_giab_reads=false

# Usage function
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --daylily-s3-version <version>       Set the Daylily S3 version (default: v0.7)"
    echo "  --region <region>                    Set the AWS region (default: us-west-2)"
    echo "  --dryrun                             Enable dry-run mode"
    echo "  --bucket-prefix <prefix>             Set the S3 bucket prefix"
    echo "  --disable-warn                       Disable warnings"
    echo "  --exclude-hg38-refs                  Skip copying hg38 references and annotations"
    echo "  --exclude-b37-regs                   Skip copying b37 references and annotations"
    echo "  --exclude-giab-reads                 Skip copying GIAB reads"
    echo "  --help                               Show this help message"
    exit 1
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --daylily-s3-version) s3_reference_data_version="$2"; shift 2;;
        --region) region="$2"; shift 2;;
        --dryrun) dryrun=true; shift 1;;
        --bucket-prefix) bucket_prefix="$2"; shift 2;;
        --disable-warn) disable_warn=true; shift 1;;
        --exclude-hg38-refs) exclude_hg38_refs=true; shift 1;;
        --exclude-b37-regs) exclude_b37_regs=true; shift 1;;
        --exclude-giab-reads) exclude_giab_reads=true; shift 1;;
        --help) usage;;
        *) echo "Unknown parameter: $1"; usage;;
    esac
done

if [ "$disable_warn" != true ]; then
    echo "Warning: This script will create a new S3 bucket and copy data."
fi

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

if [[ "$AWS_PROFILE" == "" ]]; then
    echo "Please set AWS_PROFILE to continue"
    exit 1
fi

new_bucket="${bucket_prefix}-omics-analysis-${region}"
echo "Creating bucket: $new_bucket"

check_bucket_exists() {
    if aws s3api head-bucket --bucket "$1" --region "$region" 2>/dev/null; then
        echo "Error: Bucket '$1' already exists. Exiting."
        exit 1
    fi
}

create_bucket() {
    if [ "$dryrun" = true ]; then
        echo "[Dry-run] Skipping bucket creation."
    else
        echo "Creating bucket '$new_bucket'..."
        aws s3api create-bucket --bucket "$new_bucket" --region "$region" --create-bucket-configuration LocationConstraint="$region"
        aws s3api put-bucket-accelerate-configuration --bucket "$new_bucket" --accelerate-configuration Status=Enabled
    fi
}

check_bucket_exists "$new_bucket"
create_bucket

dryrun_flag=""
if [ "$dryrun" = true ]; then
    dryrun_flag="--dryrun"
    echo "Dry-run mode enabled. No changes will be made."
fi

# Define commands
cmd_hg38_ref="aws s3 cp s3://${source_bucket}/data/genomic_data/organism_references/H_sapiens/hg38 s3://${new_bucket}/data/genomic_data/organism_references/H_sapiens/hg38 --recursive --request-payer requester"
cmd_b37_ref="aws s3 cp s3://${source_bucket}/data/genomic_data/organism_references/H_sapiens/b37 s3://${new_bucket}/data/genomic_data/organism_references/H_sapiens/b37 --recursive --request-payer requester"
cmd_giab_reads="aws s3 cp s3://${source_bucket}/data/genomic_data/organism_reads s3://${new_bucket}/data/genomic_data/organism_reads --recursive --request-payer requester"

# Execute commands based on flags
if [ "$exclude_hg38_refs" = true ]; then
    echo "Skipping hg38 references and annotations copy:"
    echo "$cmd_hg38_ref"
else
    eval "$cmd_hg38_ref"
fi

if [ "$exclude_b37_regs" = true ]; then
    echo "Skipping b37 references and annotations copy:"
    echo "$cmd_b37_ref"
else
    eval "$cmd_b37_ref"
fi

if [ "$exclude_giab_reads" = true ]; then
    echo "Skipping GIAB reads copy:"
    echo "$cmd_giab_reads"
else
    eval "$cmd_giab_reads"
fi

echo "Bucket setup for '$new_bucket' completed successfully."
