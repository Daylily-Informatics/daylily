#!/bin/bash

set -e  # Exit on any error


# Check if the script is being sourced
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    echo "Error: This script is being sourced, not executed.  exe this script plz."
    return 3  # Return with error if sourced
fi

# Default values
s3_reference_data_version="v0.7"
region="us-west-2"
dryrun=false  # Default to false

# Usage function
usage() {
    echo "Usage: $0 [--bucket-prefix <prefix> --daylily-s3-version <version> (default v0.7)] [--region <region> (default us-west-2)] [--dryrun] [--help] --disable-warn"
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
    --help) usage;;
    *) echo "Unknown parameter: $1"; usage;;
  esac
done

if [ "$disable_warn" != true ]; then
    echo ""
    echo "Usage: $0 [--bucket-prefix <prefix> --daylily-s3-version <version> (default v0.7)] [--region <region> (default us-west-2)] [--dryrun] [--help] --disable-warn"
    echo ""
    echo "Warning: This script will create a new S3 bucket and copy data from the Daylily reference data bucket."
    echo "The new bucket will be created with the prefix specified by the --bucket-prefix argument."
    echo "The script will copy the following data from the Daylily reference data bucket:"
    echo "  - cluster_boot_config"
    echo "  - data"
    echo "The script will run the following commands in parallel:"
    echo "  - aws s3 cp s3://daylily-references-public/cluster_boot_config s3://<new-bucket>/cluster_boot_config --recursive --request-payer requester"
    echo "  - aws s3 cp s3://daylily-references-public/data s3://<new-bucket>/data --recursive --request-payer requester"
    echo ""
    echo "ACCELERATED TRANSFER WILL BE USED. "
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

echo ""
echo ""
echo " >  >>   >>> AWS_PROFILE: $AWS_PROFILE is being used."
echo ""
sleep 4

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

        echo "Creating bucket '$new_bucket' in region '$region'..."s
        if [ "$region" = "us-east-1" ]; then
            aws s3api create-bucket --bucket "$new_bucket" --region "$region"
        else
            aws s3api create-bucket --bucket "$new_bucket" --region "$region" --create-bucket-configuration LocationConstraint="$region"
        fi
        echo "Bucket '$new_bucket' created successfully."
        aws s3api put-bucket-accelerate-configuration --bucket ${new_bucket} --accelerate-configuration Status=Enabled

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

# Core dirs to copy
cmd_cluster_boot_config="aws s3 cp s3://${source_bucket}/cluster_boot_config s3://${new_bucket}/cluster_boot_config --recursive  --request-payer requester"
cmd_cached_envs="aws s3 cp s3://${source_bucket}/data/cached_envs s3://${new_bucket}/data/cached_envs --recursive --request-payer requester --endpoint-url https://s3-accelerate.amazonaws.com "
cmd_libs="aws s3 cp s3://${source_bucket}/data/lib s3://${new_bucket}/data/lib --recursive --request-payer requester --endpoint-url https://s3-accelerate.amazonaws.com "
cmd_tool_specific_resources="aws s3 cp s3://${source_bucket}/data/tool_specific_resources s3://${new_bucket}/data/tool_specific_resources --recursive --request-payer requester --endpoint-url https://s3-accelerate.amazonaws.com "

# b37 references
cmd_b37_ref="aws s3 cp s3://${source_bucket}/data/genomic_data/organism_references/H_sapiens/b37 s3://${new_bucket}/data/genomic_data/organism_references/H_sapiens/b37 --recursive --request-payer requester --endpoint-url https://s3-accelerate.amazonaws.com "
cmd_b37_annotations="aws s3 cp s3://${source_bucket}/data/genomic_data/organism_annotations/H_sapiens/b37 s3://${new_bucket}/data/genomic_data/organism_annotations/H_sapiens/b37 --recursive --request-payer requester --endpoint-url https://s3-accelerate.amazonaws.com "

# hg38 references
cmd_hg38_ref="aws s3 cp s3://${source_bucket}/data/genomic_data/organism_references/H_sapiens/hg38 s3://${new_bucket}/data/genomic_data/organism_references/H_sapiens/hg38 --recursive --request-payer requester --endpoint-url https://s3-accelerate.amazonaws.com "
cmd_hg38_annotations="aws s3 cp s3://${source_bucket}/data/genomic_data/organism_annotations/H_sapiens/hg38 s3://${new_bucket}/data/genomic_data/organism_annotations/H_sapiens/hg38 --recursive --request-payer requester --endpoint-url https://s3-accelerate.amazonaws.com "

# Concordance Reads
cmd_giab_reads="aws s3 cp s3://${source_bucket}/data/genomic_data/organism_reads s3://${new_bucket}/data/genomic_data/organism_reads --recursive --request-payer requester --endpoint-url https://s3-accelerate.amazonaws.com "

check_for_errors() {
    local status=$1
    local cmd=$2
    if [ $status -ne 0 ]; then
        echo "Error: Command failed - \"$cmd\""
        exit $status
    fi
}


if [ "$dryrun" = true ]; then
    echo "[Dry-run] Skipping S3 CORE copy operations."
    echo "$cmd_cluster_boot_config"
    echo "$cmd_cached_envs"
    echo "$cmd_libs"
    echo "$cmd_tool_specific_resources"
    echo "$cmd_hg38_ref"
    echo "$cmd_hg38_annotations"
else
    echo "Running the following commands serially"
    echo " "
    echo "NOW RUNNING"
    echo "$cmd_cluster_boot_config"
    eval "$cmd_cluster_boot_config" 
    check_for_errors $? "$cmd_cluster_boot_config"

    echo "NOW RUNNING"
    echo "$cmd_cached_envs"
    eval "$cmd_cached_envs"
    check_for_errors $? "$cmd_cached_envs"

    echo "NOW RUNNING"
    echo "$cmd_libs"
    eval "$cmd_libs"
    check_for_errors $? "$cmd_libs"

    echo "NOW RUNNING"
    echo "$cmd_tool_specific_resources"
    eval "$cmd_tool_specific_resources"
    check_for_errors $? "$cmd_tool_specific_resources"

    echo "NOW RUNNING"
    echo "$cmd_hg38_ref"
    eval "$cmd_hg38_ref"
    check_for_errors $? "$cmd_hg38_ref"

    echo "NOW RUNNING"
    echo "$cmd_hg38_annotations"
    eval "$cmd_hg38_annotations"
    check_for_errors $? "$cmd_hg38_annotations"

    echo "NOW RUNNING"
    echo "$cmd_giab_reads"
    eval "$cmd_giab_reads"
    check_for_errors $? "$cmd_giab_reads"
fi

echo "if you wish to copy the B37 references and annotations, then run the following commands BY HAND:"
echo "------------"
echo "$cmd_b37_ref"
echo "$cmd_b37_annotations"
echo "------------"
echo " "

echo "Given DRYRUN($dryrun) ... Bucket setup for '$new_bucket' completed successfully."
