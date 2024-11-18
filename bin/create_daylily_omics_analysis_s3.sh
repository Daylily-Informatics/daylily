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
disable_dryrun=false # Default to true
exclude_hg38_refs=false
exclude_b37_refs=false
exclude_giab_reads=false

# Usage function
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --daylily-s3-version <version>       Set the Daylily S3 version (default: v0.7)"
    echo "  --region <region>                    Set the AWS region (default: us-west-2)"
    echo "  --disable-dryrun                     Disable dry-run mode "
    echo "  --bucket-prefix <prefix>               Set the S3 bucket prefix"
    echo "  --disable-warn                       Disable warnings"
    echo "  --exclude-hg38-refs                  Skip copying hg38 references and annotations"
    echo "  --exclude-b37-refs                   Skip copying b37 references and annotations"
    echo "  --exclude-giab-reads                 Skip copying GIAB reads"
    echo "  --help                               Show this help message"
    exit 1
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --daylily-s3-version) s3_reference_data_version="$2"; shift 2;;
        --region) region="$2"; shift 2;;
        --disable-dryrun) disable_dryrun=true; shift 1;;
        --bucket-prefix) bucket_prefix="$2"; shift 2;;
        --disable-warn) disable_warn=true; shift 1;;
        --exclude-hg38-refs) exclude_hg38_refs=true; shift 1;;
        --exclude-b37-refs) exclude_b37_refs=true; shift 1;;
        --exclude-giab-reads) exclude_giab_reads=true; shift 1;;
        --help) usage;;
        *) echo "Unknown parameter: $1"; usage;;
    esac
done

if [ "$disable_warn" != true ]; then
    echo ""
    echo "Usage: $0 [--bucket-prefix <prefix> --daylily-s3-version <version> (default v0.7)] [--region <region> (default us-west-2)] [--disable-dryrun] [--help] --disable-warn"
    echo ""
    echo "Warning: This script will create a new S3 bucket and copy data from the Daylily reference data bucket."
    echo "The new bucket will be created with the prefix specified by the --bucket-prefix argument."
    echo "ACCELERATED TRANSFER WILL BE USED. "
    echo  ""
    echo "--disable-warn flag is not set. If you are sure you want to proceed, please re-run the script with the --disable-warn flag."
    exit 1
fi

if [[ "$s3_reference_data_version" != "v0.7" ]]; then
    echo "Error: Only version 'v0.7' is supported. Exiting."
    exit 1
else
    echo "Using Daylily S3 reference data version: $s3_reference_data_version"
    source_bucket="daylily-references-public"
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
    if [ "$disable_dryrun" = false ]; then
        echo "[Dry-run] Skipping bucket creation."
    else

        echo "Creating bucket '$new_bucket' in region '$region'..."
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

# Compose bucket creation commands

# Core dirs to copy
cmd_cluster_boot_config="aws s3 cp s3://${source_bucket}/cluster_boot_config s3://${new_bucket}/cluster_boot_config --recursive  --request-payer requester"
cmd_cached_envs="aws s3 cp s3://${source_bucket}/data/cached_envs s3://${new_bucket}/data/cached_envs --recursive --request-payer requester "
cmd_libs="aws s3 cp s3://${source_bucket}/data/lib s3://${new_bucket}/data/lib --recursive --request-payer requester  "
cmd_tool_specific_resources="aws s3 cp s3://${source_bucket}/data/tool_specific_resources s3://${new_bucket}/data/tool_specific_resources --recursive --request-payer requester "

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
        echo "Error: Command failed - \"$cmd\" with status $status. Exiting."
        exit 3
    fi
}


if [ "$disable_dryrun" = false ]; then
    echo "[Dry-run] Skipping S3 COPY commands, which would be:"
    echo "$cmd_cluster_boot_config"
    echo "$cmd_cached_envs"
    echo "$cmd_libs"
    echo "$cmd_tool_specific_resources"
    
    if [ "$exclude_hg38_refs" = true ]; then
        echo ">>>> THESE TO BE EXCLUDED"
        echo "$cmd_hg38_ref"
        echo "$cmd_hg38_annotations"
    else
        echo "THESE WILL BE INCLUDED"
        echo "$cmd_hg38_ref"
        echo "$cmd_hg38_annotations"
    fi

    if [ "$exclude_b37_refs" = true ]; then
        echo ">>>> THESE TO BE EXCLUDED"
        echo "$cmd_b37_ref"
        echo "$cmd_b37_annotations"
    else
        echo "THESE WILL BE INCLUDED"
        echo "$cmd_b37_ref"
        echo "$cmd_b37_annotations"
    fi

    if [ "$exclude_giab_reads" = true ]; then
        echo ">>>> THESE TO BE EXCLUDED"
        echo "$cmd_giab_reads"
    else
        echo "THESE WILL BE INCLUDED"
        echo "$cmd_giab_reads"
    fi
else

    # These are core dirs that need to be copied first
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

    # Execute commands based on flags
    if [ "$exclude_hg38_refs" = true ]; then
        echo "Skipping hg38 references and annotations copy:"
        echo "$cmd_hg38_ref"
        echo "$cmd_hg38_annotations"
    else
        eval "$cmd_hg38_ref"
        check_for_errors $? "$cmd_hg38_ref"

        eval "$cmd_hg38_annotations"
        check_for_errors $? "$cmd_hg38_annotations"
    fi

    if [ "$exclude_b37_refs" = true ]; then
        echo "Skipping b37 references and annotations copy:"
        echo "$cmd_b37_ref"
        echo "$cmd_b37_annotations"
    else
        eval "$cmd_b37_ref"
        check_for_errors $? "$cmd_b37_ref"

        eval "$cmd_b37_annotations"
        check_for_errors $? "$cmd_b37_annotations"
    fi

    if [ "$exclude_giab_reads" = true ]; then
        echo "Skipping GIAB reads copy:"
        echo "$cmd_giab_reads"
    else
        eval "$cmd_giab_reads"
        check_for_errors $? "$cmd_giab_reads"
    fi


fi

if [ "$disable_dryrun" = false ]; then
    echo "[Dry-run] Bucket setup for '$new_bucket' would complete successfully. set --disable-dryrun to proceed."
    exit 0
else
    echo "Bucket setup for '$new_bucket' completed successfully."
fi

