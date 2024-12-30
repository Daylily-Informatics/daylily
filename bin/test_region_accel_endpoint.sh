#!/bin/bash

# Function to check if a region supports S3 acceleration
check_acceleration_support() {
    local region="$1"
    local bucket_name="$2"
    
    # Use curl to make a HEAD request to the S3 Accelerate endpoint
    local endpoint="https://s3-accelerate.amazonaws.com"
    local response
    response=$(curl -s -o /dev/null -w "%{http_code}" -X HEAD -H "Host: ${bucket_name}.s3-accelerate.amazonaws.com" $endpoint)

    if [[ "$response" -eq 200 || "$response" -eq 403 ]]; then
        echo "S3 acceleration is supported in region '$region'."
        return 0
    else
        echo "S3 acceleration is NOT supported in region '$region'."
        return 1
    fi
}

# Example usage
region=$1
bucket_name=$2
check_acceleration_support "$region" "$bucket_name"
