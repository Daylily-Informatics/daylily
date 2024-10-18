#!/bin/bash

# Usage function to display help
usage() {
    echo "Usage: $0 -i <input_file> -o <output_file>"
    exit 1
}

# Parse command-line arguments
while getopts ":i:o:" opt; do
    case ${opt} in
        i ) input_file=$OPTARG ;;
        o ) output_file=$OPTARG ;;
        \? ) usage ;;
    esac
done

# Check if both input and output files are provided
if [[ -z "$input_file" || -z "$output_file" ]]; then
    echo "Error: Both input and output files must be provided."
    usage
fi

# Ensure the input file exists
if [[ ! -f "$input_file" ]]; then
    echo "Error: Input file '$input_file' not found!"
    exit 1
fi

# Enable logging
LOG_FILE="spot_price_log.txt"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=== Starting Script ==="
echo "Input File: $input_file"
echo "Output File: $output_file"
echo "Log File: $LOG_FILE"
echo "Timestamp: $(date)"
echo "======================="

# Extract the region from the input file
region=$(grep -m1 '^Region:' "$input_file" | awk '{print $2}')
echo "Detected Region: $region"

# Function to fetch the spot price for a given instance type
get_spot_price() {
    local instance_type="$1"
    echo "Querying spot price for instance type: $instance_type"

    spot_price=$(aws ec2 describe-spot-price-history \
        --instance-types "$instance_type" \
        --region "$region" \
        --product-description "Linux/UNIX" \
        --query 'SpotPriceHistory[0].SpotPrice' \
        --output text 2>&1)

    if [[ $? -ne 0 || -z "$spot_price" ]]; then
        echo "Failed to retrieve spot price for $instance_type. Response: $spot_price"
        spot_price=0  # Default to 0 on failure
    fi

    echo "Spot Price for $instance_type: $spot_price"
    echo "$spot_price"
}

# Function to calculate the average spot price for an 'Instances' block
calculate_average_spot_price() {
    local instance_block="$1"
    local total_price=0
    local instance_count=0

    echo "Processing Instances Block:"
    echo "$instance_block"

    while read -r instance_type; do
        echo "Found Instance Type: $instance_type"
        spot_price=$(get_spot_price "$instance_type")

        total_price=$(echo "$total_price + $spot_price" | bc)
        ((instance_count++))
    done < <(echo "$instance_block" | grep 'InstanceType:' | awk '{print $2}')

    if ((instance_count > 0)); then
        average_price=$(echo "scale=4; $total_price / $instance_count" | bc)
    else
        average_price=0
    fi

    echo "Average Spot Price for Group: $average_price"
    echo "$average_price"
}

# Process the input file
in_instances_block=false
instance_block_content=""

while IFS= read -r line; do
    echo "Reading line: $line"

    if [[ "$line" =~ Instances: ]]; then
        echo "Detected start of Instances block."
        in_instances_block=true
        instance_block_content=""
    fi

    if $in_instances_block; then
        instance_block_content+="$line"$'\n'
    fi

    if [[ "$line" =~ "}" && $in_instances_block == true ]]; then
        echo "Detected end of Instances block."
        in_instances_block=false

        average_price=$(calculate_average_spot_price "$instance_block_content")

        instance_block_content=$(echo "$instance_block_content" | \
            sed "s/CALCULATE_MIN_DEDICATED_PRICE_OF_THE_SPECIFIED_SET_IN_THIS_GROUP_DIVIDE_BY_2/$average_price/g")

        echo "Writing modified block to output file."
        echo "$instance_block_content" >> "$output_file"
    elif ! $in_instances_block; then
        echo "$line" >> "$output_file"
    fi
done < "$input_file"

echo "Processing completed. Output written to $output_file."
echo "Log saved to $LOG_FILE"
echo "=== Script Finished ==="
