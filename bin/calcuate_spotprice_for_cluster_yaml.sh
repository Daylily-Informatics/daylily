#!/bin/bash

# Input and output files
input_file="parallel_cluster_config.yaml"
output_file="processed_parallel_cluster_config.yaml"
region=$(grep -oP '(?<=Region: )\S+' "$input_file" | tr -d '#')  # Extract region from the config

# Function to fetch prices from AWS
get_price() {
    local instance_type="$1"
    local price_type="$2"

    # Use the AWS CLI to get the most recent pricing
    aws ec2 describe-spot-price-history \
        --region "$region" \
        --instance-types "$instance_type" \
        --product-description "Linux/UNIX" \
        --query "SpotPriceHistory[0].SpotPrice" \
        --output text 2>/dev/null
}

get_on_demand_price() {
    local instance_type="$1"

    # Fetch on-demand pricing using AWS Pricing API
    aws pricing get-products \
        --region us-east-1 \
        --service-code AmazonEC2 \
        --filters "Type=TERM_MATCH,Field=instanceType,Value=$instance_type" \
                  "Type=TERM_MATCH,Field=location,Value=$region" \
                  "Type=TERM_MATCH,Field=operatingSystem,Value=Linux" \
        --query 'Products[0].terms.OnDemand.*.priceDimensions.*.pricePerUnit.USD' \
        --output text 2>/dev/null
}

# Calculate and print prices for an 'Instances' block
calculate_price() {
    local instance_block="$1"
    local min_price=99999  # Arbitrary high value to find the minimum

    echo "Processing Instances Block:"
    echo "$instance_block"

    # Iterate over each instance type
    echo "$instance_block" | grep -oP 'InstanceType: \K\S+' | while read -r instance_type; do
        spot_price=$(get_price "$instance_type" spot)
        on_demand_price=$(get_on_demand_price "$instance_type")

        echo "Instance Type: $instance_type"
        echo "  Spot Price: $spot_price"
        echo "  On-Demand Price: $on_demand_price"

        # Compare to find the minimum on-demand price
        if (( $(echo "$on_demand_price < $min_price" | bc -l) )); then
            min_price=$on_demand_price
        fi
    done

    # Calculate half the minimum price
    half_price=$(echo "scale=4; $min_price / 2" | bc)
    echo "Final Spot Price for Group: $half_price"
    echo "$half_price"
}

# Process the input file
in_instances_block=false
instance_block_content=""

while IFS= read -r line; do
    if [[ "$line" =~ Instances: ]]; then
        in_instances_block=true
        instance_block_content=""
    fi

    # If inside an 'Instances' block, accumulate lines
    if $in_instances_block; then
        instance_block_content+="$line"$'\n'
    fi

    # End of the 'Instances' block; perform calculation
    if [[ "$line" =~ "}" && $in_instances_block == true ]]; then
        in_instances_block=false
        half_price=$(calculate_price "$instance_block_content")

        # Replace the placeholder in the accumulated content
        instance_block_content=$(echo "$instance_block_content" | sed "s/CALCULATE_MIN_DEDICATED_PRICE_OF_THE_SPECIFIED_SET_IN_THIS_GROUP_DIVIDE_BY_2/$half_price/g")

        # Print the modified content to the output file
        echo "$instance_block_content" >> "$output_file"
    elif ! $in_instances_block; then
        # Print lines outside 'Instances' blocks directly to output
        echo "$line" >> "$output_file"
    fi
done < "$input_file"

echo "Processing completed. Output written to $output_file."
