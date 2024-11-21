#!/bin/bash

# Check if the input file and vCPU cost are provided
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <input_file> <vcpu_cost_per_min>"
    exit 1
fi

# Input file and vCPU cost per minute
input_file=$1
vcpu_cost_per_min=$2

# Output header
echo -e "Category\tTotal_Minutes\tvCPU_Min\tTotal_Cost_USD"

# Initialize the total variables
total_minutes=0
total_vcpu_min=0
total_cost=0

# Process the file to filter rows matching 'mrkdup' and calculate totals
awk -F '\t' -v cost_per_vcpu_min="$vcpu_cost_per_vcpu_min" '
    NR > 1 && $3 ~ /mrkdup/ {  # Skip header and match "mrkdup" in the rule column
        total_s += $4;  # Accumulate the "s" column (4th column)
    }
    END {
        total_minutes = total_s / 60;  # Convert total seconds to minutes
        vcpu_minutes = total_minutes * 192;  # Calculate vCPU minutes
        total_cost = vcpu_minutes * cost_per_vcpu_min;  # Calculate total cost
        printf "mrkdup\t%.2f\t%.2f\t%.2f\n", total_minutes, vcpu_minutes, total_cost;
    }
' "$input_file"

# Example: Store in environment variables if needed
export MRKDUP_TOTAL_MINUTES=$(awk -F '\t' 'NR > 1 && $3 ~ /mrkdup/ {total_s += $4} END {print total_s / 60}' "$input_file")
export MRKDUP_VCPU_MIN=$(awk "BEGIN {print $MRKDUP_TOTAL_MINUTES * 192}")
export MRKDUP_TOTAL_COST=$(awk "BEGIN {print $MRKDUP_VCPU_MIN * $vcpu_cost_per_min}")

# Output environment variables
echo
echo "MRKDUP_TOTAL_MINUTES: $MRKDUP_TOTAL_MINUTES"
echo "MRKDUP_VCPU_MIN: $MRKDUP_VCPU_MIN"
echo "MRKDUP_TOTAL_COST: $MRKDUP_TOTAL_COST"
