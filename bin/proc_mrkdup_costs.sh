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
echo -e "Category\tAvg_Minutes\tAvg_vCPU_Min\tAvg_Cost_USD"

# Process the file to filter rows where the final column is "mrkdup" and calculate averages
awk -F '\t' -v cost_per_vcpu_min="$vcpu_cost_per_min" '
    NR > 1 && $NF == "mrkdup" {  # Skip header and match rows where the final column is "mrkdup"
        total_s += $4;  # Accumulate the "s" column (4th column)
        count++;  # Count rows
    }
    END {
        if (count > 0) {
            avg_minutes = (total_s / count) / 60;  # Average runtime in minutes
            avg_vcpu_minutes = avg_minutes * 192;  # Calculate average vCPU minutes
            avg_cost = avg_vcpu_minutes * cost_per_vcpu_min;  # Calculate average cost
            printf "mrkdup\t%.2f\t%.2f\t%.2f\n", avg_minutes, avg_vcpu_minutes, avg_cost;
        } else {
            printf "mrkdup\t0.00\t0.00\t0.00\n";
        }
    }
' "$input_file"

# Export environment variables for average values with 2 significant digits
export MRKDUP_AVG_MINUTES=$(awk -F '\t' 'NR > 1 && $NF == "mrkdup" {total_s += $4; count++} END {if (count > 0) printf "%.2f", (total_s / count) / 60; else print 0}' "$input_file")
export MRKDUP_AVG_VCPU_MIN=$(awk "BEGIN {printf \"%.2f\", $MRKDUP_AVG_MINUTES * 192}")
export MRKDUP_AVG_COST=$(awk "BEGIN {printf \"%.2f\", $MRKDUP_AVG_VCPU_MIN * $vcpu_cost_per_min}")

# Output environment variables
echo
echo "MRKDUP_AVG_MINUTES: $MRKDUP_AVG_MINUTES"
echo "MRKDUP_AVG_VCPU_MIN: $MRKDUP_AVG_VCPU_MIN"
echo "MRKDUP_AVG_COST: $MRKDUP_AVG_COST"

