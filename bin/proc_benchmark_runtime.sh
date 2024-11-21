#!/bin/bash

echo "SOURCE ME"
echo ""

# Check if input arguments are provided
if [[ $# -lt 2 ]]; then
  echo "Usage: source script.sh <input_file> <spot_price_per_vcpu_min>"
  return 1
fi

# Path to the input file and spot price per vCPU minute
input_file=$1
spot_price_per_vcpu_min=$2

# Extract the 's' column, sum its values, and divide by 60 to get the total runtime in minutes
total_runtime=$(awk -F '\t' 'NR>1 {sum += $4} END {print sum / 60}' "$input_file")
export TOTAL_RUNTIME_MIN=$total_runtime

# Calculate total vCPU minutes (runtime in minutes * 192)
total_vcpu_min=$(awk -v runtime="$total_runtime" 'BEGIN {print runtime * 192}')
export TOTAL_VCPU_MIN=$total_vcpu_min

# Calculate the total cost based on the spot price
total_cost=$(awk -v vcpu_min="$total_vcpu_min" -v price="$spot_price_per_vcpu_min" 'BEGIN {print vcpu_min * price}')
export TOTAL_COST=$total_cost

# Output the results
echo "Total Runtime (min): $TOTAL_RUNTIME_MIN"
echo "Total vCPU Minutes: $TOTAL_VCPU_MIN"
echo "Total Cost (USD): $TOTAL_COST"
