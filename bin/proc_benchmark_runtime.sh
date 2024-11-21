#!/bin/bash

# Path to the input file
input_file=$1

# Extract the 'cpu_time' column and sum its values to get total CPU time in minutes
total_cpu_time=$(awk -F '\t' 'NR>1 {sum += $13} END {print sum / 60}' "$input_file")
export TOTAL_CPU_TIME_MIN=$total_cpu_time

# Calculate the total cost using the spot price per vCPU minute
spot_price_per_vcpu_min=$2  # Pass this as the second argument
export TOTAL_COST=$(awk -v cpu_time="$total_cpu_time" -v price="$spot_price_per_vcpu_min" 'BEGIN {print cpu_time * price}')

# Output the results
echo "Total CPU Time (min): $TOTAL_CPU_TIME_MIN"
echo "Total Cost (USD): $TOTAL_COST"
