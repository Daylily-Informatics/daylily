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
echo -e "Aligner\tAvg_Minutes\tvCPU_Min\tCost_USD"

# Initialize the summary string
summary=""

# Process the file to filter by 'alN' rule, calculate avg runtime, vCPU minutes, and cost
summary=$(awk -F '\t' -v cost_per_vcpu_min="$vcpu_cost_per_min" '
    NR > 1 && $3 ~ /alN/ {  # Skip header and match "alN" in the rule column
        aligner[$3] += $4;  # Accumulate the "s" column (4th column) for each aligner
        count[$3]++;  # Count occurrences for each aligner
    }
    END {
        summary_str = "";
        for (a in aligner) {
            avg_minutes = aligner[a] / count[a] / 60;  # Average runtime in minutes
            vcpu_minutes = avg_minutes * 192;  # Calculate vCPU minutes
            cost = vcpu_minutes * cost_per_vcpu_min;  # Calculate the cost
            printf "%s\t%.2f\t%.2f\t%.2f\n", a, avg_minutes, vcpu_minutes, cost;

            # Append to summary string
            summary_str = summary_str (summary_str ? ",  " : "") a ": \\$" sprintf("%.2f", cost);
        }
        print summary_str;
    }
' "$input_file")

# Split the summary into table output and the environment variable
export ALNR_SUMMARY_COST=$(echo "$summary" | tail -n 1)

# Output the results
echo "$summary" | head -n -1
echo
echo "ALNR_SUMMARY_COST: $ALNR_SUMMARY_COST"
