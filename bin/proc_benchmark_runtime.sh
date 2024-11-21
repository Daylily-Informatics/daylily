echo "SOURCE ME"
echo ""

# Path to the input file
input_file=$1

# Extract the 's' column, sum its values, and divide by 60 to get the total runtime in minutes
total_runtime=$(awk -F '\t' 'NR>1 {sum += $4} END {print sum / 60}' "$input_file")
export TOTAL_RUNTIME_MIN=$total_runtime

# Calculate total vCPU minutes (runtime in minutes * 192)
total_vcpu_min=$(awk "BEGIN {print $total_runtime * 192}")
export TOTAL_VCPU_MIN=$total_vcpu_min

# Output the results
echo "Total Runtime (min): $TOTAL_RUNTIME_MIN"
