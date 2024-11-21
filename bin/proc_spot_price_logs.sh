echo "SOURCE ME"

echo ""
echo ""

# Read and process the log file output
data=$(cat /fsx/scratch/*price.log)

# Extract unique instance types
instance_types=$(echo "$data" | awk -F ', ' '{print $3}' | sort | uniq)
export INSTANCE_TYPES="$instance_types"


# Extract instance types in alpha order on one line
instance_types_line=$(echo "$instance_types" | tr '\n' ' ' | sed 's/ $//')
export INSTANCE_TYPES_LINE=$(echo "$instance_types_line" | sed 's/Instance type://g')

# Extract all spot prices
spot_prices=$(echo "$data" | grep -oP 'Spot price: \K[0-9.]+')

# Calculate average spot price
average_spot_price=$(echo "$spot_prices" | awk '{sum += $1; count++} END {print sum / count}')
export AVERAGE_SPOT_PRICE=$(printf "%.2f" "$average_spot_price")

# Calculate median spot price
sorted_prices=$(echo "$spot_prices" | sort -n)
price_count=$(echo "$sorted_prices" | wc -l)
median_spot_price=$(echo "$sorted_prices" | awk -v n=$price_count '
  BEGIN { if (n % 2 == 0) { mid1 = n / 2; mid2 = mid1 + 1 } else { mid1 = (n + 1) / 2; mid2 = mid1 } }
  { if (NR == mid1) val1 = $1; if (NR == mid2) val2 = $1 }
  END { print (val1 + val2) / 2 }')
export MEDIAN_SPOT_PRICE=$(printf "%.2f" "$median_spot_price")

# Calculate vCPU cost per minute
vcpu_cost_per_min=$(echo "$average_spot_price" | awk '{print ($1 / 192) / 60}')
export VCPU_COST_PER_MIN=$vcpu_cost_per_min

# Output results
echo "Unique Instance Types:"
echo "$INSTANCE_TYPES_LINE"
echo "$INSTANCE_TYPES"
echo
echo "Average Spot Price (USD/hour): $AVERAGE_SPOT_PRICE"
echo "Median Spot Price (USD/hour): $MEDIAN_SPOT_PRICE"
echo "vCPU Cost (USD/min): $VCPU_COST_PER_MIN"