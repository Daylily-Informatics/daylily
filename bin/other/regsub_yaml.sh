#!/bin/bash

cluster_init_values=$1
cluster_cfg_yaml=$2

# Check if the files are provided as arguments
if [[ -z "$cluster_init_values" || -z "$cluster_cfg_yaml" ]]; then
    echo "Usage: $0 <cluster_init_values.txt> <cluster_cfg_yaml.yaml>"
    exit 1
fi

# Check if the files exist
if [[ ! -f "$cluster_init_values" ]]; then
    echo "Error: '$cluster_init_values' not found!"
    exit 1
fi

if [[ ! -f "$cluster_cfg_yaml" ]]; then
    echo "Error: '$cluster_cfg_yaml' not found!"
    exit 1
fi


# Iterate over each line in cluster_init_values.txt
while IFS='=' read -r key value; do
    # Skip empty lines or lines without key-value pairs
    [[ -z "$key" || -z "$value" ]] && continue

    # Escape special characters in the value to prevent sed issues
    escaped_value=$(printf '%s\n' "$value" | sed -e 's/[\/&]/\\&/g')

    # Substitute occurrences of the key in the cluster config YAML file
    sed -i "" "s/$key/$escaped_value/g" "$cluster_cfg_yaml"
done < "$cluster_init_values"

echo "" 
echo "Substitutions completed in $cluster_cfg_yaml."
echo ""
