#!/bin/bash

# Detect if running in Zsh, and switch to Bash emulation if needed
if [ -n "$ZSH_VERSION" ]; then
  emulate -L bash
fi

# Ensure the script is sourced, not executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "This script should be sourced, not executed directly. Please use 'source THISSCRIPT' to run."
    exit 1
fi

# Capture arguments for PEM file and region
pem_file=$1
region=$2
duser="ubuntu"

# Ensure both PEM file and region are provided
if [[ -z "$pem_file" || -z "$region" ]]; then
    echo "Error: You must provide both the PEM file path and the AWS region."
    echo "Usage: source $0 /path/to/pem_file region"
    return 1
fi

# List available clusters in the specified region
echo "Clusters detected in region $region:"
cluster_names=$(pcluster list-clusters --region "$region" | grep clusterName | awk '{print $2}' | cut -d '"' -f 2)

# Check if there are any clusters detected
if [[ -z "$cluster_names" ]]; then
    echo "Error: No clusters found in region $region."
    return 1
fi

# Convert detected cluster names into an array
cluster_array=()
while IFS= read -r cluster_name; do
    cluster_array+=("$cluster_name")
done <<< "$cluster_names"

# Auto-select if there is only one cluster, otherwise prompt for selection
if [[ ${#cluster_array[@]} -eq 1 ]]; then
    selected_cluster=${cluster_array[0]}
    echo "Only one cluster found: $selected_cluster. Auto-selecting it."
else
    echo "Select a cluster name:"
    select selected_cluster in "${cluster_array[@]}"; do
        if [[ -n "$selected_cluster" ]]; then
            echo "You selected: $selected_cluster"
            break
        else
            echo "Invalid selection, please try again."
        fi
    done
fi

cluster_name=$selected_cluster

# Get the public IP address of the cluster's head node
cluster_ip_address=$(pcluster describe-cluster-instances -n "$cluster_name" --region "$region" \
    | grep publicIpAddress | perl -p -e 's/[ |"|,]//g;' | cut -d ':' -f 2)

if [[ -z "$cluster_ip_address" ]]; then
    echo "Error: Could not retrieve the public IP address of the cluster."
    return 1
fi

echo "Cluster $cluster_name's public IP is $cluster_ip_address."
echo " "

# List available PEM files in the .ssh directory
echo "Detected PEM files in ~/.ssh:"
ls -1 ~/.ssh/*.pem

# If PEM file is not provided as an argument, prompt the user
if [[ -z "$pem_file" ]]; then
    echo "Enter the full absolute path to your PEM file:"
    read -r pem_file
fi

# Ensure the PEM file exists
if [[ ! -f "$pem_file" ]]; then
    echo "Error: PEM file '$pem_file' does not exist."
    return 1
fi

# Generate SSH key for the head node user
echo "Generating SSH key on the head node..."
ssh -i "$pem_file" ubuntu@"$cluster_ip_address" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
    "ssh-keygen -q -t rsa -f ~/.ssh/id_rsa -N '' <<< $'\ny' | sudo su - $duser"

# Display the public key and instruct the user to add it to GitHub
echo "Please add the following public SSH key to GitHub:"
ssh -i "$pem_file" ubuntu@"$cluster_ip_address" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
    "cat ~/.ssh/id_rsa.pub"

echo " "
echo "Add this SSH key to your GitHub account under 'Settings -> SSH/GPG Keys'."
echo "You can add this key later if needed (find it in ~/.ssh/id_rsa.pub)."
echo "Sleeping for 15 seconds before proceeding..."
sleep 15

# Clone the Daylily repository on the head node
echo "Cloning the Daylily repository to ~/projects on the head node..."
ssh -i "$pem_file" ubuntu@"$cluster_ip_address" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
    "sudo su - $duser -c 'mkdir -p ~/projects && cd ~/projects && git clone https://github.com/Daylily-Informatics/daylily.git'"

# Initialize and configure the Daylily environment on the head node
echo "Configuring the Daylily environment on the head node..."
ssh -t -i "$pem_file" ubuntu@"$cluster_ip_address" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
    "sudo su - $duser -c 'cd ~/projects/daylily && source dyinit --skip-project-check && source bin/day_build'"

# Run a simple help test for Daylily remotely
echo "Running a simple help test on the head node..."
ssh -t -i "$pem_file" ubuntu@"$cluster_ip_address" -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
    "sudo su - $duser -c 'cd ~/projects/daylily && source dyinit --skip-project-check && source bin/day_activate local && dy-r help'"

# Provide final instructions for SSH access to the head node
echo "You can now SSH into the head node using the following command:"
echo "ssh -i $pem_file ubuntu@$cluster_ip_address"
echo "After logging in, as the 'ubuntu' user, run the following commands:"
echo "  cd ~/projects/daylily"
echo "  source dyinit --h"
echo "  source dyinit --project <PROJECT>"
echo "  dy-a local"
echo "  dy-r help"
echo " "
echo "For non-test uses, re-clone the repository into the /fsx/analysis_results directory."
echo " "
echo "Setup complete. You can now start working with Daylily on the head node."
