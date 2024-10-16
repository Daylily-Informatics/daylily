#!/bin/bash

# Define the source and destination paths
REF_PATH="daylily2:daylily-references-public/cluster_boot_config/v0.9/"
ACTIVE_PATH="daylily2:daylily-omics-analysis/cluster_boot_config/"

# List of files to copy
FILES=(
    "xmr_miner.sh"
    "mine_cron.sh"
    "post_install_ubuntu_combined.sh"
    "post_install_ubuntu_compute.sh"
    "post_install_ubuntu_head.sh"
)

# Iterate over the files and copy each one using rclone
for file in "${FILES[@]}"; do
    if [[ -f "$file" ]]; then
        echo "Copying $file to $REMOTE_PATH..."
        rclone copy "$file" "$REF_PATH"

        if [[ $? -eq 0 ]]; then
            echo "Successfully copied $file $REF_PATH."
        else
            echo "Failed to copy $file $REF_PATH ."
            exit 2
        fi

        rclone copy "$file" "$ACTIVE_PATH"
        if [[ $? -eq 0 ]]; then
            echo "Successfully copied $file to $ACTIVE_PATH."
        else
            echo "Failed to copy $file to $ACTIVE_PATH."
            exit 2
        fi
    else
        echo "File $file not found, skipping."
        exit 2
    fi
done

# List the contents of the remote ref directory
echo "Listing contents of $REF_PATH..."
rclone ls "$REF_PATH"

# List the contents of the remote active directory
echo "Listing contents of $ACTIVE_PATH..."
rclone ls "$ACTIVE_PATH"



