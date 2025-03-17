#!/usr/bin/env bash

################################################################################
# Script Name: day_env_installer.sh
# Description: Sets up Miniconda and installs the DAY conda environment.
# Usage:       source ./day_env_installer.sh DAY
#              Provide 'DAY' as the argument to start the installation.
#              If the 'DAY' environment already exists, the script will prompt accordingly.
################################################################################


# Function to display usage information
usage() {
    echo "Usage: source $0 DAY"
    echo "This script installs Miniconda and sets up the DAY conda environment."
    echo "Provide 'DAY' as the argument to start the installation."
    exit 2
}

# Check if the correct argument is provided
DY_ENVNAME="DAY"
if [[ "$1" != "$DY_ENVNAME" ]]; then
    echo "Hello! This is the __ $DY_ENVNAME __ installation script."
    echo ""
    echo "The DAY environment installs the software needed to trigger Snakemake and run the Day (dy-) CLI."
    echo "To run and start the install, provide 'DAY' as the argument."
    echo ""
    echo "Usage: $0 DAY"
    echo ""
    echo "If you have an existing DAY install, you may need to remove it first:"
    echo "  conda env remove -n DAY"
    echo ""
    usage
fi

# Check if the shell is bash
if [[ "$SHELL" != "/bin/bash" ]]; then
    echo "Warning: This script is designed to work with bash."
    echo "Your current shell is $SHELL. Proceeding, but compatibility is not guaranteed."
    sleep 2
fi

# Set the script directory
SCRIPT_DIR=config/day/
echo "Path to environment working directory is $SCRIPT_DIR"

# Create .parallel directory if it doesn't exist
mkdir -p "$HOME/.parallel"

# Function to install Miniconda
install_miniconda() {
    echo "No conda environment detected."
    echo "Installing Miniconda to $CONDA_DIR"

    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/Miniconda3.sh
    bash /tmp/Miniconda3.sh -b -p "$CONDA_DIR"
    rm /tmp/Miniconda3.sh

    source "$CONDA_DIR/etc/profile.d/conda.sh"
    conda init bash
    source "$HOME/.bashrc"
    conda activate

    echo "Conda installation complete."
    
    echo "Adding repo tag pinning stuff"
    conda install -y -n base -c conda-forge yq || (echo 'Failed to install yq' && exit 1)
    
    mkdir -p ~/.config/daylily
    cp config/daylily_cli_global.yaml ~/.config/daylily/daylily_cli_global.yaml

    CONFIG_FILE=~/.config/daylily/daylily_cli_global.yaml
    git_tag=$(yq -r '.daylily.git_tag' "$CONFIG_FILE")
    touch ~/.config/daylily/$git_tag

    cp bin/day-clone "~/miniconda3/condabin/day-clone" || (echo 'Failed to copy day-clone script' && exit 1)
}

# Detect or install conda
if command -v conda &> /dev/null; then
    CONDA_DIR="$(dirname "$(dirname "$(which conda)")")"
    echo "Conda detected at $CONDA_DIR"
else
    CONDA_DIR="$HOME/miniconda3"
    install_miniconda
fi

# Ensure conda is initialized
source "$CONDA_DIR/etc/profile.d/conda.sh"

# Update Conda Config
conda config --add channels conda-forge
conda config --add channels bioconda

conda config --set channel_priority strict || echo 'Failed to set conda priority to strict'
conda config --set repodata_threads 10 || echo 'Failed to set repodata_threads'
conda config --set verify_threads 4 || echo 'Failed to set verify_threads'
conda config --set execute_threads 4 || echo 'Failed to set execute_threads'
conda config --set always_yes yes || echo 'Failed to set always_yes'
conda config --set default_threads 10 || echo 'Failed to set default_threads'

# Check if the DAY environment already exists
if conda env list | grep -q "^$DY_ENVNAME\s"; then
    echo ""
    echo "It appears you have a DAY environment already."
    echo "You may need to manually remove the conda env dir for DAY and try again."
    echo "To remove the environment, run:"
    echo "  conda env remove -n DAY"
    exit 0
else
    conda install -y -n base -c conda-forge yq || echo 'Failed to install yq'
    cp bin/day-clone "~/miniconda3/condabin/day-clone" || echo 'Failed to copy day-clone script to ~/minconda3/condabin'
    echo "Installing DAY environment..."
    # Create the DAY environment
    if conda env create -n "$DY_ENVNAME" -f "$SCRIPT_DIR/day.yaml"; then
        echo "DAY environment created successfully."
        echo ""
        echo "Try the following commands to get started:"
        echo "  source dyinit --project <PROJECT>"
        echo "  dy-a local"
        echo "  dy-r help"
    else
        echo "Failed to create DAY environment."
        exit 1
    fi
fi

echo ""
echo "Installation complete."
echo "Please log out and log back in, then run:"
echo "  source dyinit --project <PROJECT>"
echo "  dy-a local"
echo "  dy-r help"

exit 0
