#!/usr/bin/env python3
import yaml
import subprocess
import statistics
import argparse
from copy import deepcopy

def parse_arguments():
    """Parse command-line arguments for input and output files."""
    parser = argparse.ArgumentParser(description="Calculate and insert average spot prices into config.yaml.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input YAML configuration file.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output YAML configuration file.")
    parser.add_argument("--subnet-id", required=True, help="Subnet ID to determine the AZ (e.g., subnet-1234abcd).")
    return parser.parse_args()

def get_availability_zone(subnet_id):
    """Retrieve the availability zone for the specified subnet ID."""
    try:
        az = subprocess.check_output([
            "aws", "ec2", "describe-subnets",
            "--subnet-ids", subnet_id,
            "--query", "Subnets[0].AvailabilityZone",
            "--output", "text"
        ]).decode("utf-8").strip()
        return az
    except subprocess.CalledProcessError as e:
        print(f"Error retrieving AZ for subnet {subnet_id}: {e}")
        return None

def get_spot_price(instance_type, az):
    """Query AWS to get the current spot price of an instance type in a specific AZ."""
    region = az[:-1]  # Extract region from AZ
    try:
        result = subprocess.check_output([
            "aws", "ec2", "describe-spot-price-history",
            "--instance-types", instance_type,
            "--availability-zone", az,
            "--product-description", "Linux/UNIX",
            "--query", "SpotPriceHistory[0].SpotPrice",
            "--output", "text"
        ]).decode("utf-8").strip()
        return float(result)
    except subprocess.CalledProcessError as e:
        print(f"Error querying spot price for {instance_type} in {az}: {e}")
        return None

def process_instance_group(group, az):
    """Process each compute resource group to calculate the average spot price."""
    for resource in group.get('ComputeResources', []):
        spot_prices = []
        for instance in resource.get('Instances', []):
            instance_type = instance.get('InstanceType')
            spot_price = get_spot_price(instance_type, az)
            if spot_price is not None:
                spot_prices.append(spot_price)

        if spot_prices:
            # Calculate the average spot price
            avg_spot_price = statistics.mean(spot_prices)
            # Insert the SpotPrice tag
            resource['SpotPrice'] = round(avg_spot_price, 4)

def main():
    args = parse_arguments()

    # Detect the AZ from the provided subnet ID
    az = get_availability_zone(args.subnet_id)
    if not az:
        print("Error: Could not determine the availability zone from the subnet ID.")
        return

    # Load the input YAML configuration file
    with open(args.input, 'r') as f:
        config = yaml.safe_load(f)

    # Deepcopy the config to avoid modifying the original during iteration
    new_config = deepcopy(config)

    # Iterate over all slurm queues and process each one
    for queue in new_config.get('Scheduling', {}).get('SlurmQueues', []):
        process_instance_group(queue, az)

    # Write the updated configuration to the output file
    with open(args.output, 'w') as f:
        yaml.dump(new_config, f, sort_keys=False)

    print(f"Updated configuration saved to {args.output}.")

if __name__ == "__main__":
    main()
