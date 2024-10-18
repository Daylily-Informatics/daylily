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
    return parser.parse_args()

def get_spot_price(instance_type, region):
    """Query AWS to get the current spot price of an instance type."""
    try:
        result = subprocess.check_output([
            "aws", "ec2", "describe-spot-price-history",
            "--instance-types", instance_type,
            "--region", region,
            "--product-description", "Linux/UNIX",
            "--query", "SpotPriceHistory[0].SpotPrice",
            "--output", "text"
        ]).decode("utf-8").strip()
        return float(result)
    except subprocess.CalledProcessError as e:
        print(f"Error querying spot price for {instance_type}: {e}")
        return None

def process_instance_group(group, region):
    """Process each compute resource group to calculate the average spot price."""
    for resource in group.get('ComputeResources', []):
        spot_prices = []
        for instance in resource.get('Instances', []):
            instance_type = instance.get('InstanceType')
            spot_price = get_spot_price(instance_type, region)
            if spot_price is not None:
                spot_prices.append(spot_price)

        if spot_prices:
            # Calculate the average spot price
            avg_spot_price = statistics.mean(spot_prices)
            # Insert the SpotPrice tag
            resource['SpotPrice'] = round(avg_spot_price, 4)

def main():
    args = parse_arguments()

    # Load the input YAML configuration file
    with open(args.input, 'r') as f:
        config = yaml.safe_load(f)

    # Get the AWS region from the configuration
    region = config.get("Region", "us-west-2")

    # Deepcopy the config to avoid modifying the original during iteration
    new_config = deepcopy(config)

    # Iterate over all slurm queues and process each one
    for queue in new_config.get('Scheduling', {}).get('SlurmQueues', []):
        process_instance_group(queue, region)

    # Write the updated configuration to the output file
    with open(args.output, 'w') as f:
        yaml.dump(new_config, f, sort_keys=False)

    print(f"Updated configuration saved to {args.output}.")

if __name__ == "__main__":
    main()
