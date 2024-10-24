#!/usr/bin/env python3
import yaml
import subprocess
import statistics
import argparse
from copy import deepcopy

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Insert average spot prices into config.yaml.")
    parser.add_argument("-i", "--input", required=True, help="Input YAML configuration file path.")
    parser.add_argument("-o", "--output", required=True, help="Output YAML configuration file path.")
    parser.add_argument("--az", required=True, help="availability zone")
    parser.add_argument("--profile", help="AWS CLI profile to use.")
    return parser.parse_args()

def run_aws_command(command):
    """Run an AWS CLI command and return the result as a string."""
    try:
        return subprocess.check_output(command, text=True).strip()
    except subprocess.CalledProcessError as e:
        print(f"Error executing {' '.join(command)}: {e}")
        return None

def get_availability_zone(subnet_id, profile):
    """Retrieve the AZ for the given subnet ID."""
    command = [
        "aws", "ec2", "describe-subnets", 
        "--subnet-ids", subnet_id, 
        "--query", "Subnets[0].AvailabilityZone", 
        "--output", "text", 
        "--profile", profile,
    ]
    return run_aws_command(command)

def get_spot_price(instance_type, az, profile):
    """Retrieve the spot price for an instance type in a specific AZ."""
    region = az[:-1]  # Derive region from AZ
    command = [
        "aws", "ec2", "describe-spot-price-history", 
        "--instance-types", instance_type, 
        "--availability-zone", az, 
        "--region", region,
        "--product-description", "Linux/UNIX", 
        "--query", "SpotPriceHistory[0].SpotPrice", 
        "--output", "text",
        "--profile", profile
    ]
    result = run_aws_command(command)
    return float(result) if result else None

def calculate_average_spot_price(resource, az, profile):
    """Calculate the average spot price for all instances in the compute resource."""
    spot_prices = [
        get_spot_price(instance.get('InstanceType'), az, profile)
        for instance in resource.get('Instances', [])
        if get_spot_price(instance.get('InstanceType'), az, profile) is not None
    ]
    if spot_prices:
        resource['SpotPrice'] = round(statistics.mean(spot_prices), 4)

def process_slurm_queues(config, az, profile):
    """Process all Slurm queues to calculate and update spot prices."""
    for queue in config.get('Scheduling', {}).get('SlurmQueues', []):
        for resource in queue.get('ComputeResources', []):
            calculate_average_spot_price(resource, az, profile)

def main():
    args = parse_arguments()

    az = args.az
    profile = args.profile
    with open(args.input, 'r') as f:
        config = yaml.safe_load(f)

    new_config = deepcopy(config)
    process_slurm_queues(new_config, az, profile)

    with open(args.output, 'w') as f:
        yaml.dump(new_config, f, sort_keys=False)

    print(f"Updated configuration saved to {args.output}.")

if __name__ == "__main__":
    main()
