#!/usr/bin/env python3

import yaml
import subprocess
import statistics
import argparse
from math import isnan, fsum
from tabulate import tabulate
from colr import color

# Define neon teal, neon orange, and bright gray colors using ANSI codes
NEON_TEAL = (0, 255, 222)
NEON_ORANGE = (255, 140, 0)
BRIGHT_GRAY = (211, 211, 211)

def parse_arguments():
    """Parse command-line arguments for input YAML, output TSV, and zones. """
    parser = argparse.ArgumentParser(description="Look up spot prices for instances in the i192 partition.")
    parser.add_argument("-i", "--input", default="config/day_cluster/prod_cluster.yaml", help="Path to the input YAML configuration file. default( config/day_cluster/prod_cluster.yaml)")
    parser.add_argument("-o", "--output", required=True, help="Path to the output TSV file.")
    parser.add_argument(
        "--zones",
        default="us-west-2a,us-west-2b,us-west-2c,us-west-2d,us-east-1a,us-east-1b,us-east-1c,us-east-1d,us-east-2a,us-east-2b,us-east-2c,us-west-1a,us-west-1a,us-west-1c",
        help="Comma-separated list of availability zones (default us-west-2a,us-west-2b,us-west-2c,us-west-2d,us-east-1a,us-east-1b,us-east-1c,us-east-1d,us-east-2a,us-east-2b,us-east-2c,us-west-1a,us-west-1c). HOWEVER, this is the current list of all AZs with the 192vcpu instances we need: ap-south-1a,ap-south-1b,ap-south-1c,ap-northeast-1a,ap-northeast-1c,ap-northeast-1d,ap-northeast-2a,ap-northeast-2b,ap-northeast-2c,ap-northeast-2d,ap-southeast-1a,ap-southeast-1b,ap-southeast-1c,ap-southeast-2a,ap-southeast-2b,ap-southeast-2c,ca-central-1a,ca-\
central-1b,ca-central-1d,eu-central-1a,eu-central-1b,eu-central-1c,eu-north-1a,eu-north-1b,eu-north-1c,eu-west-1a,eu-west-1b,eu-west-1c,eu-west-2a,eu-west-2b,eu-west-2c,eu-west-3a,eu-west-3b,eu-west-3c,sa-east-1a,sa-east-1b,sa-east-1c,us-east-1a,us-east-1b,us-east-1c,us-east-1d,us-east-1f,us-east-2a,us-east-2b,us-east-2c,us-west-1a,us-west-1c,us-west-2a,us-west-2b,us-west-2c,us-west-2d")
    parser.add_argument(
        "--approx-spot-hours-per-30x-genome",
        type=float,
        default=3,
        help="Estimated spot instance hours required for processing a 30x genome. (default: 3)"
    )
    return parser.parse_args()

def get_region_from_zone(zone):
    """Extract the region from the zone by trimming the trailing letter."""
    return zone[:-1]

def get_spot_price(instance_type, zone):
    """Query AWS to get the current spot price of an instance type in a specific AZ."""
    region = get_region_from_zone(zone)
    try:
        result = subprocess.check_output([
            "aws", "ec2", "describe-spot-price-history",
            "--instance-types", instance_type,
            "--availability-zone", zone,
            "--product-description", "Linux/UNIX",
            "--query", "SpotPriceHistory[0].SpotPrice",
            "--output", "text", "--region", region
        ]).decode("utf-8").strip()
        try:
            return float(result)
        except Exception as e:
            print(f"Error converting spot price for {instance_type} in {zone}: {e}")
            return 100
    except subprocess.CalledProcessError as e:
        print(f"Error querying spot price for {instance_type} in {zone}: {e}")
        return None

def harmonic_mean(data):
    """Calculate the harmonic mean of a list of numbers."""
    if not data or any(d == 0 for d in data):
        return float('nan')
    return len(data) / fsum(1 / d for d in data)

def stability_metric(prices):
    """Calculate a simple stability metric (price variability)."""
    return max(prices) - min(prices) if prices else float('nan')

def apply_color(value, best, worst):
    """Colorize a value based on whether it's the best, worst, or identical."""
    if value == best:
        return color(f"{value:.4f}", fore=NEON_TEAL)
    elif value == worst:
        return color(f"{value:.4f}", fore=NEON_ORANGE)
    elif best == worst:
        return color(f"{value:.4f}", fore=BRIGHT_GRAY)
    else:
        return f"{value:.4f}"

def extract_i192_instances(config):
    """Extract the instance types from the i192 partition."""
    instance_types = []
    for queue in config.get('Scheduling', {}).get('SlurmQueues', []):
        if queue.get('Name') == 'i192':
            for resource in queue.get('ComputeResources', []):
                for instance in resource.get('Instances', []):
                    instance_types.append(instance.get('InstanceType'))
    return instance_types

def collect_spot_prices(instance_types, zones):
    """Collect spot prices and availability for each instance type in specified AZs."""
    spot_data = {}
    availability_data = {zone: 0 for zone in zones}

    for instance_type in instance_types:
        spot_data[instance_type] = {}
        for zone in zones:
            spot_price = get_spot_price(instance_type, zone)
            spot_data[instance_type][zone] = spot_price if spot_price is not None else float('nan')
            if spot_price is not None:
                availability_data[zone] += 1

    return spot_data, availability_data

def calculate_statistics(spot_data, zones, hours_per_genome):
    """Calculate statistics for each zone."""
    zone_stats = []
    for zone in zones:
        prices = [data.get(zone, float('nan')) for data in spot_data.values() if not isnan(data.get(zone))]
        total_with_pricing = len(prices)
        harmonic_mean_price = harmonic_mean(prices)
        stability = stability_metric(prices)
        avg_lowest_3 = statistics.mean(sorted(prices)[:3]) if len(prices) >= 3 else statistics.mean(prices)
        est_cost_per_genome = avg_lowest_3 * hours_per_genome if avg_lowest_3 else float('nan')

        if prices:
            zone_stats.append({
                'Zone': zone,
                'Total Instances Considered': len(spot_data),
                'Instances with Pricing': total_with_pricing,
                'Median Spot Price': statistics.median(prices),
                'Mean Spot Price': statistics.mean(prices),
                'Min Spot Price': min(prices),
                'Max Spot Price': max(prices),
                'Avg of Lowest 3 Prices': avg_lowest_3,
                'Harmonic Mean Price': harmonic_mean_price,
                'Estimated EC2 Cost / Genome': est_cost_per_genome,
                'Stability (Max-Min Spread)': stability
            })

    return sorted(zone_stats, key=lambda x: x['Harmonic Mean Price'])

def get_best_worst(stats, key):
    """Extract the best and worst values for a specific key."""
    values = [z[key] for z in stats]
    return min(values), max(values)

def display_statistics(zone_stats):
    """Display statistics in a colorized table with a title."""
    headers = [
        "Zone", "Total Instances Considered", "Instances with Pricing", "Median Spot Price", 
        "Mean Spot Price", "Min Spot Price", "Max Spot Price", "Avg of Lowest 3 Prices", 
        "Harmonic Mean Price", "Estimated EC2 Cost / Genome", "Stability (Max-Min Spread)"
    ]

    # Get best and worst values for each numeric column
    best_worst = {key: get_best_worst(zone_stats, key) for key in headers[1:]}

    # Prepare colorized table rows
    table = []
    for z in zone_stats:
        row = [
            color(z['Zone'], fore='blue'),
            z['Total Instances Considered'],
            z['Instances with Pricing'],
            apply_color(z['Median Spot Price'], *best_worst['Median Spot Price']),
            apply_color(z['Mean Spot Price'], *best_worst['Mean Spot Price']),
            apply_color(z['Min Spot Price'], *best_worst['Min Spot Price']),
            apply_color(z['Max Spot Price'], *best_worst['Max Spot Price']),
            apply_color(z['Avg of Lowest 3 Prices'], *best_worst['Avg of Lowest 3 Prices']),
            apply_color(z['Harmonic Mean Price'], *best_worst['Harmonic Mean Price']),
            apply_color(z['Estimated EC2 Cost / Genome'], *best_worst['Estimated EC2 Cost / Genome']),
            apply_color(z['Stability (Max-Min Spread)'], *best_worst['Stability (Max-Min Spread)'])
        ]
        table.append(row)

    print(color("Spot Price Statistics", fore='green', style='bright'))
    print(tabulate(table, headers=headers, floatfmt=".4f"))

def main():
    args = parse_arguments()

    with open(args.input, 'r') as f:
        config = yaml.safe_load(f)

    instance_types = extract_i192_instances(config)
    zones = args.zones.split(',')

    spot_data, availability_data = collect_spot_prices(instance_types, zones)

    zone_stats = calculate_statistics(spot_data, zones, args.approx_spot_hours_per_30x_genome)
    display_statistics(zone_stats)

if __name__ == "__main__":
    main()
