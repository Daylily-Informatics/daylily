#!/usr/bin/env python3

import yaml
import subprocess
import statistics
import argparse
from math import isnan, fsum
from tabulate import tabulate
from colr import color
import re
# Define ANSI colors for different modes
COLORS = {
    'dark': {
        'zone': (135, 206, 235),  
        'teal': (0, 255, 222),    
        'orange': (255, 140, 0),  
        'gray': (211, 211, 211),  
        'yellow': (255, 255, 0),  
        'highlight': (255, 255, 0)  
    },
    'light': {
        'zone': (0, 0, 139),      
        'teal': (0, 128, 128),    
        'orange': (255, 69, 0),   
        'gray': (169, 169, 169),  
        'yellow': (255, 215, 0),  
        'highlight': (255, 215, 0)
    },
    'tacky': {  
        'zone': (255, 0, 255),    
        'teal': (0, 255, 0),      
        'orange': (0, 255, 255),  
        'gray': (255, 20, 147),   
        'yellow': (255, 255, 0),  
        'highlight': (255, 255, 0)
    },
    'neon_rave': {
        'zone': (255, 0, 255),    
        'teal': (57, 255, 20),    
        'orange': (255, 0, 0),    
        'gray': (200, 200, 200),  
        'yellow': (0, 255, 255),  
        'highlight': (0, 255, 255)
    }
}

def apply_color(value, best, worst, highlight=False, mode='dark'):
    """Colorize values based on thresholds using the selected mode."""
    colors = COLORS[mode]
    # Choose the appropriate color based on value comparison
    if highlight:
        color_value = colors['highlight']
    elif value == best:
        color_value = colors['teal']
    elif value == worst:
        color_value = colors['orange']
    else:
        color_value = colors['gray']
    
    # Return the value wrapped in the color using the colr module
    return color(f"{value:.4f}", fore=color_value)

TRANSFER_RATES = {
    'internet': 0.09,  # $/GB to internet
    'within_region': 0.01,  # $/GB within region
    'cross_region': 0.02  # $/GB across regions
}

STORAGE_COSTS = {'standard': 0.023, 'glacier': 0.004}  # $/GB/month

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Calculate genomic workflow costs.")
    parser.add_argument("-i", "--input", default="config/day_cluster/prod_cluster.yaml", help="YAML config path.")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file path.")
    parser.add_argument("--profile", help="AWS CLI profile.")
    parser.add_argument("--partition", default="i192", help="Partition name (default: i192).")
    parser.add_argument("--zones", default="us-west-2a,us-west-2b,us-west-2c", help="Comma-separated zones.")
    parser.add_argument("--x-coverage", type=int, default=30.1,required=True, help="Coverage for analysis (e.g., 30x).")
    parser.add_argument("--vcpu-min-per-x", type=float, default=3.7,required=True, help="vCPU-minutes per X coverage.")
    
    # File size flags
    parser.add_argument("--bam-size-per-x", type=float, default=1.0,help="BAM size per X coverage (GB).")
    parser.add_argument("--cram-size-per-x", type=float, default=1.0,help="CRAM size per X coverage (GB).")
    parser.add_argument("--vcf-size-per-x", type=float, default=1.0,help="VCF size per X coverage (GB).")
    parser.add_argument("--gvcf-size-per-x", type=float, default=1.0,help="gVCF size per X coverage (GB).")
    parser.add_argument("--bcf-size-per-x", type=float, default=1.0,help="BCF size per X coverage (GB).")
    parser.add_argument("--gbcf-size-per-x", type=float, default=1.0,help="gBCF size per X coverage (GB).")
    parser.add_argument("--mode", choices=['dark', 'light', 'neon_rave', 'tacky'], default='dark', help="Color mode (default: dark).")
    parser.add_argument("--input-data-size-per-x", type=float, default=1.0, help="Input FASTQ size per X coverage (GB).")
    return parser.parse_args()

def calculate_vcpu_mins(args):
    """Calculate total vCPU-minutes based on coverage."""
    return args.x_coverage * args.vcpu_min_per_x


def calculate_file_sizes(args):
    bam = args.x_coverage * args.bam_size_per_x
    cram = args.x_coverage * args.cram_size_per_x
    vcf = args.x_coverage * args.vcf_size_per_x
    gvcf = args.x_coverage * args.gvcf_size_per_x
    bcf = args.x_coverage * args.bcf_size_per_x
    gbcf = args.x_coverage * args.gbcf_size_per_x
    fastq = args.x_coverage * args.input_data_size_per_x
    return bam, cram, vcf, gvcf, bcf, gbcf, fastq


def calculate_storage_costs(size_gb, storage_type):
    """Calculate monthly storage costs for a given size."""
    return size_gb * STORAGE_COSTS[storage_type]

def calculate_transfer_costs(size_gb, transfer_type):
    """Calculate transfer costs based on the type."""
    return size_gb * TRANSFER_RATES[transfer_type]

def harmonic_mean(data):
    """Calculate the harmonic mean of a list of numbers."""
    if not data or any(d == 0 for d in data):
        return float('nan')  # Avoid division by zero
    return len(data) / fsum(1 / d for d in data)

def stability_metric(prices):
    """Calculate the volatility as the difference between the max and min prices."""
    if not prices:
        return float('nan')
    return max(prices) - min(prices)

def collect_spot_prices(instance_types, zones, profile):
    """Fetch spot prices for each instance type in the given zones."""
    spot_data = {}
    for instance_type in instance_types:
        spot_data[instance_type] = {}
        for zone in zones:
            price = get_spot_price(instance_type, zone, profile)
            spot_data[instance_type][zone] = price if not isnan(price) else float('nan')
    return spot_data

def get_spot_price(instance_type, zone, profile):
    region = zone[:-1]
    try:
        result = subprocess.check_output([
            "aws", "ec2", "describe-spot-price-history",
            "--instance-types", instance_type,
            "--availability-zone", zone,
            "--product-description", "Linux/UNIX",
            "--query", "SpotPriceHistory[0].SpotPrice",
            "--output", "text", "--region", region,
            "--profile", profile,
        ]).decode().strip()
        return float(result)
    except subprocess.CalledProcessError as e:
        print(f"Failed to fetch price for {zone}: {e}")
        return float('nan')


def extract_instances(config, partition_name):
    """Extract instance types from the given partition in the YAML configuration."""
    instances = []
    for queue in config.get('Scheduling', {}).get('SlurmQueues', []):
        if queue.get('Name') == partition_name:
            for resource in queue.get('ComputeResources', []):
                for instance in resource.get('Instances', []):
                    instances.append(instance.get('InstanceType'))
    if not instances:
        raise ValueError(f"No instances found for partition: {partition_name}")
    return instances


def display_statistics(zone_stats, args):
    """Display the statistics with proper colorizing."""
    headers = [
        "Zone", "Instances", "Avg Spot Price", "Min Spot Price", 
        "Max Spot Price", "Harmonic Mean Price", "Stability (Max-Min)",
        "BAM Size (GB)", "CRAM Size (GB)", "VCF Size (GB)",
        "Cost per vCPU", "Est. Cost of Analysis"
    ]

    # Calculate best and worst values for colorizing
    best_worst = {key: get_best_worst(zone_stats, key) for key in headers[2:]}
    min_cost = min([z['Est. Cost of Analysis'] for z in zone_stats if not isnan(z['Est. Cost of Analysis'])])

    table = []
    for index, z in enumerate(zone_stats, 1):
        row = [
            color(f"{index}. {z['Zone']}", fore=COLORS[args.mode]['zone']),
            z['Instances'],
            apply_color(z['Avg Spot Price'], *best_worst['Avg Spot Price'], mode=args.mode),
            apply_color(z['Min Spot Price'], *best_worst['Min Spot Price'], mode=args.mode),
            apply_color(z['Max Spot Price'], *best_worst['Max Spot Price'], mode=args.mode),
            apply_color(z['Harmonic Mean Price'], *best_worst['Harmonic Mean Price'], highlight=True, mode=args.mode),
            apply_color(z['Stability (Max-Min)'], *reversed(best_worst['Stability (Max-Min)']), mode=args.mode),
            z['BAM Size (GB)'], z['CRAM Size (GB)'], z['VCF Size (GB)'],
            apply_color(z['Cost per vCPU'], *best_worst['Cost per vCPU'], mode=args.mode),
            apply_color(z['Est. Cost of Analysis'], min_cost, min_cost, highlight=True, mode=args.mode)
        ]
        table.append(row)

    print(tabulate(table, headers=headers, floatfmt=".4f"))

    with open(args.output, 'w') as f:
        f.write(tabulate(table, headers=headers, floatfmt=".4f", tablefmt="tsv"))


def get_best_worst(stats, key):
    """Extract the best (min) and worst (max) values for a specific numeric key."""
    values = [z.get(key, float('nan')) for z in stats if not isnan(z.get(key, float('nan')))]
    if not values:
        return float('nan'), float('nan')  # Handle case where all values are NaN
    return min(values), max(values)

def calculate_statistics(spot_data, zones, vcpu_mins, args):
    """Calculate statistics including cost per vCPU and estimated cost of analysis."""
    zone_stats = []

    for zone in zones:
        prices = [spot_data[instance].get(zone, float('nan')) for instance in spot_data]
        prices = [p for p in prices if not isnan(p)]

        if prices:
            harmonic_price = harmonic_mean(prices)
            stability = stability_metric(prices)
            avg_price = statistics.mean(prices)
            min_price = min(prices)
            max_price = max(prices)

            # Use args to calculate sizes
            bam_size = args.x_coverage * args.bam_size_per_x if args.bam_size_per_x else 0
            cram_size = args.x_coverage * args.cram_size_per_x if args.cram_size_per_x else 0
            vcf_size = args.x_coverage * args.vcf_size_per_x if args.vcf_size_per_x else 0

            # Calculate cost per vCPU and estimated cost of analysis
            cost_per_vcpu = avg_price
            est_cost = (vcpu_mins / 60) * avg_price

            # Append statistics for the current zone
            zone_stats.append({
                'Zone': zone,
                'Instances': len(spot_data),
                'Avg Spot Price': avg_price,
                'Min Spot Price': min_price,
                'Max Spot Price': max_price,
                'Harmonic Mean Price': harmonic_price,
                'Stability (Max-Min)': stability,
                'Cost per vCPU': cost_per_vcpu,
                'Est. Cost of Analysis': est_cost,
                'BAM Size (GB)': bam_size,
                'CRAM Size (GB)': cram_size,
                'VCF Size (GB)': vcf_size,
            })

    return zone_stats

def main():
    args = parse_arguments()

    # Load YAML configuration
    with open(args.input, 'r') as f:
        config = yaml.safe_load(f)

    # Extract instances from the specified partition
    try:
        instance_types = extract_instances(config, args.partition)
    except ValueError as e:
        print(e)
        return

    zones = args.zones.split(',')
    spot_data = collect_spot_prices(instance_types, zones, args.profile)
    vcpu_mins = calculate_vcpu_mins(args)

    zone_stats = calculate_statistics(spot_data, zones, vcpu_mins, args)
    display_statistics(zone_stats, args)

    # Prompt user to select a zone by number
    try:
        selected_index = int(input("\nSelect the zone by number: ")) - 1
        selected_zone = zone_stats[selected_index]['Zone']
    except (IndexError, ValueError):
        print("Invalid selection. Exiting.")
        return


    alignment_output = input("Select alignment output [bam/cram]: ").strip().lower()
    if alignment_output not in ['bam', 'cram']:
        alignment_output = 'bam'

    variant_output = input("Select SNV output [vcf/gvcf/bcf/gbcf]: ").strip().lower()
    if variant_output not in ['vcf', 'gvcf', 'bcf', 'gbcf']:
        variant_output = 'vcf'


    # Assign the appropriate size variables based on user selection
    if alignment_output == 'cram':
        alignment_size = args.x_coverage * args.cram_size_per_x if args.cram_size_per_x else 0
    else:
        alignment_size = args.x_coverage * args.bam_size_per_x if args.bam_size_per_x else 0

    if variant_output == 'vcf':
        variant_size = args.x_coverage * args.vcf_size_per_x if args.vcf_size_per_x else 0
    elif variant_output == 'gvcf':
        variant_size = args.x_coverage * args.gvcf_size_per_x if args.gvcf_size_per_x else 0
    elif variant_output == 'bcf':
        variant_size = args.x_coverage * args.bcf_size_per_x if args.bcf_size_per_x else 0
    else:
        variant_size = args.x_coverage * args.gbcf_size_per_x if args.gbcf_size_per_x else 0

    # Prompt for transfer scheme and calculate the fully burdened cost
    input_transfer = input("Input transfer scheme [internet/within_region/cross_region]: ").lower()
    output_transfer = input("Output transfer scheme [internet/within_region/cross_region]: ").lower()

    total_transfer_cost = (calculate_transfer_costs(alignment_size, input_transfer) + 
                           calculate_transfer_costs(variant_size, output_transfer))
    # Validate input data size before using it
    fastq_size = args.input_data_size_per_x or 0

    # Calculate the fully burdened cost
    total_cost = total_transfer_cost + calculate_storage_costs(fastq_size, 'standard')
    print(f"\nFully Burdened Cost for {selected_zone}: ${total_cost:.2f}")

    print(f"\nFully Burdened Cost for {selected_zone}: ${total_cost:.2f}")

    # Prompt for FASTQ retention and display S3 storage costs
    retain_fastq = input("Retain FASTQ files? (y/n): ").strip().lower()
    if retain_fastq == 'y':
        s3_standard = calculate_storage_costs(fastq_size, 'standard')
        s3_glacier = calculate_storage_costs(fastq_size, 'glacier')

        print(f"\nS3 Storage Cost (Standard): ${s3_standard:.2f}/month")
        print(f"S3 Storage Cost (Glacier): ${s3_glacier:.2f}/month")

if __name__ == "__main__":
    main()