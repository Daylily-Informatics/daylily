# app.py

import uvicorn
from fastapi import FastAPI, Request, Form
from fastapi.responses import HTMLResponse
from typing import List
import yaml
import subprocess
import statistics
from math import isnan, fsum
from tabulate import tabulate
import os

app = FastAPI()

# Define ANSI colors for different modes (simplified for web)
COLORS = {
    'dark': {
        'highlight': 'green',
        'header': 'white',
        'zone': 'skyblue',
        'dim': 'grey',
        'gray': 'lightgrey',
        'teal': 'teal',
        'orange': 'orange'
    }
}

TRANSFER_RATES = {
    'internet': 0.09,  # $/GB to internet
    'within_region': 0.01,  # $/GB within region
    'cross_region': 0.02  # $/GB across regions
}

STORAGE_COSTS = {'standard': 0.023, 'glacier': 0.004}  # $/GB/month

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

def calculate_vcpu_mins(x_coverage, vcpu_min_per_x):
    """Calculate total vCPU-minutes based on coverage."""
    return x_coverage * vcpu_min_per_x

def calculate_storage_costs(size_gb, storage_type):
    """Calculate monthly storage costs for a given size."""
    return size_gb * STORAGE_COSTS[storage_type]

def calculate_transfer_costs(size_gb, transfer_type):
    """Calculate transfer costs based on the type."""
    return size_gb * TRANSFER_RATES[transfer_type]

def get_spot_price(instance_type, zone, profile=None):
    region = zone[:-1]
    try:
        cmd = [
            "aws", "ec2", "describe-spot-price-history",
            "--instance-types", instance_type,
            "--availability-zone", zone,
            "--product-description", "Linux/UNIX",
            "--query", "SpotPriceHistory[0].SpotPrice",
            "--output", "text", "--region", region,
        ]
        if profile:
            cmd.extend(["--profile", profile])
        result = subprocess.check_output(cmd).decode().strip()
        return float(result)
    except subprocess.CalledProcessError as e:
        print(f"Failed to fetch price for {zone}: {e}")
        return float('nan')

def collect_spot_prices(instance_types, zones, profile=None):
    """Fetch spot prices for each instance type in the given zones."""
    spot_data = {}
    for instance_type in instance_types:
        spot_data[instance_type] = {}
        for zone in zones:
            price = get_spot_price(instance_type, zone, profile)
            spot_data[instance_type][zone] = price if not isnan(price) else float('nan')
    return spot_data

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

def get_best_worst(stats, attr):
    """Extract the best (min) and worst (max) values for a specific numeric attribute."""
    values = [getattr(z, attr) for z in stats if not isnan(getattr(z, attr))]
    if not values:
        return float('nan'), float('nan')  # Handle case where all values are NaN
    return min(values), max(values)

class ZoneStat:
    def __init__(self, zone_name):
        self.zone_name = zone_name
        self.instances = 0
        self.prices = []
        self.avg_price = float('nan')
        self.min_price = float('nan')
        self.max_price = float('nan')
        self.harmonic_price = float('nan')
        self.stability = float('nan')
        self.cost_per_vcpu = float('nan')
        self.est_cost = float('nan')
        self.bam_size = float('nan')
        self.cram_size = float('nan')
        self.vcf_size = float('nan')
        self.gvcf_size = float('nan')
        self.bcf_size = float('nan')
        self.gbcf_size = float('nan')
        self.fastq_size = float('nan')
        self.qc_size = float('nan')  # Added for QC data size

    def calculate_statistics(self, spot_data, vcpu_mins, args):
        # Calculate the various statistics for this zone
        self.prices = [spot_data[instance].get(self.zone_name, float('nan')) for instance in spot_data]
        self.prices = [p for p in self.prices if not isnan(p)]
        if self.prices:
            self.instances = len(self.prices)
            self.avg_price = statistics.mean(self.prices)
            self.min_price = min(self.prices)
            self.max_price = max(self.prices)
            self.harmonic_price = harmonic_mean(self.prices)
            self.stability = stability_metric(self.prices)
            self.cost_per_vcpu = self.avg_price
            self.est_cost = (vcpu_mins / 60) * self.avg_price

            # Use args to calculate sizes
            self.bam_size = args['x_coverage'] * args['bam_size_per_x'] if args['bam_size_per_x'] else 0
            self.cram_size = args['x_coverage'] * args['cram_size_per_x'] if args['cram_size_per_x'] else 0
            self.vcf_size = args['x_coverage'] * args['vcf_size_per_x'] if args['vcf_size_per_x'] else 0
            self.gvcf_size = args['x_coverage'] * args['gvcf_size_per_x'] if args['gvcf_size_per_x'] else 0
            self.bcf_size = args['x_coverage'] * args['bcf_size_per_x'] if args['bcf_size_per_x'] else 0
            self.gbcf_size = args['x_coverage'] * args['gbcf_size_per_x'] if args['gbcf_size_per_x'] else 0
            self.fastq_size = args['x_coverage'] * args['input_data_size_per_x'] if args['input_data_size_per_x'] else 0
            self.qc_size = args['x_coverage'] * args['qc_size_per_x'] if args['qc_size_per_x'] else 0  # QC data size
        else:
            # No prices, leave values as NaN or zero
            pass

@app.get("/", response_class=HTMLResponse)
async def index():
    html_content = """
    <html>
    <head>
        <title>Genomic Workflow Cost Calculator</title>
    </head>
    <body>
        <h1>Genomic Workflow Cost Calculator</h1>
        <form action="/calculate" method="post">
            <label>Availability Zones (comma-separated):</label><br>
            <input type="text" name="zones" value="us-west-2a,us-west-2b,us-west-2c"><br><br>

            <label>Coverage (X):</label><br>
            <input type="number" name="x_coverage" value="30.1" step="0.1"><br><br>

            <label>vCPU-minutes per X coverage:</label><br>
            <input type="number" name="vcpu_min_per_x" value="3.7" step="0.1"><br><br>

            <!-- Add more input fields as needed -->
            <input type="submit" value="Calculate">
        </form>
    </body>
    </html>
    """
    return html_content

@app.post("/calculate", response_class=HTMLResponse)
async def calculate(request: Request,
                    zones: str = Form(...),
                    x_coverage: float = Form(...),
                    vcpu_min_per_x: float = Form(...)):
    # Default arguments (you can expose these in the form if needed)
    args = {
        'x_coverage': x_coverage,
        'vcpu_min_per_x': vcpu_min_per_x,
        'bam_size_per_x': 1.0,
        'cram_size_per_x': 1.0,
        'vcf_size_per_x': 1.0,
        'gvcf_size_per_x': 1.0,
        'bcf_size_per_x': 1.0,
        'gbcf_size_per_x': 1.0,
        'qc_size_per_x': 0.5,
        'input_data_size_per_x': 1.0,
        'mode': 'dark',
        'partition': 'i192',
        'profile': None,
        'input': 'config/day_cluster/prod_cluster.yaml'
    }

    # Load YAML configuration
    with open(args['input'], 'r') as f:
        config = yaml.safe_load(f)

    # Extract instances from the specified partition
    try:
        instance_types = extract_instances(config, args['partition'])
    except ValueError as e:
        return f"<h3>Error: {e}</h3>"

    zone_list = [zone.strip() for zone in zones.split(',')]
    spot_data = collect_spot_prices(instance_types, zone_list, args['profile'])
    vcpu_mins = calculate_vcpu_mins(args['x_coverage'], args['vcpu_min_per_x'])

    # Calculate statistics for each zone
    zone_stats = []
    for zone in zone_list:
        zone_stat = ZoneStat(zone)
        zone_stat.calculate_statistics(spot_data, vcpu_mins, args)
        zone_stats.append(zone_stat)

    # Generate table data
    headers = [
        "Availability Zone", "Instances", "Avg Spot Price", "Min Spot Price",
        "Max Spot Price", "Harmonic Mean Price", "Stability (Max-Min)",
        "BAM (GB)", "CRAM (GB)", "VCF (GB)", "gVCF (GB)", "BCF (GB)", "gBCF (GB)",
        "Cost per vCPU", "Est. EC2 Cost"
    ]

    table = []
    for index, z in enumerate(zone_stats, 1):
        row = [
            f"{index}. {z.zone_name}",
            z.instances,
            f"{z.avg_price:.4f}",
            f"{z.min_price:.4f}",
            f"{z.max_price:.4f}",
            f"{z.harmonic_price:.4f}",
            f"{z.stability:.4f}",
            f"{z.bam_size:.4f}",
            f"{z.cram_size:.4f}",
            f"{z.vcf_size:.4f}",
            f"{z.gvcf_size:.4f}",
            f"{z.bcf_size:.4f}",
            f"{z.gbcf_size:.4f}",
            f"{z.cost_per_vcpu:.4f}",
            f"{z.est_cost:.4f}"
        ]
        table.append(row)

    # Convert table to HTML
    table_html = tabulate(table, headers=headers, tablefmt="html")

    # Generate the response HTML
    html_content = f"""
    <html>
    <head>
        <title>Calculation Results</title>
    </head>
    <body>
        <h1>Calculation Results</h1>
        <p>{args['x_coverage']}-cov genome @ {args['vcpu_min_per_x']} vCPU-min per coverage</p>
        {table_html}
        <br>
        <a href="/">Go Back</a>
    </body>
    </html>
    """
    return html_content

if __name__ == "__main__":
    # Run the app with: uvicorn app:app --reload
    uvicorn.run(app, host="0.0.0.0", port=8000)
