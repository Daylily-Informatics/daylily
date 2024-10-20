import subprocess
import sys

def get_all_regions():
    """Retrieve all available AWS regions."""
    try:
        command = ["aws", "ec2", "describe-regions", "--query", "Regions[*].RegionName", "--output", "text"]
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        return result.stdout.strip().split()
    except subprocess.CalledProcessError as e:
        print(f"Error fetching regions: {e.stderr.strip()}", file=sys.stderr)
        sys.exit(1)

def check_spot_availability(instance_type, region):
    """Check spot price history to determine if spot instances are available for the given instance type."""
    try:
        command = [
            "aws", "ec2", "describe-spot-price-history",
            "--instance-types", instance_type,
            "--product-description", "Linux/UNIX",
            "--region", region,
            "--start-time", "2024-10-20T00:00:00Z",
            "--query", "SpotPriceHistory[*].AvailabilityZone",
            "--output", "text"
        ]
        result = subprocess.run(command, check=True, capture_output=True, text=True)

        return set(result.stdout.strip().split()) if result.stdout.strip() else None
    except subprocess.CalledProcessError as e:
        print(f"Error checking spot availability for {instance_type} in {region}: {e.stderr.strip()}", file=sys.stderr)
        return None

def main():
    # List of instance types to check
    instance_types = [
        "c7i.48xlarge",
        "c7i.metal-48xl",
        "m7i.metal-48xl",
        "m7i.48xlarge",
        "r7i.48xlarge",
        "r7i.metal-48xl"
    ]

    # Fetch all regions
    regions = get_all_regions()

    # Check spot instance availability for each instance type in all regions
    for instance_type in instance_types:
        print(f"\nChecking spot availability for {instance_type}...")
        for region in regions:
            azs = check_spot_availability(instance_type, region)
            if azs:
                print(f"  Region: {region} - Available in: {', '.join(azs)}")
            else:
                print(f"  Region: {region} - Not available.")

if __name__ == "__main__":
    main()
