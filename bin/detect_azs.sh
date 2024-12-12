#!/bin/bash
echo -e "Region\tAvailabilityZone" > azs.tsv

for region in $(aws ec2 describe-regions --query "Regions[].RegionName" --output text); do
  aws ec2 describe-availability-zones --region "$region" --query "AvailabilityZones[?State=='available'].ZoneName" --output text | while read -r az; do
    echo -e "$region\t$az" >> azs.tsv
  done
done
