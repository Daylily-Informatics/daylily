#!/usr/bin/env python3

import sys
from pybedtools import BedTool

def calculate_total_unique_bases(bed_file):
    bed = BedTool(bed_file)
    merged = bed.merge()
    total_bases = sum(interval.length for interval in merged)
    return total_bases

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input.bed>")
        sys.exit(1)

    bed_file = sys.argv[1]
    total_unique_bases = calculate_total_unique_bases(bed_file)
    print(f"Total unique bases covered: {total_unique_bases}")
