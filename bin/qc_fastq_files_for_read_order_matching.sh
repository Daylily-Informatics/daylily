#!/bin/bash

# Usage: ./fastq_qc.sh <skip> <R1_1.fastq(.gz)> <R2_1.fastq(.gz)> ...
# Example: ./fastq_qc.sh 10000 R1_1.fastq.gz R2_1.fastq R1_2.fastq.gz R2_2.fastq

if [[ $# -lt 3 || $((($# - 1) % 2)) -ne 0 ]]; then
    echo "Usage: $0 <skip> <R1_1.fastq(.gz)> <R2_1.fastq(.gz)> ..."
    exit 1
fi

# Get the skip value and shift the arguments to get the file list
SKIP=$1
shift

# Function to count reads with skipping
count_reads() {
    local file=$1
    local skip=$2

    # Check if the file is gzipped or not
    if [[ "$file" == *.gz ]]; then
        reader="zcat"
    else
        reader="cat"
    fi

    # Use the appropriate reader (zcat or cat) and count reads by skipping
    $reader "$file" | awk -v skip="$skip" 'NR % (4 * skip) == 1' | wc -l
}

# Function to compare two FASTQ files
compare_fastq_pair() {
    local r1=$1
    local r2=$2

    r1_count=$(count_reads "$r1" "$SKIP")
    r2_count=$(count_reads "$r2" "$SKIP")

    if [[ "$r1_count" -ne "$r2_count" ]]; then
        echo "❌ Mismatch: $r1 ($r1_count reads) and $r2 ($r2_count reads)"
        return 1
    else
        echo "✅ Paired correctly: $r1 and $r2 ($r1_count reads)"
        return 0
    fi
}

# Export functions for GNU Parallel
export -f count_reads
export -f compare_fastq_pair
export SKIP

# Run comparisons in parallel
parallel --will-cite --jobs 4 compare_fastq_pair ::: "$@"

# Check if all comparisons passed
if [[ $? -eq 0 ]]; then
    echo "🎉 All FASTQ pairs passed QC!"
else
    echo "⚠️ Some FASTQ pairs failed QC. See above for details."
fi
