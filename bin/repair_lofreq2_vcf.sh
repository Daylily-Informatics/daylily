#!/bin/bash

# Check input arguments
if [[ $# -lt 4 ]]; then
    echo "Usage: $0 <input_calls_vcf> <output_fixed_vcf> <tmp_dir> <sample>"
    exit 1
fi

# Input arguments
calls_vcf=$1
output_vcf=$2
tmp_dir=$3
sample=$4

# Ensure tmp_dir ends with a slash
tmp_dir="${tmp_dir%/}/"

# Unique identifier for this process
unique_id=$$

# Temporary files with unique identifiers
temp_full="${tmp_dir}temp_full_${unique_id}.vcf"
temp_header="${tmp_dir}temp_header_${unique_id}.vcf"
temp_body="${tmp_dir}temp_body_${unique_id}.vcf"
temp_corrected_header="${tmp_dir}temp_corrected_header_${unique_id}.vcf"
temp_corrected_body="${tmp_dir}temp_corrected_body_${unique_id}.vcf"
temp_vcf="${tmp_dir}final_temp_${unique_id}.vcf"

# Step 1: Decompress VCF if needed and split into header and body
if [[ $calls_vcf == *.gz ]]; then
    zcat "$calls_vcf" > "$temp_full"
else
    cp "$calls_vcf" "$temp_full"
fi

grep "^#" "$temp_full" > "$temp_header"
grep -v "^#" "$temp_full" > "$temp_body"

to_add_a='##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">'
to_add_b='##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
new_title="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample"

perl -pi -e 's/\#CHROM\t.*$//g' "$temp_header"
echo -e "$to_add_a" >> "$temp_header"
echo -e "$to_add_b" >> "$temp_header"
echo -e "$new_title" >> "$temp_header"

cp "$temp_header" "$temp_corrected_header"

# Step 3: Handle body (if present)
if [[ -s $temp_body ]]; then
    echo "Adding FORMAT and sample columns to body..."
    awk 'BEGIN {OFS="\t"}
        {
            if (NF < 9) {  # If FORMAT and Sample columns are missing
                $9 = "GT:GQ";
                $10 = "1/1:35";
            } else if ($9 !~ /GQ/) {  # Add GQ to FORMAT if missing
                $9 = $9 ":GQ";
                $10 = $10 ":99";
            }
            print
        }' "$temp_body" > "$temp_corrected_body"
    # Combine header and corrected body
    cat "$temp_corrected_header" "$temp_corrected_body" > "$temp_vcf"
else
    echo "No variant calls found. Creating VCF with updated headers only..."
    cp "$temp_corrected_header" "$temp_vcf"
fi

# Step 4: Combine header and body into final VCF
cat "$temp_corrected_header" "$temp_corrected_body" > "$temp_vcf"

perl -pi -e 's/^\n//g;' "$temp_vcf"

cp "$temp_vcf" "$output_vcf"

# Optional: Compress and index the output VCF
# bgzip "$output_vcf"
# tabix "$output_vcf.gz"

# Cleanup
rm -f "$temp_full" "$temp_header" "$temp_body" \
      "$temp_corrected_header" "$temp_corrected_body" "$temp_vcf"

echo "Fixed VCF written to: $output_vcf"
