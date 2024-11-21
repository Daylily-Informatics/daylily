#!/bin/bash

# Check input arguments
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <input_calls_vcf> <output_fixed_vcf>"
    exit 1
fi

# Input arguments
calls_vcf=$1
output_vcf=$2
prefix=$3
sample=$4

# Temporary files
temp_full="temp_full.vcf"
temp_header="temp_header.vcf"
temp_body="temp_body.vcf"
temp_vcf="$prefix.vcf"

# Step 1: Decompress VCF if needed and split into header and body
if [[ $calls_vcf == *.gz ]]; then
    zcat $calls_vcf > $temp_full
else
    cp $calls_vcf $temp_full
fi

grep "^#" $temp_full > $temp_header
grep -v "^#" $temp_full > $temp_body


to_add_a='##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">'
to_add_b='##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
new_title="#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  $sample"

perl -pi -e 's/\#CHROM\t.*$//g' $temp_header
echo "$to_add_a" >> $temp_header
echo "$to_add_b" >> $temp_header
echo "$new_title" >> $temp_header


cp $temp_header temp_corrected_header.vcf

# Step 3: Add missing FORMAT and sample columns to the body
echo "Adding FORMAT and sample columns if missing..."
awk 'BEGIN {OFS="\t"}
    {
        if (NF < 9) {  # If FORMAT and Sample columns are missing
            $9 = "GT:GQ";
            $10 = "0/1:99";
        } else if ($9 !~ /GQ/) {  # Add GQ to FORMAT if missing
            $9 = $9 ":GQ";
            $10 = $10 ":99";
        }
        print
    }' $temp_body > temp_corrected_body.vcf

# Step 4: Combine header and body into final VCF
cat temp_corrected_header.vcf temp_corrected_body.vcf > $temp_vcf

perl -pi -e 's/^\n//g;' $temp_vcf

mv $temp_vcf $output_vcf

bgzip $output_vcf
tabix $output_vcf.gz

# Cleanup
rm -f $temp_full $temp_header $temp_body chrom_header.tmp metadata_header.tmp temp_corrected_header.vcf temp_corrected_body.vcf

echo "Fixed VCF written to: $output_vcf"
