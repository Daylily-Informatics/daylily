#!/usr/bin/env python3


#chmod +x generate_comprehensive_clinvar_bed.py
#pip install pandas pybedtools requests
#sudo apt install bedtools

#!/usr/bin/env python3

import pandas as pd
import requests
import pybedtools
import gzip
import os
from io import StringIO

CLINVAR_VCF = {
    'hg38': 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz',
    'b37': 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz'
}

GENCODE_GTF = {
    'hg38': 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.basic.annotation.gtf.gz',
    'b37': 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz'
}

def download_and_extract(url, outfile):
    print(f"Downloading {url}...")
    response = requests.get(url, stream=True)
    response.raise_for_status()
    with open(outfile, 'wb') as f_out:
        for chunk in response.iter_content(chunk_size=8192):
            f_out.write(chunk)
    print(f"Downloaded to {outfile}")

def parse_clinvar_vcf(vcf_gz):
    print(f"Parsing ClinVar VCF {vcf_gz}...")
    clinvar_bed = []
    with gzip.open(vcf_gz, 'rt') as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom, pos = fields[0], int(fields[1])

            # Explicitly add 'chr' prefix if missing (for hg38 compatibility)
            if not chrom.startswith('chr'):
                chrom = 'chr' + chrom

            start = max(0, pos - 501)
            end = pos + 500
            clinvar_bed.append([chrom, start, end])
    print(f"Extracted {len(clinvar_bed)} variants from ClinVar VCF.")
    return pybedtools.BedTool(clinvar_bed)

def parse_gencode_gtf(gtf_gz):
    print(f"Parsing Gencode GTF {gtf_gz}...")
    intervals = []
    with gzip.open(gtf_gz, 'rt') as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            feature_type = fields[2]
            if feature_type in ('exon', 'UTR'):
                chrom, start, end = fields[0], int(fields[3])-1, int(fields[4])
                
                # Ensure Gencode chromosome names match ClinVar ('chr' prefix standard)
                if not chrom.startswith('chr'):
                    chrom = 'chr' + chrom

                start = max(0, start - 500)
                end += 500
                intervals.append([chrom, start, end])
    print(f"Extracted {len(intervals)} exonic/UTR intervals from GTF.")
    return pybedtools.BedTool(intervals)

def create_bed(genome):
    vcf_url = CLINVAR_VCF[genome]
    gtf_url = GENCODE_GTF[genome]

    os.system("mkdir -p ./tmp")
    clinvar_vcf_file = f"./tmp/clinvar_{genome}.vcf.gz"
    gtf_file = f"./tmp/gencode_{genome}.gtf.gz"

    download_and_extract(vcf_url, clinvar_vcf_file)
    download_and_extract(gtf_url, gtf_file)

    clinvar_variants = parse_clinvar_vcf(clinvar_vcf_file)
    gene_regions = parse_gencode_gtf(gtf_file)

    # Subtract gene regions from clinvar variants to find non-coding variants
    print("Subtracting gene regions from ClinVar variants to isolate non-coding variants...")
    noncoding_variants = clinvar_variants.subtract(gene_regions)

    # Combine noncoding variants with gene regions
    print("Combining gene regions and non-coding ClinVar variants...")
    combined = gene_regions.cat(noncoding_variants, postmerge=False)

    # Merge overlapping intervals
    print("Merging overlapping intervals...")
    merged = combined.sort().merge()

    output_bed = f"comprehensive_clinvar_genes_{genome}.bed"
    merged.saveas(output_bed)
    print(f"Created BED file: {output_bed}")

    # Clean up
    #os.remove(clinvar_vcf_file)
    #os.remove(gtf_file)

if __name__ == "__main__":
    for genome in ['b37','hg38']:
        create_bed(genome)
