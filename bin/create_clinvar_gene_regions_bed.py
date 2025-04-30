#!/usr/bin/env python3
import pandas as pd
import gzip
import requests
from io import StringIO
import pybedtools
import os

CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
GFF_HG38_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz"
GFF_GRCH37_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz"

# download and parse ClinVar gene list
def get_clinvar_genes():
    response = requests.get(CLINVAR_URL, stream=True)
    with gzip.open(response.raw, 'rt') as f:
        df = pd.read_csv(f, sep='\t', low_memory=False)
    genes = set(df['GeneSymbol'].dropna().unique())
    print(f"ClinVar Genes extracted: {len(genes)}")
    return genes

# parse GFF, extract exons/UTRs by gene, extend by 500 bp, output BED
def generate_bed(gff_url, genes, assembly):
    print(f"Processing assembly: {assembly}")
    response = requests.get(gff_url, stream=True)
    with gzip.open(response.raw, 'rt') as f:
        gff_lines = [line for line in f if not line.startswith('#')]

    gff_df = pd.read_csv(StringIO(''.join(gff_lines)), sep='\t', header=None,
                         names=['chr','source','feature','start','end','score','strand','phase','attributes'],
                         low_memory=False)

    gff_df = gff_df[gff_df['feature'].isin(['exon','UTR'])]
    
    def parse_attrs(attr):
        attrs = {k:v for k,v in [x.split('=') for x in attr.split(';') if '=' in x]}
        return attrs.get('gene', attrs.get('Name',''))

    gff_df['gene'] = gff_df['attributes'].apply(parse_attrs)
    gff_df = gff_df[gff_df['gene'].isin(genes)]

    bed_lines = []
    for gene, group in gff_df.groupby('gene'):
        intervals = pybedtools.BedTool.from_dataframe(group[['chr','start','end']])
        intervals = intervals.slop(b=500, genome=assembly).merge()
        for interval in intervals:
            bed_lines.append([interval.chrom, interval.start, interval.end, gene])

    bed_df = pd.DataFrame(bed_lines, columns=['chr','start','end','gene'])
    bed_df = bed_df.sort_values(['chr','start'])

    bed_file = f"clinvar_genes_{assembly}.bed"
    bed_df.to_csv(bed_file, sep='\t', header=False, index=False)
    print(f"Saved BED file: {bed_file}")

if __name__ == "__main__":
    genes = get_clinvar_genes()
    generate_bed(GFF_HG38_URL, genes, "hg38")
    generate_bed(GFF_GRCH37_URL, genes, "hg19")
