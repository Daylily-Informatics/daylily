#!/bin/bash

huref=$1
incram=$2
out_cram=$3

# Create contigs list from reference
#samtools view -H "$incram" | grep '^@SQ' | cut -f2 | sed 's/SN://' > hg38_contigs.txt
grep '^>' $huref | sed 's/>//' > hg38_contigs.txt

# Filter input cram to specified contigs and output new cram
samtools view -@ 8 -T "$huref" -C -o "$out_cram" "$incram" $(cat hg38_contigs.txt)

# Re-index the new cram file
samtools index -@ 8 "$out_cram"



#grep '^>' hg38.fa | sed 's/>//' > hg38_contigs.txt

#samtools view -@ 8 -T /fsx/data/genomic_data/organism_references/H_sapiens/hg38/fasta_fai_minalt/GRCh38_no_alt_analysis_set.fasta -C -o 409943-NA12878-Z0025-CTCGAGATTGATGAT.filtered.cram /home/ubuntu/xroom/ug/Feb27_2025_ugsb4_data/409943-NA12878-Z0025-CTCGAGATTGATGAT.cram $(cat hg38_contigs.txt) 
