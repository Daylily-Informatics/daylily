# LoFreq -- Variant Caller
# ------------------------
#
# TO BE ASSESSED-- possible faster alternative than haplotypecaller
#
# Github: https://csb5.github.io/lofreq/installation/
#
# https://pubmed.ncbi.nlm.nih.gov/23066108/
#
# conda lofreq_v0.1.yaml

rule lofreq2:
    conda:
        "../envs/lofreq_v0.1.yaml"
    shell:
        "lofreq call-parallel --pp-threads 12 -f /fsx/data/genomic_data/organism_references/H_sapiens/b37/human_1kg_v37/human_g1k_v37.fasta -o ./lo.vcf --bed /fsx/data/genomic_data/organism_refrences/H_sapiens/b37/controls/giab/snv/v4.2.1/HG002/b37GIAB_hc/HG002.bed /fsx/SAVEME_ANA/HG002/NOVA/fresh/day/results/mod/b37/RIH0_ANA0-HG002-NOVA-19_DBC0_0/align/bwa2a/RIH0_ANA0-HG002-NOVA-19_DBC0_0.bwa2a.sort.bam;"
