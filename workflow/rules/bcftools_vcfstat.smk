######### BCFTOOLS stats
# ----------------------------
# This is a population agnostic contamination screening tool that can
# operate on single sample or multi sample BAM files.
# github: https://github.com/samtools/bcftools
# paper: http://samtools.github.io/bcftools/howtos/publications.html
# docs: https://samtools.github.io/bcftools/bcftools.html


rule bcftools_vcfstat:
    input:
        snv_vcf=(
            MDIR
            + "{sample}/align/{alnr}/snv/{snv_caller}/{sample}.{alnr}.{snv_caller}.snv.sort.vcf.gz"
        ),
        snv_vcf_tbi=(
            MDIR
            + "{sample}/align/{alnr}/snv/{snv_caller}/{sample}.{alnr}.{snv_caller}.snv.sort.vcf.gz.tbi"
        ),
    output:
        MDIR + "{sample}/align/{alnr}/snv/{snv_caller}/bcfstats/{sample}.{alnr}.{snv_caller}.bcfstats.tsv",
    log:
        MDIR + "{sample}/align/{alnr}/snv/{snv_caller}/bcfstats/logs/{sample}.{alnr}.{snv_caller}.bcfstats.log",
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.{snv_caller}.bcfstat.bench.tsv",
    conda:
        config["vanilla"]["env_yaml"]
    threads: config["bcftools_vcfstat"]["threads"]
    params:
        cluster_sample=ret_sample,
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
    shell:
        """
        bcftools stats --threads {threads} {input.snv_vcf} -F {params.huref} > {output};
        {latency_wait};
        ls {output};
        """


localrules:
    produce_bcfvcfstats,

rule produce_bcfvcfstats:  # TARGET:  jusg genvcfstats
    input:
        expand(
            MDIR + "{sample}/align/{alnr}/snv/{snv_caller}/bcfstats/{sample}.{alnr}.{snv_caller}.bcfstats.tsv",
            sample=SSAMPS,
            alnr=ALIGNERS,
            snv_caller=snv_CALLERS
        ),
