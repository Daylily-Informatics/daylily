##### RTG TOOLS ARE A GENERAL SET OF BFX UTILITIES
# ------------------------------------------------
#
# This is uysig them to calculate qc info on
# VCF files.
#
# They are best known for their tool 'vcfeval'
# which compuites the concordance between 2 VCFS
# and also attempt to define the best hard filters
# for your data.  They are also heavy early adopters to
# The distant change in VCF spec to BND format.


rule rtg_vcfstats:
    """https://github.com/RealTimeGenomics/rtg-tools"""
    input:
        svgz=MDIR
        + "{sample}/align/{alnr}/snv/{snv_caller}/{sample}.{alnr}.{snv_caller}.snv.sort.vcf.gz",
        svtbi=MDIR
        + "{sample}/align/{alnr}/snv/{snv_caller}/{sample}.{alnr}.{snv_caller}.snv.sort.vcf.gz.tbi",
    output:
        MDIR
        + "{sample}/align/{alnr}/snv/{snv_caller}/vcf_stats/{sample}.{alnr}.{snv_caller}.rtg.vcfstats.txt",
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.{snv_caller}.rtgvcfstats.bench.tsv"
    threads: config["rtg_vcfstats"]["threads"]
    resources:
        threads=config["rtg_vcfstats"]["threads"],
        partition=config["rtg_vcfstats"]["partition"],
    params:
        work_dir=MDIR + "{sample}/align/{alnr}/snv/{snv_caller}/vcf_stats/",
        cluster_sample=ret_sample,
    conda:
        config["rtg_vcfstats"]["env_yaml"]
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/{snv_caller}/vcf_stats/logs/{sample}.{alnr}.{snv_caller}.rtg.vcf.stats.log",
    shell:
        """
        mkdir -p {params.work_dir} > {log} 2>&1 ;
        rtg vcfstats --allele-lengths {input.svgz} > {output}  ;
        {latency_wait};
        ls {output};
        """


