import os

# ### Stub For Sentieon Haplotyper
# ...
#


rule sentieon_haplotypecaller:
    # Run sentieon on the the imputation bam, butg only call in the chrm active in this instance
    priority: 48
    input:
        bam=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.sort.bam",
        bambai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.sort.bam.bai",
    output:
        ovcf=MDIR + "{sample}/align/{alnr}/snv/senthc/{sample}.{alnr}.senthc.vcf",
        ovcfsort=MDIR
        + "{sample}/align/{alnr}/snv/senthc/{sample}.{alnr}.senthc.sort.vcf",
        ovcfgz=MDIR
        + "{sample}/align/{alnr}/snv/senthc/{sample}.{alnr}.senthc.sort.vcf.gz",
        ovcfgztbi=MDIR
        + "{sample}/align/{alnr}/snv/senthc/{sample}.{alnr}.senthc.sort.vcf.gz.tbi",
    benchmark:
        (MDIR + "{sample}/benchmarks/{sample}.{alnr}.senthc.bench.tsv")
    threads: config["sentieon"]["threads"]
    params:
        huref=config["supporting_files"]["files"]["glimpse"][
            f"""huref_{config["glimpse"]["genome_build"]}"""
        ]["name"],
    log:
        MDIR + "{sample}/align/{alnr}/snv/senthc/logs/{sample}.{alnr}.senthc.log",
    conda:
        config["sentieon"]["env_yaml"]
    params:
        call_interval=" "
        if int(wildcards.call_chrm) == 0
        else f" --interval   {wildcards.call_chrm} ",
        mask_regions=" "
        if "only_sites_regions" in config["sentieon"]
        else config["sentieon"]["only_sites_regions"],
        cluster_sample=ret_sample,
        #caller_rc=" --algo Haplotyper --emit_mode variant --genotype_model coalescent ",
    shell:
        " export bwt_max_mem=120G && sentieon bwa mem && sentieon driver -t {threads}  {params.call_interval}  -r {params.huref} -i {input.bam} --algo Haplotyper --emit_mode ALL --ploidy 2 {output.ovcf} > {log} 2>&1;"
        "bedtools sort -header -i {output.ovcf} > {output.ovcfsort};"
        "bgzip -@ {threads} {output.ovcfsort};"
        "tabix {output.ovcfgz};"
        "touch {output.ovcfsort}; "
        "{latency_wait};"
        "ls {output};"
