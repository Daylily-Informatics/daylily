0#### snpeff
# -------------------------------------


rule snpeff:
    input:
        vcfgz=MDIR
        + "{sample}/align/{alnr}/snv/{snv}/{sample}.{alnr}.{snv}.snv.sort.vcf.gz",
    output:
        annovcf=MDIR
        + "{sample}/align/{alnr}/snv/{snv}/snpeff/{sample}.{alnr}.{snv}.snpeff.vcf.gz",
        annovcftbi=MDIR
        + "{sample}/align/{alnr}/snv/{snv}/snpeff/{sample}.{alnr}.{snv}.snpeff.vcf.gz.tbi",        
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/{snv}/snpeff/log/{sample}.{alnr}.{snv}.snpeff.log",
    threads: config["snpeff"]["threads"]
    resources:
        vcpu=config["snpeff"]["threads"],
        partition=config["snpeff"]["partition"],
        threads=config["snpeff"]["threads"],
    params:
        cluster_sample=ret_sample,
        snpeff_genome_build=config["supporting_files"]["files"]["snpeff"]["ref_build_code"],
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        snpeff_xmx="16g" if "xmx" not in config["snpeff"] else config["snpeff"]["xmx"],
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.{snv}.snpeff.bench.tsv"
    conda:
        "../envs/snpeff_v0.1.yaml"
    shell:
        """=
        ( java -Xmx{params.snpeff_xmx} -jar  \
        $(find $CONDA_PREFIX -name snpEff.jar) \
        -v {params.snpeff_genome_build} \
        {input.vcfgz} \
        ) > {log} 2>&1;
        """


localrules:
    produce_snpeff,


rule produce_snpeff:  # TARGET: just produce snpeff results
    input:
        expand(
            MDIR
            + "{sample}/align/{alnr}/snv/{snv}/snpeff/{sample}.{alnr}.{snv}.snpeff.vcf.gz.tbi",
            sample=SSAMPS,
            alnr=ALL_ALIGNERS,
            snv=snv_CALLERS,
        ),
    output:
        "logs/snpeff_gathered.done",
    shell:
        "touch {output};"
