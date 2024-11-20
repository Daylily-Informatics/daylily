##### MULTIQC
# ----------
#
# A brilliant meta tool.  It aggregates all of the qc info from other tools
# into a single report by searching the working and sub dirs for files/patterns.
# there is a vibrant multiqc module community, and many bfx tools have modules already
# There is alsoa  meta-meta tool called megaqc, which rolls up up another layer
# multiqc can send all of it's data to a central megaqc server where batch analysis
# longitudinal analysis, etc, can be done. https://megaqc.info/docs/index.html
# Written a a massive prolific s/w tools producer

ss_ex = EX[0]
ss_ru = RU[0]


rule multiqc_for_raw_fastqs:
    """https://github.com/ewels/MultiQC"""
    input:
        #MDIR+"logs/seqfu.done",
	#        expand(MDIR + "{sample}/seqqc/fastv/{sample}.fastv.html", sample=SAMPS),
        expand(MDIR + "{sample}/seqqc/fastqc/{sample}.fastqc.done", sample=SAMPS),
        #expand(MDIR + "{sample}/seqqc/kat/{sample}.kat.done", sample=SAMPS),
        #expand(MDIR + "{sample}/seqqc/fastp/{sample}.fastp.done", sample=SAMPS),
    output:
        MDIRreportsd + "SEQQC_multiqc.html",
    benchmark:
        f"{MDIR}benchmarks/SEQQC-multiqc_.raw_fastqc.bench.tsv"
    threads: config["multiqc"]["threads"]
    params:
        odir=MDIRreportsd,
        fn=f"SEQQC_multiqc.html",
        l="{",
        r="}",
        macro_cfg=config["multiqc"]["config_yaml"].replace('compiled',''),
        micro_cfg=config["multiqc"]["seqqc"]["config_yaml"].replace('compiled',''),
        name_cfg=config["multiqc_sampname_cfg"].replace('compiled',''),
        gtag=config["gittag"],
        ghash=config["githash"],
        gbranch=config["gitbranch"],
        cluster_sample=f"{RU[0]}_{EX[0]}",
    log:
        f"{MDIR}logs/multiqc/SEQQC_multiqc.FQ.log",
    conda:
        config["multiqc"]["seqqc"]["env_yaml"]  #  ...
    shell:
        """
        multiqc -f  --config config/external_tools/multiqc_header.yaml  --config    config/external_tools/multiqc_config.yaml -p -k tsv \
        -e picard -e goleft_indexcov -e mosdepth -e qualimap -x '*/snv/*' -x '*/sv/*' -x '*.js' -x '*bench' \
        -x '*.bam' -x '*.fastq.gz'  --profile-runtime  --interactive -x '.*impute.*' -x '*alignstats*' \
        -x '*impute*' -x '*impute*glm*' -x '*multiqc*' -x '*pyc' -x '*.fastq.gz' -i 'SEQQC Multiq Report' \
        -b '{samples.RU[0]}_{samples.EX[0]} ___ {params.gbranch} {params.gtag} {params.ghash}' -n {params.fn}  \
        -o {params.odir}   $(dirname {log} )/../../   >> {log} 2>&1;
        ls {output};
        """


localrules:
    seqqc,


# These '#' TARGET tagged rules are meant to allow clear execution of the sub-dags
rule seqqc:  # TARGET : Run Just Sequence QC Rules - No Alignment or Variant Calling Performed. A Multiqc Report Is Produced in results/mod*/reports
    conda:
        config["vanilla"]["env_yaml"]
    threads: 2
    input:
        f"{MDIRreportsd}SEQQC_multiqc.html",
    output:
        f"{MDIRlogd}seqqc.done",
    shell:
        "touch {output};"
