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
        expand(MDIR + "{sample}/seqqc/fastqc/{sample}.fastqc.done", sample=SAMPS),
    output:
        MDIRreportsd + "SEQQC_multiqc.html",
    benchmark:
        f"{MDIR}benchmarks/SEQQC-multiqc_.raw_fastqc.bench.tsv"
    threads: config["multiqc"]["threads"]
    params:
        fn=f"SEQQC_multiqc.html",
        gtag=config["gittag"],
        ghash=config["githash"],
        gbranch=config["gitbranch"],
        cluster_sample=f"{RU[0]}_{EX[0]}",
    log:
        f"{MDIR}logs/multiqc/SEQQC_multiqc.FQ.log",
    container:
        "docker://daylilyinformatics/daylily_multiqc:0.2"
    shell:
        """


        dbill='$';
        cp ./config/external_tools/multiqc_header.yaml ./$(dirname {output})/rfq_multiqc_header.yaml;

        multiqc -f  \
        --config   ./$(dirname {output})/rfq_multiqc_header.yaml \
        --config  ./config/external_tools/multiqc_config.yaml  \
        --custom-css-file ./config/external_tools/multiqc.css \
        --template default \
        --filename {output} \
        -i 'FASTQ Multiqc Report' \
        -b 'https://github.com/Daylily-Informatics/daylily (BRANCH:{params.gbranch}) (TAG:{params.gtag}) (HASH:{params.ghash}) ' \
        $(dirname {input} )/../ > {log} 2>&1;
        ls -lt {output};

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
