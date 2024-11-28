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
	    f"{MDIR}other_reports/rules_benchmark_data_mqc.tsv",
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
    conda:
        config["multiqc"]["seqqc"]["env_yaml"]  #  ...
    shell:
        """


        dbill='$';
        cp config/external_tools/multiqc_header.yaml $(dirname {output})/multiqc_header.yaml;
        perl -pi -e "s/REGSUB_PROJECT/$DAY_PROJECT/g;" $(dirname {output})/multiqc_header.yaml;
        perl -pi -e "s/REGSUB_BUDGET/\\\$dbill$USED_BUDGET of \\\$dbill$TOTAL_BUDGET spent ( $PERCENT_USED\%)/g;" $(dirname {output})/multiqc_header.yaml;

        size=$(du -hs results | cut -f1);
        perl -pi -e "s/REGSUB_TOTALSIZE/$size/g;" $(dirname {output})/multiqc_header.yaml;

        sed -i "s/REGSUB_EMAIL/${{DAY_CONTACT_EMAIL//@/\\\\@}}/g" $(dirname {output})/multiqc_header.yaml

        source bin/proc_spot_price_logs.sh >> {log} 2>&1;
        perl -pi -e "s/REGSUB_SPOTCOST/median: \\\$dbill$MEDIAN_SPOT_PRICE  mean: \\\$dbill$AVERAGE_SPOT_PRICE ( avg cost per vcpu,per min: \\\$dbill$VCPU_COST_PER_MIN ) /g;" $(dirname {output})/multiqc_header.yaml;
        perl -pi -e "s/REGSUB_SPOTINSTANCES/ $INSTANCE_TYPES_LINE /g;" $(dirname {output})/multiqc_header.yaml;

        source bin/proc_aligner_costs.sh {input[1]} $VCPU_COST_PER_MIN;
        perl -pi -e "s/REGSUB_TOTALCOST/$ALNR_SUMMARY_COST/g;" $(dirname {output})/multiqc_header.yaml;
        
        source bin/proc_mrkdup_costs.sh {input[1]} $VCPU_COST_PER_MIN;  
        perl -pi -e "s/REGSUB_MRKDUPCOST/$MRKDUP_AVG_MINUTES min, costing \\\$dbill$MRKDUP_AVG_COST/g;" $(dirname {output})/multiqc_header.yaml;

        multiqc -f  \
        --config   $(dirname {output})/multiqc_header.yaml \
        --config  config/external_tools/multiqc_config.yaml  \
        --custom-css-file config/external_tools/multiqc.css \
        --template default \
        --filename {output} \
        -i 'FASTQ Multiqc Report' \
        -b 'https://github.com/Daylily-Informatics/daylily (BRANCH:{params.gbranch}) (TAG:{params.gtag}) (HASH):{params.ghash}) ' \
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
