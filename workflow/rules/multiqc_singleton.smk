import os
# This rule is actually gathering up all of the sub rules completed states
# so it can generate a final QC report, and that report will satisfy the input
# requirement of all to run

localrules:
    collect_rules_benchmark_data2,

rule collect_rules_benchmark_data2:
    input:
        f"{MDIR}logs/report_components_aggregated.done",
    output:
        f"{MDIR}other_reports/rules_benchmark_data_mqc2.tsv",
    params:
        cluster_sample="rules_benchmark_collect",
        working_file=f"{MDIR}reports/benchmarks_summary.tsv",
        ref_code=config["genome_build"],
    log:
        f"{MDIR}other_reports/logs/rules_benchmarks_summary2.log",
    container: None
    shell:
        "bin/util/benchmarks/collect_day_benchmark_data.sh {params.ref_code} > {log};"
        "python bin/util/benchmarks/split_bench_rule_col.py {params.working_file} {output} > {log};"
        "sed -i -E 's/\t$/\tNA/' {output};"


rule multiqc_final_wgs:  # TARGET: the big report
    input:
	    f"{MDIR}other_reports/rules_benchmark_data_mqc2.tsv",
    output:
        f"{MDIR}reports/multiqc_singleton.html",
    benchmark:
        f"{MDIR}benchmarks/DAY_all.final_multiqc2.bench.tsv"
    threads: config["multiqc"]["threads"]
    resources:
        threads=config["multiqc"]["threads"],
        partition=config["multiqc"]["partition"],
    priority: 50
    params:
        fnamef=f"DAY_final_multiqc.html",
        ghash=config["githash"],
        gbranch=config["gitbranch"],
        gtag=config["gittag"],
        cluster_sample=f"multiqc_final",
        cemail=config["day_contact_email"],
    log:
        f"{MDIR}reports/logs/all__mqc_fin_a2.log",
    container:
        "docker://daylilyinformatics/daylily_multiqc:0.2"
    shell:
        """
        dbill='$';
        cp config/external_tools/multiqc_header.yaml $(dirname {output})/multiqc_header.yaml >> {log} 2>&1;
        perl -pi -e "s/REGSUB_PROJECT/$DAY_PROJECT/g;" $(dirname {output})/multiqc_header.yaml >> {log} 2>&1;
        perl -pi -e "s/REGSUB_BUDGET/\\\$dbill$USED_BUDGET of \\\$dbill$TOTAL_BUDGET spent ( $PERCENT_USED\%)/g;" $(dirname {output})/multiqc_header.yaml >> {log} 2>&1;

        size=$(du -hs results | cut -f1) >> {log} 2>&1;
        perl -pi -e "s/REGSUB_TOTALSIZE/$size/g;" $(dirname {output})/multiqc_header.yaml >> {log} 2>&1;

        source bin/proc_spot_price_logs.sh >> {log} 2>&1;
        perl -pi -e "s/REGSUB_SPOTCOST/median: \\\$dbill$MEDIAN_SPOT_PRICE  mean: \\\$dbill$AVERAGE_SPOT_PRICE ( avg cost per vcpu,per min: \\\$dbill$VCPU_COST_PER_MIN ) /g;" $(dirname {output})/multiqc_header.yaml >> {log} 2>&1;
        perl -pi -e "s/REGSUB_SPOTINSTANCES/ $INSTANCE_TYPES_LINE /g;" $(dirname {output})/multiqc_header.yaml >> {log} 2>&1;

        source bin/proc_aligner_costs.sh {input[1]} $VCPU_COST_PER_MIN >> {log} 2>&1;
        perl -pi -e "s/REGSUB_TOTALCOST/$ALNR_SUMMARY_COST/g;" $(dirname {output})/multiqc_header.yaml >> {log} 2>&1;
        
        source bin/proc_mrkdup_costs.sh {input[1]} $VCPU_COST_PER_MIN  >> {log} 2>&1;
        perl -pi -e "s/REGSUB_MRKDUPCOST/$MRKDUP_AVG_MINUTES min, costing \\\$dbill$MRKDUP_AVG_COST/g;" $(dirname {output})/multiqc_header.yaml >> {log} 2>&1;

        multiqc -f  \
        --config   $(dirname {output})/multiqc_header.yaml \
        --config  config/external_tools/multiqc_config.yaml  \
        --custom-css-file config/external_tools/multiqc.css \
        --template default \
        --filename {output} \
        -i 'Final Multiqc Report' \
        -b 'https://github.com/Daylily-Informatics/daylily (BRANCH:{params.gbranch}) (TAG:{params.gtag}) (HASH:{params.ghash}) ' \
        $(dirname {input} )/../ >> {log} 2>&1;
        ls -lt {output}  >> {log} 2>&1;
        """


localrules:
    produce_multiqc_singleton,


rule produce_multiqc_singleton:  # TARGET : Generated All WGS Reports
    input:
        MDIR+ "reports/multiqc_final.html"
