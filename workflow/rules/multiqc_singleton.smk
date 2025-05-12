import os
# This rule is actually gathering up all of the sub rules completed states
# so it can generate a final QC report, and that report will satisfy the input
# requirement of all to run

RPT_TITLE=os.environ.get("RPT_TITLE", "Final")

localrules:
    collect_rules_benchmark_data2,

rule collect_rules_benchmark_data2:
    output:
        f"{MDIR}other_reports/rules_benchmark_data2_mqc.tsv",
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


localrules:
	multiqc_singleton,
	
rule multiqc_singleton:  # TARGET: the big report
    input:
	    f"{MDIR}other_reports/rules_benchmark_data2_mqc.tsv",
    output:
        f"{MDIR}reports/multiqc_singleton.html",
        f"{MDIR}reports/multiqc_header2.yaml",
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
        rtitle=RPT_TITLE,
    log:
        f"{MDIR}reports/logs/all__mqc_fin_a2.log",
    container:
        "docker://daylilyinformatics/daylily_multiqc:0.2"
    shell:
        """
        dbill='$';

        echo '''
report_header_info:
  - Project/Budget: "REGSUB_PROJECT"
  - Budget @ Runtime: "REGSUB_BUDGET"
  - Spot Instances: "REGSUB_SPOTINSTANCES"
  - Spot Costs per hr: "REGSUB_SPOTCOST"
  - FQ->BAM.sort avg Costs: "REGSUB_TOTALCOST"
  - BAM mrkdup avg Cost: "REGSUB_MRKDUPCOST"
  - Results Dir (GB): "REGSUB_TOTALSIZE"
  ''' > {output[1]};

        perl -pi -e "s/REGSUB_PROJECT/$DAY_PROJECT/g;" {output[1]} >> {log} 2>&1;
        perl -pi -e "s/REGSUB_BUDGET/\\\$dbill$USED_BUDGET of \\\$dbill$TOTAL_BUDGET spent ( $PERCENT_USED\%)/g;" {output[1]} >> {log} 2>&1;

        size=$(du -hs results | cut -f1) >> {log} 2>&1;
        perl -pi -e "s/REGSUB_TOTALSIZE/$size/g;" {output[1]} >> {log} 2>&1;

        source bin/proc_spot_price_logs.sh >> {log} 2>&1;
        perl -pi -e "s/REGSUB_SPOTCOST/median: \\\$dbill$MEDIAN_SPOT_PRICE  mean: \\\$dbill$AVERAGE_SPOT_PRICE ( avg cost per vcpu,per min: \\\$dbill$VCPU_COST_PER_MIN ) /g;" {output[1]} >> {log} 2>&1;
        perl -pi -e "s/REGSUB_SPOTINSTANCES/ $INSTANCE_TYPES_LINE /g;" {output[1]} >> {log} 2>&1;

        source bin/proc_aligner_costs.sh {input} $VCPU_COST_PER_MIN >> {log} 2>&1;
        perl -pi -e "s/REGSUB_TOTALCOST/$ALNR_SUMMARY_COST/g;" {output[1]} >> {log} 2>&1;
        
        source bin/proc_mrkdup_costs.sh {input} $VCPU_COST_PER_MIN  >> {log} 2>&1;
        perl -pi -e "s/REGSUB_MRKDUPCOST/$MRKDUP_AVG_MINUTES min, costing \\\$dbill$MRKDUP_AVG_COST/g;"  {output[1]} >> {log} 2>&1;

        multiqc -f  \
        --config  {output[1]} \
        --config  ./config/external_tools/multiqc_config.yaml  \
        --custom-css-file ./config/external_tools/multiqc.css \
        --template default \
        --ignore "*sort_metrics/*" \
        --ignore "*/norm_cov_eveness/*" \
        --filename {output[0]} \
        -i '{params.rtitle} Multiqc Report' \
        -b 'https://github.com/Daylily-Informatics/daylily (BRANCH:{params.gbranch}) (TAG:{params.gtag}) (HASH:{params.ghash}) ' \
        $(dirname {input} )/../ >> {log} 2>&1;
        ls -lt {output[0]} {output[1]}  >> {log} 2>&1;
        """


localrules:
    produce_multiqc_singleton,


rule produce_multiqc_singleton:  # TARGET : Generated All WGS Reports
    input:
        MDIR+ "reports/multiqc_singleton.html"
