import os
# This rule is actually gathering up all of the sub rules completed states
# so it can generate a final QC report, and that report will satisfy the input
# requirement of all to run


RPT_TITLE=os.environ.get("RPT_TITLE", "Final")

localrules:
    collect_rules_benchmark_data,


rule collect_rules_benchmark_data:
    input:
        f"{MDIR}logs/report_components_aggregated.done",
    output:
        f"{MDIR}other_reports/rules_benchmark_data_mqc.tsv",
    params:
        cluster_sample="rules_benchmark_collect",
        working_file=f"{MDIR}reports/benchmarks_summary.tsv",
        ref_code=config["genome_build"],
    log:
        f"{MDIR}other_reports/logs/rules_benchmarks_summary.log",
    container: None
    shell:
        "bin/util/benchmarks/collect_day_benchmark_data.sh {params.ref_code} > {log};"
        "python bin/util/benchmarks/split_bench_rule_col.py {params.working_file} {output} > {log};"
        "sed -i -E 's/\t$/\tNA/' {output};"




localrules:
    aggregate_report_components,

rule aggregate_report_components:
    input:
        f"{MDIR}other_reports/giab_concordance_mqc.tsv",
        expand(
            MDIR
            + "{sample}/align/{alnr}/snv/{snv_caller}/bcfstats/{sample}.{alnr}.{snv_caller}.bcfstats.tsv",
            sample=SSAMPS,
            alnr=ALL_ALIGNERS,
            snv_caller=snv_CALLERS,
        ),
        f"{MDIR}other_reports/normcovevenness_combo_mqc.tsv",
        expand(
            MDIR + "{sample}/align/{alnr}/alignqc/cov_calcs_complete.done",
            sample=SSAMPS,
            alnr=ALL_ALIGNERS,
        ),
        expand(
            MDIR
            + "{sample}/align/{alnr}/alignqc/picard/picard/{sample}.{alnr}.done",
            sample=SSAMPS,
            alnr=ALL_ALIGNERS,
        ),
        expand(
            MDIR
            + "{sample}/align/{alnr}/alignqc/mosdepth/{sample}.{alnr}.mosdepth.summary.sort.bed",
            sample=SSAMPS,
            alnr=ALL_ALIGNERS,
        ),
        expand(
            MDIR + "{sample}/align/{alnr}/alignqc/goleft.done",
            sample=SSAMPS,
            alnr=ALL_ALIGNERS,
        ),
        expand(
            MDIR + "{sample}/align/{alnr}/alignqc/contam/vb2/{sample}.{alnr}.vb2.tsv",
            sample=SSAMPS,
            alnr=ALL_ALIGNERS,
        ),
        expand(
            MDIR + "{sample}/align/{alnr}/alignqc/qmap/{sample}.{alnr}/{sample}.{alnr}.qmap.done",
            sample=SSAMPS,
            alnr=ALL_ALIGNERS,
        ),
        expand(MDIR + "{sample}/align/{alnr}/alignqc/samtmetrics/{sample}.{alnr}.complete",
               sample=SSAMPS,
               alnr=ALL_ALIGNERS,
        ),
        MDIR + "other_reports/samtools_metrics_gather.done",
        "logs/peddy_gathered.done",
        f"{MDIR}other_reports/alignstats_combo_mqc.tsv",
        #f"{MDIR}logs/all_svVCF_dupheld.done",
        expand(
            MDIR
            + "{sample}/align/{alnr}/snv/{snv_caller}/vcf_stats/{sample}.{alnr}.{snv_caller}.rtg.vcfstats.txt",
            sample=SSAMPS,
            alnr=ALL_ALIGNERS,
            snv_caller=snv_CALLERS,
        ),
    threads: 2
    output:
        f"{MDIR}logs/report_components_aggregated.done",
    shell:
        "touch {output};" 


rule multiqc_final_wgs:  # TARGET: the big report
    input:
        f"{MDIR}logs/report_components_aggregated.done",
	    f"{MDIR}other_reports/rules_benchmark_data_mqc.tsv",
    output:
        f"{MDIR}reports/DAY_final_multiqc.html",
        f"{MDIR}reports/multiqc_header.yaml",
    benchmark:
        f"{MDIR}benchmarks/DAY_all.final_multiqc.bench.tsv"
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
        f"{MDIR}reports/logs/all__mqc_fin_a.log",
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
  ''' > {output[1]} >> {log} 2>&1;

        #cp $PWD/config/external_tools/multiqc_header.yaml {output[1]} >> {log} 2>&1;
        perl -pi -e "s/REGSUB_PROJECT/$DAY_PROJECT/g;"  {output[1]} >> {log} 2>&1;
        perl -pi -e "s/REGSUB_BUDGET/\\\$dbill$USED_BUDGET of \\\$dbill$TOTAL_BUDGET spent ( $PERCENT_USED\%)/g;" {output[1]} >> {log} 2>&1;

        size=$(du -hs results | cut -f1) >> {log} 2>&1;
        perl -pi -e "s/REGSUB_TOTALSIZE/$size/g;" {output[1]} >> {log} 2>&1;

        source bin/proc_spot_price_logs.sh >> {log} 2>&1;
        perl -pi -e "s/REGSUB_SPOTCOST/median: \\\$dbill$MEDIAN_SPOT_PRICE  mean: \\\$dbill$AVERAGE_SPOT_PRICE ( avg cost per vcpu,per min: \\\$dbill$VCPU_COST_PER_MIN ) /g;"  {output[1]} >> {log} 2>&1;
        perl -pi -e "s/REGSUB_SPOTINSTANCES/ $INSTANCE_TYPES_LINE /g;" {output[1]} >> {log} 2>&1;

        source bin/proc_aligner_costs.sh {input[1]} $VCPU_COST_PER_MIN >> {log} 2>&1;
        perl -pi -e "s/REGSUB_TOTALCOST/$ALNR_SUMMARY_COST/g;" {output[1]} >> {log} 2>&1;
        
        source bin/proc_mrkdup_costs.sh {input[1]} $VCPU_COST_PER_MIN  >> {log} 2>&1;
        perl -pi -e "s/REGSUB_MRKDUPCOST/$MRKDUP_AVG_MINUTES min, costing \\\$dbill$MRKDUP_AVG_COST/g;" {output[1]} >> {log} 2>&1;

        multiqc -f  \
        --config   ./$(dirname {output[0]})/multiqc_header.yaml \
        --config  ./config/external_tools/multiqc_config.yaml  \
        --custom-css-file ./config/external_tools/multiqc.css \
        --ignore "*/norm_cov_eveness/*" \
	--ignore "*sort_metrics/*" \
        --template default \
        --filename {output[0]} \
        -i '{params.rtitle} Multiqc Report ' \
        -b 'https://github.com/Daylily-Informatics/daylily (BRANCH:{params.gbranch}) (TAG:{params.gtag}) (HASH:{params.ghash}) ' \
        $(dirname {input} )/../ >> {log} 2>&1;
        ls -lt {output}  >> {log} 2>&1;
        """

localrules:
    produce_multiqc_final_wgs,


def get_fin_mqc(wildcards):
    from IPython import embed
    embed()

rule produce_multiqc_final_wgs:  # TARGET : Generated All WGS Reports
    input:
        MDIR+ "reports/DAY_final_multiqc.html"

