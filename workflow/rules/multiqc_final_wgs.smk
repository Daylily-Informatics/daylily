import os
# This rule is actually gathering up all of the sub rules completed states
# so it can generate a final QC report, and that report will satisfy the input
# requirement of all to run

RU = [""]
EX = [""]
if len(samples) > 0:
    RU = samples.RU
    EX = samples.EX


localrules:
    collect_rules_benchmark_data,


rule collect_rules_benchmark_data:
    input:
        f"{MDIR}logs/report_components_aggregated.done",
    output:
        f"{MDIR}other_reports/rules_benchmarkdata_mqc.tsv",
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
        #f"{MDIR}other_reports/giabhcr_concordance_mqc.tsv",
        f"{MDIR}other_reports/norm_cov_evenness_mqc.tsv",
        expand(
            MDIR + "{sample}/align/{alnr}/alignqc/cov_calcs_complete.done",
            sample=SSAMPS,
            alnr=ALIGNERS,
            snv_caller=snv_CALLERS,
        ),
        expand(
            MDIR
            + "{sample}/align/{alnr}/snv/{snv_caller}/bcfstats/{sample}.{alnr}.{snv_caller}.bcfstats.tsv",
            sample=SSAMPS,
            alnr=ALIGNERS,
            snv_caller=snv_CALLERS,
        ),
        "logs/peddy_gathered.done",
        f"{MDIR}other_reports/alignstats_combo_mqc.tsv",
        #f"{MDIR}logs/all_svVCF_dupheld.done",
        f"{MDIRreportsd}SEQQC_multiqc.html",
        expand(
            MDIR
            + "{sample}/align/{alnr}/snv/{snv_caller}/vcf_stats/{sample}.{alnr}.{snv_caller}.rtg.vcfstats.txt",
            sample=SSAMPS,
            alnr=ALIGNERS,
            snv_caller=snv_CALLERS,
        ),
    output:
        f"{MDIR}logs/report_components_aggregated.done",
    shell:
        "touch {output};"


rule multiqc_final_wgs:  # TARGET: the big report
    input:
        f"{MDIR}other_reports/rules_benchmarkdata_mqc.tsv",
    output:
        f"{MDIR}reports/DAY_final_multiqc.html",
    benchmark:
        f"{MDIR}benchmarks/DAY_all.final_multiqc.bench.tsv"
    threads: config["multiqc"]["threads"]
    priority: 50
    params:
        odir2=MDIRreportsd,
        fnamef=f"DAY_final_multiqc.html",
        macro_cfg=config["multiqc"]["config_yaml"],
        micro_cfg=config["multiqc"]["final"]["config_yaml"],
        name_cfg=config["multiqc_sampname_cfg"],
        gtag=config["gittag"],
        ghash=config["githash"],
        gbranch=config["gitbranch"],
        jid="COMP-1" if "jid" not in config else config["jid"],
        ruu=f"{RU[0]}",
        exx=f"{EX[0]}",
        cluster_sample=f"{RU[0]}_{EX[0]}",
        ld_pre=" "
        if "ld_preload" not in config["multiqc"]
        else config["multiqc"]["ld_preload"],
        dkr="none" if "default_container" not in config else config["default_container"],
        mdir=MDIR,
        ref=config["ref_code"],
        mgroot=os.environ['DAY_ROOT'],
    log:
        f"{MDIR}reports/logs/all__mqc_fin_a.log",
    conda:
        config["multiqc"]["final"]["env_yaml"]
    shell:
        """
        multiqc -f  \
        --config config/external_tools/multiqc_header.yaml  \
        --config  config/external_tools/multiqc_config.yaml  \
        --filename {output} \
        -i 'Final Multiq Report' \
        -b 'Git Info: {params.gbranch} {params.gtag} {params.ghash}' \
        --outdir {params.odir2}  $(dirname {input} )/../ > {log} 2>&1;
        ls -lt {output};
        """


localrules:
    produce_multiqc_final_wgs,


def get_fin_mqc(wildcards):
    from IPython import embed
    embed()

rule produce_multiqc_final_wgs:  # TARGET : Generated All WGS Reports
    input:
        MDIR+ "reports/DAY_final_multiqc.html"
