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
        f"{MDIR}other_reports/rules_{{ru}}_{{ex}}_benchmarkdata_mqc.tsv",
    params:
        cluster_sample="rules_benchmark_collect",
        working_file=f"{MDIR}reports/benchmarks_summary.tsv",
        ref_code=config["genome_build"],
    log:
        f"{MDIR}other_reports/logs/{{ru}}_{{ex}}_benchmarks_summary.log",
    container: None
    shell:
        "bin/util/benchmarks/collect_day_benchmark_data.sh {params.ref_code} > {log};"
        "cp {params.working_file} {output} > {log};"


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
        f"{MDIR}other_reports/alignstats_bsummary.tsv",
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
        f"{MDIR}other_reports/rules_{{ru}}_{{ex}}_benchmarkdata_mqc.tsv",
    output:
        f"{MDIRreportsd}DAY_{{ru}}_{{ex}}_final_multiqc_hcwgs.html",
    benchmark:
        f"{MDIR}benchmarks/{{ru}}_{{ex}}_all.final_multiqc.bench.tsv"
    threads: config["multiqc"]["threads"]
    priority: 50
    params:
        odir2=MDIRreportsd,
        fnamef=f"DAY_{RU[0]}_{EX[0]}_final_multiqc_hcwgs.html",
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
        a=f"{MDIRreportsd}logs/{{ru}}_{{ex}}_mqc_fin_a.log",
        b=f"{MDIRreportsd}logs/{{ru}}_{{ex}}mqc_fin_b.log",
    conda:
        config["multiqc"]["final"]["env_yaml"]
    shell:
        """
        url_root=$(python -c "import os,sys; ap=os.path.abspath('.').replace('/fsx/data/','').replace('/','\/'); print(ap)");
        url_root2=$(python -c "import os,sys; ap2=os.path.abspath('.').replace('/fsx/data/','').replace('/','\/');print(ap2)");


        ( env bash workflow/scripts/create_final_mqc_report.sh '{params.odir2}' '{params.fnamef}' \
        '{params.mgroot}/{params.macro_cfg}' '{params.mgroot}/{params.micro_cfg}' \
        '{params.mgroot}/{params.name_cfg}' '{params.gtag}' '{params.ghash}' \
        '{params.gbranch}' '{params.jid}' '{params.ruu}' \
        '{params.exx}' '{params.cluster_sample}' \
        ' ' "$url_root" \
        '{params.mdir}'  \
        "{params.dkr}" "{params.ref}" '{log.b}' "$url_root2"  >> {log.a}  2>&1) || echo "MQC PROBLEM";
        ls -lt {output};
        """


localrules:
    produce_multiqc_final_wgs,


def get_fin_mqc(wildcards):
    from IPython import embed
    embed()

rule produce_multiqc_final_wgs:  # TARGET : Generated All WGS Reports
    input:
        expand(MDIRreportsd + "DAY_{ru}_{ex}_final_multiqc_hcwgs.html", ru=RU, ex=EX)
