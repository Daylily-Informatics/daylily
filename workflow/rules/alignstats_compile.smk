import os

# This rule is actually gathering up all of the sub rules completed states
# so it can generate a final QC report, and that report will satisfy the input
# requirement of all to run



localrules:
    alignstats_gather,

rule alignstats_gather:
    input:
        expand(
            MDIR
            + "{sample}/align/{alnr}/alignqc/alignstats/{sample}.{alnr}.alignstats.tsv",
            sample=SSAMPS,
            alnr=ALL_ALIGNERS,
            caller=sv_CALLERS,
        ),
    output:
        f"{MDIR}other_reports/alignstats_summary_gather.done",
        bm=MDIR + "benchmarks/all.alignstats_summary.bench.tsv",
    benchmark:
        MDIR + "benchmarks/all.alignstats_summary.bench.tsv"
    threads: 1
    log:
        MDIR + "logs/alignstats_summary_gather.log",
    conda:
        config["alignstats"]["env_yaml"]
    shell:
        " touch {output[0]}; touch {output.bm}"



localrules:
    alignstats_compile,


rule alignstats_compile:
    input:
        f"{MDIR}other_reports/alignstats_summary_gather.done",
    output:
        temp(f"{MDIR}other_reports/alignstats_bsummary.tsv"),
        temp(f"{MDIR}other_reports/alignstats_csummary.tsv"),
        f"{MDIR}other_reports/alignstats_combo_mqc.tsv",
        f"{MDIR}other_reports/alignstats_gs_mqc.tsv",        
    benchmark:
        MDIR + "benchmarks/all.alignstats_smmary_compile.bench.tsv"
    threads: 2
    params:
        l="{",
        r="}",
        cluster_sample="na",
    log:
        MDIR + "logs/alignstats_summary_compile.log",
    conda:
        config["alignstats"]["env_yaml"]
    shell:
        "(mkdir -p {MDIR}logs/as;"
        "echo STARTcompileAstats > {log};"
        "find {MDIR}*/align/*/alignqc/alignstats/*alignstats.tsv | head -n 1 | parallel -j 1 'head -n 1 {params.l}{params.r} > {output[0]}; echo a_{params.l}{params.r} >> {log}' >> {log} 2>&1; "
        "find {MDIR}*/align/*/alignqc/alignstats/*alignstats.tsv | parallel -j 1 'tail -n 1 {params.l}{params.r} >> {output[0]}; echo b_{params.l}{params.r}  >> {log} ' >> {log} 2>&1;  "
        "cp {output[0]} {output[1]};"
        "cp {output[0]} {output[2]};" #         "perl -pi -e 's/_DBC0_0//g;' {output};"
        ") || touch logs/ALIGNSTATSCOMPIEFAILEDw_$? ; "
        "perl -pi -e 's/ /\t/g;' {output[2]};"
        "cp {output[2]} {output[3]};"


localrules:
    produce_alignstats,


rule produce_alignstats:  # TARGET - only takes path to produce alignstats
    input:
        f"{MDIR}other_reports/alignstats_bsummary.tsv",
