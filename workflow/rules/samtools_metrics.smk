import os
# ###### SAMTOOLS METIRCS
# ------------------------
# BAM Metrics From samtools
#

rule gen_samstats:
    input:
        cram=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.cram",
        crai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.cram.crai",
    output:
        stats=MDIR + "{sample}/align/{alnr}/alignqc/samtmetrics/{sample}.{alnr}.stats.tsv",
        flagstats=MDIR + "{sample}/align/{alnr}/alignqc/samtmetrics/{sample}.{alnr}.flagstat.tsv",
        idxstats=MDIR + "{sample}/align/{alnr}/alignqc/samtmetrics/{sample}.{alnr}.idxstat.tsv",
        sent=MDIR + "{sample}/align/{alnr}/alignqc/samtmetrics/{sample}.{alnr}.complete",
    threads: 8
    conda:
        config["samtools_markdups"]["env_yaml"]
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.samt.bench.tsv"
    params:
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        cluster_sample=ret_sample,
    shell:
        "samtools stats -@  {threads} --reference {params.huref} {input.cram} > {output.stats};"
        "samtools flagstats -@  {threads} {input.cram} > {output.flagstats};"
        "samtools idxstats -@  {threads} {input.cram} > {output.idxstats};"
        "touch {output.sent};"
        "{latency_wait}; ls {output};"

localrules: produce_samtools_metrics

rule produce_samtools_metrics:  #TARGET: produce samtools BAM metrics
    input:
        expand(MDIR + "{sample}/align/{alnr}/alignqc/samtmetrics/{sample}.{alnr}.complete",
            sample=SAMPS,
            alnr=CRAM_ALIGNERS
            )
    output:
        touch(MDIR + "other_reports/samtools_metrics_gather.done"),
