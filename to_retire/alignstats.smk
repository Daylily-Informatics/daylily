##### ALIGNSTATS
# --------------
#
# It blows my mind this tool is not uniformly adopted.
# It calculates 300+ metrices relevant to monitoring NGS
# quality, and does it all in one place, often with more detail
# included than the standard tooling.  I rely on it
# for a ton of QC work. The Baylor Genome Center developed it
# github: https://github.com/jfarek/alignstat


def fetch_alnr(wildcards):
    return wildcards.alnr


rule alignstats:
    input:
        MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
    output:
        json=MDIR
        + "{sample}/align/{alnr}/alignqc/alignstats/{sample}.{alnr}.alignstats.json",
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.alignstats.bench.tsv"
    threads: config["alignstats"]["threads"]
    resources:
        attempt_n=lambda wildcards, attempt:  (attempt + 0),
        partition=config["alignstats"]["partition"],
        threads=config["alignstats"]["threads"],
        vcpu=config["alignstats"]["threads"]
    log:  MDIR + "{sample}/align/{alnr}/alignqc/alignstats/logs/{sample}.{alnr}.alignstats.log",
    params:
        P=50,
        n=config["alignstats"]["num_reads_in_mem"],
        cluster_sample=ret_sample,
        ld_preload=" "
        if "ld_preload" not in config["malloc_alt"]
        else config["malloc_alt"]["ld_preload"],
        ld_pre=" "
        if "ld_preload" not in config["alignstats"]        else config["alignstats"]["ld_preload"],
    conda:
        config["alignstats"]["env_yaml"]
    shell:
        "LD_LIBRARY_PATH=resources/lib resources/alignstats/alignstats  -C -U  -i {input} -o {output.json}  -j bam -v -P {threads} -p {threads} > {log};"


localrules:
    finish_align_stats,

rule finish_align_stats:
    input:
        json=MDIR
        + "{sample}/align/{alnr}/alignqc/alignstats/{sample}.{alnr}.alignstats.json",
    output:
        tsv=MDIR
        + "{sample}/align/{alnr}/alignqc/alignstats/{sample}.{alnr}.alignstats.tsv",
    params:
        P=50,
        n=config["alignstats"]["num_reads_in_mem"],
        cluster_sample=ret_sample,
        alnr_f=fetch_alnr,
        ld_preload=" "
        if "ld_preload" not in config["malloc_alt"]
        else config["malloc_alt"]["ld_preload"],
        ld_pre=" "
        if "ld_preload" not in config["alignstats"]
        else config["alignstats"]["ld_preload"],
    run:
        import os
        import sys
        import json

        j = json.load(open(f"{input.json}", "r"))
        aa = "sample\taligner\t" + "\t".join([str(x) for x in sorted(j)])
        bb = f"{params.cluster_sample}\t{params.alnr_f}\t" + "\t".join(
            [str(j[x]) for x in sorted(j)]
        )
        os.system(f"echo {aa} > {output.tsv} ; echo {bb} >> {output.tsv}")

