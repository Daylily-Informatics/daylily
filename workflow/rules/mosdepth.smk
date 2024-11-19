# ###### MOSDEPTH
#
# mosdepth
# github: https://github.com/brentp/mosdepth
# paper: https://academic.oup.com/bioinformatics/article/34/5/867/4583630


rule mosdepth:
    input:
        MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
    output:
        MDIR
        + "{sample}/align/{alnr}/alignqc/mosdepth/{sample}.{alnr}.mosdepth.summary.sort.bed",
    threads: config["mosdepth"]["threads"]
    resources:
        threads=config["mosdepth"]["threads"]
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.mosdepth.bench.tsv"
    log:
        a=MDIR + "{sample}/align/{alnr}/alignqc/mosdepth/{sample}.{alnr}.mosdepth.log",
        b=MDIR + "{sample}/align/{alnr}/alignqc/mosdepth/{sample}.{alnr}/{sample}.md",
    conda:
        config["mosdepth"]["env_yaml"]
    params:
        win_size=1000,
        mapq=0,
        T="0,10,20,30"
        if "depth_bins" not in config["mosdepth"]
        else config["mosdepth"]["depth_bins"],
        cluster_sample=ret_sample,
    shell:
        "rm -rf {log.b}* || echo rmlogFailedMosDepth;"
        "mosdepth --threads {threads} --by {params.win_size} --fast-mode --mapq {params.mapq} -T {params.T} {log.b} {input} > {log.a} 2>&1; "
        "touch {output};"
        "rm -rf $(dirname {log.b})/*per-base* || echo 'rm perbase failed' >> {log.a} 2>&1;"
        "{latency_wait}; ls {output};"
