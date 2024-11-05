#########  SENTIEON DEDUPING
# --------------------------

if "sent" in DDUP:

    def get_input_bams(wildcards):
        return MDIR + f"{wildcards.sample}/align/{wildcards.alnr}/{wildcards.sample}.{wildcards.alnr}.sort.bam",

    rule sentieon_markdups:
        """Runs duplicate marking on the BAM."""

        input:
            MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.sort.bam",
            MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.sort.bam.bai",
        output:
            bamo="{MDIR}{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
            baio="{MDIR}{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam.bai",
            score="{MDIR}{sample}/align/{alnr}/score.txt",
            metrics="{MDIR}{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.metrics",
        threads: config["sentieon"]["threads"]
        benchmark:
            repeat("{MDIR}{sample}/benchmarks/{sample}.{alnr}.mrkdup.bench.tsv", 0)
        conda:
            config["sentieon"]["env_yaml"]
        params:
            cluster_sample=ret_sample,
            max_mem="150G"
            if "max_mem" not in config["sentieon"]
            else config["sentieon"]["max_mem"],
            numactl=config["sentieon"]["numactl"],
        resources:
            threads=config['sentieon']['threads'],
            partition=config['sentieon']['partition'],
            vcpu=config['sentieon']['threads'],
        log:
            "{MDIR}{sample}/align/{alnr}/logs/dedupe.{sample}.{alnr}.log",
        shell:
            """
            /fsx/data/cached_envs/sentieon-genomics-202308.03/bin/sentieon driver -i {input[0]} -t {threads} \
              --algo LocusCollector {output.score} > {log}; 
            /fsx/data/cached_envs/sentieon-genomics-202308.03/bin/sentieon driver -i {input[0]} -t {threads} \
              --algo Dedup --metrics {output.metrics} \
              --bam_compression 9 \
              --optical_dup_pix_dist 100 \
              --score_info {output.score} {output.bamo} >> {log};
            samtools index -b -@ 4 -o {output.baio} {output.bamo} >> {log};
            """
