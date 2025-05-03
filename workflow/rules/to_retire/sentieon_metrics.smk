# #######  SENTIEON  METRICS TOOLS
# --------------------------------
# Sentieon Metrics Generation Utils
#


rule sentieon_qc_metrics:
    input:
        bam=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
        bam_bai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam.bai",
    output:
        wgs_mets=MDIR + "{sample}/align/{alnr}/alignqc/sentmetrics/sm/{sample}.{alnr}.mrkdup.wgsmets.tsv",
        i_mets=MDIR + "{sample}/align/{alnr}/alignqc/sentmetrics/sm/{sample}.{alnr}..mrkdup.insertmets.tsv",
        gc_summary=MDIR
        + "{sample}/align/{alnr}/alignqc/sentmetrics/sm/{sample}.{alnr}.mrkdup.gcsummary.tsv",
        gc_mets=MDIR + "{sample}/align/{alnr}/alignqc/sentmetrics/sm/{sample}.{alnr}.mrkdup.gcmets.tsv",
        alnstats=MDIR + "{sample}/align/{alnr}/alignqc/sentmetrics/sm/{sample}.{alnr}.mrkdup.alnstats.tsv",
        qualyield=MDIR
        + "{sample}/align/{alnr}/alignqc/sentmetrics/sm/{sample}.{alnr}.mrkdup.qualyield.tsv",
        complete=MDIR
        + "{sample}/align/{alnr}/alignqc/sentmetrics/sm/{sample}.{alnr}.mrkdup.metrics.complete",
    threads: 10 if "metrics_threads" not in config["sentieon"] else config["sentieon"]["metrics_threads"]
    log:
        MDIR + "{sample}/align/{alnr}/alignqc/sentmetrics/sm/logs/{sample}.{alnr}.metrics.log",
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.metrics.sm.bench.tsv"
    params:
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        cluster_sample=ret_sample,
    conda:
        config["sentieon"]["env_yaml"]
    shell:
        """

        if [ -z "$SENTIEON_LICENSE" ]; then
            echo "SENTIEON_LICENSE not set. Please set the SENTIEON_LICENSE environment variable to the license file path & make this update to your dyinit file as well.";
            exit 3;
        fi

        if [ ! -f "$SENTIEON_LICENSE" ]; then
            echo "The file referenced by SENTIEON_LICENSE ('$SENTIEON_LICENSE') does not exist. Please provide a valid file path.";
            exit 4;
        fi
        
        set +euo pipefail;
        (
        rm -rf $(dirname {output.complete} ) ;
        mkdir -p $(dirname {output.complete} ) || echo noFail;
        sentieon driver -i {input.bam} -r {params.huref} -t {threads} --algo WgsMetricsAlgo --min_base_qual  10 --min_map_qual 10 --include_unpaired true {output.wgs_mets};
        sentieon driver  -i {input.bam} -r {params.huref} -t {threads} --algo InsertSizeMetricAlgo {output.i_mets};
        sentieon driver  -i {input.bam} -r {params.huref} -t {threads} --algo QualityYield {output.qualyield} ;
        sentieon driver -t {threads} -i {input.bam} -r {params.huref} --algo GCBias --summary {output.gc_summary} {output.gc_mets} ;
        sentieon driver -i {input.bam} -t {threads} -r {params.huref} --algo AlignmentStat {output.alnstats} ;
        echo 'complete' > {output.complete};
        {latency_wait}; ls {output};) > {log} 2>&1;
        """




rule produce_sent_metrics:  # TARGET: produce just sentieon BAM metrics
    input:
        MDIR + "{sample}/align/{alnr}/alignqc/sentmetrics/{sample}.{alnr}.mrkdup.metrics.complete",
    container: None
    output:
        expand(
            MDIR + "{sample}/align/{alnr}/alignqc/sentmetrics/{sample}.{alnr}.mrkdup.metrics.complete.2",
            sample=SSAMPS,
            alnr=ALIGNERS
        ),
    shell:
        "touch {output}"
