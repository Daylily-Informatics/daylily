import os
###### Piocard QC Tools
# ----------------------------------

if os.environ.get("DAY_CRAM","") == "":

    rule picard:
        input:
            bam=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
        output:
            sent=touch(MDIR + "{sample}/align/{alnr}/alignqc/picard/picard/{sample}.{alnr}.mrkdup.sort.picard.done"),
        threads: config["picard"]["threads"]
        resources:
            vcpu=config["picard"]["threads"],
            partition=config["picard"]["partition"],
        benchmark:
            MDIR + "{sample}/benchmarks/{sample}.{alnr}.mrkdup.sort.picard.bench.tsv"
        log:
            MDIR + "{sample}/align/{alnr}/alignqc/picard/logs/{sample}.{alnr}.picard.stats.log"
        conda:
            "../envs/picard_v0.1.yaml"
        params:
            huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
            cluster_sample=ret_sample,
            metric_accumulation_level="SAMPLE" if 'metric_accumulation_level' not in config['picard'] else config['picard']['metric_accumulation_level'],
            stop_after="2000000" if 'stop_after' not in config['picard'] else config['picard']['stop_after'],
            validation_stringency="LENIENT" if 'validation_stringency' not in config['picard'] else config['picard']['validation_stringency']
        shell:
            """
            set +euo pipefail;
            (
            rm -rf $(dirname {output.sent} ) || echo rmPicardFailed;
            mkdir -p $(dirname {log} );
            pic_d="$(dirname {output.sent} )/picard/";
            mkdir -p $pic_d;
            picard CollectMultipleMetrics VALIDATION_STRINGENCY={params.validation_stringency} METRIC_ACCUMULATION_LEVEL={params.metric_accumulation_level} LEVEL={params.metric_accumulation_level} STOP_AFTER={params.stop_after} I={input.bam} O=$pic_d R={params.huref} EXT=.txt PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution PROGRAM=CollectGcBiasMetrics PROGRAM=CollectSequencingArtifactMetrics PROGRAM=CollectQualityYieldMetrics INCLUDE_UNPAIRED=true ;
            touch {output.sent}; echo empySentFilesNotBeingSeen >> {output.sent};
            {latency_wait}; ls {output.sent}; cat {output.sent} ; ) > {log} 2>&1;
            exit 0;
            """

else:

    rule picard_cram:
        input:
            cram=MDIR + "{sample}/align/{alnr}/{sample}.cram",
        output:
            sent=touch(MDIR + "{sample}/align/{alnr}/alignqc/picard/picard/{sample}.{alnr}.mrkdup.sort.picard.done"),
        threads: config["picard"]["threads"]
        resources:
            vcpu=config["picard"]["threads"],
            partition=config["picard"]["partition"],
        benchmark:
            MDIR + "{sample}/benchmarks/{sample}.{alnr}.mrkdup.sort.picard.bench.tsv"
        log:
            MDIR + "{sample}/align/{alnr}/alignqc/picard/logs/{sample}.{alnr}.picard.stats.log"
        conda:
            "../envs/picard_v0.1.yaml"
        params:
            cluster_sample=ret_sample,
            huref=config["supporting_files"]["files"]["huref"]["broad_fasta"]["name"],
            metric_accumulation_level="SAMPLE" if 'metric_accumulation_level' not in config['picard'] else config['picard']['metric_accumulation_level'],
            stop_after="2000000" if 'stop_after' not in config['picard'] else config['picard']['stop_after'],
            validation_stringency="LENIENT" if 'validation_stringency' not in config['picard'] else config['picard']['validation_stringency']
        shell:
            """
            set +euo pipefail;
            (
            rm -rf $(dirname {output.sent} ) || echo rmPicardFailed;
            mkdir -p $(dirname {log} );
            pic_d="$(dirname {output.sent} )/picard/";
            mkdir -p $pic_d;
            picard CollectMultipleMetrics VALIDATION_STRINGENCY={params.validation_stringency} METRIC_ACCUMULATION_LEVEL={params.metric_accumulation_level} LEVEL={params.metric_accumulation_level} STOP_AFTER={params.stop_after} I={input.bam} O=$pic_d R={params.huref} EXT=.txt PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution PROGRAM=CollectGcBiasMetrics PROGRAM=CollectSequencingArtifactMetrics PROGRAM=CollectQualityYieldMetrics INCLUDE_UNPAIRED=true ;
            touch {output.sent}; echo empySentFilesNotBeingSeen >> {output.sent};
            {latency_wait}; ls {output.sent}; cat {output.sent} ; ) > {log} 2>&1;
            exit 0;
            """

localrules: produce_picard,

rule produce_picard:  # $ TARGET: produce picard QC data
    input:
        expand(MDIR + "{sample}/align/{alnr}/alignqc/picard/{sample}.{alnr}.mrkdup.sort.picard.done", sample=SSAMPS,alnr=ALIGNERS)
