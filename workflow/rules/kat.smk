##### KMER ANALYSIS OF UNALIGNED READS FOR QC
# -------------------------------------------
#
# a super interesting  tool can  tell you about
# contamination, lib prep biases, etc.
# home: https://github.com/TGAC/KAT


def get_cat_attempt(wildcards, attempt):
    return int(attempt)

if config['rule_action']['kat'] in ['skip']:
    localrules: kat_skip

    rule kat_skip:
        output:
            expand(MDIR + "{sample}/seqqc/kat/{sample}.kat.done", sample=SAMPS)
        shell:
            "mkdir  -p $(dirname {output[0]} ) || echo mkdirFailed;"
            "touch {output};"

else:
    rule kat:
        input:
            #expand(MDIR + "{{sample}}/{sample_lane}.dirsetup.ready", sample_lane=SAMPLE_LANES),
            fq1=getR1s,  # method
            fq2=getR2s,  # method
        output:
            done=MDIR + "{sample}/seqqc/kat/{sample}.kat.done",  # I'm not explicitly naming files b/c I allow failing of this rule in some cases
            r1r2_stub=MDIR + "{sample}/seqqc/kat/r1r2cmp",
            ds1=temp(MDIR + "{sample}/seqqc/kat/{sample}.kat.r1.fq.gz"),
            ds2=temp(MDIR + "{sample}/seqqc/kat/{sample}.kat.r2.fq.gz"),
        benchmark:
            MDIR + "{sample}/benchmarks/{sample}.kat.bench.tsv"
        threads: config["kat"]["threads"]
        params:
            cluster_sample=ret_sample,
            kmer_len="19" if 'kmer' not in config['kat'] else config['kat']['kmer'],
            p5trim="8" if 'p5trim' not in config['kat'] else config['kat']['p5trim'],
            subsample="0.2",
            print_every="50",
        log:
            MDIR + "{sample}/seqqc/kat_logs/{sample}.kat.log",
        conda:
            config["kat"]["env_yaml"]
        resources:
            attempt_n=get_cat_attempt,  # Hacking getting the attempt number from the res block for use in shell.
        shell:
            """
            (
            touch {output};
            exit 0;

            kdir=$(dirname {output.r1r2_stub} );
            rm -rf $kdir ;
            mkdir -p $kdir ;
            touch {output.r1r2_stub};
            touch {output};
            seqkit sample --rand-seed 73 -j {threads} -p  {params.subsample} {input.fq1} > {output.ds1} &
            seqkit sample --rand-seed 73 -j {threads} -p {params.subsample} {input.fq2} > {output.ds2} &
            wait;
            if [[ {resources.attempt_n} > 1 ]]; then echo 'Kat has failed 1x, no further tries will be attempted, but since this is non-critical, we are letting the node appear to succeed.' >> {log}.multiattempt.log 2>&1; exit 0; fi;

            kat comp -v -t {threads} -h -n -p png -g -o $kdir/{params.cluster_sample}.r1r2cmp {output.ds1} {output.ds2};

            touch {output};
            {latency_wait}; ls {output}; ) >> {log} 2>&1
            """
# May be able to pulll the same trick with kat <(seqkit sample --proportion {params.subsample_pct} <(seqfu interleave -1 <(unpigz -c -q -- {input.fqr1s}) -2 <(unpigz -c -q -- {input.fqr2s}) ) )

localrules: produce_kat

rule produce_kat:  # TARGET: Produce kat
    input:
        expand(MDIR + "{sample}/seqqc/kat/{sample}.kat.done", sample=SAMPS),
    output:
        "kdone.log",
    shell:
        "touch {output};"
