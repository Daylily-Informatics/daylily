##### KMER ANALYSIS OF UNALIGNED READS FOR QC
# -------------------------------------------
#
# a super interesting  tool can  tell you about
# contamination, lib prep biases, etc.
# home: https://github.com/TGAC/KAT


def get_cat_attempt(wildcards, attempt):
    if attempt in [None, 1,"1"]:
        return '-h'
    elif attempt in [2,"2"]:
        return ''
    else:
        return '--help'




rule kat:
    input:
        #expand(MDIR + "{{sample}}/{sample_lane}.dirsetup.ready", sample_lane=SAMPLE_LANES),
        fq1=getR1s,  # method
        fq2=getR2s,  # method
    output:
        done=MDIR + "{sample}/seqqc/kat/{sample}.kat.done",  # I'm not explicitly naming files b/c I allow failing of this rule in some cases
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
    container:
        "docker://quay.io/biocontainers/kat:2.4.2--py39h7c5ebd6_3"
    resources:
        attempt_n=get_cat_attempt,  # Hacking getting the attempt number from the res block for use in shell.
    shell:
        """
    
        kdir=$(dirname {output.done})/ktmp;
        rm -rf $kdir || echo 'No kat dir to remove';
        mkdir -p $kdir ;

        kat comp -v -t {threads} {resources.attempt_n} -n -p png -m 19 {input.fq1} {input.fq2};

        touch {output};

        """
 
localrules: produce_kat

rule produce_kat:  # TARGET: Produce kat
    input:
        expand(MDIR + "{sample}/seqqc/kat/{sample}.kat.done", sample=SAMPS),
    output:
        "kdone.log",
    shell:
        "touch {output};"
