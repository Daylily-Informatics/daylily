import os

##### fastp
# ---------------------------------
#  github: https://github.com/OpenGene/fastp
#  # Warning-  fastp has wacky ways of calculating tremplate sizes (only with overlapping read pairs!)

rule fastp:
    input:
        r1=get_raw_R1s,
        r2=get_raw_R2s,
        #r1=getR1s,  # method
        #r2=getR2s,  # method
    output:
        fp1=MDIR + "{sample}/seqqc/fastp/{sample}.R1.fastp.fastq.gz",
        fp2=MDIR + "{sample}/seqqc/fastp/{sample}.R2.fastp.fastq.gz",
        up1=MDIR + "{sample}/seqqc/fastp/{sample}.unpaired.R1.fastp.fastq.gz",
        up2=MDIR + "{sample}/seqqc/fastp/{sample}.unpaired.R2.fastp.fastq.gz",
        fail=MDIR + "{sample}/seqqc/fastp/{sample}.failed.fastp.fastq.gz",
        json=MDIR + "{sample}/seqqc/fastp/{sample}.fastp.json",
        html=MDIR + "{sample}/seqqc/fastp/{sample}.fastp.html",
        done=MDIR + "{sample}/seqqc/fastp/{sample}.fastp.done",
        proceed=MDIR + "{sample}/{sample}.proc",
        bench=MDIR + "{sample}/benchmarks/{sample}.fastp.bench.tsv",
    threads: config["fastp"]["threads"]
    resources:
        threads=config["fastp"]["threads"],
        partition=config["fastp"]["partition"]
    params:
        odir=MDIR + "{sample}/seqqc/fastp/",
        qscore_filter=config["fastp"]["filter_avg_q_score"],
        cluster_end="echo done"
        if "cluster_end" not in config["fastp"]
        else config["fastp"]["cluster_end"],
        cluster_sample=ret_sample,
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.fastp.bench.tsv"
    log:
        a=MDIR + "{sample}/logs/seqqc/fastp/{sample}.fastp.log",
    conda:
        config["fastp"]["env_yaml"]
    shell:
        """mkdir -p {params.odir} >> {log.a} 2>&1;
        fastp -i <(zcat {input.r1} )  -I <(zcat {input.r2} ) -o {output.fp1} -O {output.fp2} --unpaired1 {output.up1} --unpaired2 {output.up2}  --failed_out {output.fail} -q 10 -u 60 --trim_poly_g   --detect_adapter_for_pe   --dont_overwrite    --verbose  -D   --overrepresentation_analysis  --overrepresentation_sampling 100 --low_complexity_filter --detect_adapter_for_pe  --reads_to_process 100000000 --dup_calc_accuracy 6  --trim_tail1=1 -p -j {output.json} -h {output.html} -w {threads} -z 1  >> {log.a} ;"""
        "touch {output.proceed} {output.bench} {output.done};"
        "{latency_wait};"
        "ls {output}; "


localrules:
    produce_fastp,


rule produce_fastp:  # TARGET: produce fastp qc
    input:
        expand(MDIR + "{S}.proc", S=SSI),
    output:
        expand(MDIR + "{sample}/{sample}.fq.ready", sample=SAMPS),
    shell:
        "touch {output};"
