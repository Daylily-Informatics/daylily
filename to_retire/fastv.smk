import os

##### FastV
# ---------------------------------
# github: https://github.com/OpenGene/fastv


rule fastv:
    input:
        fpqr1s=get_raw_R1s,
        fpqr2s=get_raw_R2s,
        #fpq1=getR1sS,  # MDIR + "{sample}/seqqc/fastp/{sample}.R1.fastp.fastq.gz",
        #fpq2=getR2sS,  # MDIR + "{sample}/seqqc/fastp/{sample}.R2.fastp.fastq.gz",
    output:
        fv1=MDIR + "{sample}/seqqc/fastv/{sample}.R1.fastv.fastq.gz",
        fv2=MDIR + "{sample}/seqqc/fastv/{sample}.R2.fastv.fastq.gz",
        html=MDIR + "{sample}/seqqc/fastv/{sample}.fastv.html",
        json=MDIR + "{sample}/seqqc/fastv/{sample}.fastv.json"
    threads: config["fastv"]["threads"]
    log:
        MDIR + "{sample}/seqqc/fastv/logs/{sample}.fastv.log",
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.fastv.bench.tsv"
    params:
        microbial_kmers=config["supporting_files"]["files"]["fastv"][
            "microbial_kmer_collection"
        ]["name"],
        covid_kmer=config["supporting_files"]["files"]["fastv"]["covid_kmer"]["name"],
        covid_genome=config["supporting_files"]["files"]["fastv"]["covid_genomes"][
            "name"
        ],
        jsf1="resources/fastv/js/coverage.js",
        jsf2="resources/fastv/js/plotly-1.2.0.min.js",
        odir=MDIR + "{sample}/seqqc/fastv/",
        cluster_sample=ret_sample,
    conda:
        config["fastv"]["env_yaml"]
    shell:
        "mkdir -p  {params.odir}  >> {log} 2>&1; "
        " cp {params.jsf1} {params.odir}; cp {params.jsf2} {params.odir} >> {log} 2>&1; "
        "fastv -i <(zcat {input.fpq1} ) -I <(zcat {input.fpq2} ) -o {output.fv1} -O {output.fv2} --reads_to_process 100000000 --detect_adapter_for_pe --low_complexity_filter -h {output.html} -j {output.json} -w {threads} -y -k {params.covid_kmer} -g {params.covid_genome} -c {params.microbial_kmers} >> {log} 2>&1;"
        " perl -pe 's/http\:\/\/opengene\.org\/plotly/plotly/g;' {output.html} >> {log} 2>&1; "
        " perl -pe 's/http\:\/\/opengene\.org\/fastv\///g;' {output.html} >> {log} 2>&1; "
        "{latency_wait};"
        "ls {output};"


rule produce_fastv:  # TARGET: fastv output
    input:
        expand(MDIR + "{sample}/seqqc/fastv/{sample}.fastv.json", sample=SAMPS),
    output:
        MDIR + "logs/multiqc/fastv.done",
    shell:
        "touch {output}"
