import os

# #### fastqc
# -----------
# github: https://github.com/s-andrews/FastQC


rule fastqc_subsampled:
    input:
        fqr1s=get_raw_R1s,
        fqr2s=get_raw_R2s,
    output:
        f"{MDIR}" + "{sample}/seqqc/fastqc/{sample}.fastqc.done",
    benchmark:
        f"{MDIR}" + "{sample}/benchmarks/{sample}.fastqc.bench.tsv"
    threads: config["fastqc"]["threads"]
    resources:
        vcpu=config["fastqc"]["threads"],
        partition=config['fastqc']['partition'],
    params:
        tmp=f"{MDIR}" + "{sample}/seqqc/fastqc/tmp",
        tool_dir=f"{MDIR}" + "{sample}/seqqc/fastqc",
        cluster_sample=ret_sample,
        subsample_pct="0.25" if 'subsample_pct' not in config['fastqc'] else config['fastqc']['subsample_pct'],
    log:
        f"{MDIR}" + "{sample}/logs/fastqc/{sample}.fastqc.log",
    conda:
        config["fastqc"]["env_yaml"]
    shell:
        """
        rm -rf {params.tool_dir} ;
        mkdir -p {params.tool_dir} ;
        mkdir -p {params.tmp}  ;
        #fastqc -o {params.tool_dir} -t {threads} -d {params.tmp}  <(seqkit sample --proportion {params.subsample_pct} <(seqfu interleave -1 <(unpigz -c -q -- {input.fqr1s}) -2 <(unpigz -c -q -- {input.fqr2s}) ) )  ;
        fastqc -o {params.tool_dir} -t {threads} -d {params.tmp}  {input.fqr1s}  {input.fqr2s}
        touch {output};
        touch {output}_subsampled_at_{params.subsample_pct};
        ls {output};
        """
localrules:
    just_fastqc,


rule just_fastqc:
    input:
        expand(MDIR + "{sample}/seqqc/fastqc/{sample}.fastqc.done", sample=SAMPS),
    output:
        "fqc.done",
    threads: 1
    shell:
        "touch {output}"
