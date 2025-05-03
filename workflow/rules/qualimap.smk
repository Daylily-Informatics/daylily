import os


def _get_queue():
    # Come back and have this threshold set by if the BAM is > 19G or not
    if "is_lcwgs" in config["qualimap"]:
        if config["qualimap"]["is_lcwgs"] in "lcwgs":
            return " -q dev-short "
    return " -q dev-long "

rule qualimap:
    """Run Qualimap on CRAMs"""
    input:
        cram=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.cram",
        crai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.cram.crai",
    output:
        d=MDIR + "{sample}/align/{alnr}/alignqc/qmap/{sample}.{alnr}/{sample}.{alnr}.qmap.done",
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.qmap.bench.tsv"
    resources:
        vcpu=config["qualimap"]["threads"],
        threads=config["qualimap"]["threads"],
        partition=config["qualimap"]["partition"],
    params:
        java_mem_size=config["qualimap"]["java_mem_size"],
        cluster_sample=ret_sample,
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
    conda:
        config["qualimap"]["env_yaml"]
    threads: config["qualimap"]["threads"]
    log:
        MDIR + "{sample}/align/{alnr}/alignqc/qmap/{sample}/logs/{sample}.{alnr}.qmap.log",
    shell:
        """
        set +euo pipefail;
        rm -rf $(dirname {output.d} ) || echo rmFailedQMAP ;
        mkdir -p $(echo $(dirname {output.d} )/logs ) ;
        export dn=$(dirname {output.d} );
        qualimap bamqc -bam {input.cram} -nt {threads} -c -gd {params.huref}  -outformat HTML  -outdir $dn --java-mem-size={params.java_mem_size} > {log} 2>&1   || echo QMAPERROR;
        touch {output.d};
        ls {output.d};
        """
