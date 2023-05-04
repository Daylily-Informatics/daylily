0#### ENSEMBL VEP
# -------------------------------------
# github: https://github.com/Ensembl/ensembl-vep
# docker: https://hub.docker.com/r/ensemblorg/ensembl-vep
 docker://ensemblorg/ensembl-vep

rule vep:
    input:
        vcfgz=MDIR
        + "{sample}/align/{alnr}/snv/{snv}/{sample}.{alnr}.{snv}.snv.sort.vcf.gz",
        # Generalize this for other var callers
    output:
        prefix=MDIR
        + "{sample}/align/{alnr}/snv/{snv}/peddy/{sample}.{alnr}.{snv}.peddy.",
        done=MDIR
        + "{sample}/align/{alnr}/snv/{snv}/peddy/{sample}.{alnr}.{snv}.peddy.done",
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/{snv}/peddy/log/{sample}.{alnr}.{snv}.peddy.log",
    threads: config["peddy"]["threads"]
    resources:
        vcpu=config["peddy"]["threads"],
    params:
        cluster_sample=ret_sample,
        ld_preload=config["malloc_alt"]["ld_preload"],
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.{snv}.peddy.bench.tsv"
    container:
        None
    conda:
        config["peddy"]["env_yaml"]
    shell:
        """
        mkdir -p $(dirname {log}) > {log} ;
        echo 'runnning peddy' >> {log} ;
        set +euo pipefail;
        ( ({params.ld_preload} peddy  -p {threads}  --plot --prefix {output.prefix} --loglevel DEBUG {input.vcfgz} {input.ped_f}) || (ls . && sleep 1 && echo 'peddy exited with an error.  If this is a small depth of coverate vs. the genome sample, then this may be a pca failure which can be ignored, run the command by hand to see the detailed error') ) || (sleep 1 && echo 'masking error' ) >> {log};
        set -euo pipefail;
        echo "DONE" >> {log} ;
        echo reallydone > {output.done} ;
        echo reallydone > {output.prefix} ;
        touch {output.done} ;
        touch {output.prefix} ;
        sleep 4;
        """


localrules:
    produce_peddy,


rule produce_peddy:  # TARGET: just produce peddy results
    input:
        expand(
            MDIR
            + "{sample}/align/{alnr}/snv/{snv}/peddy/{sample}.{alnr}.{snv}.peddy.done",
            sample=SSAMPS,
            alnr=ALIGNERS,
            snv=snv_CALLERS,
        ),
    output:
        "logs/peddy_gathered.done",
    shell:
        "touch {output};"
