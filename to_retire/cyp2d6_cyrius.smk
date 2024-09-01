# ###### Gauchian Cyrius
#
# Cyrius Caller
# github: see conda env for link


rule cp2d6_cyrius:
    input:
        bam=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
        bai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam.bai",
    output:
        manifest=MDIR + "{sample}/align/{alnr}/htd/cyp2d6/cyp2d6_cyrius.manifest",
    params:
        cluster_sample=ret_sample,
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.96cyriusbench.tsv"
    log:
        MDIR + "{sample}/align/{alnr}/htd/cyp2d6/logs/cyp2d6.log",
    threads: config["go_left"]["threads"]
    conda:
         "workflow/envs/cyrius_v0.1.yaml"
    shell:
        """
	mkdir -p $(dirname {output.manifest});
        echo "{input.manifest}" >  {output.manifest};
        star_caller.py -m {output.manifest}  -g 37 -o outdir $(dirname {output.manifest}) -p {sample}.cyp2d6_cyrius  -t {threads};
        ls {output};
        """

localrules: produce_cyp2d6,
rule produce_cyp2d6:
    input:
        expand(MDIR + "{sample}/align/{alnr}/htd/cyp2d6/cyp2d6_cyrius.manifest",  sample=SSAMPS, alnr=ALIGNERS)
    output:
        "./logs/gba.done"
    shell:
        """
        touch {output}
        """
