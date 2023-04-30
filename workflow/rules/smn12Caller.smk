# ###### Gauchian SMN12
#
# SMN12 Caller
# 
# github: see donda def file for link to github


rule smn12caller:
    input:
        bam=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
        bai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam.bai",
    output:
        manifest=MDIR + "{sample}/align/{alnr}/htd/smn12caller/{sample}.smn12_manifest.csv"
    params:
        cluster_sample=ret_sample,
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.smn12caller.bench.tsv"
    log:
        MDIR + "{sample}/align/{alnr}/htd/smn12caller/logs/smn12caller.log",
    threads: config["go_left"]["threads"]
    conda:
        "workflow/envs/smn12_v0.1.yaml"
    shell:
        """
        rm -rf $(dirname {output.manifest} ) || echo rmGOLfailed ;
        mkdir -p "$(dirname {output.manifest} )/logs" ; 
        echo "{input.bam}" > {output.manifest};
	smn_caller.py --manifest {output.manifest}  --genome 37 --outdir $( dirname {output.manifest}) --prefix {sample}.SMN12caller --threads {threads}      	      
        touch {output.done};;
	ls {output};
        """


localrules:

rule produce_smn12:
    input:
        expand(MDIR + "{sample}/align/{alnr}/htd/gba_gauchian/{sample}.{alnr}.gba_gauchian.manifest", sample=SSAMPS, alnr=ALIGNERS)
    output:
        "./logs/gba.done"
    shell:
        """
        touch {output};
        """
