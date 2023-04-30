# ###### Gauchian GBA
#
# GBA Caller
# gauchian
# github: https://github.com/Illumina/Gauchian.git


rule gauchian:
    input:
        bam=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
        bai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam.bai",
    output:
        manifest=MDIR + "{sample}/align/{alnr}/htd/gba_gauchian/{sample}.{alnr}.gba_gauchian.manifest",
    params:
        cluster_sample=ret_sample,
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.gba_gauchiant.bench.tsv"
    log:
        MDIR + "{sample}/align/{alnr}/htd/gba_gauchian/logs/gba_gauchian.log",
    threads: config["go_left"]["threads"]
    conda:
         "workflow/envs/gba_v0.1.yaml"
    shell:
        """
	rmdir -rf $(dirname {output.manifest});
	mkdir -p $(dirname {output.manifest});
	echo '{input.bam}' > {outpus.manifest};
	gauchian -m {output.manifest} -g 37 -o $(dirname {output.manifest}) -p {sample}.{alnr}.gbaGauchian > {log};                     
        {latency_wait};
        """


localrules: produce_gauchian

rule produce_gauchian:
    input:
        expand(MDIR + "{sample}/align/{alnr}/htd/gba_gauchian/{sample}.{alnr}.gba_gauchian.manifest", sample=SSAMPS, alnr=ALIGNERS)
    output:
        "./logs/gba.done"
    shell:
        """
        touch {output};
        """				
												      

localrules: produce_all_htd

rule produce_all_htd:
    input:
        "./logs/gba.done",
        "./logs/smn12.done",
        "./logs/cyp2d6.done"
    output:
        "./logs/all_htd_cmp_readsy"
    shell:
        "touch {output}"
