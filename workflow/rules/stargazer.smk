# ###### Stargazer -- PGX Star Allele Caller
# ------------------------------------------------------------
# tool: Stargazer
# what: Stargazer can call star alleles in 51 PGx genes (but accuracy varies!)
# docs: https://stargazer.gs.washington.edu/stargazerweb/
# paper: https://ascpt.onlinelibrary.wiley.com/doi/10.1002/cpt.1552
# github: self-hosted
# note: not clearly still supported.


rule stargazer:
    input:
        i=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
    output:
        directory(
            MDIR
            + "{sample}/align/{alnr}/xv/starg/Star_{sample}_{alnr}_.stargazer-pipeline-sges.project/"
        ),
        directory(MDIR + "{sample}/align/{alnr}/xv/starg/"),
        done=MDIR + "{sample}/align/{alnr}/xv/starg/sg.done",
    params:
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        zdbsnp=config["supporting_files"]["files"]["dbsnp"]["supersonic"]["name"],
        l="{",
        r="}",
        cluster_sample=f"wildcards.sample",
    threads: config["stargazer"]["threads"]
    log:
        MDIR + "{sample}/align/{alnr}/xv/starg/logs/{sample}.{alnr}.starg.log",
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.starg.bench.tsv"
    conda:
        config["stargazer"]["env_yaml"]
    shell:
        "mkdir -p {output[1]};"
        "python resources/Stargazer_v1.0.8/stargazer.py pipeline sges -o {output[1]}/Star_{wildcards.sample}_{wildcards.alnr}_ --assembly {params.huref}  --dbsnp {params.zdbsnp}  --gatk $PWD/resources/gatk3/GenomeAnalysisTK.jar --bam {input.i} -c egfr >> {log};"
        "blah;"
        "find {output[0]} | grep .sh  | parallel -j  15 'echo {params.l}{params.r}; bash {params.l}{params.r}|| echo ATTEMPTFAILED!'"
        "touch {output.done};"
