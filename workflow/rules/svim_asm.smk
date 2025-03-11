##### SVIM-ASM : Assembly based SV calling of discordant reads
# -------------------------------------
#
# github: https://github.com/eldariont/svim-asm
# paper: https://academic.oup.com/bioinformatics/article/36/22-23/5519/6042701?login=false


rule svim_asm:
    input:
        MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
    output:
        stub = touch(MDIR + "{sample}/align/{alnr}/sv/svasm/{sample}.{alnr}.svasm.sv"),
        vcff = touch(MDIR + "{sample}/align/{alnr}/sv/svasm/{sample}.{alnr}.svasm.sv.vcf"),
    params:
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        cluster_sample=ret_sample_alnr,
    threads: config["svim_asm"]["threads"]
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.svasm.sv.vcf.bench.tsv"
    log:
        MDIR + "{sample}/align/{alnr}/sv/svasm/logs/{sample}.{alnr}.svasm.sv.vcf.log",
    conda:
        config["svim_asm"]["env_yaml"]
    shell:
        """
        set +euo pipefail;

        ls {output} || echo passOn;
        """


rule svim_asm_sort_index:
    input:
        MDIR + "{sample}/align/{alnr}/sv/svasm/{sample}.{alnr}.svasm.sv.vcf"
    output:
        sortgz = touch(MDIR + "{sample}/align/{alnr}/sv/svasm/{sample}.{alnr}.svasm.sv.sort.vcf.gz"),
        sorttbi = touch(MDIR + "{sample}/align/{alnr}/sv/svasm/{sample}.{alnr}.svasm.sv.sort.vcf.gz.tbi"),
        tmpd=directory(MDIR+"{sample}/align/{alnr}/sv/svasm/tmp")
    threads: config["svim_asm"]["threads"]
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.svasm.sv.vcf.sort.bench.tsv"
    log:
        MDIR + "{sample}/align/{alnr}/sv/svasm/logs/{sample}.{alnr}.svasm.sv.vcf.sort.log",
    conda:
        config["svim_asm"]["env_yaml"]
    params:
        cluster_sample=ret_sample,
    shell:
        """
        set +euo pipefail;                                                                                                                      rm -rf {output.tmpd} || echo noTmpDir;
        mkdir -p {output.tmpd} || echo passOn;
        perl -pi -e 's/\\t\\*\\t/\\tN\\t/g;' {input} >> {log} || echo nomatch! ;
        bcftools sort -T {output.tmpd} -O z -o {output.sortgz} {input} || echo passOn;
        bcftools index  -f -t --threads {threads} {output.sortgz} -o {output.sorttbi} || echo passOn ;
        {latency_wait} || echo passOn;
        ls {output} || echo passOn ;
        """
