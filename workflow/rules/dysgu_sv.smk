import os

#### DYSGU - new entrant to the SV calling space ~2021
# -------------------------------------
#
# github: https://github.com/kcleal/dysgu
# paper: https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkac039/6517943#329788134
#

rule dysgu:
    input:
        bam=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
        bai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam.bai",
    output:
        vcf = MDIR + "{sample}/align/{alnr}/sv/dysgu/{sample}.{alnr}.dysgu.sv.vcf",
    params:
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        min_sv_size="30" if "dysgu" not in config else config["dysgu"]["min_sv_size"],
        cluster_sample=ret_sample_alnr,
        ld_p=config['malloc_alt']['ld_preload'] if 'ld_preload' not in config['dysgu'] else config['dysgu']['ld_preload'],
    threads: config["dysgu"]["threads"]
    resources:
        threads=config['dysgu']['threads'],
        partition=config['dysgu']['partition'],
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.dysgu.sv.vcf.bench.tsv"
    log:
        MDIR + "{sample}/align/{alnr}/sv/dysgu/logs/{sample}.{alnr}.dysgu.sv.vcf.log",
    conda:
        "../envs/dysgu_sv_v0.2.yaml"
    shell:
        """
        work_dir=$(dirname {output.vcf} );
        (rm -rf $work_dir ) || echo noFail;
        mkdir -p $(dirname {log} ) ;
        
        ( dysgu run  \
         -x -a \
        --metrics  \
        --keep-small \
        --diploid True \
        --pfix dysgu_tmp \
        -o {output.vcf} \
        -f vcf --mode pe \
        --min-support 3 \
        --min-size 30 \
        --max-cov 200  \
        {params.huref} \
        $work_dir \
        {input.bam} ) > {log} 2>&1;
        {latency_wait};
        ls {output};
        rm -rf $(dirname {output.vcf})/*bin || echo 'dysgu bins rmFailed' >> {log} 2>&1;
        rm $(dirname {output.vcf})/*_tmp.bam || echo 'rm bam failed' >> {log} 2>&1;
        """


rule dysgu_sort_index:
    input:
        MDIR + "{sample}/align/{alnr}/sv/dysgu/{sample}.{alnr}.dysgu.sv.vcf"
    output:
        sortvcf = touch(MDIR + "{sample}/align/{alnr}/sv/dysgu/{sample}.{alnr}.dysgu.sv.sort.vcf"),
        sortgz = touch(MDIR + "{sample}/align/{alnr}/sv/dysgu/{sample}.{alnr}.dysgu.sv.sort.vcf.gz"),
        sorttbi = touch(MDIR + "{sample}/align/{alnr}/sv/dysgu/{sample}.{alnr}.dysgu.sv.sort.vcf.gz.tbi"),
    priority: 35
    threads: config["dysgu_sort_index"]["threads"]
    resources:
        threads=config['dysgu_sort_index']['threads'],
        partition=config['dysgu_sort_index']['partition']
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.dysgu.sv.vcf.sort.bench.tsv"
    log:
        MDIR + "{sample}/align/{alnr}/sv/dysgu/logs/{sample}.{alnr}.dysgu.sv.vcf.sort.log",
    conda:
        "../envs/vanilla_v0.1.yaml"
    params:
        cluster_sample=ret_sample,
    shell:
        """
        (bedtools sort -header -i {input} > {output.sortvcf};
        bgzip -f -@ {threads} {output.sortvcf};
        touch {output.sortvcf};
        tabix -p vcf -f {output.sortgz};
        {latency_wait} || echo passOn;
        ls {output} || echo passOn ;
        ) > {log} 2>&1;
        {latency_wait};
        ls {output};
        """



localrules: produce_dysgu,

rule produce_dysgu:  # TARGET: Produce All Dysgu
    priority: 39
    input:
        expand(MDIR +"{sample}/align/{alnr}/sv/dysgu/{sample}.{alnr}.dysgu.sv.sort.vcf.gz.tbi", sample=SSAMPS, alnr=ALIGNERS)
