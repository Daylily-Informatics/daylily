

# ####### Extract FastQs from BAMs

rule bam_to_fastqgz:
    input:
        bam=get_bams
    output:
        fq="some.bam"
        bai="some.bam.bai"
    conda:
        "../envs/sentieon_v0.1.yaml"
    params:
        cluster_sample=ret_sample,
    threads: 4 if 'b2fq' not in config else config['b2fq']['threads']
    benchmark:
        "some.benchmark.tsv"
    shell:
        """samtools collate -n -@ {threads} -uO INPUT_BAM tmp- | samtools fastq -@ 32 \
-s >(gzip -c > single.fastq.gz) -0 >(gzip -c > unpaired.fastq.gz) \
-1 >(gzip -c > output_1.fastq.gz) -2 >(gzip -c > output_2.fastq.gz) -"""
