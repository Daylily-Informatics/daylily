##### LoFreq2
# ---------------------------
#
#
# Github: https://csb5.github.io/lofreq/installation/
#
# https://pubmed.ncbi.nlm.nih.gov/23066108/
#

import sys
import os


def get_lofreq_chrm(wildcards):
    pchr = ""
    ret_str = ""
    sl = wildcards.lfqchrm.replace('chr','').split("-")
    sl2 = wildcards.lfqchrm.replace('chr','').split("~")
    
    if len(sl2) == 2:
        ret_str = pchr + wildcards.lfqchrm
    elif len(sl) == 1:
        ret_str = pchr + sl[0]
    elif len(sl) == 2:
        start = int(sl[0])
        end = int(sl[1])
        while start <= end:
            ret_str = str(ret_str) + " " + pchr + str(start)
            start = start + 1
    else:
        raise Exception(
            "LoFreq chunks can only be one contiguous range per chunk: e.g., 1-4 with the non-numerical chromosomes assigned 23=X, 24=Y, 25=MT"
        )

    return ret_mod_chrm(ret_str)


rule lfq2_indelqual:
    input:
        bam=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
        bai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam.bai",
    output:
        bam=MDIR + "{sample}/align/{alnr}/snv/lfq2/{sample}.{alnr}.mrkdup.sort.indelqual.bam",
        bai=MDIR + "{sample}/align/{alnr}/snv/lfq2/{sample}.{alnr}.mrkdup.sort.indelqual.bam.bai",
    log:
        MDIR + "{sample}/align/{alnr}/snv/lfq2/log/{sample}.{alnr}.indelqual.log",
    conda:
        "../envs/lofreq2_v0.1.yaml"
    threads: config['lofreq2']['threads']
    resources:
        vcpu=config['lofreq2']['threads'],
        threads=config['lofreq2']['threads'],
        partition=config['lofreq2']['partition'],
        mem_mb=config['lofreq2']['mem_mb'],
    params:
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
    	cluster_sample=ret_sample, 
    shell:
        """
        touch {log};
        #lofreq indelqual --dindel -f {params.huref} -o {output.bam} {input.bam} >> {log} 2>&1;
        #samtools index {output.bam} >> {log} 2>&1;
        touch {output};
        """


rule lofreq2:
    input:
        bam=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
        bai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam.bai",
        ibam=MDIR + "{sample}/align/{alnr}/snv/lfq2/{sample}.{alnr}.mrkdup.sort.indelqual.bam",
        ibai=MDIR + "{sample}/align/{alnr}/snv/lfq2/{sample}.{alnr}.mrkdup.sort.indelqual.bam.bai",
        d=MDIR + "{sample}/align/{alnr}/snv/lfq2/vcfs/{lfqchrm}/{sample}.ready",
    output:
        vcf=temp(MDIR + "{sample}/align/{alnr}/snv/lfq2/vcfs/{lfqchrm}/{sample}.{alnr}.lfq2.{lfqchrm}.snv.vcf"),
    log:
        MDIR + "{sample}/align/{alnr}/snv/lfq2/log/{sample}.{alnr}.lfq2.{lfqchrm}.snv.log",
    threads: config['lofreq2']['threads']
    conda:
        "../envs/lofreq2_v0.1.yaml"
    priority: 45
    resources:
        vcpu=config['lofreq2']['threads'],
        threads=config['lofreq2']['threads'],
        partition=config['lofreq2']['partition'],
        mem_mb=config['lofreq2']['mem_mb'],
    benchmark:
        repeat(
            MDIR + "{sample}/benchmarks/{sample}.{alnr}.lfq2.{lfqchrm}.bench.tsv",
            0
            if "bench_repeat" not in config["lofreq2"]
            else config["lofreq2"]["bench_repeat"],
        )
    params:
        dchrm=get_lofreq_chrm,
        cluster_sample=ret_sample,
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        mdir=MDIR,
        mem_mb=config['lofreq2']['mem_mb'],
        cpre="" if "b37" == config['genome_build'] else "chr",
    shell:
        """
        touch {log};
        start_time=$(date +%s);
        dchr=$(echo {params.cpre}{params.dchrm} | sed 's/~/\:/g' | sed 's/23\:/X\:/' | sed 's/24\:/Y\:/' | sed 's/25\:/MT\:/');

        echo "DCHRM: $dchr" >> {log} 2>&1;
        
        if [[ "{params.dchrm}" == "1-24" || "{params.dchrm}" == "1-25" ]]; then
            echo "lofreq parallel" >> {log} 2>&1;
            lofreq call-parallel --pp-threads {threads}  --max-depth 10000 \
            --force-overwrite \
            -f {params.huref} \
            -o {output.vcf} {input.bam} >> {log} 2>&1;
        else
            echo "lofreq single thread" >> {log} 2>&1;
            lofreq call \
            --force-overwrite \
            -f {params.huref} \
            -r $dchr \
            -o {output.vcf} {input.bam} >> {log} 2>&1;
        fi;

        end_time=$(date +%s);
        elapsed_time=$((($end_time - $start_time) / 60));

        echo "Elapsed-Time-min:\t$elapsed_time" >> {log};
        """


rule lofreq2_sort_index_chunk_vcf:
    input:
        vcf=MDIR + "{sample}/align/{alnr}/snv/lfq2/vcfs/{lfqchrm}/{sample}.{alnr}.lfq2.{lfqchrm}.snv.vcf",
    priority: 46
    output:
        tmpvcf=temp(MDIR + "{sample}/align/{alnr}/snv/lfq2/vcfs/{lfqchrm}/{sample}.{alnr}.lfq2.{lfqchrm}.snv.tmp.vcf"),
        vcfsort=MDIR + "{sample}/align/{alnr}/snv/lfq2/vcfs/{lfqchrm}/{sample}.{alnr}.lfq2.{lfqchrm}.snv.sort.vcf",
        vcfgz=MDIR + "{sample}/align/{alnr}/snv/lfq2/vcfs/{lfqchrm}/{sample}.{alnr}.lfq2.{lfqchrm}.snv.sort.vcf.gz",
        vcftbi=MDIR + "{sample}/align/{alnr}/snv/lfq2/vcfs/{lfqchrm}/{sample}.{alnr}.lfq2.{lfqchrm}.snv.sort.vcf.gz.tbi",
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR + "{sample}/align/{alnr}/snv/lfq2/vcfs/{lfqchrm}/log/{sample}.{alnr}.lfq2.{lfqchrm}.snv.sort.vcf.gz.log",
    resources:
        vcpu=4,
        threads=config['lofreq2']['threads'],
        partition=config['lofreq2']['partition'],
    params:
        cluster_sample=ret_sample,
    threads: config['lofreq2']['threads']
    shell:
        """

        bash bin/repair_lofreq2_vcf.sh {input.vcf}  {output.tmpvcf}  $(dirname {input.vcf})/_fixvcf {params.cluster_sample} >> {log} 2>&1;
        bedtools sort -header -i {output.tmpvcf} > {output.vcfsort} 2>> {log};
        
        bgzip {output.vcfsort} >> {log} 2>&1;        
        tabix -f -p vcf {output.vcfgz} >> {log} 2>&1;

        """


rule lofreq2_concat_fofn:
    input:
        chunk_tbi=sorted(
            expand(
                MDIR + "{{sample}}/align/{{alnr}}/snv/lfq2/vcfs/{dvchm}/{{sample}}.{{alnr}}.lfq2.{dvchm}.snv.sort.vcf.gz.tbi",
                dvchm=LOFREQ_CHRMS,
            ),
            key=lambda x: float(
                str(x.replace("~", ".").replace(":", "."))
                .split("vcfs/")[1]
                .split("/")[0]
                .split("-")[0]
            ),
        ),
    output:
        fin_fofn=MDIR + "{sample}/align/{alnr}/snv/lfq2/{sample}.{alnr}.lfq2.snv.concat.vcf.gz.fofn",
        tmp_fofn=MDIR + "{sample}/align/{alnr}/snv/lfq2/{sample}.{alnr}.lfq2.snv.concat.vcf.gz.fofn.tmp",
    threads: 2
    resources:
        vcpu=2,
        threads=2,
        partition=config['lofreq2']['partition'],
    params:
        fn_stub="{sample}.{alnr}.lfq2.",
        cluster_sample=ret_sample,
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.lfq2.concat.fofn.bench.tsv"
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR + "{sample}/align/{alnr}/snv/lfq2/log/{sample}.{alnr}.lfq2.concat.fofn.log",
    shell:
        """

        for i in {input.chunk_tbi}; do
            ii=$(echo $i | perl -pe 's/\.tbi$//g');
            echo $ii >> {output.tmp_fofn};
        done;
        (workflow/scripts/sort_concat_chrm_list.py {output.tmp_fofn} {wildcards.sample}.{wildcards.alnr}.lfq2. {output.fin_fofn}) >> {log} 2>&1;
        """


rule lofreq2_concat_index_chunks:
    input:
        fofn=MDIR + "{sample}/align/{alnr}/snv/lfq2/{sample}.{alnr}.lfq2.snv.concat.vcf.gz.fofn",
        tmp_fofn=MDIR + "{sample}/align/{alnr}/snv/lfq2/{sample}.{alnr}.lfq2.snv.concat.vcf.gz.fofn.tmp",
    output:
        vcfgz=touch(MDIR + "{sample}/align/{alnr}/snv/lfq2/{sample}.{alnr}.lfq2.snv.sort.vcf.gz"),
        vcfgztbi=touch(MDIR + "{sample}/align/{alnr}/snv/lfq2/{sample}.{alnr}.lfq2.snv.sort.vcf.gz.tbi"),
    threads: 4
    resources:
        vcpu=4,
        threads=4,
        partition=config['lofreq2']['partition'],
    priority: 47
    params:
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        cluster_sample=ret_sample,
    resources:
        attempt_n=lambda wildcards, attempt: (attempt + 0)
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.lfq2.merge.bench.tsv"
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR + "{sample}/align/{alnr}/snv/lfq2/log/{sample}.{alnr}.lfq2.snv.merge.sort.gatherered.log",
    shell:
        """

        touch {log};
        mkdir -p $(dirname {log});

        bcftools concat -n --threads {threads} -f {input.fofn}  -O z -o {output.vcfgz};
        bcftools index -f -t --threads {threads} -o {output.vcfgztbi} {output.vcfgz};

        rm -rf $(dirname {output.vcfgz})/vcfs >> {log} 2>&1;
        
        """


rule produce_lofreq2_vcf:  # TARGET: lofreq2 vcfs
    input:
        vcftb=expand(
            MDIR + "{sample}/align/{alnr}/snv/lfq2/{sample}.{alnr}.lfq2.snv.sort.vcf.gz",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
        vcftbi=expand(
            MDIR + "{sample}/align/{alnr}/snv/lfq2/{sample}.{alnr}.lfq2.snv.sort.vcf.gz.tbi",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
    output:
        "gatheredall.lfq2",
    threads: 4
    priority: 48
    log:
        "gatheredall.lfq2.log",
    conda:
        "../envs/vanilla_v0.1.yaml"
    shell:
        """
        for vcf in {input.vcftb}; do
            bcf="${{vcf%.vcf.gz}}.bcf";
            bcftools view -O b -o $bcf --threads {threads} $vcf && bcftools index --threads 4 $bcf;
        done;

        touch {output};

        {latency_wait};
        ls {output} >> {log} 2>&1;

        {latency_wait}; 
        ls {output}  >> {log} 2>&1;
        """


rule prep_lofreq2_chunkdirs:
    input:
        b=MDIR + "{sample}/align/{alnr}/snv/lfq2/{sample}.{alnr}.mrkdup.sort.indelqual.bam",
    output:
        expand(
            MDIR + "{{sample}}/align/{{alnr}}/snv/lfq2/vcfs/{lfqchrm}/{{sample}}.ready",
            lfqchrm=LOFREQ_CHRMS,
        ),
    threads: 1
    log:
        MDIR + "{sample}/align/{alnr}/snv/lfq2/log/{sample}.{alnr}.chunkdirs.log",
    shell:
        """
        (echo {output};
        mkdir -p $(dirname {output});
        touch {output};
        ls {output};) > {log} 2>&1;
        """

localrules:
    prep_lofreq2_chunkdirs,
