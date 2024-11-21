
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


def get_lochrm_mod(wildcards):
    pchr=""
    if config['genome_build'] not in ['b37']:
        pchr="chr"
    ret_str = ""
    sl = wildcards.ochrm.split("-")
    sl2 = wildcards.ochrm.split("~")


    #from IPython import embed
    #embed()
    #raise

    if len(sl2) == 2:
        ret_str = pchr + wildcards.ochrm
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
            "oct chunks can only be one contiguous range per chunk : ie: 1-4 with the non numerical chrms assigned 23=X, 24=Y, 25=MT"
        )
    ret_str = ret_mod_chrm(ret_str)

    return ret_str


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
        d=MDIR + "{sample}/align/{alnr}/snv/lfq2/vcfs/{dvchrm}/{sample}.ready",
    output:
        vcf=temp(MDIR + "{sample}/align/{alnr}/snv/lfq2/vcfs/{dvchrm}/{sample}.{alnr}.lfq2.{dvchrm}.snv.vcf"),
    log:
        MDIR + "{sample}/align/{alnr}/snv/lfq2/log/{sample}.{alnr}.lfq2.{dvchrm}.snv.log",
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
            MDIR + "{sample}/benchmarks/{sample}.{alnr}.lfq2.{dvchrm}.bench.tsv",
            0
            if "bench_repeat" not in config["lofreq2"]
            else config["lofreq2"]["bench_repeat"],
        )
    params:
        dchrm=get_lochrm_mod,
        cluster_sample=ret_sample,
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        mdir=MDIR,
        mem_mb=config['lofreq2']['mem_mb'],
        cpre="" if "b37" == config['genome_build'] else "chr",
    shell:
        """
        touch {log};
        start_time=$(date +%s);
        dchr={params.cpre}{params.dchrm};

        if [[ "{params.dchrm}" == "23" ]]; then
            dchr='{params.cpre}X';
        elif [[ "{params.dchrm}" == "24" ]]; then
            dchr='{params.cpre}Y';
        elif [[ "{params.dchrm}" == "25" ]]; then
            dchr='{params.cpre}MT';
        else
            dchr={params.cpre}{params.dchrm};
        fi;

        echo 'DCHRM: $dchr' >> {log} 2>&1;

        export lochrm_mod=$(echo '{params.dchrm}' | sed 's/~/\:/g' | perl -pe 's/(^23| 23)/ X/g;' | perl -pe 's/(^24| 24)/ Y/g;' | perl -pe 's/(^25| 25)/ MT/g;');
        echo 'LOCHRM: $lochrm_mod' >> {log} 2>&1;

        lofreq call --max-depth 10000 \
        --force-overwrite \
        -f {params.huref} \
        -r $lochrm_mod \
        -o {output.vcf} {input.bam} >> {log} 2>&1;

        sleep 100000;
        end_time=$(date +%s);
        elapsed_time=$((($end_time - $start_time) / 60));

        echo "Elapsed-Time-min:\t$elapsed_time" >> {log};
        """


rule lofreq2_sort_index_chunk_vcf:
    input:
        vcf=MDIR + "{sample}/align/{alnr}/snv/lfq2/vcfs/{dvchrm}/{sample}.{alnr}.lfq2.{dvchrm}.snv.vcf",
    priority: 46
    output:
        tmpvcf=temp(MDIR + "{sample}/align/{alnr}/snv/lfq2/vcfs/{dvchrm}/{sample}.{alnr}.lfq2.{dvchrm}.snv.tmp.vcf"),
        vcfsort=MDIR + "{sample}/align/{alnr}/snv/lfq2/vcfs/{dvchrm}/{sample}.{alnr}.lfq2.{dvchrm}.snv.sort.vcf",
        vcfgz=MDIR + "{sample}/align/{alnr}/snv/lfq2/vcfs/{dvchrm}/{sample}.{alnr}.lfq2.{dvchrm}.snv.sort.vcf.gz",
        vcftbi=MDIR + "{sample}/align/{alnr}/snv/lfq2/vcfs/{dvchrm}/{sample}.{alnr}.lfq2.{dvchrm}.snv.sort.vcf.gz.tbi",
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR + "{sample}/align/{alnr}/snv/lfq2/vcfs/{dvchrm}/log/{sample}.{alnr}.lfq2.{dvchrm}.snv.sort.vcf.gz.log",
    resources:
        vcpu=4,
        threads=4,
        partition=config['lofreq2']['partition_other'],
    params:
        cluster_sample=ret_sample,
    threads: 4
    shell:
        """

        bash bin/repair_lofreq2_vcf.sh {input.vcf}  {output.tmpvcf}  $(dirname {input.vcf})/_fixvcf {params.cluster_sample} >> {log} 2>&1;
        (bedtools sort -header -i {output.tmpvcf} > {output.vcfsort} 2>> {log}) || exit 1233;
        (
        bgzip {output.vcfsort};        
        touch {output.vcfsort};
        tabix -f -p vcf {output.vcfgz};
        touch {output.vcftbi};
        ) >> {log} 2>&1;

        ls {output} >> {log} 2>&1;
        
        {latency_wait};
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
        partition="i192",
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
        (rm {output} 1> /dev/null 2> /dev/null) || echo rmFailOK >> {log} && ls ./ >> {log};

        for i in {input.chunk_tbi}; do
            ii=$(echo $i | perl -pe 's/\.tbi$//g');
            echo $ii >> {output.tmp_fofn};
        done;
        (workflow/scripts/sort_concat_chrm_list.py {output.tmp_fofn} {wildcards.sample}.{wildcards.alnr}.lfq2. {output.fin_fofn}) || echo "Python Script Error? CODE __ $? __" >> {log} && ls -lt {output.fin_fofn} {output.tmp_fofn} >> {log};
        """


rule lofreq2_concat_index_chunks:
    input:
        fofn=MDIR + "{sample}/align/{alnr}/snv/lfq2/{sample}.{alnr}.lfq2.snv.concat.vcf.gz.fofn",
        tmp_fofn=MDIR + "{sample}/align/{alnr}/snv/lfq2/{sample}.{alnr}.lfq2.snv.concat.vcf.gz.fofn.tmp",
    output:
        vcf=temp(MDIR + "{sample}/align/{alnr}/snv/lfq2/{sample}.{alnr}.lfq2.snv.sort.vcf"),
        vcfgz=touch(MDIR + "{sample}/align/{alnr}/snv/lfq2/{sample}.{alnr}.lfq2.snv.sort.vcf.gz"),
        vcfgztbi=touch(MDIR + "{sample}/align/{alnr}/snv/lfq2/{sample}.{alnr}.lfq2.snv.sort.vcf.gz.tbi"),
    threads: 4
    resources:
        vcpu=4,
        threads=4,
        partition=config['lofreq2']['partition_other'],
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
        (rm {output} 1> /dev/null 2> /dev/null) || echo rmFAIL;
        mkdir -p $(dirname {log});

        bcftools concat -a -d all --threads {threads} -f {input.fofn} -O v -o {output.vcf};
        bcftools view -O z -o {output.vcfgz} {output.vcf};
        bcftools index -f -t --threads {threads} -o {output.vcfgztbi} {output.vcfgz};
        stats_f=$(echo "{output.vcfgz}.bcf.stats");
        bcftools stats -F {params.huref} {output.vcfgz} > $stats_f;

        rm -rf $(dirname {output.vcfgz})/vcfs >> {log} 2>&1;
        {latency_wait};
        touch {log};
        """


rule produce_lofreq2_vcf:
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
            MDIR + "{{sample}}/align/{{alnr}}/snv/lfq2/vcfs/{dvchrm}/{{sample}}.ready",
            dvchrm=LOFREQ_CHRMS,
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
