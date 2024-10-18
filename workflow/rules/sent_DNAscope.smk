import sys
import os

##### sentieon DNAscope - OUR snv CALLER
# ---------------------------
#


def get_dchrm_day(wildcards):
    pchr=""
    if config['genome_build'] not in ['b37']:
        pchr="chr"


    ret_str = ""
    sl = wildcards.dchrm.replace('chr','').split("-")
    sl2 = wildcards.dchrm.replace('chr','').split("~")
    
    if len(sl2) == 2:
        ret_str = pchr + wildcards.dchrm + ':'
    elif len(sl) == 1:
        ret_str = pchr + sl[0] + ':'
    elif len(sl) == 2:
        start = int(sl[0])
        end = int(sl[1])
        while start <= end:
            ret_str = str(ret_str) + "," + pchr + str(start) + ':'
            start = start + 1
    else:
        raise Exception(
            "sentD chunks can only be one contiguous range per chunk : ie: 1-4 with the non numerical chrms assigned 23=X, 24=Y, 25=MT"
        )

    return ret_mod_chrm(ret_str).lstrip(',').replace('chr23','chrX').replace('chr24','chrY').replace('chr25','chrMT').replace('23:','X:').replace('24:','Y:').replace('25:','MT:')


rule sent_DNAscope:
    input:
        b=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
        bai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam.bai",
        d=MDIR + "{sample}/align/{alnr}/snv/sentd/vcfs/{dchrm}/{sample}.ready",
    output:
     vcf=MDIR
        + "{sample}/align/{alnr}/snv/sentd/vcfs/{dchrm}/{sample}.{alnr}.sentd.{dchrm}.snv.vcf",
     tvcf=temp(MDIR
        + "{sample}/align/{alnr}/snv/sentd/vcfs/{dchrm}/{sample}.{alnr}.sentd.{dchrm}.snv.vcf.tmp"),
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/sentd/log/vcfs/{sample}.{alnr}.sentd.{dchrm}.snv.log",
    threads: config['sentD']['threads']
    conda:
        "../envs/sentD_v0.1.yaml"
    priority: 45
    benchmark:
        repeat(
            MDIR + "{sample}/benchmarks/{sample}.{alnr}.sentd.{dchrm}.bench.tsv",
            0
            if "bench_repeat" not in config["sentD"]
            else config["sentD"]["bench_repeat"],
        )
    resources:
        attempt_n=lambda wildcards, attempt:  (attempt + 0),
        partition=config['sentD']['partition'],
        threads=config['sentD']['threads'],
        vcpu=config['sentD']['threads'],
    params:
        schrm_mod=get_dchrm_day,
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        model=config["supporting_files"]["files"]["sentieon"]["dnascope_model"]["name"],
        cluster_sample=ret_sample,
        numactl=config["sentieon"]["numactl"],
    shell:
        """
        export SENTIEON_TMPDIR='/fsx/scratch' &&   export SENTIEON_LICENSE='/fsx/SAVEME_ANA/etc/Daylily_Informatics_eval.lic'  && {params.numactl} sentieon driver --thread_count {threads} --interval {params.schrm_mod} --reference {params.huref} --input {input.b} --algo DNAscope --pcr_indel_model none --model {params.model}  {output.tvcf};
       {params.numactl} sentieon driver -t {threads} -r {params.huref} --algo DNAModelApply --model {params.model} -v {output.tvcf} {output.vcf};

        touch {output.vcf};
        """


rule sentD_sort_index_chunk_vcf:
    input:
        vcf=MDIR
        + "{sample}/align/{alnr}/snv/sentd/vcfs/{dchrm}/{sample}.{alnr}.sentd.{dchrm}.snv.vcf",
    priority: 46
    output:
        vcfsort=MDIR
        + "{sample}/align/{alnr}/snv/sentd/vcfs/{dchrm}/{sample}.{alnr}.sentd.{dchrm}.snv.sort.vcf",
        vcfgz=MDIR
        + "{sample}/align/{alnr}/snv/sentd/vcfs/{dchrm}/{sample}.{alnr}.sentd.{dchrm}.snv.sort.vcf.gz",
        vcftbi=MDIR
        + "{sample}/align/{alnr}/snv/sentd/vcfs/{dchrm}/{sample}.{alnr}.sentd.{dchrm}.snv.sort.vcf.gz.tbi",
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/sentd/vcfs/{dchrm}/log/{sample}.{alnr}.sentd.{dchrm}.snv.sort.vcf.gz.log",
    resources:
        vcpu=8,
        threads=8,
        partition="i8,i8,i64,i96,i128"
    params:
        x='y',
        cluster_sample=ret_sample,
    threads: 8 #config["config"]["sort_index_sentDna_chunk_vcf"]['threads']
    shell:
        """
        (rm {output} 1>  /dev/null  2> /dev/null )  || echo rmfailed > {log};
        (bedtools sort -header -i {input.vcf} > {output.vcfsort} 2>> {log}) || exit 1233;
        (
        bgzip {output.vcfsort};
        
        touch {output.vcfsort};
        tabix -f -p vcf {output.vcfgz};
        touch {output.vcftbi};
        {latency_wait};
        ls {output}; ) > {log} 2>&1 ;
        
        {latency_wait};
        """


localrules:
    sentD_concat_fofn,


rule sentD_concat_fofn:
    input:
        chunk_tbi=sorted(
            expand(
                MDIR
                + "{{sample}}/align/{{alnr}}/snv/sentd/vcfs/{ochm}/{{sample}}.{{alnr}}.sentd.{ochm}.snv.sort.vcf.gz.tbi",
                ochm=SENTD_CHRMS,            ),            key=lambda x: float(                str(x.replace("~", ".").replace(":", "."))               .split("vcfs/")[1]                .split("/")[0]                .split("-")[0]            ),        ),
    # This expand pattern is neat.  the escaped {} remain acting as a snakemake wildcard and expect to be derived from the dag, while th dchrm wildcard is effectively being constrained by the values in the SENTD_CHRMS array;  So you produce 1 input array of files for every sample+dchrm parir, with one list string/file name per array.  The rule will only begin when all array members are produced. It's then sorted by first sentdchrm so they can be concatenated w/out another soort as all the chunks had been sorted already.
    priority: 44
    output:
        fin_fofn=MDIR
        + "{sample}/align/{alnr}/snv/sentd/{sample}.{alnr}.sentd.snv.concat.vcf.gz.fofn",
        tmp_fofn=MDIR        + "{sample}/align/{alnr}/snv/sentd/{sample}.{alnr}.sentd.snv.concat.vcf.gz.fofn.tmp",
    threads: 2
    resources:
        threads=2
    params:
        fn_stub="{sample}.{alnr}.sentd."
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.sentd.concat.fofn.bench.tsv"
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR + "{sample}/align/{alnr}/snv/sentd/log/{sample}.{alnr}.sentd.cocncat.fofn.log",
    shell:
        """
        (rm {output} 1> /dev/null  2> /dev/null ) || echo rmFailOK >> {log} && ls ./ >> {log};
        ### export LD_LIBRARY_PATH=$PWD/resources/libs/;
        for i in {input.chunk_tbi}; do
            ii=$(echo $i | perl -pe 's/\.tbi$//g'; );
            echo $ii >> {output.tmp_fofn};
        done;
        (workflow/scripts/sort_concat_chrm_list.py {output.tmp_fofn} {wildcards.sample}.{wildcards.alnr}.sentd. {output.fin_fofn}) || echo "Python Script Error? CODE __ $? __" >> {log} && ls -lt {output} >> {log};

        """


rule sentD_concat_index_chunks:
    input:
        fofn=MDIR
        + "{sample}/align/{alnr}/snv/sentd/{sample}.{alnr}.sentd.snv.concat.vcf.gz.fofn",
    output:
        vcf=touch(
            MDIR + "{sample}/align/{alnr}/snv/sentd/{sample}.{alnr}.sentd.snv.sort.vcf"
        ),
        vcfgz=touch(
            MDIR + "{sample}/align/{alnr}/snv/sentd/{sample}.{alnr}.sentd.snv.sort.vcf.gz"
        ),
        vcfgztbi=touch(
            MDIR
            + "{sample}/align/{alnr}/snv/sentd/{sample}.{alnr}.sentd.snv.sort.vcf.gz.tbi"
        ),
    threads: 8
    resources:
        vcpu=8,
        threads=8,
        partition="i8,i8,i64,i96,i128"
    priority: 47
    params:
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        cluster_sample=ret_sample,
    resources:
        attempt_n=lambda wildcards, attempt:  (attempt + 0)
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.sentd.merge.bench.tsv"
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/sentd/log/{sample}.{alnr}.sentd.snv.merge.sort.gatherered.log",
    shell:
        """
        (rm {output} 1> /dev/null  2> /dev/null ) || echo rmFAIL;
        mkdir -p $(dirname {log});
        ### export LD_LIBRARY_PATH=$PWD/resources/libs/;
        bcftools concat -a -d all --threads {threads} -f {input.fofn}  -O v -o {output.vcf};
        bcftools view -O z -o {output.vcfgz} {output.vcf};
        bcftools index -f -t --threads {threads} -o {output.vcfgztbi} {output.vcfgz};
        stats_f=$(echo "{output.vcfgz}.bcf.stats");
        bcftools stats -F {params.huref}  {output.vcfgz} > $stats_f;
        {latency_wait}; > {log} """


localrules:
    clear_combined_sentD_vcf,


rule clear_combined_sentD_vcf:  # TARGET:  clear combined sentD vcf so the chunks can be re-evaluated if needed.
    input:
        expand(
            MDIR + "{sample}/align/{alnr}/snv/sentd/{sample}.{alnr}.sentd.snv.sort.vcf.gz",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
    priority: 42
    shell:
        "(rm {input}*  1> /dev/null  2> /dev/null ) || echo 'file not found for deletion: {input}';"


localrules:
    produce_sentD_vcf,


rule produce_sentD_vcf:  # TARGET: just gen sentD calls
    input:
        expand(
            MDIR
            + "{sample}/align/{alnr}/snv/sentd/{sample}.{alnr}.sentd.snv.sort.vcf.gz.tbi",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
    output:
        "gatheredall.sentd",
    priority: 48
    log:
        "gatheredall.sentd.log",
    shell:
        """( touch {output} ;

        {latency_wait}; ls {output} ) >> {log} 2>&1;
        """


localrules:
    prep_sentD_chunkdirs,


rule prep_sentD_chunkdirs:
    input:
        b=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.sort.bam",
    output:
        expand(
            MDIR + "{{sample}}/align/{{alnr}}/snv/sentd/vcfs/{dchrm}/{{sample}}.ready",
            dchrm=SENTD_CHRMS,
        ),
    log:
        MDIR + "{sample}/align/{alnr}/snv/sentd/logs/{sample}.{alnr}.chunkdirs.log",
    shell:
        """
        ( echo {output}  ;
        mkdir -p $(dirname {output} );
        touch {output};
        ls {output}; ) > {log} 2>&1;
        """
