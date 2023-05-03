import sys
import os

##### deepvariant
# ---------------------------
#

def get_dvchrm_day(wildcards):
    pchr=""

    ret_str = ""
    sl = wildcards.dvchrm.replace('chr','').split("-")
    sl2 = wildcards.dvchrm.replace('chr','').split("~")
    
    if len(sl2) == 2:
        ret_str = pchr + wildcards.dvchrm
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
            "deep chunks can only be one contiguous range per chunk : ie: 1-4 with the non numerical chrms assigned 23=X, 24=Y, 25=MT"
        )

    return ret_mod_chrm(ret_str)


rule deepvariant:
    input:
        bam=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
        bai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam.bai",
        d=MDIR + "{sample}/align/{alnr}/snv/deep/vcfs/{dvchrm}/{sample}.ready",
    output:
        vcf=MDIR
        + "{sample}/align/{alnr}/snv/deep/vcfs/{dvchrm}/{sample}.{alnr}.deep.{dvchrm}.snv.vcf",
        gvcf=MDIR
        + "{sample}/align/{alnr}/snv/deep/vcfs/{dvchrm}/{sample}.{alnr}.deep.{dvchrm}.snv.g.vcf",
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/deep/log/vcfs/{sample}.{alnr}.deep.{dvchrm}.snv.log",
    threads: config['deepvariant']['threads']
    container:
        "docker://daylilyinformatics/deepvariant-avx512:1.5.0"
    priority: 45
    resources:
        vcpu=config['deepvariant']['threads'],
        threads=config['deepvariant']['threads'],
        partition=config['deepvariant']['partition'],
    benchmark:
        repeat(
            MDIR + "{sample}/benchmarks/{sample}.{alnr}.deep.{dvchrm}.bench.tsv",
            0
            if "bench_repeat" not in config["deepvariant"]
            else config["deepvariant"]["bench_repeat"],
        )
    params:
        dchrm=get_dvchrm_day,
        cluster_sample=ret_sample, #
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        mdir=MDIR,
        mem_mb=config['deepvariant']['mem_mb'],
        numa="",  #config['deepvariant']['numa'],
        cpre="" if "b37" == config['genome_build'] else "chr",
    shell:
        """
        touch {log};
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

        echo 'DCHRM: $dchr';
        export LD_LIBRARY_PATH=resources/lib/;
        {params.numa} \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS --ref={params.huref} \
        --reads={input.bam} \
        --regions=$dchr \
        --output_vcf={output.vcf} \
        --output_gvcf={output.gvcf} \
        --num_shards={threads} \
        --logging_dir=$(dirname {log}) \
        --dry_run=false;
        """



rule dv_sort_index_chunk_vcf:
    input:
        vcf=MDIR
        + "{sample}/align/{alnr}/snv/deep/vcfs/{dvchrm}/{sample}.{alnr}.deep.{dvchrm}.snv.vcf",
    priority: 46
    output:
        vcfsort=MDIR
        + "{sample}/align/{alnr}/snv/deep/vcfs/{dvchrm}/{sample}.{alnr}.deep.{dvchrm}.snv.sort.vcf",
        vcfgz=MDIR
        + "{sample}/align/{alnr}/snv/deep/vcfs/{dvchrm}/{sample}.{alnr}.deep.{dvchrm}.snv.sort.vcf.gz",
        vcftbi=MDIR
        + "{sample}/align/{alnr}/snv/deep/vcfs/{dvchrm}/{sample}.{alnr}.deep.{dvchrm}.snv.sort.vcf.gz.tbi",
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/deep/vcfs/{dvchrm}/log/{sample}.{alnr}.deep.{dvchrm}.snv.sort.vcf.gz.log",
    resources:
        vcpu=8,
        threads=8,
        partition="i4-5,i16-5,i32-5,i64-5,i96-5,i128-6"
    params:
        cluster_sample=ret_sample,
    threads: 8 #config["config"]["sort_index_deepDna_chunk_vcf"]['threads']
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
    deep_concat_fofn,


rule deep_concat_fofn:
    input:
        chunk_tbi=sorted(
            expand(
                MDIR
                + "{{sample}}/align/{{alnr}}/snv/deep/vcfs/{dvchm}/{{sample}}.{{alnr}}.deep.{dvchm}.snv.sort.vcf.gz.tbi",
                dvchm=DEEPD_CHRMS,            ),            key=lambda x: float(                str(x.replace("~", ".").replace(":", "."))               .split("vcfs/")[1]                .split("/")[0]                .split("-")[0]            ),        ),
    # This expand pattern is neat.  the escaped {} remain acting as a snakemake wildcard and expect to be derived from the dag, while th dvchrm wildcard is effectively being constrained by the values in the DEEPD_CHRMS array;  So you produce 1 input array of files for every sample+dvchrm parir, with one list string/file name per array.  The rule will only begin when all array members are produced. It's then sorted by first deepdvchrm so they can be concatenated w/out another soort as all the chunks had been sorted already.
    priority: 44
    output:
        fin_fofn=MDIR
        + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.concat.vcf.gz.fofn",
        tmp_fofn=MDIR        + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.concat.vcf.gz.fofn.tmp",
    threads: 2
    resources:
        threads=2
    params:
        fn_stub="{sample}.{alnr}.deep.",
        cluster_sample=ret_sample,
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.deep.concat.fofn.bench.tsv"
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR + "{sample}/align/{alnr}/snv/deep/log/{sample}.{alnr}.deep.cocncat.fofn.log",
    shell:
        """
        (rm {output} 1> /dev/null  2> /dev/null ) || echo rmFailOK >> {log} && ls ./ >> {log};
        ### export LD_LIBRARY_PATH=$PWD/resources/libs/;
        for i in {input.chunk_tbi}; do
            ii=$(echo $i | perl -pe 's/\.tbi$//g'; );
            echo $ii >> {output.tmp_fofn};
        done;
        (workflow/scripts/sort_concat_chrm_list.py {output.tmp_fofn} {wildcards.sample}.{wildcards.alnr}.deep. {output.fin_fofn}) || echo "Python Script Error? CODE __ $? __" >> {log} && ls -lt {output} >> {log};

        """


rule deep_concat_index_chunks:
    input:
        fofn=MDIR
        + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.concat.vcf.gz.fofn",
    output:
        vcf=touch(
            MDIR + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.sort.vcf"
        ),
        vcfgz=touch(
            MDIR + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.sort.vcf.gz"
        ),
        vcfgztbi=touch(
            MDIR
            + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.sort.vcf.gz.tbi"
        ),
    threads: 8
    resources:
        vcpu=8,
        threads=8,
        partition="i4-5,i16-5,i32-5,i64-5,i96-5,i128-6"
    priority: 47
    params:
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        cluster_sample=ret_sample,
    resources:
        attempt_n=lambda wildcards, attempt:  (attempt + 0)
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.deep.merge.bench.tsv"
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/deep/log/{sample}.{alnr}.deep.snv.merge.sort.gatherered.log",
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
    clear_combined_deep_vcf,


rule clear_combined_deep_vcf:  # TARGET:  clear combined deep vcf so the chunks can be re-evaluated if needed.
    input:
        expand(
            MDIR + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.sort.vcf.gz",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
    priority: 42
    shell:
        "(rm {input}*  1> /dev/null  2> /dev/null ) || echo 'file not found for deletion: {input}';"


localrules:
    produce_deepDna_vcf,


rule produce_deepDna_vcf:  # TARGET: just gen deep calls
    input:
        expand(
            MDIR
            + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.sort.vcf.gz.tbi",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
    output:
        "gatheredall.deep",
    priority: 48
    log:
        "gatheredall.deep.log",
    shell:
        """( touch {output} ;

        {latency_wait}; ls {output} ) >> {log} 2>&1;
        """


localrules:
    prep_deep_chunkdirs,


rule prep_deep_chunkdirs:
    input:
        b=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.sort.bam",
    output:
        expand(
            MDIR + "{{sample}}/align/{{alnr}}/snv/deep/vcfs/{dvchrm}/{{sample}}.ready",
            dvchrm=DEEPD_CHRMS,
        ),
    log:
        MDIR + "{sample}/align/{alnr}/snv/deep/logs/{sample}.{alnr}.chunkdirs.log",
    shell:
        """
        ( echo {output}  ;
        mkdir -p $(dirname {output} );
        touch {output};
        ls {output}; ) > {log} 2>&1;
        """
