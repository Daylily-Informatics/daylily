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
        + "{sample}/align/{alnr}/snv/deep/log/{sample}.{alnr}.deep.{dvchrm}.snv.log",
    threads: config['deepvariant']['threads']
    container:
        "docker://daylilyinformatics/deepvariant-avx512:1.5.0"
    priority: 45
    resources:
        vcpu=config['deepvariant']['threads'],
        threads=config['deepvariant']['threads'],
        partition=config['deepvariant']['partition'],
        mem_mb=config['deepvariant']['mem_mb'],
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
        numa=config['deepvariant']['numa'],
        cpre="" if "b37" == config['genome_build'] else "chr",
    shell:
        """
        touch {log};
        TOKEN=$(curl -X PUT 'http://169.254.169.254/latest/api/token' -H 'X-aws-ec2-metadata-token-ttl-seconds: 21600');
        itype=$(curl -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-type);
        echo "INSTANCE TYPE: $itype" > {log};

        
        # Log the start time as 0 seconds
        start_time=$(date +%s);
        echo "Start-Time-sec:$itype\t0" >> {log} 2>&1;

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

        export APPTAINER_HOME=/fsx/scratch;
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
        --dry_run=false >> {log} 2>&1;

        end_time=$(date +%s);
        elapsed_time=$((($end_time - $start_time) / 60));

        # Log the elapsed time
        echo "Elapsed-Time-min:\t$itype\t$elapsed_time";
        echo "Elapsed-Time-min:\t$itype\t$elapsed_time" >> {log} 2>&1;
        """



rule dv_sort_index_chunk_vcf:
    input:
        vcf=MDIR
        + "{sample}/align/{alnr}/snv/deep/vcfs/{dvchrm}/{sample}.{alnr}.deep.{dvchrm}.snv.vcf",
        gvcf=MDIR
        + "{sample}/align/{alnr}/snv/deep/vcfs/{dvchrm}/{sample}.{alnr}.deep.{dvchrm}.snv.g.vcf",
    priority: 46
    output:
        vcfsort=MDIR
        + "{sample}/align/{alnr}/snv/deep/vcfs/{dvchrm}/{sample}.{alnr}.deep.{dvchrm}.snv.sort.vcf",
        vcfgz=MDIR
        + "{sample}/align/{alnr}/snv/deep/vcfs/{dvchrm}/{sample}.{alnr}.deep.{dvchrm}.snv.sort.vcf.gz",
        vcftbi=MDIR
        + "{sample}/align/{alnr}/snv/deep/vcfs/{dvchrm}/{sample}.{alnr}.deep.{dvchrm}.snv.sort.vcf.gz.tbi",
	gvcfsort=MDIR
        + "{sample}/align/{alnr}/snv/deep/vcfs/{dvchrm}/{sample}.{alnr}.deep.{dvchrm}.snv.g.sort.vcf",
        gvcfgz=MDIR
        + "{sample}/align/{alnr}/snv/deep/vcfs/{dvchrm}/{sample}.{alnr}.deep.{dvchrm}.snv.g.sort.vcf.gz",
        gvcftbi=MDIR
        + "{sample}/align/{alnr}/snv/deep/vcfs/{dvchrm}/{sample}.{alnr}.deep.{dvchrm}.snv.g.sort.vcf.gz.tbi",
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/deep/vcfs/{dvchrm}/log/{sample}.{alnr}.deep.{dvchrm}.snv.sort.vcf.gz.log",
    resources:
        vcpu=4,
        threads=4,
        partition=config['deepvariant']['partition_other'],
    params:
        cluster_sample=ret_sample,
    threads: 4 #config["config"]["sort_index_deepDna_chunk_vcf"]['threads']
    shell:
        """
        (rm  {output} 1>  /dev/null  2> /dev/null )  || echo rmfailed > {log};
        (bedtools sort -header -i {input.vcf} > {output.vcfsort} 2>> {log}) || exit 1233;
        (
        bgzip {output.vcfsort};        
        touch {output.vcfsort};
        tabix -f -p vcf {output.vcfgz};
        touch {output.vcftbi};
        ) >> {log} 2>&1;


        (bedtools sort -header -i {input.gvcf} > {output.gvcfsort} 2>> {log}) || exit 1233;
        (
        bgzip {output.gvcfsort};
        touch {output.gvcfsort};
        tabix -f -p vcf {output.gvcfgz};
        touch {output.gvcftbi};
        ) >> {log} 2>&1;
        
        ls {output} >> {log} 2>&1 ;
        
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
	    chunk_gtbi=sorted(
            expand(
                MDIR                + "{{sample}}/align/{{alnr}}/snv/deep/vcfs/{dvchm}/{{sample}}.{{alnr}}.deep.{dvchm}.snv.g.sort.vcf.gz.tbi",                dvchm=DEEPD_CHRMS,            ),            key=lambda x: float(                str(x.replace("~", ".").replace(":", "."))               .split("vcfs/")[1]                .split("/")[0]                .split("-")[0]            ),        ),	
    # This expand pattern is neat.  the escaped {} remain acting as a snakemake wildcard and expect to be derived from the dag, while th dvchrm wildcard is effectively being constrained by the values in the DEEPD_CHRMS array;  So you produce 1 input array of files for every sample+dvchrm parir, with one list string/file name per array.  The rule will only begin when all array members are produced. It's then sorted by first deepdvchrm so they can be concatenated w/out another soort as all the chunks had been sorted already.
    priority: 44
    output:
        fin_fofn=MDIR
        + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.concat.vcf.gz.fofn",
        tmp_fofn=MDIR        + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.concat.vcf.gz.fofn.tmp",
    	gfin_fofn=MDIR
        + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.g.concat.vcf.gz.fofn",
        gtmp_fofn=MDIR        + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.g.concat.vcf.gz.fofn.tmp",
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

        for i in {input.chunk_tbi}; do
            ii=$(echo $i | perl -pe 's/\.tbi$//g'; );
            echo $ii >> {output.tmp_fofn};
        done;
        (workflow/scripts/sort_concat_chrm_list.py {output.tmp_fofn} {wildcards.sample}.{wildcards.alnr}.deep. {output.fin_fofn}) || echo "Python Script Error? CODE __ $? __" >> {log} && ls -lt {output.fin_fofn} {output.tmp_fofn} >> {log};


        for i in {input.chunk_gtbi}; do
            ii=$(echo $i | perl -pe 's/\.tbi$//g'; );
            echo $ii >> {output.gtmp_fofn};
        done;
        (workflow/scripts/sort_concat_chrm_list.py {output.gtmp_fofn} {wildcards.sample}.{wildcards.alnr}.deep. {output.gfin_fofn}) || echo "Python Script Error? CODE __ $? __" >> {log} && ls -lt {output} >> {log};
        """


rule deep_concat_index_chunks:
    input:
        fofn=MDIR
        + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.concat.vcf.gz.fofn",
        tmp_fofn=MDIR        + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.concat.vcf.gz.fofn.tmp",
    	gfofn=MDIR
        + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.g.concat.vcf.gz.fofn",
        gtmp_fofn=MDIR        + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.g.concat.vcf.gz.fofn.tmp",
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
	gvcf=touch(
            MDIR + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.g.sort.vcf"
        ),
        gvcfgz=touch(
            MDIR + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.g.sort.vcf.gz"
        ),
        gvcfgztbi=touch(
            MDIR
            + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.g.sort.vcf.gz.tbi"
        ),
        bcf=touch(
            MDIR + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.sort.bcf"
        ),
        bcfi=touch(
            MDIR + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.sort.bcf.csi"
        ),
        gbcf=touch(
            MDIR + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.g.sort.bcf"
        ),
        gbcfi=touch(
            MDIR + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.g.sort.bcf.csi"
        ),
    threads: 4
    resources:
        vcpu=4,
        threads=4,
        partition=config['deepvariant']['partition_other'],
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

        bcftools concat -a -d all --threads {threads} -f {input.fofn}  -O v -o {output.vcf};
        bcftools view -O z -o {output.vcfgz} {output.vcf};
        bcftools index -f -t --threads {threads} -o {output.vcfgztbi} {output.vcfgz};
        stats_f=$(echo "{output.vcfgz}.bcf.stats");
        bcftools stats -F {params.huref}  {output.vcfgz} > $stats_f;

        # Convert to BCF and index it
        bcftools view -O b -o {output.bcf} --threads {threads} {output.vcfgz};
        bcftools index --threads {threads} {output.bcf};

        bcftools concat -a -d all --threads {threads} -f {input.gfofn}  -O v -o {output.gvcf};
        bcftools view -O z -o {output.gvcfgz} {output.gvcf};
        bcftools index -f -t --threads {threads} -o {output.gvcfgztbi} {output.gvcfgz};
        stats_f=$(echo "{output.gvcfgz}.bcf.stats");
        bcftools stats -F {params.huref}  {output.gvcfgz} > $stats_f;


        # Convert to BCF and index it
        bcftools view -O b -o {output.gbcf} --threads {threads} {output.gvcfgz};
        bcftools index --threads {threads} {output.gbcf};

        rm -rf $(dirname {output.gbcf})/vcfs >> {log} 2>&1;  # clean up all the crap
        {latency_wait};
        touch {log};
        """


localrules:
    clear_combined_deep_vcf,


rule clear_combined_deep_vcf:  # TARGET:  clear combined deep vcf so the chunks can be re-evaluated if needed.
    input:
        vcf=expand(
            MDIR + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.sort.vcf.gz",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
	    gvcf=expand(
            MDIR + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.g.sort.vcf.gz",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),	
    priority: 42
    shell:
        "(rm {input.vcf}* {input.gvcf}*  1> /dev/null  2> /dev/null ) || echo 'file not found for deletion: {input}';"


rule produce_deepDna_vcf:  # TARGET: just gen deep calls
    input:
        vcftb=expand(
            MDIR
            + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.sort.vcf.gz",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
        vcftbi=expand(
            MDIR
            + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.sort.vcf.gz.tbi",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
        gvcf=expand(
            MDIR
            + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.g.sort.vcf.gzi",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
	    gvcftbi=expand(
            MDIR
            + "{sample}/align/{alnr}/snv/deep/{sample}.{alnr}.deep.snv.g.sort.vcf.gz.tbi",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
    output:
        "gatheredall.deep",
    threads: 4
    priority: 48
    log:
        "gatheredall.deep.log",
    conda:
        "../envs/vanilla_v0.1.yaml"
    shell:
        """
        # Convert VCF to BCF and index it
        for vcf in {input.vcftb}; do
            bcf="${{vcf%.vcf.gz}}.bcf";
            bcftools view -O b -o $bcf --threads {threads} $vcf && bcftools index --threads 4 $bcf;
        done;

        # Convert GVCF to BCF and index it
        for gvcf in {input.gvcf}; do
            gbcf="${{gvcf%.vcf.gzi}}.bcf";
            bcftools view -O b -o $gbcf --threads {threads} $gvcf && bcftools index --threads 4 $gbcf;
        done;

        # Mark the output as completed
        touch {output};

        # Log completion and list output
        {latency_wait};
        ls {output} >> {log} 2>&1;

        #( touch {output} ;
        #echo "bcftools view -O b -o output.bcf --threads 4 input.vcf.gz && bcftools index --threads 4 output.bcf""
        {latency_wait}; 
        ls {output}  >> {log} 2>&1;
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
 
