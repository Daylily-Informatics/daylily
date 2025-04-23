import sys
import os



rule sentdhuo_snv:
    input:
        b=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
        bai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam.bai",
        cram=MDIR + "{sample}/align/{alnr}/{sample}.cram",
        crai=MDIR + "{sample}/align/{alnr}/{sample}.cram.crai",
        d=MDIR + "{sample}/align/{alnr}/snv/sentdhuo/vcfs/{dchrm}/{sample}.ready",
    output:
     vcf=temp(MDIR
        + "{sample}/align/{alnr}/snv/sentdhuo/vcfs/{dchrm}/{sample}.{alnr}.sentdhuo.{dchrm}.snv.vcf"),
     tvcf=temp(MDIR
        + "{sample}/align/{alnr}/snv/sentdhuo/vcfs/{dchrm}/{sample}.{alnr}.sentdhuo.{dchrm}.snv.vcf.tmp"),
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/sentdhuo/log/vcfs/{sample}.{alnr}.sentdhuo.{dchrm}.snv.log",
    threads: config['sentdhuo']['threads']
    conda:
        "../envs/sentdhuo_v0.2.yaml"
    priority: 45
    benchmark:
        repeat(
            MDIR + "{sample}/benchmarks/{sample}.{alnr}.sentdhuo.{dchrm}.bench.tsv",
            0
            if "bench_repeat" not in config["sentdhuo"]
            else config["sentdhuo"]["bench_repeat"],
        )
    resources:
        attempt_n=lambda wildcards, attempt:  (attempt + 0),
        partition=config['sentdhuo']['partition'],
        threads=config['sentdhuo']['threads'],
        vcpu=config['sentdhuo']['threads'],
	    mem_mb=config['sentdhuo']['mem_mb'],
    params:
        schrm_mod=get_dchrm_day,
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["namenogz"],
        model=config["sentdhuo"]["dna_scope_snv_model"],
        cluster_sample=ret_sample,
    shell:
        """



        timestamp=$(date +%Y%m%d%H%M%S);
        export TMPDIR=/fsx/scratch/sentdhuo_tmp_$timestamp;
        mkdir -p $TMPDIR;
        export APPTAINER_HOME=$TMPDIR;
        trap "rm -rf \"$TMPDIR\" || echo '$TMPDIR rm fails' >> {log} 2>&1" EXIT;
        tdir=$TMPDIR;

        if [ -z "$SENTIEON_LICENSE" ]; then
            echo "SENTIEON_LICENSE not set. Please set the SENTIEON_LICENSE environment variable to the license file path & make this update to your dyinit file as well." >> {log} 2>&1;
            exit 3;
        fi

        if [ ! -f "$SENTIEON_LICENSE" ]; then
            echo "The file referenced by SENTIEON_LICENSE ('$SENTIEON_LICENSE') does not exist. Please provide a valid file path." >> {log} 2>&1;
            exit 4;
        fi

        TOKEN=$(curl -X PUT 'http://169.254.169.254/latest/api/token' -H 'X-aws-ec2-metadata-token-ttl-seconds: 21600');
        itype=$(curl -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-type);
        echo "INSTANCE TYPE: $itype" > {log};
        echo "INSTANCE TYPE: $itype";
        start_time=$(date +%s);

        /fsx/data/cached_envs/sentieon-genomics-202503/bin/sentieon driver --thread_count {threads} --interval {params.schrm_mod} --reference {params.huref} --input {input.b} --algo DNAscope --pcr_indel_model none --model {params.model}  {output.tvcf} >> {log} 2>&1;
        /fsx/data/cached_envs/sentieon-genomics-202503/bin/sentieon driver -t {threads} -r {params.huref} --algo DNAModelApply --model {params.model} -v {output.tvcf} {output.vcf} >> {log} 2>&1;


        end_time=$(date +%s);
    	elapsed_time=$((($end_time - $start_time) / 60));
	    echo "Elapsed-Time-min:\t$itype\t$elapsed_time\n";
        echo "Elapsed-Time-min:\t$itype\t$elapsed_time" >> {log} 2>&1;

        touch {output.vcf};
        """


rule sentdhuo_sort_index_chunk_vcf:
    input:
        vcf=MDIR
        + "{sample}/align/{alnr}/snv/sentdhuo/vcfs/{dchrm}/{sample}.{alnr}.sentdhuo.{dchrm}.snv.vcf",
    priority: 46
    output:
        vcfsort=touch(MDIR
        + "{sample}/align/{alnr}/snv/sentdhuo/vcfs/{dchrm}/{sample}.{alnr}.sentdhuo.{dchrm}.snv.sort.vcf"),
        vcfgz=touch(MDIR
        + "{sample}/align/{alnr}/snv/sentdhuo/vcfs/{dchrm}/{sample}.{alnr}.sentdhuo.{dchrm}.snv.sort.vcf.gz"),
        vcftbi=touch(MDIR
        + "{sample}/align/{alnr}/snv/sentdhuo/vcfs/{dchrm}/{sample}.{alnr}.sentdhuo.{dchrm}.snv.sort.vcf.gz.tbi"),
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/sentdhuo/vcfs/{dchrm}/log/{sample}.{alnr}.sentdhuo.{dchrm}.snv.sort.vcf.gz.log",
    resources:
        vcpu=1,
        threads=1,
        partition="i192,i192mem"
    params:
        x='y',
        cluster_sample=ret_sample,
    threads: 1 #config["config"]["sort_index_sentdhuona_chunk_vcf"]['threads']
    shell:
        """
        bedtools sort -header -i {input.vcf} > {output.vcfsort} 2>> {log};
        
        bgzip {output.vcfsort} >> {log} 2>&1;
        touch {output.vcfsort};

        tabix -f -p vcf {output.vcfgz} >> {log} 2>&1;
        
        """


localrules:
    sentdhuo_concat_fofn,


rule sentdhuo_concat_fofn:
    input:
        chunk_tbi=sorted(
            expand(
                MDIR
                + "{{sample}}/align/{{alnr}}/snv/sentdhuo/vcfs/{ochm}/{{sample}}.{{alnr}}.sentdhuo.{ochm}.snv.sort.vcf.gz.tbi",
                ochm=SENTDHUO_CHRMS,            ),            key=lambda x: float(                str(x.replace("~", ".").replace(":", "."))               .split("vcfs/")[1]                .split("/")[0]                .split("-")[0]            ),        ),
    # This expand pattern is neat.  the escaped {} remain acting as a snakemake wildcard and expect to be derived from the dag, while th dchrm wildcard is effectively being constrained by the values in the sentdhuo_CHRMS array;  So you produce 1 input array of files for every sample+dchrm parir, with one list string/file name per array.  The rule will only begin when all array members are produced. It's then sorted by first sentdhuochrm so they can be concatenated w/out another soort as all the chunks had been sorted already.
    priority: 44
    output:
        fin_fofn=MDIR
        + "{sample}/align/{alnr}/snv/sentdhuo/{sample}.{alnr}.sentdhuo.snv.concat.vcf.gz.fofn",
        tmp_fofn=MDIR        + "{sample}/align/{alnr}/snv/sentdhuo/{sample}.{alnr}.sentdhuo.snv.concat.vcf.gz.fofn.tmp",
    threads: 1
    resources:
        threads=1
    params:
        fn_stub="{sample}.{alnr}.sentdhuo."
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.sentdhuo.concat.fofn.bench.tsv"
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR + "{sample}/align/{alnr}/snv/sentdhuo/log/{sample}.{alnr}.sentdhuo.cocncat.fofn.log",
    shell:
        """

        for i in {input.chunk_tbi}; do
            ii=$(echo $i | perl -pe 's/\.tbi$//g'; );
            echo $ii >> {output.tmp_fofn};
        done;
        (workflow/scripts/sort_concat_chrm_list.py {output.tmp_fofn} {wildcards.sample}.{wildcards.alnr}.sentdhuo. {output.fin_fofn}) >> {log} 2>&1;

        """


rule sentdhuo_concat_index_chunks:
    input:
        fofn=MDIR
        + "{sample}/align/{alnr}/snv/sentdhuo/{sample}.{alnr}.sentdhuo.snv.concat.vcf.gz.fofn",
    output:
        vcfgz=touch(
            MDIR + "{sample}/align/{alnr}/snv/sentdhuo/{sample}.{alnr}.sentdhuo.snv.sort.vcf.gz"
        ),
        vcfgztemp=temp(
            MDIR + "{sample}/align/{alnr}/snv/sentdhuo/{sample}.{alnr}.sentdhuo.snv.sort.temp.vcf.gz"
        ),
        vcfgztbi=touch(
            MDIR
            + "{sample}/align/{alnr}/snv/sentdhuo/{sample}.{alnr}.sentdhuo.snv.sort.vcf.gz.tbi"
        ),
    threads: 4
    resources:
        vcpu=4,
        threads=4,
        partition="i192,i192mem"
    priority: 47
    params:
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["namenogz"],
        cluster_sample=ret_sample,
    resources:
        attempt_n=lambda wildcards, attempt:  (attempt + 0)
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.sentdhuo.merge.bench.tsv"
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/sentdhuo/log/{sample}.{alnr}.sentdhuo.snv.merge.sort.gatherered.log",
    shell:
        """

        touch {log};
        mkdir -p $(dirname {log});
        # This is acceptable bc I am concatenating from the same tools output, not across tools
        touch {output.vcfgztemp};
        bcftools concat -a -d all --threads {threads} -f {input.fofn}  -O z -o {output.vcfgz};
        bcftools index -f -t --threads {threads} -o {output.vcfgztbi} {output.vcfgz};

        rm -rf $(dirname {output.vcfgz})/vcfs >> {log} 2>&1;

        """

localrules:
    clear_combined_sentdhuo_vcf,


rule clear_combined_sentdhuo_vcf:  # TARGET:  clear combined sentdhuo vcf so the chunks can be re-evaluated if needed.
    input:
        expand(
            MDIR + "{sample}/align/{alnr}/snv/sentdhuo/{sample}.{alnr}.sentdhuo.snv.sort.vcf.gz",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
    threads: 2
    priority: 42
    shell:
        """
        rm {input}*  1> /dev/null  2> /dev/null ) || echo 'file not found for deletion: {input}';
        """


localrules:
    produce_sentdhuo_vcf,


rule produce_sentdhuo_vcf:  # TARGET: sentieon dnascope vcf
    input:
        expand(
            MDIR
            + "{sample}/align/{alnr}/snv/sentdhuo/{sample}.{alnr}.sentdhuo.snv.sort.vcf.gz.tbi",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
    output:
        "gatheredall.sentdhuo",
    priority: 48
    threads: 1
    log:
        "gatheredall.sentdhuo.log",
    shell:
        """( touch {output} ;

        {latency_wait}; ls {output} ) >> {log} 2>&1;
        """


localrules:
    prep_sentdhuo_chunkdirs,


rule prep_sentdhuo_chunkdirs:
    input:
        b=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
        cram=MDIR + "{sample}/align/{alnr}/{sample}.cram",
        crai=MDIR + "{sample}/align/{alnr}/{sample}.cram.crai",
    output:
        expand(
            MDIR + "{{sample}}/align/{{alnr}}/snv/sentdhuo/vcfs/{dchrm}/{{sample}}.ready",
            dchrm=SENTDHUO_CHRMS,
        ),
    threads: 1
    log:
        MDIR + "{sample}/align/{alnr}/snv/sentdhuo/logs/{sample}.{alnr}.chunkdirs.log",
    shell:
        """
        ( echo {output}  ;
        mkdir -p $(dirname {output} );
        touch {output};
        ls {output}; ) > {log} 2>&1;
        """
