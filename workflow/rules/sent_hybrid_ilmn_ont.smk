import sys
import os

ALIGNERS_ONT = ["ont"]

rule sentdhio_snv:
    input:
        cram=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.cram",
        crai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.cram.crai",
        sr_cram=MDIR + "{sample}/align/sent/{sample}.sent.cram",
        d=MDIR + "{sample}/align/{alnr}/snv/sentdhio/vcfs/{dchrm}/{sample}.ready",
    output:
        vcf=MDIR
            + "{sample}/align/{alnr}/snv/sentdhio/vcfs/{dchrm}/{sample}.{alnr}.sentdhio.{dchrm}.snv.sort.vcf.gz",
        tbi=MDIR
            + "{sample}/align/{alnr}/snv/sentdhio/vcfs/{dchrm}/{sample}.{alnr}.sentdhio.{dchrm}.snv.sort.vcf.gz.tbi",           
    wildcard_constraints:
        alnr="|".join(ALIGNERS_ONT)
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/sentdhio/log/vcfs/{sample}.{alnr}.sentdhio.{dchrm}.snv.log",
    threads: config['sentdhio']['threads']
    conda:
        "../envs/sentieonHybrid_v0.1.yaml"
    priority: 45
    benchmark:
        repeat(
            MDIR + "{sample}/benchmarks/{sample}.{alnr}.sentdhio.{dchrm}.bench.tsv",
            0
            if "bench_repeat" not in config["sentdhio"]
            else config["sentdhio"]["bench_repeat"],
        )
    resources:
        attempt_n=lambda wildcards, attempt:  (attempt + 0),
        partition=config['sentdhio']['partition'],
        threads=config['sentdhio']['threads'],
        vcpu=config['sentdhio']['threads'],
    	mem_mb=config['sentdhio']['mem_mb'],
    params:
        schrm_mod=get_dchrm_day,
        use_threads=config["sentdhio"]["use_threads"],
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        model=config["sentdhio"]["dna_scope_snv_model"],
        cluster_sample=ret_sample,
        haploid_bed=get_haploid_bed_arg,
        diploid_bed=get_diploid_bed_arg,
        max_mem=config["sentdhio"]["max_mem"],
    shell:
        """
        export PATH=$PATH:/fsx/data/cached_envs/sentieon-genomics-202503.01.rc1/bin/
        export bwt_max_mem={params.max_mem} ;

        timestamp=$(date +%Y%m%d%H%M%S);
        export TMPDIR=/dev/shm/sentdhio_tmp_$timestamp;
        export SENTIEON_TEMP_DIR=$TMPDIR;
        export SENTIEON_TMPDIR=$TMPDIR;
        export TEMPDIR=$TMPDIR;
        export TEMP_DIR=$TMPDIR;
        export SENTIEON_TMP_DIR=$TMPDIR;
        export SENTIEON_TMPDIR=$TMPDIR;
        mkdir -p $TMPDIR;
        export APPTAINER_HOME=$TMPDIR;
        trap "rm -rf \"$TMPDIR\" || echo '$TMPDIR rm fails' >> {log} 2>&1" EXIT;

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

        ulimit -n 65536 || echo "ulimit mod failed" > {log} 2>&1;
        
        # Find the jemalloc library in the active conda environment
        jemalloc_path=$(find "$CONDA_PREFIX" -name "libjemalloc*" | grep -E '\.so|\.dylib' | head -n 1); 

        # Check if jemalloc was found and set LD_PRELOAD accordingly
        if [[ -n "$jemalloc_path" ]]; then
            LD_PRELOAD="$jemalloc_path";
            echo "LD_PRELOAD set to: $LD_PRELOAD" >> {log};
        else
            echo "libjemalloc not found in the active conda environment $CONDA_PREFIX.";
            exit 3;
        fi
        export cram_sid=$(samtools view -H {input.cram} | grep  '^@RG' | perl -pe 's/(^.*SM\:)(.*)(\w.*$)/$2/g;' | cut -d $'\t' -f 1 )
                

        LD_PRELOAD=$LD_PRELOAD sentieon-cli --verbose dnascope-hybrid \
            -t {params.use_threads} \
            -r  {params.huref} \
            --sr_aln {input.sr_cram} \
            --lr_aln {input.cram} \
            --rgsm {params.cluster_sample} \
            --sr_duplicate_marking none \
            --skip_svs \
            --skip_mosdepth \
            --skip_cnv \
            -m {params.model} \
            {params.diploid_bed}  {output.vcf} >> {log} 2>&1;


        end_time=$(date +%s);
    	elapsed_time=$((($end_time - $start_time) / 60));
	    echo "Elapsed-Time-min:\t$itype\t$elapsed_time\n";
        echo "Elapsed-Time-min:\t$itype\t$elapsed_time" >> {log} 2>&1;

        """


localrules:
    sentdhio_concat_fofn,


rule sentdhio_concat_fofn:
    input:
        chunk_tbi=sorted(
            expand(
                MDIR
                + "{{sample}}/align/{{alnr}}/snv/sentdhio/vcfs/{ochm}/{{sample}}.{{alnr}}.sentdhio.{ochm}.snv.sort.vcf.gz.tbi",
                ochm=SENTDHIO_CHRMS,            ),            key=lambda x: float(                str(x.replace("~", ".").replace(":", "."))               .split("vcfs/")[1]                .split("/")[0]                .split("-")[0]            ),        ),
    # This expand pattern is neat.  the escaped {} remain acting as a snakemake wildcard and expect to be derived from the dag, while th dchrm wildcard is effectively being constrained by the values in the sentdhio_CHRMS array;  So you produce 1 input array of files for every sample+dchrm parir, with one list string/file name per array.  The rule will only begin when all array members are produced. It's then sorted by first sentdhiochrm so they can be concatenated w/out another soort as all the chunks had been sorted already.
    priority: 44
    output:
        fin_fofn=MDIR
        + "{sample}/align/{alnr}/snv/sentdhio/{sample}.{alnr}.sentdhio.snv.concat.vcf.gz.fofn",
        tmp_fofn=MDIR        + "{sample}/align/{alnr}/snv/sentdhio/{sample}.{alnr}.sentdhio.snv.concat.vcf.gz.fofn.tmp",
    threads: 1
    wildcard_constraints:
        alnr="|".join(ALIGNERS_ONT)
    resources:
        threads=1
    params:
        fn_stub="{sample}.{alnr}.sentdhio."
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.sentdhio.concat.fofn.bench.tsv"
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR + "{sample}/align/{alnr}/snv/sentdhio/log/{sample}.{alnr}.sentdhio.cocncat.fofn.log",
    shell:
        """

        for i in {input.chunk_tbi}; do
            ii=$(echo $i | perl -pe 's/\.tbi$//g'; );
            echo $ii >> {output.tmp_fofn};
        done;
        (workflow/scripts/sort_concat_chrm_list.py {output.tmp_fofn} {wildcards.sample}.{wildcards.alnr}.sentdhio. {output.fin_fofn}) >> {log} 2>&1;

        """


rule sentdhio_concat_index_chunks:
    input:
        fofn=MDIR
        + "{sample}/align/{alnr}/snv/sentdhio/{sample}.{alnr}.sentdhio.snv.concat.vcf.gz.fofn",
    output:
        vcfgz=touch(
            MDIR + "{sample}/align/{alnr}/snv/sentdhio/{sample}.{alnr}.sentdhio.snv.sort.vcf.gz"
        ),
        vcfgztemp=temp(
            MDIR + "{sample}/align/{alnr}/snv/sentdhio/{sample}.{alnr}.sentdhio.snv.sort.temp.vcf.gz"
        ),
        vcfgztbi=touch(
            MDIR
            + "{sample}/align/{alnr}/snv/sentdhio/{sample}.{alnr}.sentdhio.snv.sort.vcf.gz.tbi"
        ),
    threads: 64
    resources:
        vcpu=64,
        threads=64,
        partition="i192,i192mem,i128"
    priority: 47
    wildcard_constraints:
        alnr="|".join(ALIGNERS_ONT)
    params:
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        cluster_sample=ret_sample,
    resources:
        attempt_n=lambda wildcards, attempt:  (attempt + 0)
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.sentdhio.merge.bench.tsv"
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/sentdhio/log/{sample}.{alnr}.sentdhio.snv.merge.sort.gatherered.log",
    shell:
        """

        
        touch {log};
        mkdir -p $(dirname {log});
        # This is acceptable bc I am concatenating from the same tools output, not across tools
        #touch {output.vcfgztemp};

        bcftools concat -a -d all --threads {threads} -f {input.fofn}  -O z -o {output.vcfgztemp} >> {log} 2>&1;

        export oldname=$(bcftools query -l {output.vcfgztemp} | head -n1) >> {log} 2>&1;
        echo -e "${{oldname}}\t{params.cluster_sample}" > {output.vcfgz}.rename.txt
        bcftools reheader -s {output.vcfgz}.rename.txt -o {output.vcfgz} {output.vcfgztemp} >> {log} 2>&1;
        bcftools index -f -t --threads {threads} -o {output.vcfgztbi} {output.vcfgz} >> {log} 2>&1;

        #rm -rf $(dirname {output.vcfgz})/vcfs >> {log} 2>&1;

        """

localrules:
    clear_combined_sentdhio_vcf,


rule clear_combined_sentdhio_vcf:  # TARGET:  clear combined sentdhio vcf so the chunks can be re-evaluated if needed.
    input:
        expand(
            MDIR + "{sample}/align/{alnr}/snv/sentdhio/{sample}.{alnr}.sentdhio.snv.sort.vcf.gz",
            sample=SSAMPS,
            alnr=ALIGNERS_ONT,
        ),
    threads: 2
    priority: 42
    shell:
        """
        rm {input}*  1> /dev/null  2> /dev/null ) || echo 'file not found for deletion: {input}';
        """


localrules:
    produce_sentdhio_vcf,


rule produce_sentdhio_vcf:  # TARGET: sentieon dnascope vcf
    input:
        expand(
            MDIR
            + "{sample}/align/{alnr}/snv/sentdhio/{sample}.{alnr}.sentdhio.snv.sort.vcf.gz.tbi",
            sample=SSAMPS,
            alnr=ALIGNERS_ONT,
        ),
    output:
        "gatheredall.sentdhio",
    priority: 48
    threads: 1
    log:
        "gatheredall.sentdhio.log",
    shell:
        """( touch {output} ;

        {latency_wait}; ls {output} ) >> {log} 2>&1;
        """


localrules:
    prep_sentdhio_chunkdirs,


rule prep_sentdhio_chunkdirs:
    input:
        DR=MDIR + "{sample}/{sample}.dirsetup.ready",
        r1=getR1s,
        r2=getR2s,
        cram=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.cram",
        crai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.cram.crai",
    output:
        expand(
            MDIR + "{{sample}}/align/{{alnr}}/snv/sentdhio/vcfs/{dchrm}/{{sample}}.ready",
            dchrm=SENTDHIO_CHRMS,
        ),
    threads: 1
    wildcard_constraints:
        alnr="|".join(ALIGNERS_ONT)
    log:
        MDIR + "{sample}/align/{alnr}/snv/sentdhio/logs/{sample}.{alnr}.chunkdirs.log",
    shell:
        """
        ( echo {output}  ;
        mkdir -p $(dirname {output} );
        touch {output};
        ls {output}; ) > {log} 2>&1;
        """
