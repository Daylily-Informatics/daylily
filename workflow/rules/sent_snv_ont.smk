import sys
import os

#
# This pipeline will use the ONT aligned cram directly and call variants
#

ALIGNERS_ONT = ["ont"]

rule sent_snv_ont:
    input:
        cram=MDIR + "{sample}/align/{alnr}/{sample}.cram",
        crai=MDIR + "{sample}/align/{alnr}/{sample}.cram.crai",
        d=MDIR + "{sample}/align/{alnr}/snv/sentdont/vcfs/{dchrm}/{sample}.ready",
    output:
        vcf=MDIR
        + "{sample}/align/{alnr}/snv/sentdont/vcfs/{dchrm}/{sample}.{alnr}.sentdont.{dchrm}.snv.vcf",
        gvcf=MDIR
        + "{sample}/align/{alnr}/snv/sentdont/vcfs/{dchrm}/{sample}.{alnr}.sentdont.{dchrm}.snv.gvcf",
        gvcfindex=MDIR
        + "{sample}/align/{alnr}/snv/sentdont/vcfs/{dchrm}/{sample}.{alnr}.sentdont.{dchrm}.snv.gvcf.idx",
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/sentdont/log/vcfs/{sample}.{alnr}.sentdont.{dchrm}.snv.log",
    threads: config['sentdont']['threads']
    conda:
        "../envs/sentieon_v0.1.yaml"
    priority: 45
    benchmark:
        repeat(
            MDIR + "{sample}/benchmarks/{sample}.{alnr}.sentdont.{dchrm}.bench.tsv",
            0
            if "bench_repeat" not in config["sentdont"]
            else config["sentdont"]["bench_repeat"],
        )
    resources:
        attempt_n=lambda wildcards, attempt:  (attempt + 0),
        partition=config['sentdont']['partition'],
        threads=config['sentdont']['threads'],
        vcpu=config['sentdont']['threads'],
	mem_mb=config['sentdont']['mem_mb'],
    params:
        schrm_mod=get_dchrm_day,
        use_threads=config['sentdont']['use_threads'],
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        model=config["sentdont"]["dna_scope_snv_model"],
	model_2=config["sentdont"]["dna_scope_apply_model"],
        cluster_sample=ret_sample,
    shell:
        """

        timestamp=$(date +%Y%m%d%H%M%S);
        export TMPDIR=/fsx/scratch/sentdont_tmp_$timestamp;
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
        
        LD_PRELOAD=$LD_PRELOAD /fsx/data/cached_envs/sentieon-genomics-202503/bin/sentieon driver -t {params.use_threads} \
            -r {params.huref} \
            -i {input.cram} \
            --interval {params.schrm_mod} \
            --algo DNAscope --model {params.model} \
            --emit_mode variant \
            {output.gvcf} >> {log} 2>&1;

         LD_PRELOAD=$LD_PRELOAD /fsx/data/cached_envs/sentieon-genomics-202503/bin/sentieon driver -t {params.use_threads} \
            -r {params.huref} \
            --algo DNAModelApply \
	    --model {params.model_2} \
            -v {output.gvcf} {output.vcf} >> {log} 2>&1;

        end_time=$(date +%s);
    	elapsed_time=$((($end_time - $start_time) / 60));
	    echo "Elapsed-Time-min:\t$itype\t$elapsed_time\n";
        echo "Elapsed-Time-min:\t$itype\t$elapsed_time" >> {log} 2>&1;

        touch {output.vcf};
        """


rule sentdont_sort_index_chunk_vcf:
    input:
        vcf=MDIR
        + "{sample}/align/{alnr}/snv/sentdont/vcfs/{dchrm}/{sample}.{alnr}.sentdont.{dchrm}.snv.vcf",
    priority: 46
    output:
        vcfsort=touch(MDIR
        + "{sample}/align/{alnr}/snv/sentdont/vcfs/{dchrm}/{sample}.{alnr}.sentdont.{dchrm}.snv.sort.vcf"),
        vcfgz=touch(MDIR
        + "{sample}/align/{alnr}/snv/sentdont/vcfs/{dchrm}/{sample}.{alnr}.sentdont.{dchrm}.snv.sort.vcf.gz"),
        vcftbi=touch(MDIR
        + "{sample}/align/{alnr}/snv/sentdont/vcfs/{dchrm}/{sample}.{alnr}.sentdont.{dchrm}.snv.sort.vcf.gz.tbi"),
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/sentdont/vcfs/{dchrm}/log/{sample}.{alnr}.sentdont.{dchrm}.snv.sort.vcf.gz.log",
    resources:
        vcpu=1,
        threads=1,
        partition="i192,i192mem"
    params:
        x='y',
        cluster_sample=ret_sample,
    threads: 64 #config["config"]["sort_index_sentdontna_chunk_vcf"]['threads']
    shell:
        """
        
        #bedtools sort -header -i {input.vcf} > {output.vcfsort} 2>> {log};
        #awk 'BEGIN{{header=1}} 
        #    header && /^#/ {{print; next}} 
        #    header && /^[^#]/ {{header=0; exit}}' {input.vcf} > {output.vcfsort} 2>> {log};
        #awk '/^[^#]/' {input.vcf} | sort --buffer-size=210G -T /fsx/scratch/ --parallel={threads} -k1,1V -k2,2n >> {output.vcfsort} 2>> {log};

        cp {input.vcf} {output.vcfsort} 2>> {log};
        touch {input.vcf};
        sleep 1;
        touch {output.vcfsort};
        bgzip  -@ {threads} {output.vcfsort} >> {log} 2>&1;
        touch {output.vcfsort};

        tabix -f -p vcf {output.vcfgz} >> {log} 2>&1;
        
        """


localrules:
    sentdont_concat_fofn,


rule sentdont_concat_fofn:
    input:
        chunk_tbi=sorted(
            expand(
                MDIR
                + "{{sample}}/align/{{alnr}}/snv/sentdont/vcfs/{ochm}/{{sample}}.{{alnr}}.sentdont.{ochm}.snv.sort.vcf.gz.tbi",
                ochm=SENTDONT_CHRMS,            ),            key=lambda x: float(                str(x.replace("~", ".").replace(":", "."))               .split("vcfs/")[1]                .split("/")[0]                .split("-")[0]            ),        ),
    # This expand pattern is neat.  the escaped {} remain acting as a snakemake wildcard and expect to be derived from the dag, while th dchrm wildcard is effectively being constrained by the values in the sentdont_CHRMS array;  So you produce 1 input array of files for every sample+dchrm parir, with one list string/file name per array.  The rule will only begin when all array members are produced. It's then sorted by first sentdontchrm so they can be concatenated w/out another soort as all the chunks had been sorted already.
    priority: 44
    output:
        fin_fofn=MDIR
        + "{sample}/align/{alnr}/snv/sentdont/{sample}.{alnr}.sentdont.snv.concat.vcf.gz.fofn",
        tmp_fofn=MDIR        + "{sample}/align/{alnr}/snv/sentdont/{sample}.{alnr}.sentdont.snv.concat.vcf.gz.fofn.tmp",
    threads: 1
    resources:
        threads=1
    params:
        fn_stub="{sample}.{alnr}.sentdont."
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.sentdont.concat.fofn.bench.tsv"
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR + "{sample}/align/{alnr}/snv/sentdont/log/{sample}.{alnr}.sentdont.cocncat.fofn.log",
    shell:
        """

        for i in {input.chunk_tbi}; do
            ii=$(echo $i | perl -pe 's/\.tbi$//g'; );
            echo $ii >> {output.tmp_fofn};
        done;
        (workflow/scripts/sort_concat_chrm_list.py {output.tmp_fofn} {wildcards.sample}.{wildcards.alnr}.sentdont. {output.fin_fofn}) >> {log} 2>&1;

        """


rule sentdont_concat_index_chunks:
    input:
        fofn=MDIR
        + "{sample}/align/{alnr}/snv/sentdont/{sample}.{alnr}.sentdont.snv.concat.vcf.gz.fofn",
    output:
        vcfgz=touch(
            MDIR + "{sample}/align/{alnr}/snv/sentdont/{sample}.{alnr}.sentdont.snv.sort.vcf.gz"
        ),
        vcfgztemp=temp(
            MDIR + "{sample}/align/{alnr}/snv/sentdont/{sample}.{alnr}.sentdont.snv.sort.temp.vcf.gz"
        ),
        vcfgztbi=touch(
            MDIR
            + "{sample}/align/{alnr}/snv/sentdont/{sample}.{alnr}.sentdont.snv.sort.vcf.gz.tbi"
        ),
    threads: 64
    resources:
        vcpu=64,
        threads=64,
        partition="i192,i192mem,i128"
    priority: 47
    params:
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        cluster_sample=ret_sample,
    resources:
        attempt_n=lambda wildcards, attempt:  (attempt + 0)
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.sentdont.merge.bench.tsv"
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/sentdont/log/{sample}.{alnr}.sentdont.snv.merge.sort.gatherered.log",
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

        rm -rf $(dirname {output.vcfgz})/vcfs >> {log} 2>&1;

        """

localrules:
    clear_combined_sentdont_vcf,


rule clear_combined_sentdont_vcf:  # TARGET:  clear combined sentdont vcf so the chunks can be re-evaluated if needed.
    input:
        expand(
            MDIR + "{sample}/align/{alnr}/snv/sentdont/{sample}.{alnr}.sentdont.snv.sort.vcf.gz",
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
    produce_sentdont_vcf,

 
rule produce_sentdont_vcf:  # TARGET: sentieon dnascope vcf
    input:
        expand(
            MDIR
            + "{sample}/align/{alnr}/snv/sentdont/{sample}.{alnr}.sentdont.snv.sort.vcf.gz.tbi",
            sample=SSAMPS,
            alnr=ALIGNERS_ONT,
        ),
    output:
        "gatheredall.sentdont",
    priority: 48
    threads: 1
    log:
        "gatheredall.sentdont.log",
    shell:
        """( touch {output} ;

        {latency_wait}; ls {output} ) >> {log} 2>&1;
        """


localrules:
    prep_sentdont_chunkdirs,


rule prep_sentdont_chunkdirs:
    input:
        cram=MDIR + "{sample}/align/{alnr}/{sample}.cram",
        crai=MDIR + "{sample}/align/{alnr}/{sample}.cram.crai",
    output:
        expand(
            MDIR + "{{sample}}/align/{{alnr}}/snv/sentdont/vcfs/{dchrm}/{{sample}}.ready",
            dchrm=SENTDONT_CHRMS,
        ),
    threads: 1
    log:
        MDIR + "{sample}/align/{alnr}/snv/sentdont/logs/{sample}.{alnr}.chunkdirs.log",
    shell:
        """
        ( echo {output}  ;
        mkdir -p $(dirname {output} );
        touch {output};
        ls {output}; ) > {log} 2>&1;
        """
