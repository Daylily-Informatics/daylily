import os

####### Sentieon
#
# Our current prod aligner
#


rule sentieon_bwa_sort:  #TARGET: sent bwa sort
    input:
        DR=MDIR + "{sample}/{sample}.dirsetup.ready",
        f1=getR1s,
        f2=getR2s,
    output:
        bamo=temp(MDIR + "{sample}/align/sent/{sample}.sent.sort.bam"),
        baio=temp(MDIR + "{sample}/align/sent/{sample}.sent.sort.bam.bai")
    log: MDIR + "{sample}/align/sent/logs/{sample}.sent.sort.log",
    threads: config["sentieon"]["threads"]
    benchmark:
        repeat(MDIR + "{sample}/benchmarks/{sample}.sent.alNsort.bench.tsv", 0)
    priority: 5
    resources:
        partition=config['sentieon']['partition'],
        vcpu=config['sentieon']['threads'],
        threads=config['sentieon']['threads'],
        mem_mb=config['sentieon']['mem_mb'],
        constraint=config['sentieon']['constraint'],
    params:
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        max_mem="130G"
        if "max_mem" not in config["sentieon"]
        else config["sentieon"]["max_mem"],
        sent_opts=config["sentieon"]["sent_opts"],
        cluster_sample=ret_sample,
        bwa_threads=config["sentieon"]["bwa_threads"],
        rgpl="presumedILLUMINA",  # ideally: passed in technology # nice to get to this point: https://support.sentieon.com/appnotes/read_groups/ :\ : note, the default sample name contains the RU_EX_SQ_Lane (0 for combined)
        rgpu="presumedCombinedLanes",  # ideally flowcell_lane(s)
        rgsm=ret_sample,  # samplename
        rgid=ret_sample,  # ideally samplename_flowcell_lane(s)_barcode  ! Imp this is unique, I add epoc seconds to the end of start of this rule
        rglb="_presumedNoAmpWGS",  # prepend with cluster sample nanme ideally samplename_libprep
        rgcn="CenterName",  # center name
        rgpg="sentieonBWAmem",  #program
        sort_thread_mem=config['sentieon']['sort_thread_mem'],
        sort_threads=config['sentieon']['sort_threads'],
        igz=config['sentieon']['igz'],
        mbuffer=config['sentieon']['mbuffer'],
        bwa_model=config['sentieon']['bwa_model'],
        subsample_head=get_subsample_head,
        subsample_tail=get_subsample_tail,
    conda:
        config["sentieon"]["env_yaml"]
    shell:
        """

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
        start_time=$(date +%s);
        export bwt_max_mem={params.max_mem} ;
        epocsec=$(date +'%s');

        ulimit -n 65536 || echo "ulimit mod failed" > {log} 2>&1;
        
        timestamp=$(date +%Y%m%d%H%M%S);
        TMPDIR=/dev/shm/sentieon_tmp_$timestamp;
        export SENTIEON_TMPDIR=$TMPDIR;

        mkdir -p $TMPDIR;
        APPTAINER_HOME=$TMPDIR;
        trap "rm -rf \"$TMPDIR\" || echo '$TMPDIR rm fails' >> {log} 2>&1" EXIT;
        tdir=$TMPDIR;

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
        LD_PRELOAD=$LD_PRELOAD /fsx/data/cached_envs/sentieon-genomics-202503.01.rc1/bin/sentieon bwa mem \
        -t {params.bwa_threads}  {params.sent_opts}  \
        -x {params.bwa_model} \
        -R "@RG\\tID:{params.cluster_sample}-$epocsec\\tSM:{params.cluster_sample}\\tLB:{params.cluster_sample}-LB-1\\tPL:ILLUMINA" \
        {params.huref} \
         {params.subsample_head} <( {params.igz} -q  {input.f1} )  {params.subsample_tail}  \
         {params.subsample_head} <( {params.igz} -q  {input.f2} )  {params.subsample_tail} {params.mbuffer} \
        | /fsx/data/cached_envs/sentieon-genomics-202503.01.rc1/bin/sentieon  util sort \
        -t  {params.sort_threads} \
        --reference {params.huref} \
        --cram_write_options version=3.0,compressor=rans,lazy_quality=true \
        --sortblock_thread_count {params.sort_threads} \
        --bam_compression 1 \
	    --temp_dir $TMPDIR \
        --intermediate_compress_level 1  \
        --block_size {params.sort_thread_mem}   \
        --sam2bam \
        -o {output.bamo} - >> {log} 2>&1;

        #samtools index -b -@ {threads} {output.baio}  >> {log} 2>&1;

        end_time=$(date +%s);
    	elapsed_time=$((($end_time - $start_time) / 60));
        echo "Elapsed-Time-min:\t$itype\t$elapsed_time" >> {log} 2>&1;
        rm -rf $tdir;
        """

localrules: produce_sentieon_bwa_sort_bam,

rule produce_sentieon_bwa_sort_bam:  # TARGET: produce_sentieon_bwa_sort_bam
     input:
         expand(MDIR + "{sample}/align/sent/{sample}.sent.sort.bam", sample=SAMPS)
