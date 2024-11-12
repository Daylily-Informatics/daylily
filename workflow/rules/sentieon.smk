import os

####### Sentieon
#
# Our current prod aligner
#


rule sentieon_bwa_sort:
    input:
        f1=getR1s,
        f2=getR2s,
    output:
        bamo=MDIR + "{sample}/align/sent/{sample}.sent.sort.bam",
        baio=MDIR + "{sample}/align/sent/{sample}.sent.sort.bam.bai"
    log:
        a=MDIR + "{sample}/align/sent/logs/{sample}.sent.sort.log",
    threads: config["sentieon"]["threads"]
    benchmark:
        repeat(MDIR + "{sample}/benchmarks/{sample}.sent.alNsort.bench.tsv", 0)
    priority: 5
    resources:
        partition=config['sentieon']['partition'],
        vcpu=config['sentieon']['threads'],
        threads=config['sentieon']['threads'],
    params:
        huref=config["supporting_files"]["files"]["huref"]["bwa_mem_sent"]["name"],
        cluster_end="echo done"
        if "cluster_end" not in config["sentieon"]
        else config["sentieon"]["cluster_end"],
        max_mem="130G"
        if "max_mem" not in config["sentieon"]
        else config["sentieon"]["max_mem"],
        K="-K 100000000" if "K" not in config["sentieon"] else config["sentieon"]["K"],
        cluster_sample=ret_sample,
        numactl=config["sentieon"]["numactl"],
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
        write_threads=config['sentieon']['write_threads'],
        igz_threads=config['sentieon']['igz_threads'],
        mbuffer_mem=config['sentieon']['mbuffer_mem'],
        bwa_model=config['sentieon']['bwa_model'],
        subsample_head=get_subsample_head,
        subsample_tail=get_subsample_tail,
    conda:
        config["sentieon"]["env_yaml"]
    shell:
        """


        TOKEN=$(curl -X PUT 'http://169.254.169.254/latest/api/token' -H 'X-aws-ec2-metadata-token-ttl-seconds: 21600');
        itype=$(curl -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-type);
        echo "INSTANCE TYPE: $itype" > {log};
        echo "INSTANCE TYPE: $itype";
        start_time=$(date +%s);
        export bwt_max_mem={params.max_mem} ;
        epocsec=$(date +'%s');

        ulimit -n 16384
        
        tdir=$(dirname {output.bamo})/tmpp;
        mkdir -p $tdir; 

        # Find the jemalloc library in the active conda environment
        jemalloc_path=$(find "$CONDA_PREFIX" -name "libjemalloc*" | grep -E '\.so|\.dylib' | head -n 1); 

        # Check if jemalloc was found and set LD_PRELOAD accordingly
        if [[ -n "$jemalloc_path" ]]; then
            export LD_PRELOAD="$jemalloc_path";
            echo "LD_PRELOAD set to: $LD_PRELOAD" >> {log.a};
        else
            echo "libjemalloc not found in the active conda environment $CONDA_PREFIX.";
            exit 3;
        fi
        LD_PRELOAD=$LD_PRELOAD /fsx/data/cached_envs/sentieon-genomics-202308.03/bin/sentieon bwa mem \
        -t {params.bwa_threads}  {params.K}  \
        -x {params.bwa_model} \
        -R '@RG\\tID:{params.rgid}_$epocsec\\tSM:{params.rgsm}\\tLB:{params.cluster_sample}{params.rglb}\\tPL:{params.rgpl}\\tPU:{params.rgpu}\\tCN:{params.rgcn}\\tPG:{params.rgpg}' \
        {params.huref} \
         {params.subsample_head} <(igzip {params.igz_threads} -q  {input.f1} )  {params.subsample_tail}  \
         {params.subsample_head} <(igzip  {params.igz_threads} -q  {input.f2} )  {params.subsample_tail}    \
        | mbuffer {params.mbuffer_mem} \
        | /fsx/data/cached_envs/sentieon-genomics-202308.03/bin/sentieon util sort \
        --thread_count {params.sort_threads} \
        --sortblock_thread_count {params.sort_threads} \
        --bam_compression 1 \
        --intermediate_compress_level 1  \
        --block_size {params.sort_thread_mem}   \
        --sam2bam \
        -o {output.bamo} - >> {log.a} 2>&1;

        samtools index -b -@ {threads} {output.bamo}  >> {log.a} 2>&1;

        end_time=$(date +%s);
    	elapsed_time=$((($end_time - $start_time) / 60));
	    echo "Elapsed-Time-min:\t$itype\t$elapsed_time\n";
        echo "Elapsed-Time-min:\t$itype\t$elapsed_time" >> {log} 2>&1;
        """
