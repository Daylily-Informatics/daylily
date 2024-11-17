####### BWA MEM2 
# An accelerated version of BWA MEM

rule bwa_mem2_sort:
    """https://github.com/bwa-mem2/bwa-mem2"""
    """    Also used here, a common bfx toolset, samtools"""
    """    https://github.com/samtools/samtools"""
    input:
        DR=MDIR + "{sample}/{sample}.dirsetup.ready",
        f1=getR1s,
        f2=getR2s,
    output:
        bami=temp(MDIR + "{sample}/align/bwa2a/{sample}.bwa2a.sort.bam.bai"),
        bamo=temp(MDIR + "{sample}/align/bwa2a/{sample}.bwa2a.sort.bam"),
    priority: 49
    log:
        MDIR + "{sample}/align/bwa2a/logs/{sample}.bwa2a_sort.log",
    resources:
        threads=config['bwa_mem2a_aln_sort']['threads'],
        mem_mb=config['bwa_mem2a_aln_sort']['mem_mb'],
        partition=config['bwa_mem2a_aln_sort']['partition'],
        vcpu=config['bwa_mem2a_aln_sort']['threads']
    threads: config['bwa_mem2a_aln_sort']['threads']
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.bwa2a.alNsort.bench.tsv"
    params:
        cluster_slots=94,
        cluster_sample=ret_sample,
        mdir=MDIR,
        samtmpd="/align/bwa2a/logs/sam_tmpd/",
        write_threads=config["bwa_mem2a_aln_sort"]["write_threads"],
        sort_threads= config["bwa_mem2a_aln_sort"]["sort_threads"],
        benchmark_runs=config["bwa_mem2a_aln_sort"]["softclip_alts"],
        softclip_alts=config["bwa_mem2a_aln_sort"]["softclip_alts"],
        bwa_mem2a_cmd=config["bwa_mem2a_aln_sort"]["cmd"],
        ldpre=config['bwa_mem2a_aln_sort']['ldpre'],
        k=get_bwa_kmer_size,  # little kay
        K=config["bwa_mem2a_aln_sort"]["K"],  # BIG KAY
        sort_thread_mem=config["bwa_mem2a_aln_sort"]["sort_thread_mem"],
        huref=config["supporting_files"]["files"]["huref"]["bwa_mem_index_vanilla"]["name"],
        rgpl="presumedILLUMINA",  # ideally: passed in technology # nice to get to this point: https://support.sentieon.com/appnotes/read_groups/  :: note, the default sample name contains the RU_EX_SQ_Lane (0 for combined)
        rgpu="presumedCombinedLanes",  # ideally flowcell_lane(s)
        rgsm='x', # ret_sample,  # samplename
        rgid='x', #ret_sample,  # ideally samplename_flowcell_lane(s)_barcode  ! Imp this is unique, I add epoc seconds to the end of start of this rule
        rglb="_presumedNoAmpWGS",  # prepend sample_name in shell block ideally samplename_libprep
        rgcn="CenterName",  # center name
        rgpg="bwamem2",  #program
        subsample_head=get_subsample_head,
        subsample_tail=get_subsample_tail,
        bwa_threads=config["bwa_mem2a_aln_sort"]["bwa_threads"],
        samp=get_samp_name,
        mbuffer_mem=config["bwa_mem2a_aln_sort"]["mbuffer_mem"],
        igz_threads=config['bwa_mem2a_aln_sort']['igz_threads']
    conda:
        config["bwa_mem2a_aln_sort"]["env_yaml"]
    shell:
        """

        TOKEN=$(curl -X PUT 'http://169.254.169.254/latest/api/token' -H 'X-aws-ec2-metadata-token-ttl-seconds: 21600');
        itype=$(curl -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-type);
        echo "INSTANCE TYPE: $itype" > {log};
        echo "INSTANCE TYPE: $itype";
        start_time=$(date +%s);
        ulimit -n 65536 || echo "ulimit mod failed";

        export tdir={params.mdir}/{params.samp}/{params.samtmpd};
        mkdir -p $tdir ;
        epocsec=$(date +'%s');
        
        # Find the jemalloc library in the active conda environment
        jemalloc_path=$(find "$CONDA_PREFIX" -name "libjemalloc*" | grep -E '\.so|\.dylib' | head -n 1); 

        # Check if jemalloc was found and set LD_PRELOAD accordingly
        if [[ -n "$jemalloc_path" ]]; then
            export LD_PRELOAD="$jemalloc_path";
            echo "LD_PRELOAD set to: $LD_PRELOAD" >> {log};
        else
            echo "libjemalloc not found in the active conda environment $CONDA_PREFIX.";
            exit 3;
        fi


         numactl --cpunodebind=0 --membind=0 {params.bwa_mem2a_cmd} mem \
         -R '@RG\\tID:{params.rgid}_$epocsec\\tSM:{params.rgsm}\\tLB:{params.samp}{params.rglb}\\tPL:{params.rgpl}\\tPU:{params.rgpu}\\tCN:{params.rgcn}\\tPG:{params.rgpg}' \
         {params.softclip_alts}  {params.K} {params.k} \
         -t {params.bwa_threads}  \
         {params.huref} \
         {params.subsample_head} <(numactl --cpunodebind=1 --membind=1 igzip {params.igz_threads} -q  {input.f1} )  {params.subsample_tail}  \
         {params.subsample_head} <(numactl --cpunodebind=1 --membind=1 igzip  {params.igz_threads} -q  {input.f2} )  {params.subsample_tail} {params.mbuffer_mem} \
        | numactl --cpunodebind=1 --membind=1 samtools sort \
        -l 1  \
        -m {params.sort_thread_mem}   \
         -@  {params.sort_threads} \
         -T $tdir \
         -O BAM  \
         --write-index \
         -o {output.bamo}##idx##{output.bami} >> {log} 2>&1;

        end_time=$(date +%s);
    	elapsed_time=$((($end_time - $start_time) / 60));
	    echo "Elapsed-Time-min:\t$itype\t$elapsed_time\n";
        echo "Elapsed-Time-min:\t$itype\t$elapsed_time" >> {log} 2>&1;
        rm -rf $tdir;
        """


localrules: produce_bwa_mem2,

rule produce_bwa_mem2:  # TARGET: only produce bwamem2a
     input:
         expand(MDIR + "{sample}/align/bwa2a/{sample}.bwa2a.sort.bam", sample=SAMPS)

	  
 