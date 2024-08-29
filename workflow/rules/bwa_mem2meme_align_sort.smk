####### BWA MEM2 -- Version meme
#
# An accelerated version of BWA MEM2


rule bwa_mem2meme_aln_sort:
    input:
        DR=MDIR + "{sample}/{sample}.dirsetup.ready",
        f1=getR1s,  # method defined in fastp.smk
        f2=getR2s,  # method defined in fastp.smk
    output:
        bami=temp(MDIR + "{sample}/align/bwa2m/{sample}.bwa2m.sort.bam.bai"),
        bamo=temp(MDIR + "{sample}/align/bwa2m/{sample}.bwa2m.sort.bam"),
    priority: 49
    resources:
        threads=config["bwa_mem2meme_aln_sort"]["threads"],
        partition=config["bwa_mem2meme_aln_sort"]["partition"],
        vcpu=config["bwa_mem2meme_aln_sort"]["threads"],
        mem_mb=config["bwa_mem2meme_aln_sort"]["mem_mb"],
    log:
        MDIR + "{sample}/align/bwa2m/logs/align_sort/{sample}.bwa2m.sort.log",
    threads: config["bwa_mem2meme_aln_sort"]["threads"]  # eval( "lambda wildcards,input,threads,  attempt:int( attempt/attempt * config['bwa_mem2meme_aln_sort']['threads'] if attempt < 2 else 254)"),
    benchmark:
        repeat(MDIR + "{sample}/benchmarks/{sample}.bwa2m.alNsort.bench.tsv", 0)
    params:
        cluster_slots=config["bwa_mem2meme_aln_sort"]["threads"], #return_aligner_slot_by_fq_size(MDIR + "{wildcards.sample}/"),
        mdir=MDIR,
        samtmpd="/align/bwa2m/logs/sam_tmpd/",
        bwa_threads=config["bwa_mem2meme_aln_sort"]["bwa_threads"],  # eval( "lambda wildcards, input,threads, attempt: int(attempt/attempt * config['bwa_mem2meme_aln_sort']['other_threads'] if attempt < 2 else 254)"),
        write_threads=config["bwa_mem2meme_aln_sort"]["write_threads"],  # eval( "lambda wildcards, input,threads, attempt: int(attempt/attempt * config['bwa_mem2meme_aln_sort']['write_threads'] if attempt < 2 else 22)"),
        sort_threads=config["bwa_mem2meme_aln_sort"]["sort_threads"],  # eval( "lambda wildcards, input,threads, attempt:  int(attempt/attempt * config['bwa_mem2meme_aln_sort']['sort_threads'] if attempt < 2 else 44)"),
        benchmark_runs=config["bwa_mem2meme_aln_sort"]["softclip_alts"],
        softclip_alts=config["bwa_mem2meme_aln_sort"]["softclip_alts"],
        bwa_mem2m_cmd=config["bwa_mem2meme_aln_sort"]["cmd"],
        k=get_bwa_kmer_size,
        K=config["bwa_mem2meme_aln_sort"]["K"],
        sort_thread_mem=config["bwa_mem2meme_aln_sort"]["sort_thread_mem"],
        huref=config["supporting_files"]["files"]["huref"]["bwa_mem_index_meme"]["name"],
        lib=" "
        if "lib" not in config["bwa_mem2meme_aln_sort"]
        else config["bwa_mem2meme_aln_sort"]["lib"],
        ulim=" ",  #"ulimit -Sn  100000 " if 'ulimit' not in config else config['ulimit'],
        rgpl="presumedILLUMINA",  # ideally: passed in technology # nice to get to this point: https://support.sentieon.com/appnotes/read_groups/  :: note, the     default sample name contains the RU_EX_SQ_Lane (Lane=0 for combined)
        rgpu="presumedCombinedLanes",  # ideally flowcell_lane(s)
        cluster_sample=ret_sample,
        rgsm=ret_sample,  # samplename
        rgid=ret_sample,  # ideally samplename_flowcell_lane(s)_barcode  ! Imp this is unique, I add epoc seconds to the end of start of this rule
        rglb="_presumedNoAmpWGS",  # prepend withs sample namne in shell ideally samplename_libprep
        rgcn="CenterName",  # center name
        rgpg="bwamem2MEME",  #program  ^^^ the Read Group Info is super important for several gatk core tools.  There is a library which will extract this info from the read names if we go full fancy.
        ld_p=" ",
        ldpre=config['bwa_mem2meme_aln_sort']['ldpre'],
        subsample_head=get_subsample_head,
        subsample_tail=get_subsample_tail,
        mbuff_mem=config["bwa_mem2meme_aln_sort"]["mbuffer_mem"],
        igz_threads=config['bwa_mem2meme_aln_sort']['igz_threads']
    conda:
        config["bwa_mem2meme_aln_sort"]["env_yaml"]
    shell: 
        """
	TOKEN=$(curl -X PUT 'http://169.254.169.254/latest/api/token' -H 'X-aws-ec2-metadata-token-ttl-seconds: 21600');
        itype=$(curl -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-type);
        echo "INSTANCE TYPE: $itype" > {log};
        echo "INSTANCE TYPE: $itype"
        start_time=$(date +%s);

        {params.lib}
        export tdir={params.mdir}/{params.cluster_sample}/{params.samtmpd};
        mkdir -p $tdir;
        epocsec=$(date +'%s');
        {params.lib}

        unset LD_PRELOAD;
        
        # mbuffer size should be determined by memory option given to samtools.
        # ex) samtools sort uses 20 threads, 1G per each thread, so mbuffer size should be 20G (= -m 1G x -@ 20)

        {params.ldpre} {params.bwa_mem2m_cmd}  mem  \
           -R '@RG\\tID:{params.rgid}_$epocsec\\tSM:{params.rgsm}\\tLB:{params.cluster_sample}{params.rglb}\\tPL:{params.rgpl}\\tPU:{params.rgpu}\\tCN:{params.rgcn}\\tPG:{params.rgpg}' \
            {params.softclip_alts}  {params.K} {params.k} -t {params.bwa_threads}  \
            -7 {params.huref} \
            {params.subsample_head} <( igzip -c -d -T {params.igz_threads} -q  {input.f1} | sed 's/ /_/')   {params.subsample_tail} \
            {params.subsample_head} <( igzip -c -d -T {params.igz_threads} -q  {input.f2} | sed 's/ /_/')  {params.subsample_tail}  \
            |  samtools sort -l 1  -m {params.sort_thread_mem}   \
            -@  {params.sort_threads} -T $tdir -O BAM --write-index -o {output.bamo}##idx##{output.bami} -  >> {log} 2>&1;
        
            end_time=$(date +%s);
            elapsed_time=$((($end_time - $start_time) / 60));
	    echo "Elapsed-Time-min:\t$itype\t$elapsed_time";
            echo "Elapsed-Time-min:\t$itype\t$elapsed_time" >> {log} 2>&1;
        """
 
