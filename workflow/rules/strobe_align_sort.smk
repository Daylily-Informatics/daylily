####### strobe aligner

rule strobe_align_sort:
    """https://github.com/ksahlin/strobealign"""
    
    input:
        DR=MDIR + "{sample}/{sample}.dirsetup.ready",
        f1=getR1s,
        f2=getR2s,
    output:
        bami=temp(MDIR + "{sample}/align/strobe/{sample}.strobe.sort.bam.bai"),
        bamo=temp(MDIR + "{sample}/align/strobe/{sample}.strobe.sort.bam")
    log:
        MDIR + "{sample}/align/strobe/logs/{sample}.strobe_sort.log",
    resources:
        threads=config['strobe_align_sort']['threads'],
        mem_mb=config['strobe_align_sort']['mem_mb'],
        partition=config['strobe_align_sort']['partition'],
        vcpu=config['strobe_align_sort']['threads']
    threads: config['strobe_align_sort']['threads']
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.strobe.alNsort.bench.tsv"
    params:
        cluster_slots=94,
        cluster_sample=ret_sample,
        mdir=MDIR,
        samtmpd="/align/strobe/logs/sam_tmpd/",
        write_threads=config["strobe_align_sort"]["write_threads"],
        sort_threads= config["strobe_align_sort"]["sort_threads"],
        benchmark_runs=config["strobe_align_sort"]["softclip_alts"],
        softclip_alts=config["strobe_align_sort"]["softclip_alts"],
        strobe_cmd=config["strobe_align_sort"]["cmd"],
        ldpre=config['strobe_align_sort']['ldpre'],
        k=config["strobe_align_sort"]["k"],  # little kay
        K=config["strobe_align_sort"]["K"],  # BIG KAY
        sort_thread_mem=config["strobe_align_sort"]["sort_thread_mem"],
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        rgpl="presumedILLUMINA",  # ideally: passed in technology # nice to get to this point: https://support.sentieon.com/appnotes/read_groups/  :: note, the default sample name contains the RU_EX_SQ_Lane (0 for combined)
        rgpu="presumedCombinedLanes",  # ideally flowcell_lane(s)
        rgsm='x', # ret_sample,  # samplename
        rgid='x', #ret_sample,  # ideally samplename_flowcell_lane(s)_barcode  ! Imp this is unique, I add epoc seconds to the end of start of this rule
        rglb="_presumedNoAmpWGS",  # prepend sample_name in shell block ideally samplename_libprep
        rgcn="CenterName",  # center name
        subsample_head=get_subsample_head,
        subsample_tail=get_subsample_tail,
        strobe_threads=config["strobe_align_sort"]["strobe_threads"],
        samp=get_samp_name,
        mbuff_mem=config["strobe_align_sort"]["mbuffer_mem"],
        rgpg="strobealigner",
        numa=config["strobe_align_sort"]["numa"],
        igz_threads=config['strobe_align_sort']['igz_threads']
    conda:
        config["strobe_align_sort"]["env_yaml"]
    shell:
        """
        TOKEN=$(curl -X PUT 'http://169.254.169.254/latest/api/token' -H 'X-aws-ec2-metadata-token-ttl-seconds: 21600');
        itype=$(curl -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-type);
        echo "INSTANCE TYPE: $itype" > {log};
        echo "INSTANCE TYPE: $itype"
        start_time=$(date +%s);

        export tdir={params.mdir}/{params.samp}/{params.samtmpd};
        mkdir -p $tdir ;
        epocsec=$(date +'%s');
        
        echo 'WARNING! SPACES IN FASTQ READ NAMES ARE REPLACED WITH  \ ' >> {log} 2>&1;

        {params.numa} \
        {params.strobe_cmd} \
        -t {params.strobe_threads} \
        --rg-id="{params.rgid}_$epocsec" \
        --rg="SM:{params.rgsm}" \
        --rg=LB:"{params.samp}{params.rglb}" \
        --rg=PL:"{params.rgpl}" \
        --rg=PU:"{params.rgpu}" \
        --rg=CN:"{params.rgcn}" \
        --rg=PG:"{params.rgpg}" \
        --use-index {params.huref}  \
        {params.subsample_head} <(igzip -c -d -T {params.igz_threads} -q  {input.f1} )  {params.subsample_tail} \
        {params.subsample_head}  <(igzip -c -d -T {params.igz_threads} -q {input.f2} )  {params.subsample_tail} \
        | samtools sort -l 1  -m {params.sort_thread_mem}   \
        -@  {params.sort_threads} -T $tdir -O BAM --write-index -o {output.bamo}##idx##{output.bami} - >> {log} 2>&1;

        end_time=$(date +%s);
        elapsed_time=$((($end_time - $start_time) / 60));
        echo "Elapsed-Time-min:\t$itype\t$elapsed_time";
        echo "Elapsed-Time-min:\t$itype\t$elapsed_time" >> {log} 2>&1;
        """


localrules: produce_strobe_align,

rule produce_strobe_align:  # TARGET: only produce strobe align
     input:
         expand(MDIR + "{sample}/align/strobe/{sample}.strobe.sort.bam", sample=SAMPS)
 
 
