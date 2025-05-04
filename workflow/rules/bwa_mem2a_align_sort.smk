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
        vcpu=config['bwa_mem2a_aln_sort']['threads'],
        constraint=config['bwa_mem2a_aln_sort']['constraint'],
    threads: config['bwa_mem2a_aln_sort']['threads']
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.bwa2a.alNsort.bench.tsv"
    params:
        cluster_sample=ret_sample,
        sort_threads= config["bwa_mem2a_aln_sort"]["sort_threads"],
        bwa_mem2a_cmd=config["bwa_mem2a_aln_sort"]["cmd"],
        bwa_opts=config["bwa_mem2a_aln_sort"]["bwa_opts"],
        sort_thread_mem=config["bwa_mem2a_aln_sort"]["sort_thread_mem"],
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        rgpl="ILLUMINA",  # ideally: passed in technology # nice to get to this point: https://support.sentieon.com/appnotes/read_groups/  :: note, the default sample name contains the RU_EX_SQ_Lane (0 for combined)
        rgpu="presumedCombinedLanes",  # ideally flowcell_lane(s)
        rgsm=ret_sample, # ret_sample,  # samplename
        rgid=ret_sample, #ret_sample,  # ideally samplename_flowcell_lane(s)_barcode  ! Imp this is unique, I add epoc seconds to the end of start of this rule
        rglb="_presumedNoAmpWGS",  # prepend sample_name in shell block ideally samplename_libprep
        rgcn="CenterName",  # center name
        rgpg="bwamem2",  #program
        subsample_head=get_subsample_head,
        subsample_tail=get_subsample_tail,
        bwa_threads=config["bwa_mem2a_aln_sort"]["bwa_threads"],
        samp=get_samp_name,
        mbuffer=config["bwa_mem2a_aln_sort"]["mbuffer"],
        igz=config['bwa_mem2a_aln_sort']['igz']
    conda:
        config["bwa_mem2a_aln_sort"]["env_yaml"]
    shell:
        """

        TOKEN=$(curl -X PUT 'http://169.254.169.254/latest/api/token' -H 'X-aws-ec2-metadata-token-ttl-seconds: 21600');
        itype=$(curl -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-type);
        echo "INSTANCE TYPE: $itype" > {log};
        start_time=$(date +%s);
        ulimit -n 65536 || echo "ulimit mod failed";


        timestamp=$(date +%Y%m%d%H%M%S);
        TMPDIR=/dev/shm/bwa2a_tmp_$timestamp;
        mkdir -p $TMPDIR;
        APPTAINER_HOME=$TMPDIR;
        trap "rm -rf \"$TMPDIR\" || echo '$TMPDIR rm fails' >> {log} 2>&1" EXIT;

        tdir=$TMPDIR;

        epocsec=$(date +'%s');
        

        {params.bwa_mem2a_cmd} mem \
        -R '@RG\\tID:{params.cluster_sample}-$epocsec\\tSM:{params.cluster_sample}\\tLB:{params.cluster_sample}-LB-1\\tPL:{params.rgpl}\\tPU:{params.rgpu}\\tCN:{params.rgcn}\\tPG:{params.rgpg}' \
        {params.bwa_opts}  -t {params.bwa_threads}  \
        {params.huref} \
        {params.subsample_head} <( {params.igz} -q  {input.f1} )  {params.subsample_tail}  \
        {params.subsample_head} <( {params.igz} -q  {input.f2} )  {params.subsample_tail} {params.mbuffer} \
        |  samtools sort \
        -l 1  \
        -m {params.sort_thread_mem}   \
        -@  {params.sort_threads} \
        -T $tdir \
        -O BAM  \
        --write-index \
        -o {output.bamo}##idx##{output.bami} >> {log} 2>&1;

        end_time=$(date +%s);
    	elapsed_time=$((($end_time - $start_time) / 60));
        echo "Elapsed-Time-min:\t$itype\t$elapsed_time" >> {log} 2>&1;

        """


localrules: produce_bwa_mem2_sort_bam,

rule produce_bwa_mem2_sort_bam:  # TARGET: produce_bwa_mem2_sort_bam
     input:
         expand(MDIR + "{sample}/align/bwa2a/{sample}.bwa2a.sort.bam", sample=SAMPS)
