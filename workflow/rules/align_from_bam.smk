import os

# #### RULES and METHODS for normalizing the inputs DAY accepts.
# This includes: locally stored: fasta, fastq, single sample BAM(easily extended to multisample BAMS),
#  and BCLs (which have a rule 'run_bcl2fastq' which produces the samplesheet and fastqs to begin from here).
#  Remotely stored, mod can work with BAMS (and with a few tweaks to the BAM rule, remote fastq/fasta)-- these
#  can be stored anyplace rclone can access, whhich stands presently at 14 S3 based providers, and a total of 28 different
#  non-S3 storage systems or protocols, including arbitrary unionFS views of these - most useful likely to be:
#  (local drive, S3 , HTTP, S/FTP, WebDAV, Dropbox, Gdrive, Box, Dropbox).  This requires preliminary configuration of each
# remote source via rclone, but then can be reused as long as the credntials are valid.
# ---------



# dont want to fix this unecessary if quite yet
if True:


    def get_rclone_handle(wildcards):
        rc_handle =samples.loc[(wildcards.sample), "rclone_handle"][0]
        return rc_handle


    # Special case, if the sample sheet contains remote data sources, go fetch the files.  Possibly even align while
    #   receiving.
    # NOTE-  rclone supports the definition of a 'remote' data source which is the local file system, so any local data, ie:
    #   BAM files, may be processed via this method now.
    # Analysis manifests should only contain all of 'remote bams'(local bams are processed by treating them as remote),
    #   'remote fastq', or 'local fastq'. The implementation anticipates, and should work with mixed types, but this has not
    #   been tested yet, so consider it experimental.
    if "remote_bam_file_path" in samples.columns:
        BAM_READ_GROUPS=[]

        if 'rclone_conf_file' not in config:
            raise Exception("\n\n\tERROR: You must specify a rclone conf file via --config rclone_conf_file=/path/to/rclone.conf, if you wush to use a manifest which has BAM files to be fetched via rclone (even if using a remote backend handle locally)'\n\n")

        else:
            # Will be doing this redundantly, but that is fine for now.
            pass


        localrules: build_bam_header_files

        rule build_bam_header_files:
            # Store the distinct header info from each input bam, wich may be combined
            # from severarl upstream BAMs and have multiple entries.  This info should
            # be stored with the reads individually as the sam header/readgroup handling
            # by all tools, except GATK(but only very recently) is ... inconsistent
            # and worse, edits, removes, remixes. Happily, I believe of the 40 tools used
            # in DAY, a total of 1 has any interest in the read group info, and that
            # is picard-  which is on the short list to be cut for reliability and performance
            # reasons largely.
            input:
                MDIR + "{sample}/{sample}.dirsetup.ready",
            output:
                MDIR + "{sample}/align/{sample}.bamheader.txt"
            run:
                # I do not advise use of 'run' unless absolutely necessary.  It just does not play nicely with most advanced snakemake features.   I would move this to a script and call from a shell block if making an changes here in the future.
                for k in samples.iterrows():
                    rclone_handle = k[1][17] + k[1][16]
                    resp = os.popen(f"rclone --config {config['rclone_conf_file']} -L  cat {rclone_handle} | samtools view -H | head -n 10000 | grep '@RG'")
                    outputf=f"{output[0]}"
                    print(outputf,file=sys.stderr)
                    os.system(f"echo yyy; rm -rf {output[0]} || sleep .1 ; mkdir -p $(dirname {output[0]} ) || sleep .1 ;")
                    print(outputf,file=sys.stderr)
                    of=open(outputf,'w')
                    for ii in resp:
                        iii = ii.rstrip()
                        of.write(iii)

                of.close()
                os._exit(0)


        def get_remote_bam(wildcards):
            bam = samples.loc[(wildcards.sample), "remote_bam_file_path"][0]
            return bam

        def get_rclone_handle(wildcards):
            rc_handle =samples.loc[(wildcards.sample), "rclone_handle"][0]
            return rc_handle


        # In the interest of time, I am just making copies of the paired reads from the
        # orig BAM, I would weave in using them directly with process substitution, but
        # am to short on time presently
        rule produce_fastqs_from_bams:
            input:
                MDIR + "{sample}/align/{sample}.bamheader.txt",
            output:
                pr1=MDIR + "{sample}/{sample_lane}.R1.fastq.gz",
                pr2=MDIR + "{sample}/{sample_lane}.R2.fastq.gz",
                ur1=MDIR + "{sample}/{sample_lane}.unpair.R1.fastq.gz",
                ur2=MDIR + "{sample}/{sample_lane}.unpair.R2.fastq.gz",
                singletons=MDIR + "{sample}/{sample_lane}.singletons.fastq.gz",
            conda:
                "../envs/biobambam2_v0.1.yaml" if 'conda_env' not in config['produce_fastqs_from_bams'] else config['produce_fastqs_from_bams']['conda_env']
            benchmark:
                MDIR + "{sample}/benchmarks/{sample_lane}.bench.tsv"
            params:
                cluster_sample=ret_sample,
                rc_config=config['rclone_conf_file'],
                rclone_handle=get_rclone_handle,
                rclone_remote_bam=get_remote_bam,
                samtools_view_threads=config['produce_fastqs_from_bams']['samtools_view_threads']
            threads: config['produce_fastqs_from_bams']['threads']
            log:
                MDIR+"{sample}/logs/{sample_lane}.biobambam2_to_fastq.log"
            shell:
                """
                ( rm -rf {output} || echo okToProceed ;
                touch {output.ur1} {output.ur2} {output.singletons};   # because they may legitimately not be present ### export LD_LIBRARY_PATH=$PWD/resources/ && \
                rclone -L --use-mmap --buffer-size=96K --transfers 8 --checkers 8 --config {params.rc_config} cat {params.rclone_handle}{params.rclone_remote_bam} | samtools view \
                --output-fmt=SAM -h -@ {params.samtools_view_threads} - | \
                bamtofastq  \
                collate=1   S={output.singletons}  O={output.ur1}  02={output.ur2} F={output.pr1} F2={output.pr2} \
                combs=1 gz=1 level=2 inputformat=sam inputbuffersize=192M \
                inputbuffersize=10000000000 ;
                echo ok;
                {latency_wait};
                ls {output}; ) > {log};
                """

        localrules: produce_fastqs_from_all_bams

        rule produce_fastqs_from_all_bams:  # TARGET just produce bam fastq files
            input:
                ret_all_local_R1_R2_lane_fqs


        #  Has not been fully tested, finish when more cases arise && do a few more benchmarks
        #  to establish there is a real benefit, not just cool factor.
        #  It is just too big, would need to really be a time saver-  off for the time being anyway.
        if 1 in ['2']:
            rule direct_bwamem2a:
                input:
                    MDIR + "{sample}/align/{sample}.bamheader.txt",
                    bamo=MDIR+"{sample}/align/bwa2af/{sample}.bwa2af.sort.bam",
                    bami=MDIR + "{sample}/align/bwa2af/{sample}.bwa2af.sort.bam.bai",
                    gzur1=MDIR + "{sample}/{sample}.unpair.R1.fastq.gz",
                    gzur2=MDIR + "{sample}/{sample}.unpair.R2.fastq.gz",
                    gzsingletons=MDIR + "{sample}/{sample}.singletons.fastq.gz",
                    gzinterleaved=MDIR + "{sample}/{sample}.interleavedpairs.fastq.gz",
                    ur1=MDIR + "{sample}/{sample}.unpair.R1.fastq",
                    ur2=MDIR + "{sample}/{sample}.unpair.R2.fastq",
                    singletons=MDIR + "{sample}/{sample}.singletons.fastq",
                    interleaved=MDIR + "{sample}/{sample}.interleavedpairs.fastq"
                log:
                    MDIR + "{sample}/align/bwa2af/logs/{sample}.bwa2af_sort.log",
                threads: config["bwa_mem2a_aln_sort"]["threads"]
                benchmark:
                    MDIR + "{sample}/benchmarks/{sample}.bwa2af.sort.bench.tsv"
                params:
                    mdir=MDIR,
                    rc_config=config['rclone_conf_file'],
                    rclone_handle=get_rclone_handle,
                    rclone_remote_bam=get_remote_bam,
                    samtmpd="/align/bwa2af/logs/sam_tmpd/",
                    other_threads=config["bwa_mem2a_aln_sort"]["other_threads"],
                    write_threads=config["bwa_mem2a_aln_sort"]["write_threads"],
                    sort_threads=config["bwa_mem2a_aln_sort"]["sort_threads"],
                    benchmark_runs=config["bwa_mem2a_aln_sort"]["softclip_alts"],
                    softclip_alts=config["bwa_mem2a_aln_sort"]["softclip_alts"],
                    bwa_mem2a_cmd=config["bwa_mem2a_aln_sort"]["cmd"],
                    k=config["bwa_mem2a_aln_sort"]["k"],
                    K=config["bwa_mem2a_aln_sort"]["K"],
                    thread_mem=config["bwa_mem2a_aln_sort"]["thread_mem"],
                    bwv=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
                    ri=ri,
                    ulim=" ",
                    cluster_sample=ret_sample,
                    rgpl="notset", # ideally: passed in technology # nice to get to this point: https://support.sentieon.com/appnotes/read_groups/ :: note, the default sample name contains the RU_EX_SQ_Lane (0 for combined)
                    rgpu="presumedCombinedLanes", # ideally flowcell_lane(s)
                    rgsm=ret_sample,  # samplename
                    rgid=ret_sample,  # ideally samplename_flowcell_lane(s)_barcode  ! Imp this is unique, I add epoc seconds to the end of start of this rule
                    rglb="_notset",  # prepend sample_name in shell block ideally samplename_libprep
                    rgcn="External",  # center name
                    rgpg="bwamem2",  #program
                    #subsample_head=get_subsample_head,   # Highly experimental feature
                    #subsample_tail=get_subsample_tail,
                conda:
                    "../envs/biobambam2_v0.2.yaml"
                shell:
                    """
                    ( rm -rf {output} || echo okToProceed;
                    export tdir={params.mdir}/{params.cluster_sample}/{params.samtmpd};
                    mkdir -p $tdir ;
                    echo "STARTING NOW ~~~" 2>&1;
                    touch {output.ur1} {output.ur2} {output.singletons};  # because they may legitimately not be present
                    ### export LD_LIBRARY_PATH=$PWD/resources/ && \
                     {params.bwa_mem2a_cmd} mem \
                    -R '@RG\\tID:{params.rgid}\\tSM:{params.rgsm}\\tLB:{params.rglb}\\tPL:{params.rgpl\\tPU:{params.rgpu}' \
                    -a {params.softclip_alts}  {params.K} {params.k} \
                    -t {params.other_threads} -p  {params.bwv} \
                    <(rclone -L --use-mmap --buffer-size 32M --bwlimit 2G --config {params.rc_config} cat {params.rclone_handle}{params.rclone_remote_bam} \
                    | samtools view --output-fmt=SAM -h -@ 1 - \
                    | bamtofastq collate=1 S={output.singletons}  \
                    O={output.ur1}  \
                    02={output.ur2} \
                    combs=1 \
                    gz=0 \
                    inputformat=sam \
                    inputbuffersize=192M - \
                    | tee {output.interleavedfq} ) \
                    | samtools sort -l 0  -m {params.thread_mem} \
                    -@  {params.sort_threads} \
                    -T $tdir -O SAM  - \
                    | samtools view -b -1  \
                    -h -@ {params.write_threads} \
                    -O BAM --write-index \
                    -o {output.bamo}##idx##{output.bami} - ;
                    {latency_wait};
                    (pigz -5 {output.singletons}; touch {output.singletons}:   ") & # these can go off to the background to run locally in parallel, with the wait command holding progress until they all return.
                    (pigz -5 {output.ur1}; touch {output.ur1};   ") &
                    (pigz -5 {output.ur2}; touch {output.ur2}:   ") &
                    (pigz -5 {output.interleaved}; touch {output.interleaved}; ") &
                    wait;
                    {latency_wait}; ) > {log} 2>&1
                    """

                    # The manual command I first experimented with
                    # What it looks like to align/sort & dedup a file from S3 as it transfers.bwa-mem2 mem  -Y -K 150000000000 -k 19 -t 254 -p  /l/data/external_data/research_experiments/PHYTO__RESOURCES/genomic_data/organism_refrences/human/human_g1k_v37_modified.fasta/bwa-mem2vanilla/human_g1k_v37_modified.fasta <(rclone cat rcloneHandle:/l/home/jmajor/public_html/EXPERIMENTS/20220302/day_experiments/main/PMGRC-108-108-0.bam | samtools view --input-fmt=BAM --output-fmt=BAM -u -h -@ 5 | bamtofastq  S=s O=o 02=o2   ) | bamsort level=1 SO=coordinate blockmb=60000 index=1 indexfilename=testO.bam.bai inputformat=sam outputformat=bam inputthreads=5 outputthreads=5 O=testO.bam fixmates=1 calmdn=1 calmdnmreference=/l/data/external_data/research_experiments/PHYTO__RESOURCES/genomic_data/organism_refrences/human/human_g1k_v37_modified.fasta/bwa-mem2vanilla/human_g1k_v37_modified.fasta calmdnmrecompindetonly=1 adddupmarksupport=1 markduplicates=1 sortthreads=250

    elif 'remote_fastq' in config:
        # These are rclone described fastqs.  Here, we simply fetch them and make copies to where we usually link to local fastqs

        localrules: generate_remote_fastq_dir_sentinels

        rule generate_remote_fastq_dir_sentinels:
            input:
                MDIR + "{sample}/{sample}.dirsetup.ready",
            output:
                MDIR + "{sample}/align/{sample}.remote.fq.ready"
            log:
                MDIR + "{sample}/align/logs/{sample}.remote.fq.ready.log"
            shell:
                "touch {output} > {log} 2>&1;"


        def get_remote_r1_fq_path(wildcards):
            r1 = samples.loc[(wildcards.sample), "remote_r1_fq_file_path"][0]
            return r1

        def get_remote_r2_fq_path(wildcards):
            r2 = samples.loc[(wildcards.sample), "remote_r2_fq_file_path"][0]
            return r2

        rule pre_prep_remote_raw_fq:
            input:
                MDIR + "{sample}/align/{sample}.remote.fq.ready"
            output:
                or1=MDIR + "{sample}/{sample_lane}.R1.fastq.gz",
                or2=MDIR + "{sample}/{sample_lane}.R2.fastq.gz",
            conda:
                "../envs/biobambam2_v0.1.yaml"
            params:
                cluster_sample=ret_sample,
                rc_config=config['rclone_conf_file'],
                rclone_handle=get_rclone_handle,
                rclone_remote_r1_fq_path=get_remote_r1_fq_path,
                rclone_remote_r2_fq_path=get_remote_r2_fq_path,
            shell:
                """(
                rclone -L --use-mmap  --config {params.rc_config} copyto {params.rclone_handle}{params.rclone_remote_r1_fq_path} {output.or1};
                rclone -L --use-mmap  --config {params.rc_config} copyto {params.rclone_handle}{params.rclone_remote_r2_fq_path} {output.or2};
                ) > {log};
                """


    else:
        # LOCAL FASTQS
        # Fetch and stage our input data, but only as links. deal with fastqs seperately until BAM creation

        def get_fqs(wildcards):
            r1=os.path.abspath(samples[samples['sample_lane'] == wildcards.sample]['r1_path'][0]) #os.path.abspath(samples.loc[(wildcards.sample, wildcards.sample_lane), "r1_path"])
            r2=os.path.abspath(samples[samples['sample_lane'] == wildcards.sample]['r2_path'][0])
            #r2=os.path.abspath(samples.loc[(wildcards.sample, wildcards.sample_lane), "r2_path"])
            return [r1, r2]

        localrules:
            pre_prep_raw_fq,

        rule pre_prep_raw_fq:
            input:
                get_fqs,
            output:
                or1=MDIR + "{sample}/{sample_lane}.R1.fastq.gz",
                or2=MDIR + "{sample}/{sample_lane}.R2.fastq.gz",
            params:
                c=config["prep_input_sample_files"]["source_read_method"],
            shell:
                "{params.c} {input[0]} {output.or1};"
                "{params.c} {input[1]} {output.or2};"

    localrules: prep_inputs,

    rule prep_inputs:  # TARGET: Just Pre
        input:
            ret_all_local_R1_R2_lane_fqs
        output:
            "staged",
        shell:
            "touch  {output}"
