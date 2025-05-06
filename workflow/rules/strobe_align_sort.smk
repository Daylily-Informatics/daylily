import os
####### strobe aligner


if os.environ.get("DAY_STROBE_TOGGLE","") == "":


    rule strobe_align_sort_bam:
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
            vcpu=config['strobe_align_sort']['threads'],
            constraint=config['strobe_align_sort']['constraint'],
        threads: config['strobe_align_sort']['threads']
        benchmark:
            MDIR + "{sample}/benchmarks/{sample}.strobe.alNsort.bench.tsv"
        params:
            cluster_sample=ret_sample,
            sort_threads= config["strobe_align_sort"]["sort_threads"],
            strobe_cmd=config["strobe_align_sort"]["cmd"],
            strobe_opts=config["strobe_align_sort"]["strobe_opts"],  # BIG KAY
            sort_thread_mem=config["strobe_align_sort"]["sort_thread_mem"],
            huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
            rgpl="ILLUMINA",  # ideally: passed in technology # nice to get to this point: https://support.sentieon.com/appnotes/read_groups/  :: note, the default sample name contains the RU_EX_SQ_Lane (0 for combined)
            rgpu="presumedCombinedLanes",  # ideally flowcell_lane(s)
            rgsm=ret_sample, # ret_sample,  # samplename
            rgid=ret_sample, #ret_sample,  # ideally samplename_flowcell_lane(s)_barcode  ! Imp this is unique, I add epoc seconds to the end of start of this rule
            rglb="_presumedNoAmpWGS",  # prepend sample_name in shell block ideally samplename_libprep
            rgcn="CenterName",  # center name
            subsample_head=get_subsample_head,
            subsample_tail=get_subsample_tail,
            strobe_threads=config["strobe_align_sort"]["strobe_threads"],
            samp=get_samp_name,
            mbuffer=config["strobe_align_sort"]["mbuffer"],
            rgpg="strobealigner",
            tofq_threads=config["strobe_align_sort"]["tofq_threads"],
            igz=config['strobe_align_sort']['igz']
        conda:
            config["strobe_align_sort"]["env_yaml"]
        shell:
            """
            TOKEN=$(curl -X PUT 'http://169.254.169.254/latest/api/token' -H 'X-aws-ec2-metadata-token-ttl-seconds: 21600');
            itype=$(curl -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-type);
            echo "INSTANCE TYPE: $itype" > {log};
            start_time=$(date +%s);

            timestamp=$(date +%Y%m%d%H%M%S);
            TMPDIR=/dev/shm/strobe_tmp_$timestamp;
            mkdir -p $TMPDIR;
            APPTAINER_HOME=$TMPDIR;
            trap "rm -rf \"$TMPDIR\" || echo '$TMPDIR rm fails' >> {log} 2>&1" EXIT;
            tdir=$TMPDIR;

            epocsec=$(date +'%s');

            ulimit -n 65536 || echo "ulimit mod failed";


            {params.strobe_cmd} \
            -t {params.strobe_threads} {params.strobe_opts} \
            --rg-id="{params.cluster_sample}-$epocsec" \
            --rg="SM:{params.cluster_sample}" \
            --rg=LB:"{params.cluster_sample}-LB-1" \
            --rg=PL:"{params.rgpl}" \
            --rg=PU:"{params.rgpu}" \
            --rg=CN:"{params.rgcn}" \
            --rg=PG:"{params.rgpg}" \
            --use-index {params.huref}  \
            {params.subsample_head} <(  {params.igz} -q  {input.f1} )  {params.subsample_tail} \
            {params.subsample_head}  <( {params.igz} -q {input.f2} )  {params.subsample_tail}  {params.mbuffer} \
            |  samtools sort \
            -l 1  \
            -m {params.sort_thread_mem}   \
            -@  {params.sort_threads} \
            -T $tdir \
            -O BAM \
            --write-index \
            -o {output.bamo}##idx##{output.bami} - >> {log} 2>&1;

            end_time=$(date +%s);
            elapsed_time=$((($end_time - $start_time) / 60));
            echo "Elapsed-Time-min:\t$itype\t$elapsed_time" >> {log} 2>&1;
            """

else:
        
    rule strobe_align_sort_cram:
        """https://github.com/ksahlin/strobealign"""
        
        input:
            cram=MDIR+"{sample}/align/ug/{sample}.ug.cram",
            crai=MDIR+"{sample}/align/ug/{sample}.ug.cram.crai",
        output:
            bami=temp(MDIR + "{sample}/align/strobe/{sample}.strobe.sort.bam.bai"),
            bamo=temp(MDIR + "{sample}/align/strobe/{sample}.strobe.sort.bam")
        log:
            MDIR + "{sample}/align/strobe/logs/{sample}.strobe_sort.log",
        resources:
            threads=config['strobe_align_sort']['threads'],
            mem_mb=config['strobe_align_sort']['mem_mb'],
            partition=config['strobe_align_sort']['partition'],
            vcpu=config['strobe_align_sort']['threads'],
            constraint=config['strobe_align_sort']['constraint'],
        threads: config['strobe_align_sort']['threads']
        benchmark:
            MDIR + "{sample}/benchmarks/{sample}.strobe.alNsort.bench.tsv"
        params:
            cluster_sample=ret_sample,
            sort_threads= config["strobe_align_sort"]["sort_threads"],
            strobe_cmd=config["strobe_align_sort"]["cmd"],
            strobe_opts=config["strobe_align_sort"]["strobe_opts"],  # BIG KAY
            sort_thread_mem=config["strobe_align_sort"]["sort_thread_mem"],
            huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
            rgpl="ILLUMINA",  # ideally: passed in technology # nice to get to this point: https://support.sentieon.com/appnotes/read_groups/  :: note, the default sample name contains the RU_EX_SQ_Lane (0 for combined)
            rgpu="presumedCombinedLanes",  # ideally flowcell_lane(s)
            rgsm=ret_sample, # ret_sample,  # samplename
            rgid=ret_sample, #ret_sample,  # ideally samplename_flowcell_lane(s)_barcode  ! Imp this is unique, I add epoc seconds to the end of start of this rule
            rglb="_presumedNoAmpWGS",  # prepend sample_name in shell block ideally samplename_libprep
            rgcn="CenterName",  # center name
            subsample_head=get_subsample_head,
            subsample_tail=get_subsample_tail,
            strobe_threads=config["strobe_align_sort"]["strobe_threads"],
            samp=get_samp_name,
            mbuffer=config["strobe_align_sort"]["mbuffer"],
            rgpg="strobealigner",
            tofq_threads=config["strobe_align_sort"]["tofq_threads"],
            igz=config['strobe_align_sort']['igz']
        conda:
            config["strobe_align_sort"]["env_yaml"]
        shell:
            """
            TOKEN=$(curl -X PUT 'http://169.254.169.254/latest/api/token' -H 'X-aws-ec2-metadata-token-ttl-seconds: 21600');
            itype=$(curl -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-type);
            echo "INSTANCE TYPE: $itype" > {log};
            start_time=$(date +%s);

            timestamp=$(date +%Y%m%d%H%M%S);
            TMPDIR=/fsx/scratch/strobe_tmp_$timestamp;
            mkdir -p $TMPDIR;
            APPTAINER_HOME=$TMPDIR;
            trap "rm -rf \"$TMPDIR\" || echo '$TMPDIR rm fails' >> {log} 2>&1" EXIT;
            tdir=$TMPDIR;

            epocsec=$(date +'%s');

            ulimit -n 65536 || echo "ulimit mod failed";


            samtools fastq -n -T {params.huref} -@ {params.sort_threads}  {input.cram} \
            {params.mbuffer} | {params.strobe_cmd} \
            -t {params.strobe_threads} {params.strobe_opts} \
            --rg-id="{params.cluster_sample}-$epocsec" \
            --rg="SM:{params.cluster_sample}" \
            --rg=LB:"{params.cluster_sample}-LB-1" \
            --rg=PL:"{params.rgpl}" \
            --rg=PU:"{params.rgpu}" \
            --rg=CN:"{params.rgcn}" \
            --rg=PG:"{params.rgpg}" \
            --use-index {params.huref}  {params.mbuffer} \
            |  samtools sort \
            -l 1  \
            -m {params.sort_thread_mem}   \
            -@  {params.sort_threads} \
            -T $tdir \
            -O BAM \
            --write-index \
            -o {output.bamo}##idx##{output.bami} - >> {log} 2>&1;

            end_time=$(date +%s);
            elapsed_time=$((($end_time - $start_time) / 60));
            echo "Elapsed-Time-min:\t$itype\t$elapsed_time" >> {log} 2>&1;
            """


localrules: produce_strobe_align_sort_bam,

rule produce_strobe_align_sort_bam:  # TARGET: strobe aligner sorted bam
     input:
         expand(MDIR + "{sample}/align/strobe/{sample}.strobe.sort.bam", sample=SAMPS)
 
 
