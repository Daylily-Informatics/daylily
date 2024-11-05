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
        baio=MDIR + "{sample}/align/sent/{sample}.sent.sort.bam.bai",
        samo=temp(MDIR + "{sample}/align/sent/{sample}.sent.sort.sam"),
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
        igz_threads=config['sentieon']['igz_threads']
    conda:
        config["sentieon"]["env_yaml"]
    shell:
        """
        export bwt_max_mem={params.max_mem} ;
        epocsec=$(date +'%s');

        ulimit -n 16384
        
        touch {output.samo};
        tdir="/fsx/scratch/";

        /fsx/data/cached_envs/sentieon-genomics-202308.03/bin/sentieon bwa mem \
        -t {threads}  {params.K}  \
        -R '@RG\\tID:{params.rgid}_$epocsec\\tSM:{params.rgsm}\\tLB:{params.cluster_sample}{params.rglb}\\tPL:{params.rgpl}\\tPU:{params.rgpu}\\tCN:{params.rgcn}\\tPG:{params.rgpg}' \
        {params.K} -t {threads} {params.huref} \
         <(igzip -c -d -T  {params.igz_threads} -q  {input.f1} ) \
         <(igzip -c -d -T  {params.igz_threads} -q  {input.f2} )  \
        |  /fsx/data/cached_envs/sentieon-genomics-202308.03/bin/sentieon util sort \
        --thread_count {threads} \
        --sortblock_thread_count {params.sort_threads} \
        --bam_compression 9 \
        --intermediate_compress_level 9  \
        --block_size {params.sort_thread_mem}   \
        --sam2bam \
        -o {output.bamo} - >> {log.a} ;

        samtools index -b -@ {threads} {output.bamo};  > {log} 2>&1;

        """
