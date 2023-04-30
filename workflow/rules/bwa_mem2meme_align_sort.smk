####### BWA MEM2c - meme m-- Version radix-tree
#
# An accelerated version of BWA MEM
# bwa-mem being the defactor standard mapper/aligner
#
# - bwa mem2 orig : https://github.com/bwa-mem2/bwa-mem2
# - 2meme paper : https://www.biorxiv.org/content/10.1101/2021.09.01.457579v1.full
# - 2meme github: https://github.com/kaist-ina/BWA-MEME/


rule bwa_mem2c_meme_aln_sort:
    # bwamem2 branch using ML to improve seeding speeds
    # claims 33% speed increase
    input:
        DR=MDIR + "{sample}/{sample}.dirsetup.ready",
        f1=getR1s,  # method defined in fastp.smk
        f2=getR2s,  # method defined in fastp.smk
    output:
        bamo=MDIR + "{sample}/align/bwa2c/{sample}.bwa2c.sort.bam",
        bami=MDIR + "{sample}/align/bwa2c/{sample}.bwa2c.sort.bam.bai",
        samo=MDIR + "{sample}/align/bwa2c/{sample}.bwa2c.sam",
    log:
        MDIR + "{sample}/align/bwa2c/logs/align_sort/{sample}.bwa2c.sort.log",
    threads: config["bwa_mem2c_aln_sort"]["threads"]  #eval( "lambda wildcards: attempt, input,threads, attempt: int( int(attempt/attempt * config['bwa_mem2c_aln_sort']['threads'] if attempt < 2 else 254))"),
    priority: 5
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.bwa2c.alNsort.bench.tsv"
    params:
        cluster_slots=254,  #return_aligner_slot_by_fq_size(MDIR + "{sample}/"),
        mdir=MDIR,
        samtmpd="/align/bwa2c/logs/sam_tmpd/",
        other_threads=config["bwa_mem2c_aln_sort"]["other_threads"],
        write_threads=config["bwa_mem2c_aln_sort"]["write_threads"],
        sort_threads=config["bwa_mem2c_aln_sort"]["sort_threads"],
        benchmark_runs=config["bwa_mem2c_aln_sort"]["softclip_alts"],
        softclip_alts=config["bwa_mem2c_aln_sort"]["softclip_alts"],
        bwa_mem2c_cmd=config["bwa_mem2c_aln_sort"]["cmd"],
        k=config["bwa_mem2c_aln_sort"]["k"],
        K=config["bwa_mem2c_aln_sort"]["K"],
        thread_mem=config["bwa_mem2c_aln_sort"]["thread_mem"],
        bwv=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        lib=" "
        if "lib" not in config["bwa_mem2c_aln_sort"]
        else config["bwa_mem2c_aln_sort"]["lib"],
        ulim="ulimit -Sn  100000 " if "ulimit" not in config else config["ulimit"],
        cluster_sample=ret_sample,
        ld_preload=config["malloc_alt"]["ld_preload"],
        ld_pre=config["malloc_alt"]["ld_preload"],
        rgpl="presumedILLUMINA",  # ideally: passed in technology # nice to get to this point: https://support.sentieon.com/appnotes/read_groups/  :: note, the    default sample name contains the RU_EX_SQ_Lane (0 for combined)
        rgpu="presumedCombinedLanes",  # ideally flowcell_lane(s)
        rgsm=ret_sample,  # samplename
        rgid=ret_sample,  # ideally samplename_flowcell_lane(s)_barcode  ! Imp this is unique, I add epoc seconds to the end of start of this rule
        rglb="_presumedNoAmpWGS",  # prepend with cluster name in shell ideally samplename_libprep
        rgcn="CenterName",  # center name
        rgpg="bwamem2meme",  #program
        piped_align=config["bwa_mem2c_aln_sort"]['piped_align'],
    container:
        None
    conda:
        config["bwa_mem2c_aln_sort"]["env_yaml"]
    shell:
        """
        export tdir={params.mdir}{params.cluster_sample}/{params.samtmpd};
        mkdir -p $tdir;
        epocsec=$(date +'%s');
        if [[ "{params.piped_align}" == "yes" ]]; then
            {params.numactl} {params.ld_preload} {params.bwa_mem2c_cmd} mem -7 \
            {params.softclip_alts} {params.K} {params.k} \
            -R '@RG\\tID:{params.rgid}_$epocsec\\tSM:{params.rgsm}\\tLB:{params.cluster_sample}{params.rglb}\\tPL:{params.rgpl}\\tPU:{params.rgpu}\\tCN:{params.rgcn}\\tPG:{params.rgpg}' \
            -t {params.other_threads} {params.bwv}    \
            <(unpigz -c -q -- {input.f1} )  \
            <(unpigz -c  -q -- {input.f2} ) \
            | {params.ld_preload} samtools sort -l 9  -m {params.thread_mem} -@  {params.sort_threads} -T $tdir -O BAM  - \
            | {params.ld_preload} samtools view -b -1  -h -@ {params.write_threads} -O BAM --write-index -o {output.bamo}##idx##{output.bami} - ;
        else
            {params.numactl} {params.ld_preload} {params.bwa_mem2c_cmd} mem -7 \
            {params.softclip_alts} -o  {output.samo}  {params.K} {params.k} \
            -R '@RG\\tID:{params.rgid}_$epocsec\\tSM:{params.rgsm}\\tLB:{params.cluster_sample}{params.rglb}\\tPL:{params.rgpl}\\tPU:{params.rgpu}\\tCN:{params.rgcn}\\tPG:{params.rgpg}' \
            -t {params.other_threads} {params.bwv}    \
            <(unpigz -c -q -- {input.f1} )  \
            <(unpigz -c  -q -- {input.f2} ) >> {log} ;
            {params.ld_preload} samtools sort -l 9  -m {params.thread_mem} -@  {params.sort_threads} -T $tdir -O BAM  {output.samo} \
            | {params.ld_preload} samtools view -b   -h -@ {params.write_threads} -O BAM --write-index -o {output.bamo}##idx##{output.bami} - >> {log};
        fi;


        echo "File deleted after completing or writing bam and bai. But File this needs to exist for snakemake to be happy.... I find the snkemake temp() to be buggy, so I'm leaving this message instead of a huge sam file :-)" > {output.samo};
        {latency_wait}; ls -lt {output};

        """
