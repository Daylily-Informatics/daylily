####### strobe aligner

rule strobe_align_sort:
    """https://github.com/ksahlin/strobealign"""
    
    input:
        DR=MDIR + "{sample}/{sample}.dirsetup.ready",
        f1=getR1s,
        f2=getR2s,
    output:
        bami=MDIR + "{sample}/align/strobe/{sample}.strobe.sort.bam.bai",
        bamo=MDIR + "{sample}/align/strobe/{sample}.strobe.sort.bam",
    log:
        MDIR + "{sample}/align/strobe/logs/{sample}.strobe_sort.log",
    resources:
        threads=config['strobe_align_sort']['threads'],
        mem_mb=60000 if 'mem_mb' not in config['strobe_align_sort'] else config['strobe_align_sort']['mem_mb'],
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
	numa=config["strobe_align_sort"]["mbuffer_mem"]
    conda:
        config["strobe_align_sort"]["env_yaml"]
    shell:
        """
        export tdir={params.mdir}/{params.samp}/{params.samtmpd};
        mkdir -p $tdir ;

        epocsec=$(date +'%s');
        
        
        {params.numa} {params.strobe_cmd} -v \
         --rg '@RG\\tID:{params.rgid}_$epocsec\\tSM:{params.rgsm}\\tLB:{params.samp}{params.rglb}\\tPL:{params.rgpl}\\tPU:{params.rgpu}\\tCN:{params.rgcn}\\tPG:{params.rgpg}' \
          -t {params.strobe_threads}  \
	  --use-index {params.huref} \
         {params.subsample_head} <(unpigz -t 4 -c  -q -- {input.f1} )  {params.subsample_tail}  \
         {params.subsample_head} <(unpigz -t 4 -c  -q -- {input.f2} )  {params.subsample_tail}    \
        |   samtools sort -l 0  -m {params.sort_thread_mem}   \
         -@  {params.sort_threads} -T $tdir  -O SAM - \
        |  samtools view -b -@ {params.write_threads} -O BAM --write-index -o {output.bamo}##idx##{output.bami} -  >> {log};
        """


localrules: produce_strobe_align,

rule produce_strobe_align:  # TARGET: only produce strobe align
     input:
         expand(MDIR + "{sample}/align/strobe/{sample}.strobe.sort.bam", sample=SAMPS)
	 
