####### BWA MEM2 -- 2 versions!
# Alert: Process Substitution Ahead: https://en.wikipedia.org/wiki/Process_substitution
#
# An accelerated version of BWA MEM
# bwa-mem being the defactor standard mapper/aligner
#
# OG BWAMEM: https://github.com/lh3/bwa
#    A bit more: https://informatics.fas.harvard.edu/short-introduction-to-bwa.html
#       bwa requires large pre-computed indexes for the reference genome it works with
#        these indexes have been pre-generated and are part of the suporting data.
#

# if a subsample_pct column is present in the sample sheet, return the back end of the process substitution
# whcih will do the subsampling
def get_subsample_head_tail(sample_id):
    ss_head = ""
    ss_tail = ""
    raisable = False
    try:
        ss_pct = samples.loc[(sample_id), "subsample_pct"][0]
        if ss_pct in ["", "na", 0, "0", None, "None"]:
            pass  # no subsampling requested
        else:
            ss_pct_float = 1000.1
            try:
                ss_pct_float = float(ss_pct)
            except Exception as e:
                raisable = True
                raise (e)

            if ss_pct_float > 1.0 or ss_pct_float < 0.0:
                raisable = True
                raise Exception(
                    "ERROR:::: NO SUBSAMPLING WILL BE EXECUTED::: you must specify a float from 0.0-1.0"
                )
            else:
                ss_head = f" <( seqkit sample --line-width=0 --quiet --rand-seed=7  --seq-type=dna --proportion={ss_pct_float}  "
                ss_tail = " ) "

    except Exception as e:
        if raisable:
            print(
                "Samplesheet Error with subsample_pct column ---- \n\n", file=sys.stderr
            )
            raise (e)
        else:
            pass

    return (ss_head, ss_tail)


def get_subsample_head(wildcards):
    return get_subsample_head_tail(wildcards.sample)[0]


def get_subsample_tail(wildcards):
    return get_subsample_head_tail(wildcards.sample)[1]

def get_samp_name(wildcards):
    return wildcards.sample

# used if you wish to restart and not have bwa rerun b/c snakemake wants to for some reason

rule bwa_mem2_sort:
    """https://github.com/bwa-mem2/bwa-mem2"""
    """    Also used here, a common bfx toolset, samtools"""
    """    https://github.com/samtools/samtools"""
    input:
        DR=MDIR + "{sample}/{sample}.dirsetup.ready",
        f1=getR1s,
        f2=getR2s,
    output:
        bami=MDIR + "{sample}/align/bwa2a/{sample}.bwa2a.sort.bam.bai",
        bamo=MDIR + "{sample}/align/bwa2a/{sample}.bwa2a.sort.bam",
    priority: 49
    log:
        MDIR + "{sample}/align/bwa2a/logs/{sample}.bwa2a_sort.log",
    resources:
        threads=config['bwa_mem2a_aln_sort']['threads'],
        mem_mb=60000 if 'mem_mb' not in config['bwa_mem2a_aln_sort'] else config['bwa_mem2a_aln_sort']['mem_mb'],
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
        piped_align=config["bwa_mem2a_aln_sort"]['piped_align'],
        samp=get_samp_name,
    conda:
        config["bwa_mem2a_aln_sort"]["env_yaml"]
    shell:
        """
        export tdir={params.mdir}/{params.samp}/{params.samtmpd};
        mkdir -p $tdir ;
        epocsec=$(date +'%s');
        
        {params.ldpre} {params.bwa_mem2a_cmd} mem \
         -R '@RG\\tID:{params.rgid}_$epocsec\\tSM:{params.rgsm}\\tLB:{params.samp}{params.rglb}\\tPL:{params.rgpl}\\tPU:{params.rgpu}\\tCN:{params.rgcn}\\tPG:{params.rgpg}' \
         {params.softclip_alts}  {params.K} {params.k} -t {params.bwa_threads}  {params.huref} \
         {params.subsample_head} <(unpigz -c  -q -- {input.f1} )  {params.subsample_tail}  \
         {params.subsample_head} <(unpigz -c  -q -- {input.f2} )  {params.subsample_tail}    \
        |  samtools sort -l 0  -m {params.sort_thread_mem}   \
         -@  {params.sort_threads} -T $tdir -O SAM - \
        |  samtools view -b -@ {params.write_threads} -O BAM --write-index -o {output.bamo}##idx##{output.bami} -  >> {log};
        """


localrules: produce_bwa_mem2,

rule produce_bwa_mem2:  # TARGET: only produce bwamem2a
     input:
         expand(MDIR + "{sample}/align/bwa2a/{sample}.bwa2a.sort.bam", sample=SAMPS)
	 
