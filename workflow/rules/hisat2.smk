# ###### Hisat2 - Specialized Aligner and Caller for Many HTDs
# ------------------------------------------------------------
# tool: Hisat2
# what: graph based aligner and var caller!
# docs: https://daehwankimlab.github.io/
# paper: https://genome.cshlp.org/content/31/7/1290
# github: https://github.com/DaehwanKimLab/hisat2
#




rule hisat2_align_sort:
    input:
        f1=getR1s,
        f2=getR2s,
    output:
        bamo=MDIR + "{sample}/align/hi2/{sample}.hi2.sort.bam",
        bami=MDIR + "{sample}/align/hi2/{sample}.hi2.sort.bam.bai",
        samo=temp(MDIR + "{sample}/align/hi2/{sample}.hi2.sam"),
    log:
        MDIR + "{sample}/align/hs2/logs/{sample}.hi2.sort.log",
    threads: config["hisat2"]["threads"]
    priority: 5
    benchmark:
        repeat(
            MDIR + "{sample}/benchmarks/{sample}.hi2.aln.bench.tsv",0
        )
    resources:
        threads=config["hisat2"]["threads"],
        vcpu=config["hisat2"]["threads"],
        partition="i192,i192mem"
    params:
        K=400000400,
        huref=config["supporting_files"]["files"]["huref"]["hisat2"]["name"],
        tmpd=MDIR + "{sample}/align/hi2/tmp",
        cluster_sample="{wildcards.sample}",
        extra="",
        smem=config["hisat2"]["smem"],
        wthreads=config["hisat2"]["wthreads"],
        sthreads=config["hisat2"]["sthreads"],
        tmpdir="/fsx/scratch",
        ldpre=config["hisat2"]["ldpre"],
        rgpl="presumedILLUMINA",  # ideally: passed in technology # nice to get to this point: https://support.sentieon.com/appnotes/read_groups/  :: note, the default sample name contains the RU_EX_SQ_Lane (0 for combined)
        rgpu="presumedCombinedLanes",  # ideally flowcell_lane(s)
        rgsm=ret_sample,  # samplename
        rgid=ret_sample,  # ideally samplename_flowcell_lane(s)_barcode  ! Imp this is unique, I add epoc seconds to the end of start of this rule
        rglb="_presumedNoAmpWGS",  # prepend sample_name in shell block ideally samplename_libprep
        rgcn="CenterName",  # center name
        rgpg="hisat2",  #program
    conda:
        config["hisat2"]["env_yaml"]
    shell:
        """
        export epocsec=$(date +'%s');
        {params.ldpre} hisat2 -p {threads}  -x {params.huref} \
        --rg-id 'ID:{params.rgid}_$epocsec' \
        --rg 'SM:{params.rgsm}' \
        --rg 'LB:{params.rgsm}{params.rglb}' \
        --rg 'PL:{params.rgpl}' \
        --rg 'PU:{params.rgpu}' \
        --rg 'CN:{params.rgcn}' \
        --rg 'PG:{params.rgpg}' \
        -1  {input.f1}    \
        -2  {input.f2}    \
        |  samtools sort -l 0  -m {params.smem}   \
        -@  {params.sthreads} -T {params.tmpdir} -O SAM  - \
        |  samtools view -b -1  -@ {params.wthreads} -O BAM --write-index -o {output.bamo}##idx##{output.bami} -  >> {log};
        touch {output.samo};
        """


localrules: produce_hisat2

rule produce_hisat2:  # TARGET: hisat2 only
    input:
        expand(MDIR + "{sample}/align/hi2/{sample}.hi2.sort.bam", sample=SSAMPS)
