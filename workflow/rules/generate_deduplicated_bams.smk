# #### Terminal Rule to end processing at generating DDUPED BAMS



localrules:
    produce_deduplicated_bams,


rule produce_deduplicated_bams:  # TARGET : Generate Just BAMs with Dups Marked .
    input:
        expand(
            MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
    output:
        expand(MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.ddupgen.complete",sample=SAMPS, alnr=ALIGNERS)
    threads: 1
    shell:
        "touch {output};"
