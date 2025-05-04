# #### Terminal Rule to end processing at generating DDUPED BAMS



localrules:
    produce_deduplicated_crams,


rule produce_deduplicated_crams:  # TARGET : Generate Just BAMs with Dups Marked .
    input:
        expand(
            MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.cram",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
    output:
        expand(MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.ddupgen.complete",sample=SAMPS, alnr=ALIGNERS)
    threads: 1
    shell:
        "touch {output};"
