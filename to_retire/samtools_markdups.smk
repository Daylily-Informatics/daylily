##### SAMTOOLS MRKDUP
# ---------------------------------------
# Old reliable.  SAMBAMBA has a random
# memory leak
#


rule samtools_markdups:
    input:
        "{MDIR}{sample}/align/{alnr}/{sample}.{alnr}.sort.bam",
    output:
        mdbam=("{MDIR}{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam"),
        mdbai=("{MDIR}{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam.bai"),
        colbam=temp("{MDIR}{sample}/align/{alnr}/{sample}.{alnr}.collate.bam"),
        fixbam=temp("{MDIR}{sample}/align/{alnr}/{sample}.{alnr}.fixmate.bam"),
        posbam=temp("{MDIR}{sample}/align/{alnr}/{sample}.{alnr}.posn.bam"),
        tmpd=(directory("{MDIR}{sample}/align/{alnr}/tmpmkdup/")),
        stats=("{MDIR}{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.stats"),
    priority: 35
    benchmark:
        "{MDIR}{sample}/benchmarks/{sample}.{alnr}.mrkdup.bench.tsv"
    log:
        "{MDIR}{sample}/logs/{sample}.{alnr}.mrkdup.log",
    threads: config["samtools_markdups"]["threads"]
    params:
        cluster_sample=f"wildcards.sample",
    conda:
        config["samtools_markdups"]["env_yaml"]
    wildcard_constraints:
        alnr="bwa2.*",
    shell:
        "mkdir -p {output.tmpd};"
        "ulimit -Sn 100000 &> {log};"
        "echo START &>> {log};"
        "samtools collate -n 64 -u -@ {threads} -o {output.colbam} {input} &>> {log};"
        "echo COLLATE &>> {log};"
        "samtools fixmate  -u -@ {threads} -m {output.colbam} {output.fixbam} &>> {log};"
        "echo FIXMATE &>> {log};"
        "samtools sort -@ {threads} -o {output.posbam} {output.fixbam} &>> {log};"
        "echo POSN &>> {log};"
        "samtools markdup -T {output.tmpd}/$JOB_ID. -s -f {output.stats} -c -d 2500 -m s -S -@ {threads} {output.posbam} {output.mdbam} &>> {log};"
        "echo DONE &>> {log};"
        "samtools index -b -@ {threads} {output.mdbam} &> {log};"
        "echo VDONE &>> {log}; "
        # The bai index is made from the bam full name
