#########  biobambam2
# --------------------------
# name==biobambam2
# comment="pretty fast."  
# ref==https://gitlab.com/german.tischler/biobambam2
#



def get_dppl_inputs(wildcards):
    ret_files = []
    for samp in samples[samples['samp'] == wildcards.sx ]['sample_lane']:
        ret_files.append( MDIR + f"{samp}/align/{wildcards.alnr}/{samp}.{wildcards.alnr}.sort.bam")
        ret_files.append( MDIR + f"{samp}/align/{wildcards.alnr}/{samp}.{wildcards.alnr}.sort.bam.bai")

    return ret_files


if "bbb2" in DDUP:


    rule biobambam2_mkdup:
        """Runs duplicate marking on the raw BAM."""
        input:
            bam=MDIR + "{sx}/align/{alnr}/{sx}.{alnr}.sort.bam",
            bai=MDIR + "{sx}/align/{alnr}/{sx}.{alnr}.sort.bam.bai",
        priority: 3
        output:
            bamo="{MDIR}{sx}/align/{alnr}/{sx}.{alnr}.mrkdup.sort.bam",
            bami="{MDIR}{sx}/align/{alnr}/{sx}.{alnr}.mrkdup.sort.bam.bai",
        threads: config["doppelmark_markdups"]["threads"]
        benchmark:
            repeat("{MDIR}{sx}/benchmarks/{sx}.{alnr}.dppl.mrkdup.bench.tsv", 0)
        container:
            "docker://daylilyinformatics/biobambam2:2.0.0"
        resources:
            threads=config["doppelmark_markdups"]["threads"],
            partition=config["doppelmark_markdups"]["partition"],
            vcpu=config["doppelmark_markdups"]["threads"],
        params:
            na=1,
            cluster_sample=ret_sx, #
        log:
            "{MDIR}{sx}/align/{alnr}/logs/dedupe.{sx}.{alnr}.log",
        shell:
            """
            bammarkduplicates2 -@ {threads} \
            I={input.bam} \
            O={output.bamo} \
            colhashbits=22 \
            collistsize=2147483648 \
            fragbufsize=4294967296 \
            inputbuffersize=262144 \
            optminpixeldif=1000 \
            index=1;
            {latency_wait};
            """
