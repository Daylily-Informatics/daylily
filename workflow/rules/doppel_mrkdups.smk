#########  doppelmark
# --------------------------
# name=doppelmark
# comment="v fast dedup marking."
# note="seems fussy about some aligners output bams, ie: bwamem2-ert fails/...."  
# git=https://github.com/grailbio/doppelmark
#



if "dppl" in DDUP:


    rule doppelmark_dups:
        """Runs duplicate marking on the merged(mg) BAM."""
        input:
            bam=MDIR + "{sx}/align/{alnr}/{sx}.{alnr}.sort.bam",
            bai=MDIR + "{sx}/align/{alnr}/{sx}.{alnr}.sort.bam.bai",
        priority: 3
        output:
            bamo="{MDIR}{sx}/align/{alnr}/{sx}.{alnr}.mrkdup.sort.bam",
            bami="{MDIR}{sx}/align/{alnr}/{sx}.{alnr}.mrkdup.sort.bam.bai",
        threads: config["doppelmark"]["threads"]
        benchmark:
            repeat("{MDIR}{sx}/benchmarks/{sx}.{alnr}.mrkdup.bench.tsv", 0)
        conda:
            "../envs/vanilla_v0.1.yaml"
        resources:
            threads=config["doppelmark"]["threads"],
            partition=config["doppelmark"]["partition"],
            vcpu=config["doppelmark"]["threads"],
        params:
            cluster_sample='na',
            shard_size=config['doppelmark']['shard_size'],
            clip_padding=config['doppelmark']['clip_padding'],
            min_bases=config['doppelmark']['min_bases'],
            queue_length=config['doppelmark']['queue_length'],
        log:
            "{MDIR}{sx}/align/{alnr}/logs/dedupe.{sx}.{alnr}.log",
        shell:
            """

            resources/DOPPLEMARK/doppelmark \
             -parallelism {threads} \
             -bam {input.bam} \
             -clip-padding {params.clip_padding} \
             -output {output.bamo} \
             -logtostderr \
             -min-bases {params.min_bases} \
             -queue-length {params.queue_length} \
             -shard-size {params.shard_size} > {log};

            samtools index -b {output.bamo} >> {log};
            
            # using biobambabm2
            #LD_LIBRARY_PATH=resources/lib/  resources/biobambam2/bammarkduplicates I={input.bam} O={output.bamo} index=1 markthreads={threads}  M={output.bamo}.metrics fragbufsize=500331648000 colhashbits=23 collistsize=93554432   indexfilename={output.bami} ;
            touch {output};
            """