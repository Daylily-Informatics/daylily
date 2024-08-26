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
	    instance_log="{MDIR}{sx}/align/{alnr}/{sx}.{alnr}.mrkdup.sort.instance.log",
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
	    numa=config['doppelmark']['numa'],
            shard_size=config['doppelmark']['shard_size'],
            clip_padding=config['doppelmark']['clip_padding'],
            min_bases=config['doppelmark']['min_bases'],
            queue_length=config['doppelmark']['queue_length'],
        log:
            "{MDIR}{sx}/align/{alnr}/logs/dedupe.{sx}.{alnr}.log",
        shell:
            """
            touch {log} {output.instance_log};
            TOKEN=$(curl -X PUT 'http://169.254.169.254/latest/api/token' -H 'X-aws-ec2-metadata-token-ttl-seconds: 21600');
            itype=$(curl -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-type);
            echo "INSTANCE TYPE: $itype" >> {output.instance_log};
            start_time=$(date +%s);

            {params.numa} resources/DOPPLEMARK/doppelmark \
             -parallelism {threads} \
             -bam {input.bam} \
             -clip-padding {params.clip_padding} \
             -output {output.bamo} \
             -logtostderr \
             -min-bases {params.min_bases} \
             -queue-length {params.queue_length} \
             -shard-size {params.shard_size} > {log};

            samtools index -b {output.bamo} >> {log};
            
            touch {output};

            end_time=$(date +%s);
            elapsed_time=$((end_time - start_time));
            echo "Elapsed-Time-sec:\t$itype\t$elapsed_time >> {output.instance_log} 2>&1";
            """ 
