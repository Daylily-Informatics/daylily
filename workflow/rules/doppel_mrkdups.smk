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
            mem_mb=config["doppelmark"]["mem_mb"],
        params:
            cluster_sample='na',
	    numa=config['doppelmark']['numa'],
            shard_size=config['doppelmark']['shard_size'],
            clip_padding=config['doppelmark']['clip_padding'],
            min_bases=config['doppelmark']['min_bases'],
            queue_length=config['doppelmark']['queue_length'],
	    read_buffer_size=config['doppelmark']['read_buffer_size'],
	    huref_fasta=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        log:
            "{MDIR}{sx}/align/{alnr}/logs/dedupe.{sx}.{alnr}.log",
        shell:
            """
            touch {log};
            TOKEN=$(curl -X PUT 'http://169.254.169.254/latest/api/token' -H 'X-aws-ec2-metadata-token-ttl-seconds: 21600');
            itype=$(curl -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-type);
            echo "INSTANCE TYPE: $itype" > {log};
	    echo "INSTANCE TYPE: $itype";
            start_time=$(date +%s);

            {params.numa} resources/DOPPLEMARK/doppelmark \
             -parallelism {threads} \
             -bam {input.bam} \
             -clip-padding {params.clip_padding} \
             -output {output.bamo} \
             -logtostderr \
	     -disk-mate-shards 0 \
	     -max-depth 300000 \
	     -scratch-dir $TMPDIR \
             -min-bases {params.min_bases} \
             -queue-length {params.queue_length} \
             -shard-size {params.shard_size}  >> {log};

            samtools index -b {output.bamo} >> {log};
            
            touch {output};

            end_time=$(date +%s);
	    elapsed_time=$((($end_time - $start_time) / 60));
	    echo "Elapsed-Time-min:\t$itype\t$elapsed_time";
            echo "Elapsed-Time-min:\t$itype\t$elapsed_time" >> {log} 2>&1;

	    cram_cmd="samtools view -@ {threads} -m 2G  -C -T {params.huref_fasta}   --write-index  -o  {output.bamo}.cram  {output.bamo}";
	    echo "$cram_cmd";
	    echo "$cram_cmd" >> {log};

            """ 
