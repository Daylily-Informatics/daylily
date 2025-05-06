#########  doppelmark
# --------------------------
# name=doppelmark
# comment="v fast dedup marking."
# git=https://github.com/grailbio/doppelmark
#



if "dppl" in DDUP:


    rule doppelmark_dups:
        """Runs duplicate marking on the merged(mg) CRAM."""
        input:
            bam=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.sort.bam",
            bai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.sort.bam.bai",
        priority: 3
        output:
            cram="{MDIR}{sample}/align/{alnr}/{sample}.{alnr}.cram",
            crai="{MDIR}{sample}/align/{alnr}/{sample}.{alnr}.cram.crai",
        wildcard_constraints:
            alnr="|".join(OG_ALIGNERS)
        threads: config["doppelmark"]["threads"]
        benchmark:
            repeat("{MDIR}{sample}/benchmarks/{sample}.{alnr}.mrkdup.bench.tsv", 0)
        conda:
            "../envs/doppelmark_v0.1.yaml"
        resources:
            threads=config["doppelmark"]["threads"],
            partition=config["doppelmark"]["partition"],
            vcpu=config["doppelmark"]["threads"],
            mem_mb=config["doppelmark"]["mem_mb"],
            constraint=config["doppelmark"]["constraint"],
        params:
            cluster_sample=ret_sample,
	        numa=config['doppelmark']['numa'],
            doppelmark_threads=config['doppelmark']['doppelmark_threads'],
            cram_compression=config["doppelmark"]["cram_compression"],
            shard_size=config['doppelmark']['shard_size'],
            clip_padding=config['doppelmark']['clip_padding'],
            min_bases=config['doppelmark']['min_bases'],
            queue_length=config['doppelmark']['queue_length'],
	        huref_fasta=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
	        compress_mem=config["doppelmark"]["compress_mem"],
	        compress_threads=config["doppelmark"]["compress_threads"],
            view_threads=config["doppelmark"]["view_threads"],
            view_mem=config["doppelmark"]["view_mem"],
            mbuffer_mem=config["doppelmark"]["mbuffer_mem"],
        log:
            "{MDIR}{sample}/align/{alnr}/logs/dedupe.{sample}.{alnr}.log",
        shell:
            """
            touch {log};
            TOKEN=$(curl -X PUT 'http://169.254.169.254/latest/api/token' -H 'X-aws-ec2-metadata-token-ttl-seconds: 21600');
            itype=$(curl -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-type);
            echo "INSTANCE TYPE: $itype" > {log};
	        echo "INSTANCE TYPE: $itype";
            start_time=$(date +%s);
            
            timestamp=$(date +%Y%m%d%H%M%S);

            TMPDIR=/fsx/scratch/doppel_tmp_$timestamp;
            #mkdir -p $TMPDIR;
            APPTAINER_HOME=$TMPDIR;
            #trap "rm -rf \"$TMPDIR\" || echo '$TMPDIR rm fails' >> {log} 2>&1" EXIT;
            tdir=$TMPDIR; 

            {params.numa} resources/DOPPLEMARK/doppelmark \
             -parallelism {params.doppelmark_threads} \
             -bam {input.bam} \
             -clip-padding {params.clip_padding} \
             -logtostderr \
    	     -disk-mate-shards 0 \
	         -scratch-dir $tdir \
             -min-bases {params.min_bases} \
             -queue-length {params.queue_length} \
             -shard-size {params.shard_size}   \
            | mbuffer -m {params.mbuffer_mem} \
            | samtools view -@ {params.view_threads} -m {params.view_mem} --output-fmt-option level={params.cram_compression} -C -T {params.huref_fasta}   --write-index  -o  {output.cram} - >> {log} 2>&1; 


            end_time=$(date +%s);
    	    elapsed_time=$((($end_time - $start_time) / 60));
	        echo "Elapsed-Time-min:\t$itype\t$elapsed_time";
            echo "Elapsed-Time-min:\t$itype\t$elapsed_time" >> {log} 2>&1;
            
            """ 
