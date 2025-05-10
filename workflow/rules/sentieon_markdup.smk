#########  SENTIEON DEDUPING
# --------------------------

if "sent" in DDUP:

    def get_input_bams(wildcards):
        return MDIR + f"{wildcards.sample}/align/{wildcards.alnr}/{wildcards.sample}.{wildcards.alnr}.sort.bam",

    rule sentieon_markdups:
        """Runs duplicate marking on the BAM."""
        input:
            bam=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.sort.bam",
            bai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.sort.bam.bai",
        priority: 3
        output:
            cram="{MDIR}{sample}/align/{alnr}/{sample}.{alnr}.cram",
            crai="{MDIR}{sample}/align/{alnr}/{sample}.{alnr}.cram.crai",
            score="{MDIR}{sample}/align/{alnr}/{sample}.{alnr}.cram.score.txt",
            metrics="{MDIR}{sample}/align/{alnr}/{sample}.{alnr}.cram.mrkdup.metrics",
        threads: config["sentieon"]["threads"]
        benchmark:
            repeat("{MDIR}{sample}/benchmarks/{sample}.{alnr}.mrkdup.bench.tsv", 0)
        conda:
            config["sentieon"]["env_yaml"]
        params:
            cluster_sample=ret_sample,
	        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        resources:
            threads=config['sentieon']['threads'],
            partition=config['sentieon']['partition'],
            vcpu=config['sentieon']['threads'],
            mem_mb="160G",
        log:
            "{MDIR}{sample}/align/{alnr}/logs/dedupe.{sample}.{alnr}.log",
        shell:
            """     
            export bwt_max_mem={params.max_mem} ;

            if [ -z "$SENTIEON_LICENSE" ]; then
                echo "SENTIEON_LICENSE not set. Please set the SENTIEON_LICENSE environment variable to the license file path & make this update to your dyinit file as well." >> {log} 2>&1;
                exit 3;
            fi

            if [ ! -f "$SENTIEON_LICENSE" ]; then
                echo "The file referenced by SENTIEON_LICENSE ('$SENTIEON_LICENSE') does not exist. Please provide a valid file path." >> {log} 2>&1;
                exit 4;
            fi

            TOKEN=$(curl -X PUT 'http://169.254.169.254/latest/api/token' -H 'X-aws-ec2-metadata-token-ttl-seconds: 21600');
            itype=$(curl -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-type);
            echo "INSTANCE TYPE: $itype" > {log};
            start_time=$(date +%s);

            ulimit -n 65536 || echo "ulimit mod failed" > {log} 2>&1;
            
            timestamp=$(date +%Y%m%d%H%M%S);
            TMPDIR=/dev/shm/sentieon_mkduo_tmp_$timestamp;
            export SENTIEON_TMPDIR=$TMPDIR;
            mkdir -p $TMPDIR;
            APPTAINER_HOME=$TMPDIR;
            trap "rm -rf \"$TMPDIR\" || echo '$TMPDIR rm fails' >> {log} 2>&1" EXIT;
            tdir=$TMPDIR;

            # Find the jemalloc library in the active conda environment
            jemalloc_path=$(find "$CONDA_PREFIX" -name "libjemalloc*" | grep -E '\.so|\.dylib' | head -n 1); 

            # Check if jemalloc was found and set LD_PRELOAD accordingly
            if [[ -n "$jemalloc_path" ]]; then
                LD_PRELOAD="$jemalloc_path";
                echo "LD_PRELOAD set to: $LD_PRELOAD" >> {log};
            else
                echo "libjemalloc not found in the active conda environment $CONDA_PREFIX.";
                exit 3;
            fi
            LD_PRELOAD=$LD_PRELOAD  /fsx/data/cached_envs/sentieon-genomics-202503.01.rc1/bin/libexec/driver \
             --input {input.bam} \
             --reference {params.huref} \
             --thread_count {threads} \
             --interval_padding 0 \
             --algo Dedup \
             --score_info {output.score} \
             --cram_write_options version=3.0,compressor=rans \
             --metrics {output.metrics} \
            {output.cram} >> {log} 2>&1;


            """
