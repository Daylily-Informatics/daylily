import sys
import os

##### Clair3
# ---------------------------

def get_clair_model_path(wildcards):
    # using the aligner is a hack, should move this to config and pulled from seq technology
    instrument = samples[samples["samp"] == wildcards.sample]["instrument"][0].lower()
    
    model_path = "/opt/models/ilmn"
    if instrument == "ont":
        model_path = "/opt/models/ont"
    elif instrument == "pacbio":
        model_path = "/opt/models/pacbio"
    elif instrument in ["illumina","novaseq","ilmn"]:
        model_path = "/opt/models/ilmn"
    elif instrument == "ultima":
         model_path = "/opt/models/ultima"
    else:
        model_path = f"/opt/models/{instrument}"
        print(
            f"Unknown instrument {instrument}, going to take a shot and set the clair3 model path to: {model_path} !!!",
            file=sys.stderr,
        )
    return model_path

def get_clair_platform(wildcards):
    # using the aligner is a hack, should move this to config and pulled from seq technology
    instrument = samples[samples["samp"] == wildcards.sample]["instrument"][0].lower()
    
    platform = "ilmn"

    if instrument == "ont":
        platform = "ont"
    elif instrument in ["pacbio","hifi"]:
        platform = "hifi"
    elif instrument in ["illumina","novaseq","ilmn"]:
        platform = "ilmn"
    elif instrument == "ultima":
        platform = "ilmn"
    else:
        platform = f"{instrument}"
        print(
            f"Unknown instrument {instrument}, going to take a shot and set the clair3 model path to: {platform} !!!",
            file=sys.stderr,
        )
    return platform


def get_clair3_chrom(wildcards):
    pchr="" #prefix handled already
    ret_str = ""
    sl = wildcards.clairchrm.replace('chr','').split("-")
    sl2 = wildcards.clairchrm.replace('chr','').split("~")
    
    if len(sl2) == 2:
        ret_str = pchr + wildcards.clairchrm
    elif len(sl) == 1:
        ret_str = pchr + sl[0]
    elif len(sl) == 2:
        start = int(sl[0])
        end = int(sl[1])
        while start <= end:
            ret_str = str(ret_str) + " " + pchr + str(start)
            start = start + 1
    else:
        raise Exception(
            "Clair3 chunks can only be one contiguous range per chunk: e.g., 1-4 with non-numerical chromosomes assigned 23=X, 24=Y, 25=MT, or 25=M in hg38"
        )

    return ret_mod_chrm(ret_str)

rule clair3:
    input:
        cram=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.cram",
        crai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.cram.crai",
        d=MDIR + "{sample}/align/{alnr}/snv/clair3/vcfs/{clairchrm}/{sample}.ready",
    output:
        vcf=MDIR
        + "{sample}/align/{alnr}/snv/clair3/vcfs/{clairchrm}/{sample}.{alnr}.clair3.{clairchrm}.snv.vcf.gz"
    log:
        MDIR + "{sample}/align/{alnr}/snv/clair3/log/{sample}.{alnr}.clair3.{clairchrm}.snv.log",
    threads: config['clair3']['threads']
    container:
        "docker://hkubal/clair3:latest"
    priority: 45
    resources:
        vcpu=config['clair3']['threads'],
        threads=config['clair3']['threads'],
        partition=config['clair3']['partition'],
        mem_mb=config['clair3']['mem_mb'],
    benchmark:
        repeat(
            MDIR + "{sample}/benchmarks/{sample}.{alnr}.clair3.{clairchrm}.bench.tsv",
            0
            if "bench_repeat" not in config["clair3"]
            else config["clair3"]["bench_repeat"],
        )
    params:
        cchrm=get_clair3_chrom,
        cluster_sample=ret_sample,
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        mdir=MDIR,
        samview_threads=config['clair3']['samview_threads'],
        mem_mb=config['clair3']['mem_mb'],
        numa=config['clair3']['numa'],
        clair3_threads=config['clair3']['clair3_threads'],  
        cpre="" if "b37" == config['genome_build'] else "chr",
        mito_code="MT" if "b37" == config['genome_build'] else "M",
        model_path=get_clair_model_path,
        platform=get_clair_platform,
    shell:
        """
        touch {log};
        TOKEN=$(curl -X PUT 'http://169.254.169.254/latest/api/token' -H 'X-aws-ec2-metadata-token-ttl-seconds: 21600');
        itype=$(curl -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-type);
        echo "INSTANCE TYPE: $itype" > {log};



        timestamp=$(date +%Y%m%d%H%M%S);
        TMPDIR=/fsx/scratch/clair3_tmp_$timestamp;
        mkdir -p $TMPDIR;
        APPTAINER_HOME=$TMPDIR;
        trap "rm -rf $TMPDIR" EXIT;
        tdir=$TMPDIR;

        # Log the start time as 0 seconds
        start_time=$(date +%s);
        echo "Start-Time-sec:$itype\t0" >> {log} 2>&1;

        cchr=$(echo {params.cpre}{params.cchrm} | sed 's/~/\:/g' | sed 's/23\:/X\:/' | sed 's/24\:/Y\:/' | sed 's/25\:/{params.mito_code}\:/');
        ##--bam_fn=<(samtools view -@ {params.samview_threads} -b -T {params.huref} {input.cram}) 
        echo "CCHRM: $cchr" >> {log} 2>&1;
        {params.numa}   /opt/bin/run_clair3.sh \
        --bam_fn={input.cram} \
        --ref_fn={params.huref} \
        --threads={params.clair3_threads} \
        --platform={params.platform} \
        --model_path={params.model_path}  \
        --sample_name={params.cluster_sample} \
        --ctg_name=$cchr \
        --output=$(dirname {input.d})  >> {log} 2>&1;

        ls -lth $(dirname {input.d})  >> {log} 2>&1;
        echo "CCHRM: $cchr" >> {log} 2>&1;
        
        mv $(dirname {input.d})/merge_output.vcf.gz {output.vcf};
        mv $(dirname {input.d})/merge_output.vcf.gz.tbi {output.vcf}.tbi;
        end_time=$(date +%s);
        elapsed_time=$((($end_time - $start_time) / 60));

        # Log the elapsed time
        echo "Elapsed-Time-min:\t$itype\t$elapsed_time" >> {log} 2>&1;
        """

rule clair3_sort_index_chunk_vcf:
    input:
        vcf=MDIR + "{sample}/align/{alnr}/snv/clair3/vcfs/{clairchrm}/{sample}.{alnr}.clair3.{clairchrm}.snv.vcf.gz"
    priority: 46
    output:
        vcfsort=MDIR + "{sample}/align/{alnr}/snv/clair3/vcfs/{clairchrm}/{sample}.{alnr}.clair3.{clairchrm}.snv.sort.vcf",
        vcfgz=MDIR + "{sample}/align/{alnr}/snv/clair3/vcfs/{clairchrm}/{sample}.{alnr}.clair3.{clairchrm}.snv.sort.vcf.gz",
        vcftbi=MDIR + "{sample}/align/{alnr}/snv/clair3/vcfs/{clairchrm}/{sample}.{alnr}.clair3.{clairchrm}.snv.sort.vcf.gz.tbi",
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR + "{sample}/align/{alnr}/snv/clair3/vcfs/{clairchrm}/log/{sample}.{alnr}.clair3.{clairchrm}.snv.sort.vcf.gz.log",
    resources:
        vcpu=4,
        threads=4,
        partition=config['clair3']['partition_other'],
    params:
        cluster_sample=ret_sample,
    threads: 4
    shell:
        """
        (bcftools sort -O v -o {output.vcfsort} {input.vcf} 2>> {log}) || exit 1233;
        (
        bgzip {output.vcfsort};        
        touch {output.vcfsort};
        tabix -f -p vcf {output.vcfgz};
        ) >> {log} 2>&1;
        
        ls {output} >> {log} 2>&1 ;
        
        {latency_wait};
        """

localrules:
    clair3_concat_fofn,

rule clair3_concat_fofn:
    input:
        chunk_tbi=sorted(
            expand(
                MDIR
                + "{{sample}}/align/{{alnr}}/snv/clair3/vcfs/{clairchm}/{{sample}}.{{alnr}}.clair3.{clairchm}.snv.sort.vcf.gz.tbi",
                clairchm=CLAIR3_CHRMS,
            ),
            key=lambda x: float(
                str(x.replace("~", ".").replace(":", "."))
               .split("vcfs/")[1]
                .split("/")[0]
                .split("-")[0]
            ),
        ),
    priority: 44
    output:
        fin_fofn=MDIR
        + "{sample}/align/{alnr}/snv/clair3/{sample}.{alnr}.clair3.snv.concat.vcf.gz.fofn",
        tmp_fofn=MDIR        + "{sample}/align/{alnr}/snv/clair3/{sample}.{alnr}.clair3.snv.concat.vcf.gz.fofn.tmp",
    threads: 2
    resources:
        threads=2,
        vcpu=2,
    params:
        fn_stub="{sample}.{alnr}.clair3.",
        cluster_sample=ret_sample,
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.clair3.concat.fofn.bench.tsv"
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR + "{sample}/align/{alnr}/snv/clair3/log/{sample}.{alnr}.clair3.concat.fofn.log",
    shell:
        """
        (rm {output} 1> /dev/null  2> /dev/null ) || echo rmFailOK >> {log} && ls ./ >> {log};

        for i in {input.chunk_tbi}; do
            ii=$(echo $i | perl -pe 's/\.tbi$//g'; );
            echo $ii >> {output.tmp_fofn};
        done;
        (workflow/scripts/sort_concat_chrm_list.py {output.tmp_fofn} {wildcards.sample}.{wildcards.alnr}.clair3. {output.fin_fofn})  >> {log} 2>&1;
        """

rule clair3_concat_index_chunks:
    input:
        fofn=MDIR
        + "{sample}/align/{alnr}/snv/clair3/{sample}.{alnr}.clair3.snv.concat.vcf.gz.fofn",
        tmp_fofn=MDIR        + "{sample}/align/{alnr}/snv/clair3/{sample}.{alnr}.clair3.snv.concat.vcf.gz.fofn.tmp",
    output:
        vcfgz=touch(
            MDIR + "{sample}/align/{alnr}/snv/clair3/{sample}.{alnr}.clair3.snv.sort.vcf.gz"
        ),
        vcfgztemp=temp(
            MDIR + "{sample}/align/{alnr}/snv/clair3/{sample}.{alnr}.clair3.snv.sort.temp.vcf.gz"
        ),
        vcfgztbi=touch(
            MDIR
            + "{sample}/align/{alnr}/snv/clair3/{sample}.{alnr}.clair3.snv.sort.vcf.gz.tbi"
        ),
    threads: 4
    resources:
        vcpu=4,
        threads=4,
        partition=config['clair3']['partition_other'],
    priority: 47
    params:
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        cluster_sample=ret_sample,
    resources:
        attempt_n=lambda wildcards, attempt:  (attempt + 0)
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.clair3.merge.bench.tsv"
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/clair3/log/{sample}.{alnr}.clair3.snv.merge.sort.gathered.log",
    shell:
        """

        touch {log};
        mkdir -p $(dirname {log});
        # This is acceptable bc I am concatenating from the same tools output, not across tools
        #touch {output.vcfgztemp};

        bcftools concat -a -d all --threads {threads} -f {input.fofn}  -O z -o {output.vcfgztemp} >> {log} 2>&1;

        export oldname=$(bcftools query -l {output.vcfgztemp} | head -n1) >> {log} 2>&1;
        echo -e "${{oldname}}\t{params.cluster_sample}" > {output.vcfgz}.rename.txt
        bcftools reheader -s {output.vcfgz}.rename.txt -o {output.vcfgz} {output.vcfgztemp} >> {log} 2>&1;
        bcftools index -f -t --threads {threads} -o {output.vcfgztbi} {output.vcfgz} >> {log} 2>&1;

        rm -rf $(dirname {output.vcfgz})/vcfs >> {log} 2>&1;

        """

localrules:
    clear_combined_clair3_vcf,

rule clear_combined_clair3_vcf:
    input:
        vcf=expand(
            MDIR + "{sample}/align/{alnr}/snv/clair3/{sample}.{alnr}.clair3.snv.sort.vcf.gz",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
    priority: 42
    shell:
        "(rm {input.vcf}*  1> /dev/null  2> /dev/null ) || echo 'file not found for deletion: {input}';"

rule produce_clair3_vcf:  # TARGET: clair3 vcf
    input:
        vcftb=expand(
            MDIR
            + "{sample}/align/{alnr}/snv/clair3/{sample}.{alnr}.clair3.snv.sort.vcf.gz",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
        vcftbi=expand(
            MDIR
            + "{sample}/align/{alnr}/snv/clair3/{sample}.{alnr}.clair3.snv.sort.vcf.gz.tbi",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
    output:
        "gatheredall.clair3",
    threads: 4
    priority: 48
    log:
        "gatheredall.clair3.log",
    conda:
        "../envs/vanilla_v0.1.yaml"
    shell:
        """
        # Convert VCF to BCF and index it
        for vcf in {input.vcftb}; do
            bcf="${{vcf%.vcf.gz}}.bcf";
            bcftools view -O b -o $bcf --threads {threads} $vcf && bcftools index --threads 4 $bcf;
        done;

        # Mark the output as completed
        touch {output};

        # Log completion and list output
        {latency_wait};
        ls {output} >> {log} 2>&1;

        {latency_wait}; 
        ls {output}  >> {log} 2>&1;
        """

localrules:
    prep_clair3_chunkdirs,

rule prep_clair3_chunkdirs:
    input:
        b=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.cram",
        ba=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.cram.crai",
    output:
        expand(
            MDIR + "{{sample}}/align/{{alnr}}/snv/clair3/vcfs/{clairchrm}/{{sample}}.ready",
            clairchrm=CLAIR3_CHRMS,
        ),
    threads: 2
    log:
        MDIR + "{sample}/align/{alnr}/snv/clair3/log/{sample}.{alnr}.chunkdirs.log",
    shell:
        """
        ( echo {output}  ;
        mkdir -p $(dirname {output} );
        touch {output};
        ls {output}; ) > {log} 2>&1;
        """
