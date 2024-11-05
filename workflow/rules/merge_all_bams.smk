import os

# #### merge_bam
# --------------
# github: 

  
def get_merge_samp(wildcards):
    ret_files = []
    for sample in samples[samples['samp'] == wildcards.sample ]['sample_lane']:
        ret_files.append( MDIR + f"{sample}/align/{wildcards.alnr}/{sample}.{wildcards.alnr}.sort.bam")
    return ret_files

 
def get_lane_from_samp(wildcards):
    if len( samples[samples['samp'] == wildcards.sample]) in [0,'0']:
        raise Exception(f"LEN ZERO 0 :: {samples[samples['samp'] == wildcards.sample]}")
        
    return len(samples[samples['samp'] == wildcards.sample])
    
    
rule merge_bam:
    input:
        get_merge_samp
    output:
        bamo=MDIR+ "{sample}/align/{alnr}/{sample}.{alnr}.mg.sort.bam",
        baio=MDIR+ "{sample}/align/{alnr}/{sample}.{alnr}.mg.sort.bam.bai",
        bamm=MDIR+ "{sample}/align/{alnr}/{sample}.{alnr}.mg.sort.bam.merge.bam",
    log:
        MDIR+ "{sample}/align/{alnr}/log_bammerge/{sample}.{alnr}.mg.sort.bam.log"
    benchmark:
        MDIR+ "{sample}/benchmark/{sample}.{alnr}.mg.sort.bam.log"
    threads: config['merge_bam']['threads']
    params:
        glfs=get_lane_from_samp,
        smem=config['merge_bam']['smem'],
        mmem=config['merge_bam']['mmem'],
        cluster_sample=ret_sample,
        mthreads=config['merge_bam']['mthreads'],
        sthreads=config['merge_bam']['sthreads'],
    resources:
        threads=config['merge_bam']['threads'],
        vcpu=config['merge_bam']['threads'],
        partition=config['merge_bam']['partition'],
    conda:
        "../envs/sambamba_v0.1.yaml"
    shell:
        """

        if [[ '{params.glfs}' == '1' ]]; then
           cp {input} {output.bamo};
           cp {input}.bai {output.baio};
        else
            sambamba merge -t {params.mthreads} -p  {output.bamm} {input} ;
            sambamba sort -t {params.sthreads} -m {params.smem} -p -o {output.bamo} -l 9 --tmpdir /fsx/scratch  {output.bamm}  ;
            sambamba index -t {threads} -p {output.bamo} {output.baio};

        fi;

        touch {input};
        sleep 0.2;
        touch {output};

        """

localrules:
    merge_all_bams,


rule merge_all_bams:  # TARGET merge_all_bams
    input:
        b=expand(MDIR+ "{sample}/align/{alnr}/{sample}.{alnr}.mg.sort.bam", sample=SSAMPS.keys(), alnr=ALIGNERS)
    output:
        "merge.done",
    shell:
        "touch {output}"
