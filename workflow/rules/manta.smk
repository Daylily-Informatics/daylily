##### MANTA A VENERABLE sv CALLER
# -------------------------------
# The defacto standard SV caller, if there can be said to be one
# it hits f-scores of .6 +/- .2 depending how the truth set is framed
# which is better tha most.  It locked in python 2.6, and is a huge
# pain to get running even with conda helping.
# github: https://github.com/Illumina/manta
# a paper: https://doi.org/10.1093/bioinformatics/btv710

rule manta_get_centos_env:
    priority: 30
    conda:
        "../envs/vanilla_v0.1.yaml"  #  "../envs/manta_uge_v0.2.yaml"
    shell:
        "echo got it"


rule manta:
    """https://github.com/Illumina/manta"""
    input:
        bamo=f"{MDIR}" + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
        bami=f"{MDIR}" + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam.bai",
    output:
        vcf=f"{MDIR}" + "{sample}/align/{alnr}/sv/manta/{sample}.{alnr}.manta.sv.vcf",
    threads: config["manta"]["threads"]
    priority: 36
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.manta.bench.tsv"
    log:
        MDIR + "{sample}/align/{alnr}/sv/manta/logs/{sample}.{alnr}.manta.log",
    resources:
        vcpu=config["manta"]["threads"],
        partition=config["manta"]["partition"],
    params:
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        work_dir=MDIR + "{sample}/align/{alnr}/sv/manta/manta_work/",
        mdir=MDIR,
        log=MDIR + "{sample}/align/{alnr}/sv/manta/logs/{sample}.{alnr}.manta.log",
        tb=os.popen("which tabix").readline().rstrip(),
        bg=os.popen("which bgzip").readline().rstrip(),
        cluster_sample=ret_sample_alnr,
    conda:
        "../envs/manta_v0.1.yaml"  #config["manta"]["env_yaml"]
    shell:
        """

        
        timestamp=$(date +%Y%m%d%H%M%S);
        export TMPDIR=/fsx/scratch/manta_tmp_$timestamp;
        mkdir -p $TMPDIR;
        export APPTAINER_HOME=$TMPDIR;
        trap "rm -rf \"$TMPDIR\" || echo '$TMPDIR rm fails' >> {log} 2>&1" EXIT;

        set +euo pipefail;
        (rm -rf {params.work_dir} || echo rmFail;
        mkdir -p {params.work_dir}|| echo mkdirERROR ;
        configManta.py --bam {input.bamo} --reference {params.huref} --runDir {params.work_dir}  ;
        python {params.work_dir}/runWorkflow.py  -j {threads} || echo noFail;) > {log} 2>&1;
        python workflow/scripts/manta_uniter.py {params.work_dir}/results/variants/diploidSV.vcf.gz {params.work_dir}/results/variants/candidateSV.vcf.gz > {output.vcf} ;
        {latency_wait}; ls {output};
	    rm -rf {params.work_dir}/workspace;
        """


rule manta_sort_and_index:
    input:
        vcf=MDIR + "{sample}/align/{alnr}/sv/manta/{sample}.{alnr}.manta.sv.vcf",
    priority:37
    output:
        vcfsort=MDIR + "{sample}/align/{alnr}/sv/manta/{sample}.{alnr}.manta.sv.sort.vcf",
        vcfgz=MDIR + "{sample}/align/{alnr}/sv/manta/{sample}.{alnr}.manta.sv.sort.vcf.gz",
        vcfgztbi=MDIR + "{sample}/align/{alnr}/sv/manta/{sample}.{alnr}.manta.sv.sort.vcf.gz.tbi",
    log:
        MDIR+ "{sample}/align/{alnr}/sv/manta/logs/{sample}.{alnr}.manta.sv.sort.vcf.log",
    conda:
        "../envs/vanilla_v0.1.yaml"
    params:
        cluster_sample=ret_sample
    threads: 16
    resources:
        vcpu=16
    shell:
        """( (rm {output}) || echo noRm > {log};
        echo 0Index >> {log};
        bedtools sort -header -i {input.vcf} > {output.vcfsort};
        bgzip -f {output.vcfsort} ;
        tabix -f -p vcf {output.vcfgz};
        touch {output.vcfsort}; ) >> {log} 2>&1;
        {latency_wait}; ls {output};"""

localrules: produce_manta,

rule produce_manta:   # TARGET: just produce manta vcfs
    priority: 38
    input:
        expand(MDIR + "{sample}/align/{alnr}/sv/manta/{sample}.{alnr}.manta.sv.sort.vcf.gz.tbi",sample=SSAMPS,alnr=ALIGNERS)
