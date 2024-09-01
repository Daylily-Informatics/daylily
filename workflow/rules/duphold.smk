##### DUPHOLD USESS snv PLOIDY RATIOS TO FINE TUNE sv CALLS
# --------------------------------------------------------
#
# Effectively a SV call filter which uses variant call zygosity
# to support or downgrade SV calls.
# github: https://github.com/brentp/duphold
# paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6479422/

rule duphold:
    """https://github.com/brentp/duphold"""
    """     Static binary release https://github.com/brentp/duphold/releases/tag/v0.2.3"""
    input:
        bam=(MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam"),
        bai=(MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam.bai"),
        snv_vcf=(
            MDIR
            + "{sample}/align/{alnr}/snv/{snv_caller}/{sample}.{alnr}.{snv_caller}.snv.sort.vcf.gz"
        ),
        snv_vcftbi=(
            MDIR
            + "{sample}/align/{alnr}/snv/{snv_caller}/{sample}.{alnr}.{snv_caller}.snv.sort.vcf.gz.tbi"
        ),
        sv_vcf=(
            MDIR
            + "{sample}/align/{alnr}/sv/{s_v_caller}/{sample}.{alnr}.{s_v_caller}.sv.sort.vcf.gz"
        ),
        sv_vcftbi=(            MDIR            + "{sample}/align/{alnr}/sv/{s_v_caller}/{sample}.{alnr}.{s_v_caller}.sv.sort.vcf.gz.tbi"        ),
    output:
        vcf= MDIR + "{sample}/align/{alnr}/sv/{s_v_caller}/{sample}.{alnr}.{snv_caller}.{s_v_caller}-d.vcf",
    priority: 48
    benchmark:
        (
            MDIR
            + "{sample}/benchmarks/{sample}.{alnr}.{snv_caller}.{s_v_caller}-duphold.bench.tsv"
        )
    log: MDIR + "{sample}/align/{alnr}/sv/DUPHOLD/{s_v_caller}_{snv_caller}/duphold.log"
    threads: config["duphold"]["threads"]
    conda:
        config["vanilla"]["env_yaml"]
    resources:
        vcpu=config["duphold"]["threads"],
        threads=config['duphold']['threads'],
    params:
        duphold_bin="resources/duphold/duphold",
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        ld_p=config['dysgu']['ld_preload'] if 'ld_preload' not in config['dysgu'] else config['dysgu']['ld_preload'],
        cluster_sample=ret_sample_sv,
    shell:
        """
        ( (rm -rf {output}) || echo rmFailedDUPHOLD;
        mkdir -p $( dirname {output.vcf} ) ;
        {params.duphold_bin} -s {input.snv_vcf} -t {threads} -v {input.sv_vcf}  -b {input.bam} -f {params.huref} -o {output.vcf} ;
        {latency_wait};
        ls {output}; ) > {log};
        """



# I've broken out the sort and indexing because for a while there was a difficult to debug bug
# obscured by sort/bgzip. They would all return, but not necessarry
rule duphold_sort_index:
    input:
        vcf= MDIR + "{sample}/align/{alnr}/sv/{s_v_caller}/{sample}.{alnr}.{snv_caller}.{s_v_caller}-d.vcf",
    output:
        vcfsort =  MDIR + "{sample}/align/{alnr}/sv/{s_v_caller}/{sample}.{alnr}.{snv_caller}.{s_v_caller}-d.sort.vcf",
        vcfgz=   MDIR + "{sample}/align/{alnr}/sv/{s_v_caller}/{sample}.{alnr}.{snv_caller}.{s_v_caller}-d.sort.vcf.gz",
        vcfgztbi=MDIR + "{sample}/align/{alnr}/sv/{s_v_caller}/{sample}.{alnr}.{snv_caller}.{s_v_caller}-d.sort.vcf.gz.tbi",
    benchmark:
        (
            MDIR
            + "{sample}/benchmarks/{sample}.{alnr}.{snv_caller}.{s_v_caller}-d.sort.bench.tsv"
        )
    threads: config["duphold"]["threads"]
    priority: 49
    resources:
        vcpu=config["duphold"]["threads"],
        threads=config['duphold']['threads'],
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR + "{sample}/align/{alnr}/sv/{s_v_caller}/logs/{sample}.{alnr}.{snv_caller}.{s_v_caller}-d.sort.log",
    params:
        cluster_sample=ret_sample_sv,
    shell:
        """
        set +euo pipefail;
        (
        (rm {output}) || echo rmFailedOK;
        bedtools sort -header -i {input.vcf} > {output.vcfsort}  ;
        bgzip -f  {output.vcfsort};
        touch {output.vcfsort};
        tabix  -p vcf -f {output.vcfgz} ;
        ls {output.vcfgztbi} ;
        {latency_wait}; ) > {log} 2>&1;
        """


localrules:
    produce_all_svs,


rule produce_all_svs:   # TARGET: gather cnvs calls and duphold
    input:
        expand(
            MDIR
            + "{sample}/align/{alnr}/sv/{s_v_caller}/{sample}.{alnr}.{snv_caller}.{s_v_caller}-d.sort.vcf.gz.tbi",
            sample=SSAMPS,
            alnr=ALIGNERS,
            s_v_caller=sv_CALLERS,
            snv_caller=snv_CALLERS,
        ),
    output:
        MDIR + "logs/all_svVCF_dupheld.done",
    threads: 1
    conda:
        config["vanilla"]["env_yaml"]
    shell:
        "echo {input}; "
        "touch {output};  "
