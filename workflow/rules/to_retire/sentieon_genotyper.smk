#### Sentieon Genotyper
# --------------------------
#
# Sentieons Implementation of the GATK Genotyper
#


rule sentieon_genotyper:
    input:
        bam=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
    output:
        vcf=MDIR + "{sample}/align/{alnr}/snv/sgt/{sample}.{alnr}.sgt.snv.vcf",
        sort_vcf_gz=MDIR
        + "{sample}/align/{alnr}/snv/sgt/{sample}.{alnr}.sgt.snv.sort.vcf.gz",
    log:
        MDIR + "{sample}/align/{alnr}/snv/sgt/log/{sample}.{alnr}.sgt.snv.log",
    threads: config["octopus"]["threads"]
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.sgt.bench.tsv"
    conda:
        config["sentieon"]["env_yaml"]
    params:
        huref=config["supporting_files"]["files"]["glimpse"]["GLIMPSE_huref"]["name"],  #make this dynamic
        genome_build=config["glimpse"]["genome_build"],
        genome_file=config["supporting_files"]["files"]["huref"]["genome_file"]["name"],
        skr=config["supporting_files"]["files"]["octopus"]["skip_regions_f"]["name"],
        mdir=MDIR,
        caller_rc=" --algo Genotyper --var_type both --call_conf 15 --emit_conf  15 --emit_mode variant --genotype_model multinomial ",
        cluster_sample=ret_sample,
        jemalloc=os.popen(f"head -n 1 {config['jem_dot']} || echo '';")
        .readline()
        .rstrip(),
    shell:
        "/fsx/data/cached_envs/sentieon-genomics-202308.03/bin/sentieon driver -t {threads} -r {params.huref} -i {input.bam} {params.caller_rc}  --ploidy 2 {output.vcf} > {log} 2>&1;"
        "bcftools sort -m 100G -O z -o {output.sort_vcf_gz}  {output.vcf};"
        "bcftools index --threads {threads} -t {output.sort_vcf_gz};"
        "{latency_wait};"
        "ls {output};"
        "touch {params.mdir}{wildcards.sample}/benchmarks/{wildcards.sample}.{wildcards.alnr}.sgt.bench.tsv ;"
        "ls {params.mdir}{wildcards.sample}/benchmarks/{wildcards.sample}.{wildcards.alnr}.sgt.bench.tsv ;"


rule index_sent_genotype_snv:
    input:
        MDIR + "{sample}/align/{alnr}/snv/sgt/{sample}.{alnr}.sgt.snv.sort.vcf.gz",
    output:
        MDIR + "{sample}/align/{alnr}/snv/sgt/{sample}.{alnr}.sgt.snv.sort.vcf.gz.tbi",
    conda:
        config["sentieon"]["env_yaml"]
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/sgt/log/{sample}.{alnr}.sgt.snv.sort.vcf.gz.tbi.log",
    params:
        huref=config["supporting_files"]["files"]["glimpse"]["GLIMPSE_huref"]["name"],  #make this dynamic
        genome_build=config["glimpse"]["genome_build"],
        genome_file=config["supporting_files"]["files"]["huref"]["genome_file"]["name"],
        caller_rc=" --algo Genotyper --emit_mode variant --genotype_model multinomial ",
        constrain_to_chrm=(" " if "constrain_to_chrm" not in config["glimpse"]  else config["glimpse"]["constrain_to_chrm"]        ),
    shell:
        "sentieon driver -t {threads} {params.constrain_to_chrm} {params.constrain_to_hc_regions} -r {params.huref} -i {input.bam} {params.caller_rc} --given {input.sites_ref_vcfgz} --ploidy 2 {output.vcf} > {log.a} 2>&1;"
        "bcftools sort -m 100G -O z -o {output.sort_vcf_gz}  {output.vcf};"
        "bcftools index --threads {threads} -t {output.sort_vcf_gz};"
        "python3 workflow/scripts/add_custom_annotations2.py {input.sites_ref_vcfgz} {output.sort_vcf_gz} {wildcards.gli\mpse_sample} {output.fin_vcfgz} {output.norm_sort_vcf} {output.norm_sort_vcf_gz} {log.a} {threads} {log.b} {params.genome_file} {params.cstart} {params.cend} >> {log.a} 2>&1;"
        "echo COMPLETEDphase_2;"



"tabix -f {input};"
        "{latency_wait};"
        "ls {output};"


localrules:
    produce_sent_unified_gtyper,


rule produce_sent_unified_gtyper:  # TARGET: collect all sent results
    input:
        expand(
            MDIR
            + "{sample}/align/{alnr}/snv/sgt/{sample}.{alnr}.sgt.snv.sort.vcf.gz.tbi",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
    output:
        "gatheredall.sgt",
    shell:
        "touch {output};"

