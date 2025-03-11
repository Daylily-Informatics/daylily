##### SVABA SV Detection By Local Assembly
# -------------------------------------
#
# github: https://github.com/walaj/svaba
# paper:  https://github.com/walaj/svaba/blob/master/gitfig_schematic.png
#
# Microbial PanDB:  https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0399-2
#

#
#Author: Jeremiah Wala <jwala@broadinstitute.org>
#Date:   Wed Nov 11 15:57:01 2020 -0800
#
#    issue 92 patch
#
# Vanilla opts suggested somatic as:
# wget "https://data.broadinstitute.org/snowman/dbsnp_indel.vcf" ## get a DBSNP known indel file
# DBSNP=dbsnp_indel.vcf
#   svaba run -t $TUM_BAM -n $NORM_BAM -p $CORES -D $DBSNP -a somatic_run -G $REF
##
#  germline
### Set -I to not do mate-region lookup if mates are mapped to different chromosome.
##   This is appropriate for germline-analysis, where we don't have a built-in control
##   to against mapping artifacts, and we don't want to get bogged down with mate-pair
##   lookups.
##  Set -L to 6 which means that 6 or more mate reads must be clustered to
##   trigger a mate lookup. This also reduces spurious lookups as above, and is more
##   appropriate the expected ALT counts found in a germline sample
##   (as opposed to impure, subclonal events in cancer that may have few discordant reads).
##
#
#  A flag offered in the program: --germline Sets recommended settings for case-only analysis (eg germline). (-I, -L5, assembles NM >= 3 reads)  Variant filtering and classification
#
# svaba run -t $GERMLINE_BAM -p $CORES -L 5 -I -a PREFIX -G $REF
#
#    I am going to add -k=ucsc mappable regions bed of sam style regios(blargh), -p # treads(35 suggested, thugh no mem control) -x 30000 -L 5
#
# SVABA optimizations taken from here
# https://github.com/walaj/svaba/issues/38
# https://github.com/walaj/svaba/issues/99
# https://github.com/walaj/svaba/issues/38

rule svaba:
    input:
        bam=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
        bai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam.bai",
    priority: 38
    output:
        stub =touch( MDIR + "{sample}/align/{alnr}/sv/svaba/{sample}.{alnr}"),
        vcff = touch(MDIR + "{sample}/align/{alnr}/sv/svaba/{sample}.{alnr}.svaba.sv.vcf"),
        vcfindel= touch(MDIR+"{sample}/align/{alnr}/sv/svaba/{sample}.{alnr}.svaba.indel.vcf"),
    params:
        ld_p=config['malloc_alt']['ld_preload'] if 'ld_preload' not in config['svaba'] else config['svaba']['ld_preload'],
        dbsnp=config['supporting_files']['files']['dbsnp']['broad_snowman_indelvcf']['name'],
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        svaba_cmd="./resources/svaba/svaba ",
        svaba_blacklist_regions=config["supporting_files"]["files"]["ucsc"]["problem_regions_bed"]['name'] ,
        cluster_sample=ret_sample_alnr,
        bwa_mem2a_cmd=config["bwa_mem2a_aln_sort"]["cmd"],  #Try and speed up by aliasing bwamem call to bwamem2
        problem_regions_bed=config["supporting_files"]["files"]["ucsc"]["problem_regions_bed"]["name"],
        simple_repeats_bed=config["supporting_files"]["files"]["ucsc"]["simple_repeats_bed"]["name"],   ## ! From UCSC
        bwv=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        svaba_calling_regions=config["supporting_files"]["files"]["giab_sv"]["v0.6_NA24385_ties_1_2"]['name'],
        cmd=config['svaba']['cmd'],
        #microbial_genomes=config['supporting_files']['files']['microbial']['reprDB']['name'],
    threads: config['svaba']['threads']  #{threads}  # eval( "lambda wildcards, input, threads, attempt: attempt, 254 if int(attempt) > 1  else 255 ")
    resources:
        vcpu=config['svaba']['threads'],
        partition=config['svaba']['partition'],
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.svb.svabavcf.bench.tsv"
    log:
        MDIR + "{sample}/align/{alnr}/sv/svaba/logs/{sample}.{alnr}.svaba.sv.vcf.log",
    conda:
        "../envs/svaba_v0.1.yaml"
    shell:
        """
        export DBSNP={params.dbsnp};
        rm -rf $(dirname {output.stub} ) || echo rmFailed;
        mkdir -p $(dirname {output.stub} );
        mkdir -p $(dirname {log});
        touch {output.stub};
        {params.ld_p} svaba run  -D {params.dbsnp} \
        --germline  -t {input.bam} -p {threads} \
        -k {params.svaba_calling_regions} \
        -R {params.simple_repeats_bed} \
        -x 30000 -L 5 -I -a {output.stub} \
        -G {params.huref} >> {log} 2>&1;

        perl -pi -e 's/(^.CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t)(.*\/)(.+)(\..*\.mrkdup\..*)/$1$3/g;' {output.vcff} >> {log} 2>&1 ;
        # Should put a fix for the header to define the PL tag as numeric=G..... ugh vcf.
        perl -pi -e 's/(^.CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t)(.*\/)(.+)(\..*\.mrkdup\..*)/$1$3/g;' {output.vcfindel}; >> {log} 2>&1 ;
        sleep 1;
        """

rule svaba_sort_index:
    input:
        vcf = MDIR + "{sample}/align/{alnr}/sv/svaba/{sample}.{alnr}.svaba.sv.vcf",
        vcfindel= MDIR+"{sample}/align/{alnr}/sv/svaba/{sample}.{alnr}.svaba.indel.vcf",
    output:
        sortvcf=MDIR + "{sample}/align/{alnr}/sv/svaba/{sample}.{alnr}.svaba.sv.sort.vcf",
        sortgz = MDIR + "{sample}/align/{alnr}/sv/svaba/{sample}.{alnr}.svaba.sv.sort.vcf.gz",
        sorttbi = MDIR + "{sample}/align/{alnr}/sv/svaba/{sample}.{alnr}.svaba.sv.sort.vcf.gz.tbi",
        isortvcf = MDIR + "{sample}/align/{alnr}/sv/svaba/{sample}.{alnr}.svaba.indel.sort.vcf",
        isortgz = MDIR + "{sample}/align/{alnr}/sv/svaba/{sample}.{alnr}.svaba.indel.sort.vcf.gz",
        isorttbi = MDIR + "{sample}/align/{alnr}/sv/svaba/{sample}.{alnr}.svaba.indel.sort.vcf.gz.tbi"
    threads: 16 if 'svaba' not in config else int(config["svaba"]["threads"])
    priority: 39
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.svaba.sv.vcf.sort.bench.tsv"
    log:
        MDIR + "{sample}/align/{alnr}/sv/svaba/logs/{sample}.{alnr}.svaba.sv.vcf.sort.log",
    conda:
        "../envs/vanilla_v0.1.yaml"
    resources:
        vcpu=16
    params:
        cluster_sample=ret_sample,
        ld_p=config['malloc_alt']['ld_preload'] if 'ld_preload' not in config['svaba'] else config['svaba']['ld_preload'],
    shell:
        """
        set +euo pipefail ;
        (  (rm -rf {output}) || echo rmFailok ; \
        {params.ld_p} bedtools sort -header -i {input.vcf} > {output.sortvcf}  ; \
        bgzip -f -@ {threads} {output.sortvcf}; \
        touch  {output.sortvcf}; \
        tabix -f -p vcf {output.sortgz} ; \
        {params.ld_p} bedtools sort -header -i {input.vcfindel} > {output.isortvcf}  ; \
        bgzip -f -@ {threads} {output.isortvcf}; \
        touch {output.isortvcf}; \
        tabix -f -p vcf {output.isortgz} ; \
        {latency_wait} ) > {log};
        """

localrules: produce_svaba,

rule produce_svaba:  # TARGET: just generate svaba vcfs
    priority: 39
    input:
        expand(MDIR + "{sample}/align/{alnr}/sv/svaba/{sample}.{alnr}.svaba.sv.sort.vcf.gz.tbi",sample=SSAMPS,alnr=ALIGNERS)
