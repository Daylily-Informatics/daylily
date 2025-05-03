import sys
import os

##### OCTOPUS- OUR snv CALLER
# ---------------------------
#
# One of the most important parts of the pipeline
# also one of the most complicated.
# Great docs can be found here: https://luntergroup.github.io/octopus/
# paper: https://www.nature.com/articles/s41587-021-00861-3

config["snv_pos_samps"] = {}



def get_ploidy(wildcards):

    ploidy_str = " "
    try:
        if config["sample_info"][wildcards.sample]["biological_sex"].lower() in [
            "male",
            "m",
        ]:
            ploidy_str = ploidy_str + " --contig-ploidies X=1 Y=1 "
        elif config["sample_info"][wildcards.sample]["biological_sex"].lower() in [
            "female",
            "f",
        ]:
            ploidy_str = ploidy_str + " --contig-ploidies X=2 Y=0 "
        else:
            ploidy_str = ploidy_str + " --contig-ploidies X=2 Y=1 "
    except Exception as e:
        print(e, file=sys.stderr)

    return " "  # ploidy_str




def get_ochrm_mod(wildcards):
    pchr=GENOME_CHR_PREFIX
    ret_str = ""
    sl = wildcards.ochrm.split("-")
    sl2 = wildcards.ochrm.split("~")


    #from IPython import embed
    #embed()
    #raise

    if len(sl2) == 2:
        ret_str = pchr + wildcards.ochrm
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
            "oct chunks can only be one contiguous range per chunk : ie: 1-4 with the non numerical chrms assigned 23=X, 24=Y, 25=MT"
        )
    ret_str = ret_mod_chrm(ret_str)

    return ret_str



def ret_oct_clust_samp(wildcards):
    return wildcards.sample + "_" + wildcards.ochrm


rule octopus:
    """https://github.com/luntergroup/octopus"""
    input:
        c=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.cram",
        crai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.cram.crai",
        d=MDIR + "{sample}/align/{alnr}/snv/oct/vcfs/{ochrm}/{sample}.ready",
    output:
        #oct_tmpd=directory(MDIR + "{sample}/align/{alnr}/snv/oct/vcfs/{ochrm}/oct_tmp"),
        vcf=MDIR
        + "{sample}/align/{alnr}/snv/oct/vcfs/{ochrm}/{sample}.{alnr}.oct.{ochrm}.snv.vcf",
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/oct/log/vcfs/{sample}.{alnr}.oct.{ochrm}.snv.log",
    threads: config['octopus']['threads']
    container:
        "docker://daylilyinformatics/octopus-skylake:0.7.4" 
    priority: 45
    benchmark:
            MDIR + "{sample}/benchmarks/{sample}.{alnr}.oct.{ochrm}.bench.tsv"
    resources:
        vcpu= config['octopus']['threads'],
        attempt_n=lambda wildcards, attempt:  (attempt + 0),
        partition="i192,i192mem",
        threads= config['octopus']['threads']
    params:
        cluster_sample=ret_sample, 
        ochrm_mod=get_ochrm_mod,
        anno=config["octopus"]["anno"],
        addl_options=" " if "addl_options" not in config["octopus"] else config["octopus"]["addl_options"],
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        skr=config['supporting_files']['files']['ucsc']['build_gaps']['name'],
        ld_pre=config['octopus']['ld_pre'],
        mdir=MDIR,
        mito_code="MT" if "b37" == config['genome_build'] else "M",
        chrm_prefix=GENOME_CHR_PREFIX,
    shell:
        """ 
        touch {output.vcf};
        timestamp=$(date +%Y%m%d%H%M%S);
        TMPDIR=./octo_tmp_$timestamp;
        APPTAINER_HOME=$TMPDIR;
        trap "sleep 2 && rm -rf \"$TMPDIR\" || echo '$TMPDIR rm fails' >> {log} 2>&1" EXIT;
        
        oochrm_mod=$(echo '{params.ochrm_mod}' | sed 's/~/\:/g' | perl -pe 's/(^23| 23)/ X/g;' | perl -pe 's/(^24| 24)/ Y/g;' | perl -pe 's/(^25| 25)/ {params.mito_code}/g;');

        {params.ld_pre} /opt/octopus/bin/octopus -T $oochrm_mod --threads {threads}    \
        --reference {params.huref}  \
        --temp-directory $TMPDIR \
        --reads {input.c}   \
        --annotations {params.anno}   \
        --skip-regions-file {params.skr}  {params.addl_options} > {output.vcf};
        
        """


rule oct_sort_index_chunk_vcf:
    input:
        vcf=MDIR
        + "{sample}/align/{alnr}/snv/oct/vcfs/{ochrm}/{sample}.{alnr}.oct.{ochrm}.snv.vcf",
    priority: 46
    output:
        vcfsort=temp(MDIR
        + "{sample}/align/{alnr}/snv/oct/vcfs/{ochrm}/{sample}.{alnr}.oct.{ochrm}.snv.sort.vcf"),
        vcfgz=MDIR
        + "{sample}/align/{alnr}/snv/oct/vcfs/{ochrm}/{sample}.{alnr}.oct.{ochrm}.snv.sort.vcf.gz",
        vcftbi=MDIR
        + "{sample}/align/{alnr}/snv/oct/vcfs/{ochrm}/{sample}.{alnr}.oct.{ochrm}.snv.sort.vcf.gz.tbi",
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/oct/vcfs/{ochrm}/log/{sample}.{alnr}.oct.{ochrm}.snv.sort.vcf.gz.log",
    resources:
        vcpu=8,
        threads=8,
        partition="i192,i192mem"
    params:
        cluster_sample=ret_sample,
    threads: 8 #config["config"]["sort_index_oct_chunk_vcf"]['threads']
    shell:
        """
        bedtools sort -header -i {input.vcf} > {output.vcfsort} 2>> {log};
        
        bgzip {output.vcfsort} >> {log} 2>&1;
        touch {output.vcfsort};

        tabix -f -p vcf {output.vcfgz} >> {log} 2>&1;
        
        """


localrules:
    oct_concat_fofn,


rule oct_concat_fofn:
    input:
        chunk_tbi=sorted(
            expand(
                MDIR
                + "{{sample}}/align/{{alnr}}/snv/oct/vcfs/{ochm}/{{sample}}.{{alnr}}.oct.{ochm}.snv.sort.vcf.gz.tbi",
                ochm=OCTO_CHRMS,            ),            key=lambda x: float(                str(x.replace("~", ".").replace(":", "."))               .split("vcfs/")[1]                .split("/")[0]                .split("-")[0]            ),        ),
    # This expand pattern is neat.  the escaped {} remain acting as a snakemake wildcard and expect to be derived from the dag, while th ochrm wildcard is effectively being constrained by the values in the OCTO_CHRMS array;  So you produce 1 input array of files for every sample+ochrm parir, with one list string/file name per array.  The rule will only begin when all array members are produced. It's then sorted by first octochrm so they can be concatenated w/out another soort as all the chunks had been sorted already.
    priority: 44
    output:
        fin_fofn=MDIR
        + "{sample}/align/{alnr}/snv/oct/{sample}.{alnr}.oct.snv.concat.vcf.gz.fofn",
        tmp_fofn=MDIR        + "{sample}/align/{alnr}/snv/oct/{sample}.{alnr}.oct.snv.concat.vcf.gz.fofn.tmp",
    threads: 2
    resources:
        threads=2
    params:
        fn_stub="{sample}.{alnr}.oct.",
        cluster_sample=ret_sample,
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.oct.concat.fofn.bench.tsv"
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR + "{sample}/align/{alnr}/snv/oct/log/{sample}.{alnr}.oct.cocncat.fofn.log",
    shell:
        """

        for i in {input.chunk_tbi}; do
            ii=$(echo $i | perl -pe 's/\.tbi$//g'; );
            echo $ii >> {output.tmp_fofn};
        done;
        (workflow/scripts/sort_concat_chrm_list.py {output.tmp_fofn} {wildcards.sample}.{wildcards.alnr}.oct. {output.fin_fofn}) >> {log} 2>&1 ;

        """


rule oct_concat_index_chunks:
    input:
        fofn=MDIR
        + "{sample}/align/{alnr}/snv/oct/{sample}.{alnr}.oct.snv.concat.vcf.gz.fofn",
    output:
        vcfgz=touch(
            MDIR + "{sample}/align/{alnr}/snv/oct/{sample}.{alnr}.oct.snv.sort.vcf.gz"
        ),
        vcfgztemp=temp(
            MDIR + "{sample}/align/{alnr}/snv/oct/{sample}.{alnr}.oct.snv.sort.temp.vcf.gz"
        ),
        vcfgztbi=touch(
            MDIR
            + "{sample}/align/{alnr}/snv/oct/{sample}.{alnr}.oct.snv.sort.vcf.gz.tbi"
        ),
    threads: 8
    resources:
        vcpu=8,
        threads=8,
        partition="i192,i192mem"
    priority: 47
    params:
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        cluster_sample=ret_sample,
    resources:
        attempt_n=lambda wildcards, attempt:  (attempt + 0)
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.oct.merge.bench.tsv"
    conda:
        config["octopus"]["oct_gather_env"]
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/oct/log/{sample}.{alnr}.oct.snv.merge.sort.gatherered.log",
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
    clear_combined_octovcf,


rule clear_combined_octovcf:  # TARGET:  clear combined octo vcf so the chunks can be re-evaluated if needed.
    input:
        expand(
            MDIR + "{sample}/align/{alnr}/snv/oct/{sample}.{alnr}.oct.snv.sort.vcf.gz",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
    priority: 42
    threads: 2
    shell:
        "(rm {input}*  1> /dev/null  2> /dev/null ) || echo 'file not found for deletion: {input}';"


localrules:
    produce_oct_vcf,


rule produce_oct_vcf:  # TARGET: octopus vcf
    input:
        expand(
            MDIR
            + "{sample}/align/{alnr}/snv/oct/{sample}.{alnr}.oct.snv.sort.vcf.gz.tbi",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
    output:
        "gatheredall.oct",
    priority: 48
    threads: 2
    log:
        "gatheredall.oct.log",
    shell:
        """( touch {output} ;

        {latency_wait}; ls {output} ) >> {log} 2>&1;
        """


localrules:
    oct_prep_chunkdirs,


rule oct_prep_chunkdirs:
    input:
        c=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.cram",
        i=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.cram.crai",
    output:
        expand(
            MDIR + "{{sample}}/align/{{alnr}}/snv/oct/vcfs/{ochrm}/{{sample}}.ready",
            ochrm=OCTO_CHRMS,
        ),
    threads: 2
    log:
        MDIR + "{sample}/align/{alnr}/snv/oct/logs/{sample}.{alnr}.chunkdirs.log",
    shell:
        """
        ( echo {output}  ;
        mkdir -p $(dirname {output} );
        touch {output};
        ls {output}; ) > {log} 2>&1;
        """
 
