#### Paragraph - SV caller  using novel graph assembly  centric approach
# -------------------------------------
#  They claim in their paper : " Paragraph achieves a recall of 0.86 and a precision of 0.91."
#   Which is phnomenal performance, depending on the details of how they frame the calulation of these numbers.
#
# github: https://github.com/Illumina/paragraph
# paper: https://www.biorxiv.org/content/10.1101/635011v2.full
#    https://github.com/Illumina/paragraph/blob/master/doc/graph-tools.md#Genotyper

def get_bam_depth(wildcards):
    ## TODO, at time of bam generation, save the aligner-BAM-medcov  to samplesheet for use downstream, ie: here
    samp_alnr_depth=1
    #samp_alnr_depth = samples.loc("wildcards.sample", "wildcards.alnr-BAM-medcov")[0]
    return samp_alnr_depth

# this will end up playing the role of SV caller as well as like duphpld a sv refiner
#  this is very crude calling
rule paragraph:
    input:
        vcf=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
    output:
        vcf= MDIR + "{sample}/align/{alnr}/sv/para/{sample}.{alnr}.para.sv.vcf",
        vcfsort = MDIR + "{sample}/align/{alnr}/sv/para/{sample}.{alnr}.para.sv.sort.vcf",
        vcfsortgz = MDIR + "{sample}/align/{alnr}/sv/para/{sample}.{alnr}.para.sv.sort.vcf.gz",
        vcfsortgztbi = MDIR + "{sample}/align/{alnr}/sv/para/{sample}.{alnr}.para.sv.sort.vcf.gz.tbi",
        idxjson=MDIR + "{sample}/align/{alnr}/sv/para/{sample}.{alnr}.para.sv.idx.json",
        manifest=MDIR + "{sample}/align/{alnr}/sv/para/{sample}.{alnr}.para.sv.manifest",
    params:
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        cluster_sample=ret_sample_alnr,
        bam_depth=get_bam_depth,
        read_len="150",
        cmd="paragraph",
    threads: 44 if 'paragraph' not in config else config["paragraph"]["threads"]
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.para.sv.vcf.bench.tsv"
    log:
        MDIR + "{sample}/align/{alnr}/sv/para/logs/{sample}.{alnr}.para.sv.vcf.log",
    conda:
        config['paragraph']['env_yaml']
    shell:
        """
        (rm -rf $(dirname {output.vcf} ) || echo noFail;  # This is a quick example of a simple case may wish to explore fine tuning of VCFs capabilities
        work_dir="$(dirname {output.vcf} )/";
        mkdir -p $(dirname {log} );
        idxdepth -b {input.bam} -r {params.huref} -o {output.idxjson};
        echo "id\tpath\tread length\tdepth" > {output.manifest};
        echo "{params.cluster_sample}\t{input.bam}\t{params.read_len}\t{params.bam_mi_depth}" >> {output.manifest};
        multigrmpy.py -i {input.vcfgz} -m {output.manifest} -r {params.huref} -o $work_dir;
        tabix -p vcf -f {output.work_dir}/genotypes.vcf.gz;
        tabix -p vcf -f {output.work_dir}/variants.vcf.gz;
        (grep ^"#" <(unpigz -q -c -- input.vcf); grep -v ^"#" <(unpigz -q -c -- input.vcf ) | sort -k1,1 -k2,2n) > {output.vcfsort};
        bgzip -f {output.vcfsort};
        touch {output.vcfsort}
        tabix -p vcf -f {output.vcfsortgz};
        {latency_wait} ;
        ls {output};  >> {log} ;
        """


#rule para_sort_index:
#    input:
#        MDIR + "{sample}/align/{alnr}/sv/para/{sample}.{alnr}.para.sv.vcf"
#    output:
#        sortvcf = touch(MDIR + "{sample}/align/{alnr}/sv/para/{sample}.{alnr}.para.sv.sort.vcf"),
#        sortgz = touch(MDIR + "{sample}/align/{alnr}/sv/para/{sample}.{alnr}.para.sv.sort.vcf.gz"),
#        sorttbi = touch(MDIR + "{sample}/align/{alnr}/sv/para/{sample}.{alnr}.para.sv.sort.vcf.gz.tbi"),
#    threads: config["tiddit"]["threads"] if 'paragraph' not in config else config['para']['threads']
#    benchmark:
#        MDIR + "{sample}/benchmarks/{sample}.{alnr}.para.sv.vcf.sort.bench.tsv"
#    log:
#        MDIR + "{sample}/align/{alnr}/sv/para/logs/{sample}.{alnr}.para.sv.vcf.sort.log",
#    conda:
#        "../envs/bedtools_v0.1.yaml"
#    params:
#        cluster_sample=ret_sample,
#        bam_depth=grab_bam_depth,
# #   shell:
#        """
#        bedtools sort -header -i {input} > {output.sortvcf};
#        bgzip -f -@ {threads} {output.sortvcf};
#        touch {output.sortvcf};
#        tabix -p vcf -f {output.sortgz};
#        {latency_wait} || echo passOn;
#        ls {output} || echo passOn ;
#        """



localrules: produce_paragraph,

rule produce_paragraph:  # TARGET: Produce All Tiddit
    input:
        expand(MDIR +"{sample}/align/{alnr}/sv/para/{sample}.{alnr}.para.sv.sort.vcf.gz.tbi", sample=SSAMPS, alnr=ALIGNERS)
