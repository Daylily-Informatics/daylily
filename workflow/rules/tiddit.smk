0#### TIDDIT IS A MORE DAYERN sv CALLER
# -------------------------------------
# It leverages snv calls and some population
# freq data to do it's work. It is a
# nice complement to Manta
#
# github: https://github.com/SciLifeLab/TIDDIT
# paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5521161/


rule tiddit:
    input:
        bamo=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
        bai=MDIR +  "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam.bai",
    output:
        stub = MDIR + "{sample}/align/{alnr}/sv/tiddit/{sample}.{alnr}.tiddit.sv",
        vcf = MDIR + "{sample}/align/{alnr}/sv/tiddit/{sample}.{alnr}.tiddit.sv.vcf",
        #vcfgz = touch(MDIR+ "{sample}/align/{alnr}/sv/tiddit/{sample}.{alnr}.tiddit.sv.sort.vcf.gz"),
    params:
        huref=config["supporting_files"]["files"]["huref"]["ref"]["name"],
        min_sv_size=config["tiddit"]["min_sv_size"],
        cluster_sample=ret_sample_alnr,
        ld_p=config['malloc_alt']['ld_preload'] if 'ld_preload' not in config['tiddit'] else config['tiddit']['ld_preload'],
    threads: 16
    resources:
        vcpu=1,
        partition="i16-5,i32-5,i64-5,i96-5",
        threads=16
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.tiddit.sv.vcf.bench.tsv"
    log:
        MDIR + "{sample}/align/{alnr}/sv/tiddit/logs/{sample}.{alnr}.tiddit.sv.vcf.log",
    container:
        "docker://daylilyinformatics/tiddit:2.12.0"
    shell:
        """
        export APPTAINER_HOME=/fsx/scratch;
        set +euo pipefail;
        rm -rf {output.stub}* || echo rmFailedTiddit;  # verging on overkill cleanup for restarts
        mkdir -p "$( dirname {output.vcf} )/logs"  ;
        touch {output.stub};
        TIDDIT.py --help >> {log};
        echo TheFileWasCreated > {output.stub};
        TIDDIT.py  --sv --bam {input.bamo} -z {params.min_sv_size} -o {output.stub} --ref {params.huref} >> {log} ;
        touch {output};
        ls {output};  >> {log} ;
        """


rule tiddit_sort_index:
    input:
        MDIR + "{sample}/align/{alnr}/sv/tiddit/{sample}.{alnr}.tiddit.sv.vcf"
    output:
        sortvcf = touch(MDIR + "{sample}/align/{alnr}/sv/tiddit/{sample}.{alnr}.tiddit.sv.sort.vcf"),
        sortgz = touch(MDIR + "{sample}/align/{alnr}/sv/tiddit/{sample}.{alnr}.tiddit.sv.sort.vcf.gz"),
        sorttbi = touch(MDIR + "{sample}/align/{alnr}/sv/tiddit/{sample}.{alnr}.tiddit.sv.sort.vcf.gz.tbi"),
    threads: config["tiddit"]["threads"]
    priority: 8
    resources:
        vcpu=8,
        partition="i4-5,i16-5,i32-5,i64-5,i96-5",
        threads=8
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.tiddit.sv.vcf.sort.bench.tsv"
    log:
        MDIR + "{sample}/align/{alnr}/sv/tiddit/logs/{sample}.{alnr}.tiddit.sv.vcf.sort.log",
    conda:
        "../envs/vanilla_v0.1.yaml"
    params:
        cluster_sample=ret_sample,
    shell:
        """
        set +euo pipefail;
        perl -pi -e 's/\\t\\*\\t/\\tN\\t/g;' {input} >> {log} 2>&1 || echo 'nomatch!' ;
        bedtools sort -header -i {input} > {output.sortvcf};
        bgzip -f -@ {threads} {output.sortvcf};
        touch {output.sortvcf};
        tabix -p vcf -f {output.sortgz};
        {latency_wait} || echo passOn;
        ls {output} || echo passOn ;
        """



localrules: produce_tiddit,

rule produce_tiddit:  # TARGET: Produce All Tiddit
    priority: 45
    input:
        expand(MDIR +"{sample}/align/{alnr}/sv/tiddit/{sample}.{alnr}.tiddit.sv.sort.vcf.gz.tbi", sample=SSAMPS, alnr=ALIGNERS)
