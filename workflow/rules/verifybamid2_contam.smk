######### CONTAMINATION SCREEN
# ----------------------------
# This is a population agnostic contamination screening tool that can
# operate on single sample or multi sample BAM files.
# github: http://griffan.github.io/VerifyBamID/
# paper: https://doi.org/10.1101/gr.246934.118
#   2020. “Ancestry-agnostic estimation of DNA sample
#   contamination from sequence reads.” Genome Research.


rule verifybamid2_contam:
    input:
        b=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
    output:
        vb_prefix=MDIR + "{sample}/align/{alnr}/alignqc/contam/vb2/{sample}.{alnr}.vb2",
        vb_tsv=MDIR + "{sample}/align/{alnr}/alignqc/contam/vb2/{sample}.{alnr}.vb2.tsv",
    log:
        MDIR + "{sample}/align/{alnr}/alignqc/contam/vb2/logs/{sample}.{alnr}.vb2.log",
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.vb2.bench.tsv"
    conda:
        config["verifybamid2_contam"]["env_yaml"]
    threads: config["verifybamid2_contam"]["threads"]
    resources:
        vcpu=config["verifybamid2_contam"]["threads"]
    params:
        cluster_sample=ret_sample,
        huref=config["supporting_files"]["files"]["octopus"]["huref"]["name"],
        db_prefix=config["supporting_files"]["files"]["verifybam2"]["dat_files"]["name"],
    shell:
        """
        set +euo pipefail;
        rm -rf $(dirname {output.vb_tsv} ) || echo rmVerifyBAMfailed;
        mkdir -p $(dirname {output.vb_tsv} )/logs;
        verifybamid2 --BamFile {input} --Output {output.vb_prefix} --DisableSanityCheck --SVDPrefix {params.db_prefix} --Reference {params.huref} ;
        touch  {output.vb_prefix}.selfSM {output.vb_tsv};
        cp {output.vb_prefix}.selfSM {output.vb_tsv};
        touch {output.vb_prefix};
        export samplename=$(echo $(basename {output.vb_tsv}) | cut -d'.' -f 1);
        perl -pi -e 's/^x\t/$samplename\t/g;' {output.vb_prefix}.selfSM;
        perl -pi -e 's/^x\t/$samplename\t/g;' {output.vb_tsv};
        {latency_wait};
        ls {output};
        """

localrules:
    produce_contam_estimate,


rule produce_contam_estimate:  # TARGET:  jusg gen contam
    input:
        expand(
            MDIR + "{sample}/align/{alnr}/alignqc/contam/vb2/{sample}.{alnr}.vb2.tsv",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
