import os

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
        cram=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.cram",
    output:
        vb_prefix=MDIR + "{sample}/align/{alnr}/alignqc/contam/vb2/{sample}.{alnr}.vb2",
        vb_tsv=MDIR + "{sample}/align/{alnr}/alignqc/contam/vb2/{sample}.{alnr}.vb2.tsv",
        tmppile=temp(MDIR + "{sample}/align/{alnr}/alignqc/contam/vb2/{sample}.{alnr}.tmp.pileups.table"),
        contam=MDIR + "{sample}/align/{alnr}/alignqc/contam/vb2/{sample}.{alnr}.contam.tsv",
        selfSM=MDIR + "{sample}/align/{alnr}/alignqc/contam/vb2/{sample}.{alnr}.vb2.selfSM",
        mqc=MDIR + "{sample}/align/{alnr}/alignqc/contam/vb2/{sample}.{alnr}.vb2_mqc.tsv",
    log:
        MDIR + "{sample}/align/{alnr}/alignqc/contam/vb2/logs/{sample}.{alnr}.vb2.log",
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.vb2.bench.tsv"
    conda:
        config["verifybamid2_contam"]["env_yaml"]
    threads: config["verifybamid2_contam"]["threads"]
    resources:
        vcpu=config["verifybamid2_contam"]["threads"],
        partition=config["verifybamid2_contam"]["partition"],
    params:
        cluster_sample=ret_sample,
        alnr=get_alnr,
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        db_prefix=config["supporting_files"]["files"]["verifybam2"]["dat_files"]["name"],
        chrm_prefix=GENOME_CHR_PREFIX,
    shell:
        """
        set +euo pipefail;
        rm -rf $(dirname {output.vb_tsv} ) || echo rmVerifyBAMfailed;
        mkdir -p $(dirname {output.vb_tsv} )/logs;

        #verifybamid2 --WithinAncestry --BamFile {input.cram} --Output {output.vb_prefix} --DisableSanityCheck \
        #    --SVDPrefix {params.db_prefix}  --NumThread  {threads} --Reference {params.huref}  >> {log} 2>&1 ;
                    
        for chr in {params.chrm_prefix}{{1..22}} {params.chrm_prefix}X; do echo $chr; done | parallel -j {threads} '
            gatk GetPileupSummaries \
            -I {input.cram} \
            -V {params.db_prefix} \
            --reference {params.huref} \
            -O  {output.tmppile}.{{}}.pileups.table.tmp \
            --interval-set-rule INTERSECTION \
            -L {params.db_prefix} \
            -L {{}} > {log}.{{}} 2>&1 ;
        ';


        touch {output.tmppile};
        head -n 2 $(ls $(dirname {output.tmppile} )/*.tmp | head -n 1) > {output.tmppile};
        perl -pi -e 's/(^.*SAMPLE\=)(.*$)/$1{params.cluster_sample}.{params.alnr}/g;' {output.tmppile};
        tail -n +2 -q your_sample.{params.chrm_prefix}*.pileups.table >> {output.tmppile};

        gatk CalculateContamination \
        -I {output.tmppile} \
        -O {output.contam}

        export contam=$(awk 'NR==2 {{print $2}}' {output.contam});

        echo -e "SEQ_ID\tRG\tCHIP_ID\t#SNPS\t#READS\tAVG_DP\tFREEMIX\tFREELK1\tFREELK0\tFREE_RH\tFREE_RA\tCHIPMIX\tCHIPLK1\tCHIPLK0\tCHIP_RH\tCHIP_RA\tDPREF\tRDPHET\tRDPALT" > {output.selfSM};
        echo -e "{params.cluster_sample}.{params.alnr}\tNA\tNA\tNA\tNA\tNA\t$contam\t-1\t-1\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> {output.selfSM};
        cp {output.selfSM} {output.vb_tsv};
        cp {output.selfSM} {output.mqc};
        touch {output.vb_prefix};

        rm  $(dirname {output.tmppile} )/*.tmp

        """

localrules:
    produce_contam_estimate,

rule produce_contam_estimate:  # TARGET:  jusg gen contam
    input:
        expand(
            MDIR + "{sample}/align/{alnr}/alignqc/contam/vb2/{sample}.{alnr}.vb2.tsv",
            sample=SSAMPS,
            alnr=ALL_ALIGNERS,
        ),
