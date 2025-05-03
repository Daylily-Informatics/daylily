import os

# ##### PEDDY - Pedigree Tools
# ----------------------------
# ** ported from SUPERSONIC
# We are mostly using it for it's gender and ethnicity prediction
# capabilities
# github: https://github.com/brentp/peddy
# paper: http://dx.doi.org/10.1016/j.ajhg.2017.01.017

# ped file:  "family_id individual_id paternal_id maternal_id bio_sex phenotype"
def gen_ped_file(wildcards):
    bio_sex = samples[samples["samp"] == wildcards.sample]["biological_sex"][
        0
    ]  #  sample_info is created in rule_common.smk
    ped_sex = 0
    if bio_sex in ["female"]:
        ped_sex = 2
    elif bio_sex in ["male"]:
        ped_sex = 1
    else:
        ped_sex = 0
    ped_f = f"{MDIR}{wildcards.sample}/align/{wildcards.alnr}/snv/{wildcards.snv}/peddy/{wildcards.sample}.{wildcards.alnr}.{wildcards.snv}.peddy.ped"
    os.system(f"mkdir -p $(dirname {ped_f});")
    ped_fh = open(ped_f, "w")
    ped_fh.write(f"{wildcards.sample}\t{wildcards.sample}\t0\t0\t{ped_sex}\t0\n")
    ped_fh.close()
    return ped_f


rule peddy:
    input:
        vcfgz=MDIR
        + "{sample}/align/{alnr}/snv/{snv}/{sample}.{alnr}.{snv}.snv.sort.vcf.gz",
        # Generalize this for other var callers
        ped_f=gen_ped_file,
    output:
        prefix=MDIR
        + "{sample}/align/{alnr}/snv/{snv}/peddy/{sample}.{alnr}.{snv}.peddy.",
        done=MDIR
        + "{sample}/align/{alnr}/snv/{snv}/peddy/{sample}.{alnr}.{snv}.peddy.done",
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/{snv}/peddy/log/{sample}.{alnr}.{snv}.peddy.log",
    threads: config["peddy"]["threads"]
    resources:
        vcpu=config["peddy"]["threads"],
        partition=config["peddy"]["partition"],
    params:
        cluster_sample=ret_sample,
        ld_preload=config["malloc_alt"]["ld_preload"],
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.{snv}.peddy.bench.tsv"
    container:
        None
    conda:
        config["peddy"]["env_yaml"]
    shell:
        """
        mkdir -p $(dirname {log}) > {log} ;
        echo 'runnning peddy' >> {log} ;
        set +euo pipefail;
        ( ({params.ld_preload} peddy  -p {threads}  --plot --prefix {output.prefix} --loglevel DEBUG {input.vcfgz} {input.ped_f}) || (ls . && sleep 1 && echo 'peddy exited with an error.  If this is a small depth of coverate vs. the genome sample, then this may be a pca failure which can be ignored, run the command by hand to see the detailed error') ) || (sleep 1 && echo 'masking error' ) >> {log};
        set -euo pipefail;
        echo "DONE" >> {log} ;
        echo reallydone > {output.done} ;
        echo reallydone > {output.prefix} ;
        touch {output.done} ;
        touch {output.prefix} ;
        """


localrules:
    produce_peddy,

rule produce_peddy:  # TARGET: just produce peddy results
    input:
        expand(
            MDIR
            + "{sample}/align/{alnr}/snv/{snv}/peddy/{sample}.{alnr}.{snv}.peddy.done",
            sample=SSAMPS,
            alnr=ALL_ALIGNERS,
            snv=snv_CALLERS,
        ),
    output:
        "logs/peddy_gathered.done",
    shell:
        "touch {output};"
