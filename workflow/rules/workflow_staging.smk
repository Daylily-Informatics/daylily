###### Workflow Staging
# ---------------------
# Bookeeping tasks to get done before anyhting really start to run.
# Setting up needed ephemeral directories and temp/log dirs.
# By makling this part of your workflow, if you re-un a w/f
# several times this rule will likely have succeeded and not run again

if "this_is_the_template" in config:
    if int(config["this_is_the_template"]) == 1:
        raise Exception(
            "\n\n\tThe config file template sitll has the safeguard key present 'this_is_the_template' -- please delete it if you have modified the sample sheet.\n\n"
        )


localrules:
    workflow_staging,
    yield_ref,


rule workflow_staging:
    input:
        "logs/supporting_data_staging.done",
    output:
        "logs/workflow_staging.done",
    threads: 1
    params:
        cluster_sample=f"{RU[0]}_{EX[0]}",
        res=config["dirs"]["operational"]["results"],
        ops_log=config["dirs"]["operational"]["logs"],
        tmp=config["dirs"]["operational"]["tmp"],
        mrd=config["dirs"]["day"]["reports"],
        mdir=MDIR,
        oreport=MDIR+"other_reports",
    log:
        "logs/workflow_staging.log",
    conda:
        config["vanilla"]["env_yaml"]
    shell:
        """(mkdir -p {params.res} && mkdir -p {params.ops_log} && mkdir -p {params.tmp} && mkdir -p {params.mrd} && mkdir -p {params.oreport} ) > {log} 2>&1;"""
        """mkdir -p {params.mdir}other_reports;"""
        """touch {output}"""


rule yield_ref:
    input:
        "logs/workflow_staging.done",
    output:
        "human_g1k_v37_modified.fasta",
    conda:
        config["vanilla"]["env_yaml"]
    threads: 1
    shell:
        "touch {output};"

