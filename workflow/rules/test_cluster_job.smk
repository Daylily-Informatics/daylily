import datetime
import os

# ###### Test rule to confirm various cluster features work
# ---------------------------------------------------------
#


def test_func(frm, log):
    os.system(f"echo {frm} >> {log}")


rule i_am_a_test:
    input:
        "iamatest_{middle}_iamatest",
    output:
        touch("iamatest_{middle}.complete"),
    conda:
        config["vanilla"]["env_yaml"]
    log:
        a="LOGG/{middle}/iamatest_-_{middle}_.log",
        b=directory("LOGG/{middle}/cluster_logs/"),
    params:
        A="SOMETHING-{middle}",
        C="ATHIRDTHING-{wildcards.middle}",
        cluster_sample=f"wildcards.sample",
    shell:
        "echo 'A:{params.A}' &> {log};"
        "echo 'C:{params.C}' &>> {log};"
        "echo 'D:{params.D}' &>> {log};"
        "echo 'R:{rule}' &>> {log};"
        '(echo "SGE:$HOSTNAME $SGE_TASK_ID $JOB_ID - SM X &>> {log} || echo Z &>> {log});'


localrules:
    the_end,


rule the_end:  # TARGET : interal testing
    input:
        expand("iamatest_{a}.complete", a=["XXX", "QQQ", "RRR"]),
