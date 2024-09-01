# ###### GOLEFT
#
# coverage metrics calculation tools.
# goleft
# github: https://github.com/brentp/goleft
# paper:  https://academic.oup.com/gigascience/article/6/11/gix090/4160383


rule goleft:
    input:
        bam=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
        bai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam.bai",
    output:
        done=touch(MDIR + "{sample}/align/{alnr}/alignqc/goleft.done"),
        donetwo=touch(MDIR + "{sample}/align/{alnr}/alignqc/goleft/golefttwo.done"),
    params:
        cluster_sample=ret_sample,
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.goleft.bench.tsv"
    log:
        MDIR + "{sample}/align/{alnr}/alignqc/goleft/logs/goleft.log",
    threads: config["go_left"]["threads"]
    conda:
        config["go_left"]["env_yaml"]
    shell:
        """
        set +euo pipefail;                                                                                                                                 
        rm -rf $(dirname {output.donetwo} ) || echo rmGOLfailed ;
        mkdir -p $(dirname {output.donetwo} )/logs ; 
        echo setlog > {log};
        echo 'goingleft' >> {log} 2>&1;
        export gl=$(dirname {output.donetwo} ) ;
        goleft indexcov --directory $gl {input.bam} >> {log} 2>&1;
        touch {output.done}; touch {output.donetwo};
        {latency_wait};  ls {output};
        """
