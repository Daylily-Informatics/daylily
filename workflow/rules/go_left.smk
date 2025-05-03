import os
# ###### GOLEFT
#
# coverage metrics calculation tools.
# goleft
# github: https://github.com/brentp/goleft
# paper:  https://academic.oup.com/gigascience/article/6/11/gix090/4160383


if os.environ.get('DAY_CRAM','') =='':

    rule goleft:
        input:
            bam=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
            bai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam.bai",
        output:
            done=touch(MDIR + "{sample}/align/{alnr}/alignqc/goleft.done"),
            donetwo=touch(MDIR + "{sample}/align/{alnr}/alignqc/goleft/golefttwo.done"),
        params:
            cluster_sample=ret_sample,
            sexchrms="chrX,chrY" if os.environ.get('DAY_GENOME_BUILD','') == 'hg38' else "X,Y",
        benchmark:
            MDIR + "{sample}/benchmarks/{sample}.{alnr}.goleft.bench.tsv"
        resources:
            vcpu=config["go_left"]["threads"],
            partition=config["go_left"]["partition"],
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
            goleft indexcov --sex {params.sexchrms} --directory $gl {input.bam} >> {log} 2>&1;
            touch {output.done}; touch {output.donetwo};
            {latency_wait};  ls {output};
            """

else:
    rule goleft_cram:
        input:
            cram=MDIR + "{sample}/align/{alnr}/{sample}.cram",
            crai=MDIR + "{sample}/align/{alnr}/{sample}.cram.crai",
        output:
            done=touch(MDIR + "{sample}/align/{alnr}/alignqc/goleft.done"),
            donetwo=touch(MDIR + "{sample}/align/{alnr}/alignqc/goleft/golefttwo.done"),
        params:
            cluster_sample=ret_sample,
            sexchrms="X,Y" if os.environ.get('DAY_GENOME_BUILD','') == 'b37' else "chrX,chrY",
            huref=config["supporting_files"]["files"]["huref"]["broad_fasta"]["name"],
        benchmark:
            MDIR + "{sample}/benchmarks/{sample}.{alnr}.goleft.bench.tsv"
        resources:
            vcpu=config["go_left"]["threads"],
            partition=config["go_left"]["partition"],
        log:
            MDIR + "{sample}/align/{alnr}/alignqc/goleft/logs/goleft.log",
        threads: config["go_left"]["threads"]
        conda:
            config["go_left"]["env_yaml"]
        shell:
            """
            #set +euo pipefail;                                                                                                                                 
            
            rm -rf $(dirname {output.donetwo} ) || echo rmGOLfailed ;
            mkdir -p $(dirname {output.donetwo} )/logs ; 
            
            export gl=$(dirname {output.donetwo} ) ;
            export REF_PATH={params.huref};
            goleft indexcov --directory $gl --sex {params.sexchrms} --fai {params.huref}.fai {input.crai} >> {log} 2>&1;
            
            touch {output.done}; touch {output.donetwo};
            
            """
