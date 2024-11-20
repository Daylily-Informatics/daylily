import os

# #### seqfu
# -----------
#   Super efficent fasta/q manipulations
#
# github:  https://github.com/telatin/seqfu2
# paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8148589/


rule seqfu:
    input:
        #DR=MDIR + "{sample}/{sample}.dirsetup.ready",
        #fqs=get_raw_fastqs,
        #f1=getR1sS,  # method defined in fastp.smk
        #f2=getR2sS,  # method defined in fastp.smk
        f1=get_raw_R1s,
        f2=get_raw_R2s,
    output:
        mqc_r1=f"{MDIR}" + "{sample}/seqqc/seqfu/{sample}.seqfuR1.mqc.tsv",
        mqc_r2=f"{MDIR}" + "{sample}/seqqc/seqfu/{sample}.seqfuR2.mqc.tsv",
        sent=MDIR + "{sample}/seqqc/seqfu/{sample}.seqfu.done",
    benchmark:
        f"{MDIR}" + "{sample}/benchmarks/{sample}.seqfu.bench.tsv"
    threads: config["seqfu"]["threads"]
    params:
        cluster_sample=ret_sample,
        ld_preload=" "
        if "ld_preload" not in config["malloc_alt"]
        else config["malloc_alt"]["ld_preload"],
        ld_pre=" "
        if "ld_preload" not in config["seqfu"]
        else config["seqfu"]["ld_preload"],
    log:
        f"{MDIR}" + "{sample}/seqqc/seqfu/{sample}.seqfu.log",
    conda:
        config["fastqc"]["env_yaml"]
    shell:
        """
        ( cat  <(zcat {input.f1} ) | env {params.ld_preload} seqfu stats --nice -b  --verbose --multiqc ./{output.mqc_r1} - &
        cat <(zcat {input.f2} ) | env {params.ld_preload} seqfu stats --nice -b  --verbose --multiqc ./{output.mqc_r2} - &
        wait;
        touch {output.sent};) > {log}
        {latency_wait}; ls {output};
        """


localrules:
    compile_seqfu,


rule compile_seqfu:
    input:
        expand(MDIR + "{sample}/seqqc/seqfu/{sample}.seqfu.done", sample=SAMPS),
    container:
        None
    output:
        mqc1=MDIR + "other_reports/seqfu1.mqc.tsv",
        mqc2=MDIR + "other_reports/seqfu2.mqc.tsv",
        d=MDIR + "logs/seqfu.done",
    shell:
        """mkdir -p {MDIR}other_reports;
        single_file=$( find results | grep seqfuR1.mqc.tsv | head -n 1);
        if [[ "$single_file" == "" ]]; then
            echo "NO DATA FOUND" > {output.mqc1};
        else
            head -n 35 $single_file > {output.mqc1};
            find results | grep seqfuR1.mqc.tsv | parallel 'tail -n 1 {{}} >> {output.mqc1}';
        fi;

        single_file2=$( find results | grep seqfuR2.mqc.tsv  | head -n 1);
        if [[ "$single_file2" == "" ]]; then
            echo "NO DATA FOUND" > {output.mqc2};
        else
            head -n 35 $single_file2 > {output.mqc2};
            find results | grep seqfuR2.mqc.tsv  | parallel 'tail -n 1 {{}} >> {output.mqc2}';
        fi;
        touch {output.d};

        {latency_wait};
        """


localrules:
    produce_seqfu,


rule produce_seqfu:  # TARGET: seqfu output
    input:
        MDIR + "logs/seqfu.done",
    container:
        None
