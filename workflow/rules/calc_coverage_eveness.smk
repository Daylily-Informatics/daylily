#  Approximating coverage eveness by norm cov stddev and coofvar
#  Early notes: https://docs.google.com/presentation/d/1b4_gStw5D8MDR89L95wvEuAr9I7Fm1Qn_LVREwNWrcM/edit#slide=id.gcd39261d0a_0_44
#  Math proof asserting NC stddev and/or coovar are good Evenness approximations:
#  paper: Evaluation of the evenness score in next-generation sequencing | Journal of Human Genetics
#     "the general evaluation presented in this paper reveals that in most circumstances the evenness score E of a NGS output can be predicted quite well by the standard deviation σ of the normalized data"
#

rule calc_coverage_evenness:
    input:
        bam=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
        bai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam.bai",
    output:
        # tsv=MDIR  + "{sample}/align/{alnr}/norm_cov_eveness/{sample}.{alnr}.norm_cov_eveness.mqc.tsv",
        mos_pre=MDIR   + "{sample}/align/{alnr}/alignqc/norm_cov_eveness/{sample}.{alnr}.md",
    container: None
    conda:
        "../envs/mosdepth_v0.1.yaml"
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.norm_cov_eveness.bench.tsv"
    threads: config['calc_coverage_evenness']['threads']
    log:
        MDIR + "{sample}/align/{alnr}/alignqc/norm_cov_eveness/logs/norm_cov_eveness.mqc.log",
    params:
        cluster_sample=ret_sample,
        chrm_regions="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22",
        l="{",
        r="}",
        chr="" if config['genome_build'] in ['b37'] else "chr",
    shell:
        """
        set +euo pipefail;
        rm -rf {output.mos_pre}* || echo nothingToDel;
        mkdir -p $( dirname {output.mos_pre} )/logs;
        touch {output.mos_pre};
        echo "Sample\tCHRM\tmeanRawCov\tmedianRawCov\tstdevRawCov\tRawCovCoefofvar\tNCmean\tNCmedian\tstdevNC\tNCcoefofvar\tpctEQ0\tpctLT5\tpctLT10" > {output.mos_pre}.norm_cov_eveness.mqc.tsv;
        for i in {params.l}1..22{params.r};
        do
            echo "Processing {params.cluster_sample} Chrm:{params.chr}$i";
            mosdepth -x  -Q 1 -T 0 -m --by 50 -c {params.chr}$i --threads 20 {output.mos_pre}.{params.chr}$i {input.bam};
            touch {output.mos_pre}.{params.chr}$i.regions.bed.gz;
            Rscript workflow/scripts/calc_norm_cov_sd.R {output.mos_pre}.{params.chr}$i.regions.bed.gz  {params.cluster_sample} {params.chr}$i | sed 's/\\"//g;' >> {output.mos_pre}.norm_cov_eveness.mqc.tsv ;
        done;
        {latency_wait};
        ls {output};
        rm $(dirname {output.mos_pre})/*per-base* || echo 'rm perbase failed';
        """

localrules: produce_cov_uniformity

rule produce_cov_uniformity:  # TARGET: Produce cov eveness calcs, swapping out sambamba for mosdepth
    input:
        expand(MDIR       + "{sample}/align/{alnr}/alignqc/norm_cov_eveness/{sample}.{alnr}.md", sample=SSAMPS, alnr=ALIGNERS)
    container: None
    threads: 8
    output:
        mqc=MDIR+"other_reports/norm_cov_evenness_combo_mqc.tsv",
    shell:
        """
        mkdir -p $(dirname {output});
        single_file=$( find results | grep norm_cov_eveness.mqc.tsv | head -n 1);
        if [[ "$single_file" == "" ]]; then
            echo "NO DATA FOUND" > {output.mqc};
        else
            head -n 1 $single_file > {output.mqc};
            find results | grep .norm_cov_eveness.mqc.tsv | parallel 'tail -n +2 {{}} >> {output.mqc}';
        fi;
        {latency_wait};
        ls {input};
        """
