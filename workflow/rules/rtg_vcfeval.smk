import os
import sys



def get_samp_concordance_truth_dir(wildcards):
    cntrl_dir = samples[samples['samp'] == wildcards.sample]["concordance_control_path"][0]
    return cntrl_dir

def get_alt_sample_name(wildcards):
    return samples[samples['samp'] == wildcards.sample]['external_sample_id'][0]

def get_sampn(wildcards):
    return wildcards.sample


def get_snv_caller(wildcards):
    return wildcards.snv

def get_cdir(wildcards):
    ret_d =  MDIR+f"{wildcards.sample}/align/{wildcards.alnr}/snv/{wildcards.snv}/concordance/"

    if ret_d.startswith('/marigo'):
        ret_d = f"results"+ret_d
    return ret_d


def get_in_rtg_vcf(wildcards):
    return f"{MDIR}{wildcards.sample}/align/{wildcards.alnr}/snv/{wildcards.snv}/{wildcards.sample}.{wildcards.alnr}.{wildcards.snv}.snv.sort.vcf.gz"


def get_in_rtg_tbi(wildcards):
    return f"{MDIR}{wildcards.sample}/align/{wildcards.alnr}/snv/{wildcards.snv}/{wildcards.sample}.{wildcards.alnr}.{wildcards.snv}.snv.sort.vcf.gz.tbi"


if len(CONCORDANCE_SAMPLES.keys()) > 0:

    rule prep_for_concordance_check:
        input:
            cvcf=get_in_rtg_vcf,
            ctbi=get_in_rtg_tbi,
        priority: 48
        output:
            s=touch(MDIR + "{sample}/align/{alnr}/snv/{snv}/concordance/concordance.done"),
            fofn=touch(MDIR + "{sample}/align/{alnr}/snv/{snv}/concordance/concordance.fofn"),
            fin_cmds=touch( MDIR + "{sample}/align/{alnr}/snv/{snv}/concordance/concordance.fin.cmds"),
        log:
            MDIR
            + "{sample}/align/{alnr}/snv/{snv}/concordance/logs/{sample}.{alnr}.{snv}.concordance.log",
        benchmark:  MDIR+ "{sample}/benchmarks/{sample}.{alnr}.{snv}.concordance.bench.tsv",
        threads: config['rtg_vcfeval']['threads']
        resources:
            vcpu=config['rtg_vcfeval']['threads'],
            threads=config['rtg_vcfeval']['threads'],
            partition=config['rtg_vcfeval']['partition_other']
        conda:
            config["rtg_vcfeval"]["env_yaml"]
        params:
            tdir=get_samp_concordance_truth_dir,
            alnr=get_alnr,
            snv_caller=get_snv_caller,
            conc_dir=get_cdir,
            cluster_sample=get_sampn,
            ld_p=config['malloc_alt']['ld_preload'] if 'ld_preload' not in config['rtg_vcfeval'] else config['rtg_vcfeval']['ld_preload'],
            l="{",
            r="}",
            sdf=config['supporting_files']['files']['huref']['rtg_tools_genome']['name'],
            sub_threads="7" if 'sub_threads' not in config['rtg_vcfeval'] else config['rtg_vcfeval']['sub_threads'],
            alt_name=get_alt_sample_name, 
        shell:
            """
            # This is a *very* old shameful mess :-) I apologize in advance if you need to debug this before I've killed it off.

            export TMPDIR="/fsx/scratch/";
            mkdir -p $HOME/.parallel;
            parallel --record-env;

            # this needs to be disassembled and re-imagined.  Talk about rube goldbergian....
            set +euo pipefail;
            echo "" > {output.fofn};
            echo "" > {output.fin_cmds};
            ( 
            #IFS=',' read -a vararray <<< $(env LD_PRELOAD=./resources/lib/libgsl.so.25  bcftools query -f "[%DP,]" {input.cvcf} );
            var_sum_dp=na;  # $(IFS=+; echo "$((${params.l}vararray[*]{params.r}))");
            allvar_mean_depth=na; # $(( var_sum_dp  / $(echo ${params.l}#vararray[@]{params.r} ) )) ;
            export allvar_mean_dp=na  #$allvar_mean_depth;
            if [[ "$allvar_mean_depth" == "" ]]; then
                echo THEVARIANT_depth_CALCULATIONHASFAILED;
                export allvar_mean_dp=-1;
                echo sleeping~~~~~~~~~~~;
            fi;

            export aligner={params.alnr};
            export snv={params.snv_caller};

            if [[ "$aligner" == "" ]]; then
                echo "ERROR::: ALIGNER NOT SET, setting to na";
                export aligner="na";
            fi;

            if [[ "$snv" == "" ]]; then
                echo "ERROR::: snv NOT SET, settign to na";
                export snv="na";
            fi;

            rm -rf {output} || echo nothingToDel;
            mkdir -p $( dirname {output.fofn} ) || echo MKDIRfailed ;
            mkdir -p {params.conc_dir}/logs;
            export alt_name={params.alt_name}  ### $(dirname {params.tdir}/{params.alt_name}/. | perl -pe 's/^.*\///g;' );


            # Hack... if the sample entry does not have a concordance dir set, as approximated by the string length of the field being < 6, fake the output files and \
touch a sentinel noting they have been faked
            aconcdir={params.conc_dir}
            if (( ${params.l}#aconcdir{params.r} <= 6 )); then
                echo 'WARNING: concordance is not being run for sample {params.cluster_sample}.' 2>&1;
                touch {output} {output.s}.SKIPPED;
                exit 0;
            fi;

            (ls -d {params.tdir}/* | parallel --halt never --env _ -k --jobs 1 'export vcf="{params.l}{params.r}/$alt_name.vcf.gz";
            export bed="{params.l}{params.r}/$alt_name.bed";
            echo "BEDBRDBED: _$bed _";
            # Retrieve the second-to-last element (-2 in Bash)
            export subd=$(basename $(dirname "$bed"));

            echo "XXXXXXXX-----> _$subd _";
            if [[ "$subd" == "" ]]; then
                exit 33;
            fi;

            echo "BBBDFBDFBDFGDF"; #--use-all-records
            rm -rf {params.conc_dir}/_$subd || sleep 1;  export cmd="`which rtg` vcfeval --decompose --squash-ploidy  --ref-overlap -e $bed -b $vcf -c {input.cvcf} -o {params.conc_dir}/_$subd -t {params.sdf} --threads {params.sub_threads}";
            export  fin_cmd="env python workflow/scripts/parse-vcfeval-summary.py {params.conc_dir}/_$subd/summary.txt {params.cluster_sample} $bed $subd  $alt_name {params.conc_dir}_$subd/{params.cluster_sample}_$subdb_summary.txt $allvar_mean_dp $aligner $snv ";
            ccmd="$cmd >> {params.conc_dir}_a.err; $fin_cmd >> {params.conc_dir}_b.err ; ";
            echo "$ccmd" >> {output.fofn} 2>&1; '
            ) >> {log} 2>&1;
            echo  AlomstDone;
                 if [[ "$allvar_mean_dp" == "" || "$snv" == "" || "$aligner" == "" ]]; then
                echo "WARNING:::: 1 or more of the depth/snv/aligner annotations is null! ... $allvar_mean_dp ... $snv .. $aligner ...\n\n";
            fi;
            ) >> {log} ;  # || echo commandFAILSrtgvcfeval >> {log} 2>&1;
            cat {output.fofn} | env bash ;

            touch {output};
            {latency_wait};
            ls {output};
            """


else:

    localrules: no_concordance_data,

    rule no_concordance_data:
        input:
            MDIR + "{sample}/align/{alnr}/snv/{snv}/{sample}.{alnr}.{snv}.snv.sort.vcf.gz.tbi",
        output:
            MDIR + "{sample}/align/{alnr}/snv/{snv}/concordance/concordance.done",
        shell:
            "touch {output};"


localrules: produce_snv_concordances
rule produce_snv_concordances:  # TARGET:  produce snv concordances
    input:
        expand(
            MDIR + "{sample}/align/{alnr}/snv/{snv}/concordance/concordance.done",
            sample=SSAMPS,
            alnr=ALL_ALIGNERS,
            snv=snv_CALLERS
        )
    priority: 48
    conda:
        "../envs/vanilla_v0.1.yaml"
    params:
        cluster_sample="aggregate",
        mdir=MDIR,
        genome_build=config['genome_build'],
        pc=print_wildcards_etc,
    output:
        touch(MDIR+"other_reports/giab_concordance_mqc.tsv")
    shell:
        """
        set +euo pipefail;
        export wcv=$(find  results/ | grep concord | grep fofn | wc -l);

        (find results/day/{params.genome_build}/*/align/*/snv/*/concordance/ | grep  concordance.mqc  | head -n 1 | parallel 'head -n 1 {{}} > {output}';) || echo 'GetHeaderFAILS' 1>&2;
        (find {params.mdir}*/align/*/snv/*/concordance/ | grep  .mqc | parallel ' tail -n +2 {{}} >> {output}';) || echo "GETCONCORDANCECALLSfails"  1>&2;

        (perl -pi -e 's/^(.+?)(\t)(.+?)(\t)(.+)$/$3\t$1\t$5/g;' {output} ) || echo "perl regsub failed"  1>&2;
	perl -pi -e 's/^([^\t]+?)-None\t([^\t]+)/$1-$2\t$2/g;' {output}  1>&2;

        """
