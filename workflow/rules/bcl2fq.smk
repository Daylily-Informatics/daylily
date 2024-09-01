
import os
import sys
import resource

EXX=[""]
RUU=[""]

#B2FD = config["dirs"]["bcl2fq"]["root"]

LANES = []

##
# Holy guacamole, this got crazy with the setup (it grew out of doing this manually once in a while & I did not enjoy it :-)
#
#The setup stuff can all be chucked overboard, but the per-lane bits are largely taken and expanded to allow for more flags.  And most important (to me at least) is generating the mqc and interop reports to track over time.  those commands should justerally cut and paste (well... depending what our naming scheme is.##


if 'run_b2fq' in config:
    if 'b2fq_ruid' not in config.keys() or 'b2fq_exid' not in config.keys() or "source_run_path" not in config.keys() or "b2fq_lanes" not in config.keys():
        print('IF running BCL2FQ ---\n to run BCL2FQ, you must specify, in addition to the target name : " --config run_b2fq=true b2fq_ruid=RU### b2fq_exid=EX#### b2fq_lanes=1+2+3+4 source_run_path=/path/to/run/dir/  (-optionally-) b2fq_use_bases_mask=y#n#,i#,i#,y#n# b2fq_nmismatch=0|1 (and) b2fq_addl_flags= "')
        raise Exception(f'{config.keys()} .. ERROR --- \n\all three config keys must be set:  b2fq_ruid|b2fq_exid|source_run_path|b2fq_lanes \n ')
else:
    config['b2fq_b2fq_exid'] =""
    config['b2fq_ruid'] = ""
    config['source_run_path'] = ""
    config['b2fq_lanes'] = ""

if 'b2fq_lanes' in config:
    if config['b2fq_lanes'] != "":
        for l in config['b2fq_lanes'].split('+'):
            LANES.append(l)
if 'source_run_path' not in config:
    config['source_run_path']=''
if 'b2fq_ruid' not in config:
    config['b2fq_ruid'] = ''
    RUU=[""]
if 'b2fq_exid' not in config:
    config['b2fq_exid'] = ''
    EXX=[""]

config['b2fq_disambiguated_samplesheet_name'] = f"{config['b2fq_ruid']}_{config['b2fq_exid']}_B2FQ_SampleSheet.csv"
LANE_TILES =  {"L001": "s_1", "L002": "s_2", "L003": "s_3", "L004": "s_4"}


if 'run_b2fq' in config:
    if 'b2fq_use_bases_mask' in  config:
        config['b2fq_use_bases_mask_flag']= f" --use-bases-mask {config['b2fq_use_bases_mask']} "
    else:
        print('''\n  ! WARNING ----\n   " --config b2fq_use_bases_mask='...' not set.  \n    y150n1,i8,i8,y will be used by default''')
        config['b2fq_use_bases_mask_flag']=" --use-bases-mask y150n1,i8,i8,y150n1 "

    if 'b2fq_addl_flags' not in config :
        print('   WARNING --- \n    " --config b2fq_addl_flags="" not set, no addl bcltofq flags will be set\n\n')
        config['b2fq_addl_flags']= " "   #common options :  --fastq-compression-level 1  --create-fastq-for-index-reads  --ignore-missing-bcls --adapter-stringency=0.9 --mask-short-adapter-reads 22  --minimum-trimmed-read-length 35  ... specified on the command line like "b2fq_addl_flags":"++ignore-missing-bcls ++adapter-stringency=0.6" ... where the ++ will be replaced by '--', which is -- is used, confuses the SMK config parser
    else:
        config['b2fq_addl_flags']= config['b2fq_addl_flags'].replace('++','--')
else:
    pass



runn=config['b2fq_ruid']
if runn in [None,""," "]:
    runn=""
B2FD = f"results/bcl2fq/{runn}/"
config['bcl2fq_output_path']=B2FD
config['bcl2fq_out_dir'] = f"{B2FD}output"
##### RUN BCL2FQ
## -------------
#
# minimalistic BCL2FQ pipeline that offers no user input validation
# but does formalize running BCL2FQ out of band, or with unusual
# settings.  All inputs must be defined in the rules.yaml profile
# config file for this to run
# https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf



LANE_TILES = {"L001": "s_1", "L002": "s_2", "L003": "s_3", "L004": "s_4"}

if "ru_path" not in config:
    config["ru_path"] = config["source_run_path"]
    config["ru_uid"] = config["b2fq_ruid"]
    config["bcl2fq_output_path"] = f"{B2FD}"


if 'b2fq_nmismatch' in config:
    config["bcl2fq"]["bcl2fq_cmd"]["barcode_mismatch"] = config['b2fq_nmismatch']
else:
    pass # taken from default config


localrules:
    config_check,


rule config_check:
    output:
        f"{B2FD}all.clear",
    threads: 1
    conda:
        config["bcl2fq"]["env_yaml"]
    shell:
        "touch {output};"
        "{latency_wait};"

def gen_lanedirs(roots,ru):
    #for r in roots:
    for lane in LANES:
        rl_d = f"{B2FD}output/{lane}"  # rl_d = f"{B2FD}output/{r}/{lane}"
        cmd = f"mkdir -p {rl_d} > /dev/null 2>&1"
        os.system(cmd)


localrules:
    prep_bcl2fq_lane_dirs,


rule prep_bcl2fq_lane_dirs:
    """Setup the space to run bcl2fq"""
    input:
        f"{B2FD}all.clear",
    output:
          run_xml=f"{B2FD}RunInfo.xml",
    conda:
        config["bcl2fq"]["env_yaml"]
    log:
        B2FD + "logs/bcl2fq.setup.log",
    threads: 1
    params:
        lbcl_reports_d=gen_lanedirs(["Reports", "Stats", "InterOp"], runn),
        run_path=config["source_run_path"],
        cp_sample_sheet=config["b2fq_disambiguated_samplesheet_name"],
        orig_sample_sheet=config["bcl2fq"]["sample_sheet"]
    shell:
        """
        mkdir -p $(dirname {log} ) > /dev/null 2>&1 ;
        (ln -s {params.run_path}Data/ {B2FD}Data &>> {log}) || echo failedToCreate-filealreadycp ;
        (cp {params.run_path}/RunInfo.xml {output} &>> {log}) || echo "cpFileExists-filealreadycp" ;
        (cp {params.orig_sample_sheet} {params.cp_sample_sheet} &>> {log} ) || echo "cpFailed-filealradycp";
        perl -pi -e 's/(.*\,)(XL.*\-SQ)(\d+)(-.*)/$1SQ$3$4/g;' {params.cp_sample_sheet};
        perl -pi -e 's/(.*)(XL.*\-SQ)(\d+)(-.*)/$1SQ$3$4/g;' {params.cp_sample_sheet}  ;
        {latency_wait};
        """


localrules:
    spit_out_lanes,

rule spit_out_lanes:
    input:
        run_xml=f"{B2FD}RunInfo.xml",
    output:
        expand(B2FD  + "output/Reports/{lane}/html/lane.ready", lane=LANES),
    shell:
        "touch {output};"

##if "b2fqrb" not in config:
##    config["b2fqrb"] = 1
##else:
##    config["b2fqrb"] = int(config["b2fqrb"])


def get_machine_ulimit(wildcards):
     flim = resource.getrlimit( resource.RLIMIT_OFILE )
     return int(str(int(float(flim[1])/2)))-1

rule bcl2fq_by_lane:
    """RunBCL2FQ prep lanes"""
    input:
        B2FD +"output/Reports/{lane}/html/lane.ready",
        run_work_dir=f"{B2FD}",
        run_xml=f"{B2FD}RunInfo.xml",
    output:
        B2FD +"output/Reports/{lane}/html/index.html",  #  run=RU, lane=LANES),
    threads: 8 if "threads" not in config["bcl2fq"] else config["bcl2fq"]["threads"]
    benchmark:
        B2FD +"benchmarks/B2FQ_"+runn+"_{lane}.benchmark.txt"
    conda:
        config["bcl2fq"]["env_yaml"]
    params:
        bc_mm=config["bcl2fq"]["bcl2fq_cmd"]["barcode_mismatch"],
        sample_sheet=config['b2fq_disambiguated_samplesheet_name'],
        threads_loading=config["bcl2fq"]["bcl2fq_cmd"]["threads_loading"],
        threads_writing=config["bcl2fq"]["bcl2fq_cmd"]["threads_writing"],
        threads_processing=config["bcl2fq"]["bcl2fq_cmd"]["threads_processing"],
        #threads_demux=config["bcl2fq"]["bcl2fq_cmd"]["threads_demux"],
        ulimit_detected=get_machine_ulimit,
        cluster_sample=lambda wildcards:  f"b2fq_{wildcards.lane}",
        ld_preload=" ",  # if 'ld_preload' not in config['malloc_alt'] else config['malloc_alt']['ld_preload'],ld_pre=" " if 'ld_preload' not in config['bcl2fq'] else config['bcl2fq']['ld_preload'],
        ruid=config['b2fq_ruid'],
        exid=config['b2fq_exid'],
        ru_source_path=config['source_run_path'],
        addl_flags=" " if 'b2fq_addl_flags' not in config else config['b2fq_addl_flags'],
        use_bases_mask_flag=" " if 'b2fq_use_bases_mask_flag' not in config else config['b2fq_use_bases_mask_flag']
    log:
        B2FD+"logs/{lane}.bcl2fq.log",
    shell:
        """
        if [[ -f {log} ]]; then
            echo EXISTS!!! &>> {log};
        else
            export tile=s_0;
            if [[ "{wildcards.lane}" == "L001" ]]; then export tile=s_1; fi  &>> {log} ;
            if [[ "{wildcards.lane}" == "L002" ]]; then export tile=s_2; fi  &>> {log} ;
            if [[ "{wildcards.lane}" == "L003" ]]; then export tile=s_3; fi  &>> {log} ;
            if [[ "{wildcards.lane}" == "L004" ]]; then export tile=s_4; fi  &>> {log} ;
            echo $tile &>> {log} ;
            ((ulimit -Sn  {params.ulimit_detected}) || echo ulimitFAIL ) &>> {log};
            ( {params.ld_preload} bcl2fastq         \
             --loading-threads {params.threads_loading}         \
             --processing-threads {params.threads_processing}      \
             --writing-threads {params.threads_writing}         \
             --barcode-mismatches {params.bc_mm}        \
             --output-dir {input.run_work_dir}/output       \
             --runfolder-dir {params.ru_source_path}      \
             {params.use_bases_mask_flag}  {params.addl_flags}       --sample-sheet {params.sample_sheet}  \
             --interop-dir {input.run_work_dir}/output/InterOp/{wildcards.lane}     \
             --stats-dir {input.run_work_dir}/output/Stats/{wildcards.lane}      \
             --reports-dir {input.run_work_dir}/output/Reports/{wildcards.lane} \
             --tiles $tile ) &>> {log} ;
        fi ;
        """


localrules:
    multiqc_bcl2fq_ready,

rule multiqc_bcl2fq_ready:
    """call_multiqc"""
    input:
        expand(
            f"{B2FD}"+"output/Reports/{lane}/html/index.html",
            lane=LANES,
        ),
    output:
        ready=f"{B2FD}output/{runn}_postbcl2fq.ready",
        bcl2fq_go="results/bcl2fq_reports/go.run",
        #manifest="config/analysis_manifest.csv"
    threads: 1
    log:
        f"{B2FD}snk_log/{runn}_bcl2fq_multiqc.ready.log",
    conda:
        config["multiqc"]["bcl2fq"]["env_yaml"]
    shell:
        """
        touch {output};
        {latency_wait};
        """

rule bcl2fastq_run:  # TARGET : Target to run BCL2FQ
    input:
        f"{B2FD}output/"+runn+"_postbcl2fq.ready",  #
        f"results/bcl2fq_reports/bcl2fq_{config['b2fq_ruid']}.multiqc.html",
    output:
        f"{B2FD}" + f"output/bcl2fq.done",
    log:
        f"{B2FD}" + f"output/bcl2fq_finalrrule.log",
    params:
        bcl_outdir=os.path.abspath(f"{B2FD}"+"output/"),
        ru=runn,
        ex=config['b2fq_exid'],
        cluster_sample="bcl2fq",
    benchmark:
        B2FD + f"output/bcl2fq_finalrrule.bench.tsv",
    conda:
        config["multiqc"]["bcl2fq"]["env_yaml"]
    shell:
        """epocsec=$(date +'%s');
        cmd=$( echo "at2mgsamp_sheet.py {params.bcl_outdir} {params.ru} {params.ex} merge $DAY_ROOT/bcl2fq_links_$epocsec/ $DAY_ROOT/config/analysis_manifest.csv");
        echo "$cmd";
        $cmd;
        colr '  ANALYSIS SAMPLESHEET READY IN config/analysis_manifest.csv' 'yellow' 'green' 'b' 2>&1;
        touch {output};
        {latency_wait};
        """
