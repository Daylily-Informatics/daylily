import os

# ##### Multiqc For BCL2FQ and Interop
#
#
def get_ru_d(wildcards):
    ret_dirs = []
    b2fq_ruid = config['b2fq_ruid']
    b2fq_path = config['source_run_path']
    ret_d = ""
    print_err = False
    if "source_run_path" in config:
        ret_dirs.append(config["source_run_path"])
    else:
        print_err = True
    if "bcl2fq_output_path" in config:

        if config["bcl2fq_output_path"].endswith(".ready"):
            b2fq_path = os.path.dirname(config["bcl2fq_output_path"]) + "/"
        elif config["bcl2fq_output_path"].endswith("output/"):
            b2fq_path = config["bcl2fq_output_path"]
        elif config["bcl2fq_output_path"].endswith("output"):
            b2fq_path = config["bcl2fq_output_path"] + "/"
        if b2fq_path.endswith("output/"):
            ret_d = b2fq_path
        else:
            raise Exception("bcl2fq_output_path must end in output/")

    else:
        print_err = True

    if "b2fq_ruid" in config:
        b2fq_ruid = config["b2fq_ruid"]
    else:
        print_err = True
        b2fq_ruid = RU  # hack, we could theoretically allow mised runs...

    if print_err:
        os.system(
            f"colr '  +++ W A R N I N G +++ you are running target multiqc_bcl2fq w/out setting one or all of: --config ru_path=path/to/ru_dir  b2fq_ruid=RU123456 bcl2fq_output_path=/path/to/bcl2fq/output/  ~~or~~ this is being called via bcl2fq and these are not being set there-  if not, it is advisable to correct tis issue before proceeding. ' 'black' 'yellow' 'b' 2>&1; "
        )

        os._exit(1)

    return ret_d


def get_ru_dir(wildcards):
    ret_dirs = []
    b2fq_ruid = "NA"
    b2fq_path = ""
    print_err = False
    if "ru_path" in config:
        ret_dirs.append(config["ru_path"])
    else:
        print_err = True
    if "bcl2fq_output_path" in config:
        ret_dirs.append(config["bcl2fq_output_path"])
    else:
        print_err = True

    if "b2fq_ruid" in config:
        b2fq_ruid = config["b2fq_ruid"]
    else:
        print_err = True

    bcl2fq_go="results/bcl2fq_reports/go.run"
    if 'singleton_mqc_bcl2fq' in config:
        os.system(f"touch {bc2fq_go}")
    ret_dirs.append(bcl2fq_go)

    if print_err:
        os.system(
            f"colr '  ERROR -- you are running target multiqc_bcl2fq w/out setting one or all of: --config ru_path=path/to/ru_dir  b2fq_ruid=RU123456 bcl2fq_output_path=/path/to/bcl2fq/output/  ~~or~~ this is being called via bcl2fq and these are not being set there. ' 'black' 'yellow' 'b' 2>&1;"
        )
        os._exit(1)

    return ret_dirs


# !!!!! this is an old rule with a lot of crap in it.....

rule multiqc_bcl2fq:
    input:
        f"{config['bcl2fq_output_path']}",
        f"{config['source_run_path']}",
        "results/bcl2fq_reports/go.run"
    output:
        html=f"results/bcl2fq_reports/bcl2fq_{config['b2fq_ruid']}.multiqc.html",
        interop1=f"results/bcl2fq_reports/interop_data/{config['b2fq_ruid']}_interop_summary.csv",
        interop2=f"results/bcl2fq_reports/interop_data/{config['b2fq_ruid']}_interop_index-summary.csv",
    threads: config["multiqc"]["threads"]
    log:
        f"results/logs/bcl2fq_{config['b2fq_ruid']}_multiqc.log",
    benchmark:
        f"results/benchmarks/bcl2fq_{config['b2fq_ruid']}_multiqc.bench.tsv"
    params:
        rud= f"{config['source_run_path']}",
        fn=f"bcl2fq_{config['b2fq_ruid']}.multiqc.html",
        odir="results/bcl2fq_reports/",
        odir2="results/bcl2fq_reports/".replace("/", "\/"),
        interop_d="results/bcl2fq_reports/interop_data/",
        macro_cfg=config["multiqc"]["config_yaml"],
        gtag=config["gittag"],
        ghash=config["githash"],
        gbranch=config["gitbranch"],
        ru=config["b2fq_ruid"],
        mdir="results/",
        cluster_sample="bcl2fqMQC_{config['b2fq_ruid']}",
        iop_d="results/bcl2fq_reports/"+config["b2fq_ruid"],
    container:
        "docker://daylilyinformatics/daylily_multiqc:0.2"
    shell:
       """
       (ln -s {params.rud} results/bcl2fq_reports/{params.ru}) ;
       (mkdir -p $(dirname {output.interop1}) ) ;
       ### export LD_LIBRARY_PATH=./resources/lib/ &&
       mkdir -p {params.interop_d};
       ### export LD_LIBRARY_PATH=./resources/lib/ && interop_summary {params.iop_d} --csv=1 > {output.interop1};
        ### export LD_LIBRARY_PATH=./resources/lib/ && interop_index-summary {params.iop_d} --csv=1 > {output.interop2};
        ### export LD_LIBRARY_PATH=./resources/lib/ && interop_plot_flowcell {params.iop_d}  | perl -pe "s/(output \')(.*)(\.png)/\$1{params.odir2}\$2\$3_mqc\.png/g"  | gnuplot ;
        ### export LD_LIBRARY_PATH=./resources/lib/ && interop_plot_by_cycle  {params.iop_d} | perl -pe "s/(output \')(.*)(\.png)/\$1{params.odir2}\$2\$3_mqc\.png/g" | gnuplot;
        ### export LD_LIBRARY_PATH=./resources/lib/ && interop_plot_by_lane  {params.iop_d} | perl -pe "s/(output \')(.*)(\.png)/\$1{params.odir2}\$2\$3_mqc\.png/g" | gnuplot;
        ### export LD_LIBRARY_PATH=./resources/lib/ && interop_plot_qscore_histogram  {params.iop_d} | perl -pe "s/(output \')(.*)(\.png)/\$1{params.odir2}\$2\$3_mqc\.png/g" | gnuplot;
       ### export LD_LIBRARY_PATH=./resources/lib/ && interop_plot_qscore_heatmap  {params.iop_d}| perl -pe "s/(output \')(.*)(\.png)/\$1{params.odir2}\$2\$3_mqc\.png/g" | gnuplot;
       multiqc --interactive -m custom_content  -m interop -m bcl2fastq -x '*/mod/*' -x '*.js' -x '*.bam' -x '*.fastq.gz' -x '*multiqc*' -x '*pyc' -x '*.fastq.gz' -v -x '*/impute/*'  -i '{params.ru} BCL2FQ/Interop REPORT' -p  -b '{params.ru} ___ {params.gbranch} {params.gtag} {params.ghash}' -n {params.fn} -o {params.odir} --profile-runtime  -c {params.macro_cfg} -f results/ ;
        touch {params.mdir}benchmarks/bcl2fq_{params.ru}_multiqc.bench.tsv ;
        {latency_wait};
        ls {output}; """
        # in lieu of the command line latency wait
        # fail hard if anyting is missing, don't trigger the wait for file cycle as it's buggy. latency wait built in abve
