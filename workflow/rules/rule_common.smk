from snakemake.utils import min_version
from snakemake.utils import validate
import os
import pandas as pd
import sys
import yaml
import multiprocessing
import random
import shutil
import datetime as dtm

## TODO: sweep through this thing and prune out the no longer needed/experimental stuff

# ec2=os.popen('echo "EC2instance:$(bin/helpers/get_ec2_type.sh)"').readline().rstrip()
# print(ec2, file=sys.stderr)

EX = [""]
RU = [""]
# ##### these are globally avail, but  make the linters freak out b/c it's hidden from them by snakemake, quieting the complaints.
config = config  # noqa   ### Just needed to quiet linters
cluster_config = cluster_config  # noqa   ### Just needed to quiet linters


# ##### Generate the day top level directories
if "conda_prefix" not in config:
    config["conda_prefix"] = "etc/"


config["jem_dot"] = os.environ["DAY_ROOT"] + "/.jemalloc_loc"
os.system(f"touch {config['jem_dot']}")


# ##### A place to track failed samples so they can be tallied up at the end.
# NOT REALLY IMPLEMENTED YET
config["failed_samples"] = {}  ### USE SAMPLE SHEET FOR THIS


# ##### Safety sort of the yaml defined crms
config["glimpse"] = {"impute_chrms": "2"}


# ####  PARSE SUPPORTING DATA FILES/DIRS INTO CONFIG
# -----------------------------------
# suppporting dirs

config["supporting_files"] = f"config/supporting_files/{config['genome_build']}_supporting_files.yaml"

a_yaml_file = open(config["supporting_dirs"])
parsed_yaml_file = yaml.load(a_yaml_file, Loader=yaml.FullLoader)
config["supporting_dirs_obj"] = parsed_yaml_file

files_yaml_file = open(config["supporting_data"]["supporting_files_conf"])
files_yaml_file = yaml.load(files_yaml_file, Loader=yaml.FullLoader)

# supporting files
config["supporting_files"] = {
    "files": files_yaml_file["supporting_files"]["files"],
    "root": files_yaml_file["supporting_files"]["root"],
}

genome_build_chrm_prefix_map = {
    "b37": "",
    "hg38": "chr", 
    "hg38_broad": "chr",
} 

GENOME_CHR_PREFIX="na"
if os.environ.get("DAY_GENOME_BUILD","na") in genome_build_chrm_prefix_map:
    GENOME_CHR_PREFIX = genome_build_chrm_prefix_map[os.environ.get("DAY_GENOME_BUILD")]
    print(
        f"INFO::: The genome build {os.environ.get('DAY_GENOME_BUILD')} is supported.  The genome build prefix is '{GENOME_CHR_PREFIX}''.",
        file=sys.stderr,
    )
else:
    err_msg=f"ERROR::: The genome build {os.environ.get('DAY_GENOME_BUILD')} is not supported.  Please check the config file and try again. No genome build prefix was set."
    print(
        err_msg,
        file=sys.stderr,
    )
    raise Exception(
        err_msg
    )

# SNV caller chunk arrays
SENTD_CHRMS = config["sentD"][f"{config['genome_build']}_sentD_chrms"].split(",")
DEEPD_CHRMS = config["deepvariant"][f"{config['genome_build']}_deep_chrms"].split(",")
OCTO_CHRMS = config["octopus"][f"{config['genome_build']}_octo_chrms"].split(",")
CLAIR3_CHRMS = config["clair3"][f"{config['genome_build']}_clair3_chrms"].split(",")
LOFREQ_CHRMS = config["lofreq2"][f"{config['genome_build']}_lofreq_chrms"].split(",")
SENTDUG_CHRMS = config["sentdug"][f"{config['genome_build']}_sentdug_chrms"].split(",")
SENTDONT_CHRMS = config["sentdont"][f"{config['genome_build']}_sentdont_chrms"].split(",")
SENTDHUO_CHRMS = config["sentdhuo"][f"{config['genome_build']}_sentdhuo_chrms"].split(",")
SENTDHIO_CHRMS = config["sentdhio"][f"{config['genome_build']}_sentdhio_chrms"].split(",")
SENTDPB_CHRMS = config["sentdpb"][f"{config['genome_build']}_sentdpb_chrms"].split(",")
SENTDPBR_CHRMS = config["sentdpbr"][f"{config['genome_build']}_sentdpbr_chrms"].split(",")
SENTDONTR_CHRMS = config["sentdontr"][f"{config['genome_build']}_sentdontr_chrms"].split(",")
SENTDUGR_CHRMS = config["sentdugr"][f"{config['genome_build']}_sentdugr_chrms"].split(",")
SENTDHIP_CHRMS = config["sentdhip"][f"{config['genome_build']}_sentdhip_chrms"].split(",")

# ##### Setting the allowed aligners to run and to which deduper to use.
# presently, 1+ aligners may run, but all must use the same deduper


# Handle aligners


ALIGNERS = []
if 'aligners' not in config:
    
    os.system(
        f'''colr "...WARNING: No aligners set in the config." "$DY_WT1" "$DY_WB1" "$DY_WS1" 1>&2'''
    )
else:
    ALIGNERS = sorted(set([] if config.get('aligners') is None else config["aligners"]))
    # PRINT INFO
    os.system(
        f"""colr 'aligners: {ALIGNERS}' "$DY_WT1" "$DY_B1" "$DY_WS1" 1>&2;"""
    )

# Handle dedupers
DDUP = []
if 'dedupers' not in config:
    os.system(
        f'''colr "...WARNING: No dedupers set in the config." "$DY_WT1" "$DY_WB1" "$DY_WS1" 1>&2'''
    )
else:
    DDUP = sorted(set([] if config.get('dedupers') is None else config["dedupers"]))
    # PRINT INFO
    os.system(
        f"""colr 'deduper: {DDUP}' "$DY_WT1" "$DY_B1" "$DY_WS1" 1>&2;"""
    )

snv_CALLERS = []
if 'snv_callers' not in config:
     os.system(
        f'''colr "...WARNING: No snv_callers set in the config." "$DY_WT1" "$DY_WB1" "$DY_WS1" 1>&2'''
     )
else:
    snv_CALLERS = sorted(set([] if 'snv_callers' not in config or config['snv_callers'] == None else config["snv_callers"]))
    ## PRINT INFO
    os.system(
        f"""colr 'SNV Callers:{snv_CALLERS}' "$DY_WT1" "$DY_B1" "$DY_WS1" 1>&2;"""
    )

sv_CALLERS = []
if 'sv_callers' not in config:

     os.system(
        f'''colr "... WARNING: No sv_callers set in the config." "$DY_WT1" "$DY_WB1" "$DY_WS1" 1>&2'''
     )
else:
    sv_CALLERS = sorted(set([] if 'sv_callers' not in config or config['sv_callers'] == None else config["sv_callers"]))
    ## PRINT INFO
    os.system(
        f"""colr 'SV Callers:{sv_CALLERS}' "$DY_WT1" "$DY_B1" "$DY_WS1" 1>&2;"""
    )   
    


# ####  LOAD SAMPLE DATA FROM ANALYSIS MANIFEST
# ---------------------
# Load sample sheet into pandas dataframe, then validate & if pass
# expose amples as a global  Including to any other rule or called script from a rule  included from the
# main snakefile



if "analysis_manifest" in config:
    os.system(
        """colr  '     _____ COMMAND LINE ANALYSIS MANIFEST SET. This will be copied to config/analysis_manifest.csv, and this copy used' "$DY_WT0" "$DY_WB0" "$DY_WS1"  """
    )
    if os.path.exists("config/analysis_manifest.csv"):
        raise Exception(
            "\n\n A file exists in config/analysis_manifest.csv, you must remove it to use the command line specified manifest.\n\n"
        )
    else:
        user_analysis_manifest = config["analysis_manifest"]
        if os.path.exists(user_analysis_manifest):
            os.system(f"cp {user_analysis_manifest} config/analysis_manifest.csv")
            os.system(f"echo '{str(dtm.datetime.now())}\t{user_analysis_manifest}' >> config/analysis_manifest.log")
        else:
            raise Exception(
                f"\n\nERROR::: The user specified analysis manifest file was not found at {user_analysis_manifest}.  Please check the path and try again."
            )
else:
    default_analysis_manifest = config[f"{config['genome_build']}_analysis_manifest"]
    if os.path.exists("config/analysis_manifest.csv"):
        os.system(
            f"""colr '     _____ EXISTING ANALYSIS MANIFEST FILE DETECTED: config/analysis_manifest.csv --  this will be used' "$DY_WT0" "$DY_WB0" "$DY_WS1" >&2 """
        )
        os.system('sleep 1')
    elif os.path.exists(default_analysis_manifest):
        os.system(f"cp {default_analysis_manifest} config/analysis_manifest.csv")
        os.system(f"echo '{str(dtm.datetime.now())}\t{default_analysis_manifest}' >> config/analysis_manifest.log")

    else:
        raise Exception(
            f"\n\nERROR::: The default analysis manifest file was not found at {default_analysis_manifest}.  Please check the path and try again, was the genome_build specified?"
        )

config["analysis_manifest"] = "config/analysis_manifest.csv"
analysis_manifest = config["analysis_manifest"]

# Alert is analysis manifest is not found
os.system(
    f"""(colr "A    N   A   L  Y S I S    MANIFEST FILE DETECTED FOR USE ::: {analysis_manifest}" "$DY_ET1" "$DY_EB1" "$DY_ES1" 1>&2; sleep {config["warn_err_sleep"]}) || (colr "~~~ANALYSIS MANIFEST DETECTIONFAILED!!~~~~ " "$DY_ET0" "$DY_EB0" "$DY_ES0" 1>&2; sleep {config["warn_err_sleep"]}; exit 33;)"""
)

# IMPORTANT: initialize the samples dataframe from the analysis_manifest.csv
samples = pd.read_table(analysis_manifest, ",").set_index(
    ["sample", "sample_lane"], drop=False
)

# validate the analysis_manifest.csv with the appropriate yaml schema
validate(samples, schema="../schemas/analysis_manifest.schema.yaml")


# TODO: remove the RU/EX usage throughout, unecessary
if len(samples) > 0:
    RU = list(samples.RU)
    EX = list(samples.EX)

# Strip out Unclassified unless directed not to
if "keep_undetermined" in config:
    pass
else:
    new_samples = samples.query("SQ != 'Undetermined'")
    samples = new_samples


# scan the input files and calculate the per sample total fastq size, for use in scaling compute.

pos_samp_vcf_bed = {}

config["samp_fq_size"] = {}

for s in samples["sample"]:
    if s not in config["samp_fq_size"]:
        config["samp_fq_size"][s] = 0.0
        for r1f in samples.loc[s, ("r1_path")]:
            if r1f in ["","na", None]:
                pass
            else:
                try:
                    config["samp_fq_size"][s] += float(
                        float(os.path.getsize(r1f)) / 1000000000.0
                    )  # Norm to GB
                except Exception as e:
                    print(e, file=sys.stderr)

# added to cluster name to uniquely identify sets of jobs from the same snakemake execution
ri2 = random.randint(0, 999)


sub_user = config["sub_user"]
config["sub_host"] = (
    str(os.environ.get("HOSTNAME"))
    if os.environ.get("HOST") in ["", "None", None]
    else str(os.environ.get("HOST"))
)
config["sub_root"] = str(os.environ.get("PWD"))
config["calling_cli"] = str(os.environ.get("DY_CALLING_CLI"))

#pop out to ipython if set in command line config
if "ipython" in config:
    from IPython import embed

    embed()
    #    Snakefile can be reloaded, to keep this from activating again, remove it
    del config["ipython"]


# TODO: DEPRECATE
# Set the BATCH id for the cluster
try:
    cluster_config[
        "batch_id"
    ] = f"{samples['Run'][0]}_{samples['EX'][0]}_{config['githash']}"
except Exception as e:
    del e
cluster_config["profile_name"] = config["profile_name"]
cluster_config["cwd"] = os.environ.get("PWD", "./")
cluster_config["sub_user"] = os.environ.get("USER", "na")


### TO DO: CHANGE THIS TO USE COLUMN HEADERS AND NOT IMPLICIT ORDER!!!!

# Set info for concordance calculations
CONCORDANCE_SAMPLES = {}
CRAM_ALIGNERS = []
sample_info = {}
sample_lane_info = {}
for i in samples.iterrows():
    samp = i[0]
    sample = i[1][0]
    dat = dict(i[1])
    sq_i= i[1][3]
    ru_i= i[1][4]
    ex_i= i[1][5]
    lane_i= i[1][6]
    sample_lane = i[1][1]
    merge_single = i[1][14]
    if samp in sample_info and merge_single in ["single"]:
        raise (
            "\n\nMANIFEST ERROR:: "
            + samp + " ... " + sample_lane
            + f"appears 2+ times in the sample sheet. This should only occur if 'merge' has been specified. {merge_single} has been set. Also.. column order is sadly important: samp,sample,sample_lane,SQ,RU,EX,LANE,r1_path,r2_path,biological_sex,iddna_uid,concordance_control_path,is_positive_control,is_negative_control,sample_type,merge_single,external_sample_id,instrument,lib_prep,bwa_kmer"
            )

    if merge_single in ["single"]:
        raise Exception(
            f"\n\nMANIFEST ERROR '{sample}, {sample_lane}': This feature was implemented, then unusued and has not been vetted to work properly again, so if you wish to run per lane, create a manifest with the sample and sample_lane column having the same id.  merge will create symlinks for each sample_lane pair of fastqs, then use these via process substitution to appear as one file for those tools expecting 1 R1 and 1 R2. Also.. column order is sadly important: samp,sample,sample_lane,SQ,RU,EX,LANE,r1_path,r2_path,biological_sex,iddna_uid,concordance_control_path,is_positive_control,is_negative_control,sample_type,merge_single,external_sample_id,instrument,lib_prep,bwa_kmer"
        )

    if len(sq_i.split(".")) > 1 or len(sq_i.split("_")) > 1 or len(ru_i.split(".")) > 1 or len(ru_i.split("_")) > 1 or len(ex_i.split(".")) > 1 or len(ex_i.split("_")) > 1 or len(str(lane_i).split(".")) > 1 or len(str(lane_i).split("_")) > 1:
        raise Exception(
            f"\n\nMANIFEST ERROR {sample} ... {sample_lane}: The SQ & RU & EX & LANE cols may not contain a period or '-' in the name.  Please check the sample name and try again. Also.. column order is sadly important: samp,sample,sample_lane,SQ,RU,EX,LANE,r1_path,r2_path,biological_sex,iddna_uid,concordance_control_path,is_positive_control,is_negative_control,sample_type,merge_single,external_sample_id,instrument,lib_prep,bwa_kmer"
        )

    if sample_lane in sample_info or len(sample.split(".")) > 1:
        raise Exception(
            f'\n\nMANIFEST ERROR  {sample} ... {sample_lane}: It is not allowed for the sample entry to be present in the sample_lane column and vice versa.  sample can be any string, excluding containing a period or "-", and does not need to be unique w/in the sample column. it may not end in "_[1-9]+", howevr, it may end in "_0"'
        )

    
    

    samp = samp[0]
    if samp not in sample_info:
        sample_info[samp] = {}

    for iix in dat:
        val = dat[iix]
        if iix in ["biological_sex"]:
            val = val.lower()

            if val.startswith("m"):
                val = "male"
            elif val.startswith("f"):
                val = "female"
            else:
                val = "na"
            sample_info[samp][iix] = val
        elif iix in ["ultima_cram"]:
            sample_info[samp][iix] = val
        elif iix in ["ont_cram"]:
            sample_info[samp][iix] = val
        elif iix in ["ultima_cram_aligner"]:
            sample_info[samp][iix] = val
            if val not in CRAM_ALIGNERS:
                if val in ['','na',None,'None','hyb']:
                    pass
                else: 
                    CRAM_ALIGNERS.append(val)
        elif iix in ["ont_cram_aligner"]:
            sample_info[samp][iix] = val
            if val not in CRAM_ALIGNERS:
                if val in ['','na',None,'None']:
                    pass
                else: 
                    CRAM_ALIGNERS.append(val)
        elif iix in ["pb_bam_aligner"]:
            sample_info[samp][iix] = val
            if val not in CRAM_ALIGNERS:
                if val in ['','na',None,'None']:
                    pass
                else: 
                    CRAM_ALIGNERS.append(val)
        elif iix in ["deep_model"]:
            """supported models
                WGS: Whole-genome sequencing data, typically for human genomes sequenced to around 30x coverage.
                WES: Whole-exome sequencing data.
                PACBIO: Long-read PacBio data.
                HYBRID_PACBIO_ILLUMINA: Hybrid PacBio and Illumina data.
                ONT_R104: Oxford Nanopore Technologies data (model optimized for R10.4 chemistry).
                ONT_R941: Oxford Nanopore Technologies data (model optimized for R9.4.1 chemistry).
            """
            deep_models = [
                "WGS",
                "WES",
                "PACBIO",
                "HYBRID_PACBIO_ILLUMINA",
                "ONT_R104",
                "ONT_R941",
            ]
            if val not in deep_models:
                print(
                    f"\n\n\tWARNING::: The model {val} is not in the supported set of {deep_models} .  WGS will be selected!!!\n\n\n",file=sys.stderr
                )
                sample_info[samp][iix] = "WGS"
            else:
                sample_info[samp][iix] = val
        elif iix in ["ultima_cram_snv_caller"]:
            sample_info[samp][iix] = val
        elif iix in ["ont_cram_snv_caller"]:
            sample_info[samp][iix] = val
        elif iix in ["concordance_control_path"]:
            sample_info[samp][iix] = val
            if val not in ["na", "NA", "", None, "None"]:
                CONCORDANCE_SAMPLES[samp] = val
        elif iix in ["sample_type"]:
            sample_info[samp][iix] = val
        elif iix in ["is_positive_control"]:
            sample_info[samp][iix] = val
        elif iix in ["is_negative_control"]:
            sample_info[samp][iix] = val
        elif iix in ["external_sample_id"]:
            sample_info[samp][iix] = val
        elif iix in ["bwa_kmer"]:
            sample_info[samp][iix] = val
        elif iix in ["iddna_uid"]:
            val = val.split(":")
            sample_info[samp][iix] = val
        elif iix in ["instrument"]:
            sample_info[samp][iix] = val
        elif iix in ["lib_prep"]:
            sample_info[samp][iix] = val
        elif iix in ["sample_lane"]:
            if val in sample_lane_info:
                raise Exception(
                    "The sample+lane entry, 'sample_lane' must be unique in the sample sheet, but this entry was seen 2x in the column {val}.  You may find sample sheet specs in the docs, and can run 'marion the sample sheets. gold-run yo_samps' to just run the yaml validators.",
                    val,
                )
            else:
                sample_lane_info[val] = True

config["sample_info"] = sample_info

CRAM_ALIGNERS = list(set(CRAM_ALIGNERS))
# Aspirationally hoping to adopt PEPs...
# http://pep.databio.org/en/latest/

instanceid = 'na'
CLUSTER_PRI=""

# #### UTILITY METHODS
# -------------------
# Return just the R1 file for a given sampleid
# note the sort of strange invocation of this method in the rule below
# you specify simply 'get_fastq_rq', without(), and the wildcards
# special variable is inserted.  It holds all of the current pattern
# matching values for the instance of the rule being run.  So for each
# rule being run calling 'get_fastq_r1' (even many at once), the {sample}
# pattern match w/in the rule is passed along in wildcards, which here
# I'm using to get just the R1 read for this sample.
def get_fastq_r1(wildcards):
    file = ""
    for f in samples["r1_path"]:
        if len(f.split(wildcards.sample)) > 1:
            file = os.path.abspath(f)
    return file


# extract the current pattern matching samplePlusR
# from the calling rule wildcards obj passed to the
# method. Since this is matching my sample name AND
# the R1/R2 bit of the file name as well, the
# samplePlusR will be unique with the R1/2 added
# we just return one file path here


def get_raw_R1s(wildcards):
    r1s = []
    for i in samples[samples["sample"] == wildcards.sample][
        "r1_path"
    ]:  # .loc[wildcards.sample, "sample_lane"]['sample_lane']:
        r1s.append(i)
    return sorted(r1s)


def get_raw_R2s(wildcards):
    r2s = []
    for i in samples[samples["sample"] == wildcards.sample][
        "r2_path"
    ]:  # .loc[wildcards.sample, "sample_lane"]:
        r2s.append(i)
    return sorted(r2s)


def get_fastq_r1_r2(wildcards):
    # generate filepaths corresponding to {sample} : {RR} combos
    if wildcards.RR == "R1":
        return None if samples.loc[wildcards.sample, "r1_path"] in ["na","",None] else os.path.abspath(samples.loc[wildcards.sample, "r1_path"])
    elif wildcards.RR == "R2":
        return  None if samples.loc[wildcards.sample, "r1_path"] in ["na","",None] else os.path.abspath(samples.loc[wildcards.sample, "r2_path"])
    else:
        raise ValueError(f"invalid value: {wildcards.RR}")


# Helper method to read the pandas dataframe for the 'sample' IDs
# from the samples.csv spreadsheet.   Used for pattern matching all
# of the expected sample IDs. No wildcards used.
def get_samp_ids():
    samps = {}
    for ii in samples["sample"]:
        if "just_this_sample" in config:
            if config["just_this_sample"] == ii:
                samps[ii] = True
        else:
            samps[ii] = True
    if len(samps.keys()) == 0 and "run_b2fq" not in config:
        raise Exception(
            "NO SAMPLES HAVE BEEN LOADED TO THE SAMPS ARRAY -or- if running Bcl2FQ, the --config run_b2fq=true and assoc params are not set."
        )
    return samps.keys()


SAMP_SAMPI = []  #  deprecate
for ix in samples.index:
    sample_x = f"{ix[0]}/{ix[1]}"
    SAMP_SAMPI.append(sample_x)


for si in samples.iterrows():
    stup = (si[1][0], si[1][1])

SSI = SAMP_SAMPI  # deprecate
SAMPS = list(get_samp_ids())

# TODO: Revisit if this is still in use
if "remove_samples" in config:
    for rs in config["remove_samples"]:
        SAMPS.remove(rs)
        print(f"SAMPLE REMOVED: {rs}")

SSAMPS = {}
for sample in list(get_samp_ids()):
    ssamp = samples[samples["sample"] == sample]["samp"][0]
    if ssamp in SSAMPS:
        SSAMPS[ssamp].append(sample)
    else:
        SSAMPS[ssamp] = [sample]


SAMP_SAMPI_INDEX = list(samples.index)  # deprecate
RR = ["R1", "R2"]


def getR2s(wildcards):
    fr2s = []
    for sample_lane in samples.loc[wildcards.sample, "sample_lane"]:
        r2 = f"{MDIR}{wildcards.sample}/{sample_lane}.R2.fastq.gz"
        fr2s.append(r2)
    return sorted(fr2s)


def getR1s(wildcards):
    fr1s = []
    for sample_lane in samples.loc[wildcards.sample, "sample_lane"]:
        r1 = f"{MDIR}{wildcards.sample}/{sample_lane}.R1.fastq.gz"
        fr1s.append(r1)
    return sorted(fr1s)

def getCRAMs(wildcards):
    crams = []
    for ss_cram in samples.loc[wildcards.sample, "cram"]:
        if ss_cram not in [None, "None", "", "na"]:
            cram = os.path.abspath(cram)
        else:
            pass
    return sorted(crams)

def getR2sS(wildcards):
    fr2s = []
    for r2 in samples[samples["samp"] == wildcards.sample]["r2_path"]:
        fr2s.append(r2)
    return sorted(fr2s)


def getR1sS(wildcards):
    fr1s = []
    for r1 in samples[samples["samp"] == wildcards.sample]["r1_path"]:
        fr1s.append(r1)
    return sorted(fr1s)


# Call from params block to get sample ID back, without() wildcards (and others ) are added automatically if no () is included.
def ret_sample(wildcards):
    if "sample" in wildcards.keys():
        return wildcards.sample
    elif "sx" in wildcards.keys():
        return wildcards.sx
    else:
        return "get sample ERROR"

def ret_sample_sentD(wildcards):
    return wildcards.sample


def ret_sample_sv(wildcards):
    return f"{wildcards.sample}_{wildcards.s_v_caller}"


def ret_sample_snv(wildcards):
    return f"{wildcards.sample}_{wildcards.snv}"


def ret_sample_alnr(wildcards):
    return f"{wildcards.sample}_{wildcards.alnr}"


def get_bwa_kmer_size(wildcards):
    ret_k = None

    if wildcards.sample in config["sample_info"]:
        if "bwa_kmer" in config["sample_info"][wildcards.sample]:
            ret_k = int(config["sample_info"][wildcards.sample]["bwa_kmer"])
    if ret_k in ["None", None]:
        ret_k = config["bwa_mem2a_aln_sort"]["k"]
    return f" -k {ret_k} "


def ret_mod_chrm(ret_str):

    if ret_str[0] in ['c']:
        if ret_str[3:5] in ['23']:
            ret_str = ret_str[0:3] + "X" + ret_str[5:]
        elif ret_str[3:5] in ['24']:
            ret_str = ret_str[0:3] + "Y" + ret_str[5:]
        elif ret_str[3:5] in ['25']:
            ret_str = ret_str[0:3] + "MT" + ret_str[5:]
    else:
        if ret_str[0:2] in ['23']:
            ret_str =  "X" + ret_str[2:]
        elif ret_str[0:2] in ['24']:
            ret_str =  "Y" + ret_str[2:]
        elif ret_str[0:2] in ['25']:
            ret_str =  "MT" + ret_str[2:]

    return ret_str


## Method to wrap the bwa fastq reader with seqtk to subsample
# if a subsample_pct column is present in the sample sheet, return the back end of the process substitution
# whcih will do the subsampling
def get_subsample_head_tail(sample_id):
    ss_head = ""
    ss_tail = ""
    raisable = False
    try:
        ss_pct = samples.loc[(sample_id), "subsample_pct"][0]
        if ss_pct in ["", "na", 0,"0.0",100,100.0,"100","100.0", "0", None, "None"]:
            pass  # no subsampling requested
        else:
            ss_pct_float = 1000.1
            try:
                ss_pct_float = float(ss_pct)
            except Exception as e:
                raisable = True
                raise (e)

            if ss_pct_float > 1.0 or ss_pct_float < 0.0:
                raisable = True
                raise Exception(
                    "ERROR:::: NO SUBSAMPLING WILL BE EXECUTED::: you must specify a float from 0.0-1.0"
                )
            else:
                ss_head = f" <( seqkit sample --line-width=0 --quiet --rand-seed=7  --seq-type=dna --proportion={ss_pct_float}  "
                ss_tail = " ) "

    except Exception as e:
        if raisable:
            print(
                "Samplesheet Error with subsample_pct column ---- \n\n", file=sys.stderr
            )
            raise (e)
        else:
            pass

    return (ss_head, ss_tail)


def get_subsample_head(wildcards):
    return get_subsample_head_tail(wildcards.sample)[0]

def get_subsample_tail(wildcards):
    return get_subsample_head_tail(wildcards.sample)[1]

def get_samp_name(wildcards):
    return wildcards.sample


def get_instrument(wildcards):
    return samples[samples["samp"] == wildcards.sample]["instrument"][0]


def get_diploid_bed_arg(wildcards):

    diploid_bed = ""
    try:
        sample_bsex = samples[samples["samp"] == wildcards.sample]["biological_sex"][0].lower()

        if "male" == sample_bsex:
            diploid_bed = f' -b {config["supporting_files"]["files"]["huref"]["male_diploid"]["name"]} '
        elif "female" == sample_bsex:
            diploid_bed = f' -b {config["supporting_files"]["files"]["huref"]["female_diploid"]["name"]} '
        else:
            diploid_bed = f' -b {config["supporting_files"]["files"]["huref"]["core_bed"]["name"]} '
    except Exception as e:
        print(
            f"ERROR:::  Unable to get biological_sex from samples dataframe for sample {wildcards.sample} for diploid bed-- {e}",
            file=sys.stderr,
        )
        diploid_bed = f' -b {config["supporting_files"]["files"]["huref"]["bed"]["name"]} '

    return f" {diploid_bed} "


def get_haploid_bed_arg(wildcards):
    haploid_bed = " "

    try:
        sample_bsex = samples[samples["samp"] == wildcards.sample]["biological_sex"][0].lower()

        if "male" == sample_bsex:
            haploid_bed = f' --haploid_bed {config["supporting_files"]["files"]["huref"]["male_haploid"]["name"]} '
        elif "female" == sample_bsex:
            haploid_bed = ' '
    
    except Exception as e:
        print(
            f"ERROR:::  Unable to get biological_sex from samples dataframe for sample {wildcards.sample} for haploid bed-- {e}",
            file=sys.stderr,
        )
        haploid_bed = ' '

    return f" {haploid_bed} "


def print_wildcards_etc(wildcards):
    print("All Wildcards: ", wildcards,  " all snv_CALLERS: ", snv_CALLERS, " all ALIGNERS: ", ALIGNERS, " all CRAM_ALIGNERS: ", CRAM_ALIGNERS, " all samples: ", SSAMPS, " all concordance samples: ", CONCORDANCE_SAMPLES.keys(), file=sys.stderr)

def get_alnr(wildcards):
    return wildcards.alnr

def get_dchrm_day(wildcards):
    pchr= GENOME_CHR_PREFIX

    ret_str = ""
    sl = wildcards.dchrm.replace('chr','').split("-")
    sl2 = wildcards.dchrm.replace('chr','').split("~")
    
    if len(sl2) == 2:
        ret_str = pchr + wildcards.dchrm + ':'
    elif len(sl) == 1:
        ret_str = pchr + sl[0] + ':'
    elif len(sl) == 2:
        start = int(sl[0])
        end = int(sl[1])
        while start <= end:
            ret_str = str(ret_str) + "," + pchr + str(start) + ':'
            start = start + 1
    else:
        raise Exception(
            "sentD chunks can only be one contiguous range per chunk : ie: 1-4 with the non numerical chrms assigned 23=X, 24=Y, 25=MT"
        )
    mito_code="MT" if "b37" == config['genome_build'] else "M",

    return ret_mod_chrm(ret_str).lstrip(',').replace('chr23','chrX').replace('chr24','chrY').replace('chr25','chrMT').replace('23:','X:').replace('24:','Y:').replace('25:',f'{mito_code}:')


OG_ALIGNERS=list(set(ALIGNERS)-set(CRAM_ALIGNERS))
ALL_ALIGNERS=list(set(ALIGNERS+CRAM_ALIGNERS))