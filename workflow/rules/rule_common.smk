from snakemake.utils import min_version
from snakemake.utils import validate
import os
import pandas as pd
import sys
import yaml
import multiprocessing
import random
import shutil

## TODO: sweep through this thing and prune out the no longer needed/experimental stuff

ec2=os.popen('echo "EC2instance:$(bin/helpers/get_ec2_type.sh)"').readline().rstrip()
print(ec2, file=sys.stderr)

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


# SNV caller chunk arrays
SENTD_CHRMS = config["sentD"]["sentD_chrms"].split(",")
DEEPD_CHRMS = config["deepvariant"]["deep_chrms"].split(",")
OCTO_CHRMS = config["octopus"]["octo_chrms"].split(",")


# ##### Setting the allowed aligners to run and to which deduper to use.
# presently, 1+ aligners may run, but all must use the same deduper

snv_CALLERS = config["active_snv_callers"]
sv_CALLERS = config["active_s_v_callers"]


ALIGNERS_D = {}
ALIGNERS = []
DDUP = {}

# set day_i as the active alginer if in remote_mode(experimental!)
if "active_aligners" not in config and "day_i" in config:
    config["active_aligners"] = ["mgi"]
# config[''] = []

# Handle the cases where more than 1 aligner is specified.
for aligner in config["active_aligners"]:
    for ii in aligner:
        if ii in ALIGNERS_D:
            raise Exception(
                "\n\nError:  You have specified this aligner >1x in the rules.yaml file. Please correct the file in the config/day_profiles/"
                + os.environ.get("DAY_PROFILE", "na")
                + "/rules.yaml file. The problem caller is: "
                + ii
                + "\n\n"
            )
        else:
            ALIGNERS_D[ii] = True
        ALIGNERS.append(ii)

        if aligner[ii]["mkdup"] in ["sambamba", "samtools"]:
            DDUP["other"] = True
        elif aligner[ii]["mkdup"] in ["sent"]:
            DDUP["sent"] = True
        elif aligner[ii]["mkdup"] in ["dppl"]:
            DDUP["dppl"] = True
        elif aligner[ii]["mkdup"] in ["bbb2"]:
            DDUP["bbb2"] = True
        else:
            raise Exception(
                "\n\n\t ERROR ERROR:  The dedupe option set in the config file is unknown",
                aligner[ii]["mkdup"],
                "\n\n",
            )
        ALIGNERS = list(ALIGNERS_D.keys())
if len(DDUP) != 1:
    raise Exception("Only 1 dduper can run at a time", DDUP, "<---theese were set")


# DELETE THE FOLLOWING CRAP
# For multiqc sample naming, if multiple aligners or SNV callers produce data, then samples will have entries for multipletools, and the sample name must include more info so as not to collide.  This flags the most common uses case, 1 or each -or- flags if more than 1 of each so the multiqc regex is set correctly.
if len(ALIGNERS) == 1 and len(snv_CALLERS) == 1:
    # For multiqc, where to place '..' which indicates the end of the sample name, which for bakeoffs includes the aligner+snv code
    config["multiqc_sampname_cfg"] = "config/external_tools/multiqc_blank.yaml"
    os.system(
        'colr "____There is 1 ALIGNER and 1 SNV_CALLER SET, MQC NAMES WILL NOT BE EXPANDED", "$DY_IT2" "$DY_IB2" "$DY_IS2" 1>&2; sleep 0.15;'
    )

else:
    # Multiple of alignment or var calling or some other duplication in the pipe bfx analysis is happening such that sample names will be represented >1x per type of analysis duplicated. We wish to preserve this info so do not truncate the file names for ease of display in multiqc. Instead,they are left long to different treatments can be compared.
    os.system(
        'colr "____There are >1 ALIGNERS OR SNV_CALLERS SET, MQC NAMES WILL BE EXPANDED" "$DY_IT2" "$DY_IB2" "$DY_IS2" 1>&2; sleep 0.15;'
    )
    config["multiqc_sampname_cfg"] = "config/external_tools/multiqc_blank.yaml"

os.system(
    f"""colr 'aligners:{ALIGNERS} , DDUP:{DDUP}' "$DY_WT1" "$DY_WB1" "$DY_WS1" 1>&2 ; sleep {config['warn_err_sleep']} 1>&2;"""
)

os.system(
    f"""colr 'SNV Callers:{snv_CALLERS} ,  SV Callers:{sv_CALLERS}' "$DY_WT1" "$DY_WB1" "$DY_WS1" 1>&2 ; sleep {config['warn_err_sleep']} 1>&2;"""
)

# ####  LOAD SAMPLE DATA FROM ANALYSIS MANIFEST
# ---------------------
# Load sample sheet into pandas dataframe, then validate & if pass
# expose amples as a global  Including to any other rule or called script from a rule  included from the
# main snakefile

if "override_analysis_manifest" in config:
    os.system(
        """colr  '     _____ COMMAND LINE ANALYSIS MANIFEST SET. This will be copied to config/analysis_manifest.csv, and this copy used' "$DY_WT0" "$DY_WB0" "$DY_WS1"  """
    )
    if os.path.exists("config/analysis_manifest.csv"):
        raise Exception(
            "\n\n A file exists in config/analysis_manifest.csv, you must remove it to use the command line specified manifest-  be careful.\n\n"
        )
    else:
        os.system(
            f"cp {config['override_analysis_manifest']} config/analysis_manifest.csv"
        )
analysis_manifest = (
    "config/analysis_manifest.csv"
    if os.path.exists("config/analysis_manifest.csv")
    else config["analysis_manifest"]
)

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
    os.system(
        f'(colr " --config keep_undetermined=1 is set. Undetermined reads will be included." "$DY_IT2" "$DY_IB2" "$DY_IS2" 1>&2; sleep {config["warn_err_sleep"]})'
    )
else:
    os.system(
        f'(colr "NOTE! NOTE !! NOTE !!! ---- The Undetermined Sample Is Excluded. Set --config keep_undetermined=1 to process it." "$DY_IT1" "$DY_IB1" "$DY_IS1" 1>&2; sleep {config["warn_err_sleep"]})'
    )
    new_samples = samples.query("SQ != 'Undetermined'")
    samples = new_samples


# scan the input files and calculate the per sample total fastq size, for use in scaling compute.

pos_samp_vcf_bed = {}

config["samp_fq_size"] = {}

for s in samples["sample"]:
    if s not in config["samp_fq_size"]:
        config["samp_fq_size"][s] = 0.0
        for r1f in samples.loc[s, ("r1_path")]:
            try:
                config["samp_fq_size"][s] += float(
                    float(os.path.getsize(r1f)) / 1000000000.0
                )  # Norm to GB
            except Exception as e:
                print(e, file=sys.stderr)


# Set some informational items in the config dict
config["gittag"] = os.popen("git tag | tail -n 1 2>&1").readline().rstrip()
config["githash"] = (
    os.popen("git  --git-dir .git rev-parse --short HEAD").readline().rstrip()
)
config["gitbranch"] = os.popen("git branch | grep '*'").readline().rstrip()
cwd = os.path.abspath(".")
config["cwd"] = os.path.abspath(".")
config["sub_user"] = (
    "n/a"
    if str(os.environ.get("USER")) in ["", None, "None"]
    else str(os.environ.get("USER"))
)

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


# Set info for concordance calculations
CONCORDANCE_SAMPLES = {}
sample_info = {}
sample_lane_info = {}
for i in samples.iterrows():
    samp = i[0]
    sample = i[1][0]
    dat = dict(i[1])
    sample_lane = i[1][1]
    merge_single = i[1][14]
    if samp in sample_info and merge_single in ["single"]:
        raise (
            "\n\nMANIFEST ERROR:: "
            + samp
            + f"appears 2+ times in the sample sheet. This should only occur if 'merge' has been specified. {merge_single} has been set."
        )

    if merge_single in ["single"]:
        raise Exception(
            "\n\nMANIFEST ERROR: This feature was implemented, then unusued and has not been vetted to work properly again, so if you wish to run per lane, create a manifest with the sample and sample_lane column having the same id.  merge will create symlinks for each sample_lane pair of fastqs, then use these via process substitution to appear as one file for those tools expecting 1 R1 and 1 R2."
        )

    if sample_lane in sample_info or len(sample.split(".")) > 1:
        raise Exception(
            '\n\nMANIFEST ERROR: It is not allowed for the sample entry to be present in the sample_lane column and vice versa.  sample can be any string, excluding containing a period and ideally not having more than 4 "_", and does not need to be unique w/in the sample column. it may not end in "_[1-9]+", howevr, it may end in "_0"'
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
        return os.path.abspath(samples.loc[wildcards.sample, "r1_path"])
    elif wildcards.RR == "R2":
        return os.path.abspath(samples.loc[wildcards.sample, "r2_path"])
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


def ret_sx(wildcards):
    return wildcards.sx


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
        if ss_pct in ["", "na", 0, "0", None, "None"]:
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