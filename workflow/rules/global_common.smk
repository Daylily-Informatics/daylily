# noqa   # noqa
import snakemake  # noqa
from snakemake.utils import min_version
import datetime
import os
import random

config = config  # noqa   ## Unnecessary, this is merely to keep the linters happy
config['tmpsub'] = random.randint(1,10000000)

# #### GLOBAL COMMON.SMK
# -----------------------
# 


# #### MIN VERSION CHECK Snakemake
# --------------------------------
min_version("6.0.0")


## DEPRECATE BIOME STUFF
mg_biome = os.environ.get("DAY_BIOME", "AWSPC")
cfg_f = "config/global.yaml"
if mg_biome in ["AWSPC"]:
    cfg_f = "config/global_day.yaml"
elif mg_biome in ["OTHER"]:
    pass
else:
    raise Exception(
        f"\n\n\t\tERROR::: 'DAY_BIOME' was not set to one of LOCAL, found {mg_biome} instead.  This should be set correctly when running day-activate [PROFILE].  From the MGroot:$DAY_ROOT try running day-deactivate, then day_activate [PROFILE] .   echo $DAY_BIOME should then return an expected value, or there is a bug.\n\n"
    )


# Colr, if missing... add it
try:
    if int(os.popen("colr --version | wc -l").readline().rstrip()) != 1:
        os.system("pip install colr docopt")
except Exception as e:
    raise (
        str(e)
        + "\n\nERROR ERROR.  \n\t\tCan Not Proceed: packages colr and docopt are missing.  Please install with 'pip install colr docopt.'"
    )

# Check the genome_build is an allowable reference
PERMITTED_REFCODES = ["b37",  "hg38"]

if config['genome_build'] not in PERMITTED_REFCODES:
   print(f"The entered config.conf GENOME BUILD CODE {config['genome_build']} is not permitted: {config['genome_build']}, acceptable values are: {PERMITTED_REFCODES} .", file=sys.stdout)
   raise(f"The entered config.conf GENOME BUILD CODE {config['genome_build']} is not permitted: {config['genome_build']}, acceptable values are: {PERMITTED_REFCODES} .")
config["ref_code"] = config["genome_build"]


if config["ref_code"] not in PERMITTED_REFCODES:
    raise Exception(
        f" The provided code does not match {PERMITTED_REFCODES} .... should not get here now as we default to b37 at the moment."
    )
else:
    os.environ["DY_REF_CODE"] = config["ref_code"]
config["supporting_data"]["supporting_files_conf"] = config["supporting_data"]["supporting_files_conf"].replace("XREF_CODEX", config["ref_code"])


# ###################
# Prep necessary dirs (spring cleaning target)
# ---------------------

if "dirs" not in config:
    config["dirs"] = {}
if "operational" not in config["dirs"]:
    config["dirs"]["operational"] = {"bcl2fq": {}, "day": {}}
if "day" not in config["dirs"]:
    config["dirs"]["day"] = {}

config["dirs"]["operational"]["logs"] = "logs/"
config["dirs"]["operational"]["tmp"] = "tmp/"
config["dirs"]["operational"]["results"] = "results/"

config["dirs"]["operational"]["bcl2fq"]["root"] = "results/bcl2fq/"
config["dirs"]["operational"]["bcl2fq"]["logs"] = "results/bcl2fq/logs/"

config["dirs"]["day"]["root"] = f"results/day/{config['ref_code']}/"
config["dirs"]["day"]["logs"] = f"results/day/{config['ref_code']}/logs/"
config["dirs"]["day"][
    "reports"
] = f"results/day/{config['ref_code']}/reports/"
config["dirs"]["day"][
    "compiled_impute_results"
] = f"results/day/{config['ref_code']}/compiled_impute_results/"

# ##### create core dirs needed
for dirn in config["dirs"]["day"]:
    dirtomake = config["dirs"]["day"][dirn]
    os.system(f"mkdir -p {dirtomake}")

# ##### Expose common dirs as global vars, it's a little easier and less wordy.
MDIR = (
    config["dirs"]["day"]["root"].rstrip("/") + "/"
    if "mgsuffix" not in config  # oh, not sure i like this anymore
    else f"{config['mgsuffix']}/"
)  # if we are running duplicates in the day dir, make the dedicated directory
# dirs:

config["B2FQD"] = "results/bcl2fq"
os.system(f"mkdir -p {config['B2FQD']}")

config["MDIR"] = MDIR
MDIRlogd = MDIR + "logs/"  # config["dirs"]["day"]["logs"]
MDIRreportsd = MDIR + "reports/"  #  config["dirs"]["day"]["reports"]
MDIRlog = config["logs"]["day"]["top"]
IMPUTED_RESULTSD = (
    MDIR + "compiled_impute_results"
)  #  config["dirs"]["day"]["compiled_impute_results"]
config["mgrppt"] = os.environ.get("DAY_ROOT", "day_root_env_var_not_set")
os.system(f"mkdir -p {IMPUTED_RESULTSD}")


# Create some expected Dirs
for dirn in config["dirs"]["operational"]:
    os.system(f"mkdir -p {dirn}")

RDIR = config["dirs"]["operational"]["results"]
global_log = config["logs"]["operational"]["top"]

# runtime info message residence time
config["warn_err_sleep"] = 0.02

if "jid" not in config:
    config["jid"] = ""


# Pull in the other global smk files
include: "rule_common.smk"  # noqa
include: "supporting_data.smk"  # noqa


# ##### LOG some meta-meta info for the command exe
the_time_is = str(datetime.datetime.now()).split(".")[0]

os.system(
    f"""echo "THETIMEIS:{the_time_is}\nGitTag:{config['gittag']}\nGitHash:{config['githash']}\nSUBuser:{config['sub_user']}\nSUBhost:{config['sub_host']}\nSUBroot:{config['sub_root']}\nCallingCLI:{config['calling_cli']}\n" >> {config['logs']['operational']['about']}"""
)


# END OF GLOBAL SMK
