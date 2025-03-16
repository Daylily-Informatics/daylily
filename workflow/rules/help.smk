import os
import sys

####### HELP INFO
# ---------------
#
# Help info re: daylily
#

cdsx = {"2": "20", "1": "128", "17": "33"} 


# I had fun with some terminal color and appearance hacking...

os.environ["TERM"] = "xterm-256color"

# Return the terminal width of the active term session
def get_term_width():
    try:
        return os.get_terminal_size().columns
    except Exception as e:
        del e
        return 1
	
# get the colr command
def get_ccmd():
    ccmd = "echo"
    res = os.popen("which colr").readlines()
    if len(res) == 1:
        ccmd = "colr"
    return ccmd

# method to standardize the colr message printing
def print_colr_msg(
    msg="", txt_colr="orchid4", bg_colr="chartreuse", style="n", sleep=0
):
    # Print a colorized ascii string
    #    msg can be any string w/out troublesome characters to quote
    #    txt_color and bg_colr can be a name, RGB or hex
    #     style 1=normal, 7=flash, 1=bright, 4=underline
    os.system(
        f"""(colr "{msg.replace('"','-')}" "{txt_colr}" "{bg_colr}" "{style}" 1>&2; sleep {sleep};) || exit 0 """
    )


# Pull out the target rules.
targetsa = os.popen(
    """grep -s rule workflow/Snakefile | grep '# TARGET' | cut -d ' ' -f 2 | cut -d : -f 1"""
).readlines()
targetsb = os.popen(
    """grep -s rule workflow/rules/*smk | grep '# TARGET' | cut -d ' ' -f 2 | cut -d : -f 1"""
).readlines()

targetsc = targetsa + targetsb

tgts = """

"""

# Caclc pad to add to colr target list so it spans term width
def _ret_pad_line(tgtl):
    w = get_term_width()
    b = len(tgtl) - 1
    if b < 1:
        b = 1
    ret_str = ""
    pad = ""
    try:
        pad_mult = float(abs(w)) * (
            1.0 - float("0." + str(float(b) / float(abs(w))).split(".")[-1])
        )
        if pad_mult > w or pad_mult < 1:

            pass
        else:
            ret_str = tgtl + " " * int(pad_mult)
    except Exception as e:
        print("WARNING, " + str(e), file=sys.stderr)
        ret_str = tgtl
    return ret_str

# compose the new list
for i in sorted(targetsc):
    if len(i) > 2:
        tgts = (
            tgts
            + _ret_pad_line(
                """<> """
                + i.rstrip()
                + os.popen(
                    f"""grep -s rule workflow/Snakefile workflow/rules/*smk | grep -s {i.rstrip()}: | cut -d '#' -f 2"""
                )
                .readline()
                .rstrip()
                .replace("TARGET", "")
            )
            + """
"""
        )

tgts = (
    tgts
    + _ret_pad_line("""<> NO TARGET SPECIFIED == all <-- not advised """)
    + _ret_pad_line(" ")
    + """
 """
    + "~" * int(get_term_width())
    + _ret_pad_line(
        """
It is highly reccomended to run your snakemake command with the n flag set prior to running hot. This sets dry run mode. This will confirm the steps you expect to occur are happening & maybe more importantly, it is able to detect most bugs in the DAG this way, saving you a headache if a long running job has a glitch well into running. <-- though, snakemake recovers from failures rather gracefully in this unfortunate scenario."""
    )
)

# clean up junk in this string before displaying
tgts.replace("(", "").replace(")", "").replace("-", "")


# Test what failures look like. This rule does not produce the expected
# output file when executed
localrules:
    util_test_fail,

rule util_test_fail:  # TARGET: rule which will always fail
    threads: 1
    output:
        "logs/this_file_is_never_created.txt",
    shell:
        "echo 'Running Rule Which Will Fail To Produce {output}';"


# This rule prints all rules tagged with TARGET
localrules:
    util_print_targets,

rule util_print_targets:  # TARGET : Return a List of All Targets
    threads: 1
    params:
        tgts=tgts,
        c=get_ccmd(),
    shell:
        """
        set +euo pipefail; echo AAAAAA$?;
        echo 'The Following Targets Are Available: ' 1>&2;
        {params.c} '''           {params.tgts}''' $DY_IT1 $DY_IB2 $DY_IS2  1>&2 ;
        {params.c} '   .....you may run (source dyinit && mod-activate PROFILE && mod-run [TARGETNAME, ...]  -p  -n) or just  (mod-run [TARGETNAME,]  -p  -n) if you are already init+active. This will show you what rules will be run with the analysis manifest detected, if none specified then the test data will be used.' $DY_IT1 $DY_IB1 $DY_IS1 1>&2;
        """


# method to demonstrate how to dynamically inject values for use in rules
def get_slots(wildcards):
    return 2


rule util_env_check:  # TARGET : Will create a few files in the results dir, useul in debuggig container issues, and submitq a job to the relevant eneironment scheduler.
    params:
        cluster_sample="na",
        cluster_slots=get_slots,
    resources:
        partition="i8,i192,i192mem",  # if partition is not set in rule, default from profile config is used
        threads=get_slots,  # if threads is not set in rule, default from profile config is used
    conda:
        "../envs/vanilla_v0.1.yaml"  # oddity- the conda envs are declared relative to the rules directory, not like everything else which is from the execution directory.
    log:
        MDIR + "logs/env_check.log",
    threads: 8
    output: 
        "logs/env_check.out",
    shell:
        """echo 'running!' 2>&1;
        echo 'blahblah' > {output};
	echo 'postsleep' 1>&2;
        touch {output};"""


# Pop out after initializing to an ipython shell, for debugging
localrules:
    util_ipython_shell,

rule util_ipython_shell:  # TARGET: opens an ipython embed session in the run block of the rule (this is the only run block in mg). you may open a global context ipython shell by setting --config ipython=yes.
    run:
        from IPython import embed
        embed()


# TODO: deprecate this
localrules:
    initm,

rule initm:
    """ Internal use. runs after mod-activate to force creation of template files for a new profile. """
    shell:
        "echo profile:$DAY_PROFILE.initialized;"


# TODO: deprecate this
localrules:
    util_profile_config_reset,

lb = "{"
rb = "}"

# This is a rare use of the 'run' directive.  They are not encouraged if you wish to eventually run in containers. See SMK docs for reasoning
rule util_profile_config_reset:  # TARGET: clears all profile config files so they can be regenerated from templates on next mod-run o r -activate
    """ Clears all of the initialized profile config files """
    """ The active profile config will be regenerated from templates w next 'mod-run'. """
    run:
        os.system(
            f"""echo $PWD; ((rm config/day_profiles/*/config.yaml && rm config/day_profiles/*/rule_config.yaml && rm config/day_profiles/*/cluster.yaml ) || colr 'no config to remove' "$DY_WT1" "$DY_WB1"); colr 'CONFIG FILES REMOVED' "$DY_WT0" "$DY_WB0" "$DY_WB0";"""
        )
        os.system(
            f"echo 'CHECKING CONFIG DIRS::::  Don not be alarmed when this crashes, it should if successful.' 1>&2;ls config/day_profiles/*/ 1>&1;"
        )
        os._exit(3)


# help returns all TARGET annotated rules in the rules dir


rule test_benchmark:
    container: None
    threads: 2
    params:
        cluster_sample="test_benchmark"
    resources:
        vcpu=2,
        threads=2,
        partition="i8,i128,i192,i192mem"
    benchmark:
        "logs/test_bench.tsv"
    shell:
        """
	sleep 3;
	"""

rule help:
    container: None
    threads: 1
    params:
        c=get_ccmd(),
        cluster_sample="help"
    resources:
        vcpu=1,
        threads=1,
        partition="i8,i192,i192mem",
    shell:
        """

        {params.c} '##### DY-CLI HELP                                                    ' "$DY_IT2" "$DY_IB2" "$DY_IS2" 1>&2;
        {params.c} '/-----------------------------------------------------------------------------' "$DY_IT1" "$DY_IB1" "$DY_IS1" 1>&2;
        {params.c} 'Welcome to daylily. 

        More complete docs can be found here: 
              https://github.com/Daylily-Informatics/daylily 
     

        ========================================================================
        The dy- cli , first source it:

        . dyinit
        . dyinit -h


        Day CLI initialized for project $PROJECT in region $region.
        The Daylily CLI is now available.
        The following commands are available:\n\t**TAB COMPLETION(single and double tabbing) IS ENABLED FOR ALL CLI COMMANDS AND FLAGS, tab following a space or following a - to see completions**

          (day-help / dy-h)         - Display this help message.
          (day-activate / dy-a) (slurm | local) - Activate a Slurm or local environment.
          (day-run / dy-r)           - Run a command in the current environment.
          (day-deactivate / dy-d)    - Deactivate the current environment.
                             Use dy-d reset to hard reset the environment (may kill your terminal!).
          (day-build / dy-b)         - Build a new environment. **Only needs to be run once per head node user.**


        ========================================================================
        For use with dy-r [rules, try produce_<tab><tab>] [flags, try --con<tab><tab>]
        + Important snakemake arguments (all optional):
          `-n` == dryrun, generates an execution plan, but runs no commands.  
              Good practice to run with -n before pulling the trigger.
          `-p` == print shell commands, pretty much always include this one
          `--rerun-triggers mtime` == this will calculate an execution plan 
              only re-running rules which have parent files timestamps are 
              newer than children.  w/out this, snakemake cohsiders not only 
              timestamps, but edits to config files, source data changes, code 
              changes, environment changes... and will re-run agressively.  
              ! Always run with -n before you attempt to re-run, or add, 
              analysis to an existing results directory.  You can easily 
              erase a lot of work if not careful.
          `-h` == snakemake help... `--help` == snakemake exhaustive help
          `--rerun-incomplete` == if rerunning a workflow which failed 
              ungracefully, snakemake will not proceed unless this flag is 
              set so it can cleanup  the incomplete results subdirs.
          `--jobs n` == limit the snakemake scheduler to only running 
              `n` jobs at once.
          `--report FILENAME.html` == produce a summary report of the workflow
          `--dag|--rulegraph|--filegraph` == use as follows to produce 
              various graph visualizations of the TARGET being asked for: 
              `dy-r produce_deduplicated_bams --rulegraph | dot -Tpng > rg.png`
         ========================================================================
                + TARGETS (which are snakemake rules that are a terminal point
                      in the larger DAG of rules). A few common ones:
                   - produce_deduplicated_bams == will run the specified 
                      aligner(s) and the specified dup marker (only one ddup 
                      presently can run at a time). aligners, ddupers, 
                      SNV&SV(multiple of each of these allowed) callers are 
                      all configured in the: 
                       config/day_profles/{{local,slurm}}/rule_config.yaml
                   - produce_snv_concordances == runs all specified 
                      aligners/SV callers and produces a concordance report 
                      for each aligner+SVcaller+sample (in the other_reports 
                      dir)
                   - produce_final_multiqc_wgs == runs everything and produces 
                      the big final MQC report in the `reports` dir. ' ;
                      
       {params.c} '/----------------------------------------------------------\
-------------------' "$DY_IT2" "$DY_IB2" "$DY_IS2"  1>&2;

    """


# Method to update jira by creating tickets with the following info
def update_jira(jid, targets, nsamps, details, url, c1, c2, c3, ex, ru, state):
    try:
        jcmd = f"""bin/register-with-jira.py '{jid}' '{details}' '{ex}' '{ru}' '{targets}' '{state}' '{url}' '{nsamps}' '{c1}' '{c2}' '{c3}' 1>&2"""
        print(jcmd, file=sys.stderr)
        os.system(jcmd)
    except Exception as e:
        print(e)
