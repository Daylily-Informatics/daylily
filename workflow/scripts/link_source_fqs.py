import os
import sys

##  TESTING

# The script function only accepts a path to a file, no log redirects or whatnot.
# This is a workaround to capture everything w/in the 'with' block to log
with open(snakemake.config["logs"]["global_log"], "w") as f:
    sys.stderr = sys.stdout = f

    # note- the object 'snakemake' is automatically imported and available with all of the
    # variables from the rule the script directive called from.  Presumably you must use the python interpreter
    # snakemake is running in as these tend not to have shebangs at the top. So:
    print(f"this is the input fqgz: {snakemake.input.fq}", file=sys.stderr)
    print(f"this is the string matching the sample wildcard: {snakemake.wildcards.sample} ", file=sys.stderr)
    print(f"number of threads requested: {snakemake.threads}", file=sys.stderr)

    # One main reason to use 'script' to invoke python code, is if the chunk of work won't fit neatly into a few lines in the snakefile. But, more importantly, the conda virtual environment management will not function via the 'run' directives, but will for python called via the scripts directive.  So, unless it's really simple and will not need it's own environment, scripts are likely the way to go for most substantial python.

    # This is all we actually came here for:
    cmd = f"ln -s {snakemake.input.fq} {snakemake.output.gz} > {snakemake.log} 2>&1"
    print(cmd, file=sys.stderr)
    os.system(cmd)
