# The script function only accepts a path to a file, no log redirects or whatnot.
# This is a workaround to capture everything w/in the 'with' block to log

with open(snakemake.config["logs"]["global_log"], "w") as f:
    sys.stderr = sys.stdout = f
