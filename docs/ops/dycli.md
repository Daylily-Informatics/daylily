# DY-CLI

- The basic interaction looks like:
```bash

# init the dy-cli
source dyinit -h
source dyinit  --project <PROJECT>

# activate the local or slurm running environment
dy-a local # or slurm  # tab complete for avail envs

# run commands
dy-r produce_snv_concordances -p # tab w/a preceding space to get all available analysis targets, tab immediately following a - to get snakemake command line flags

# clear dy-cli
dy-d reset

```


