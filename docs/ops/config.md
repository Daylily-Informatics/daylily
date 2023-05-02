# Daylily Config

## Rulse vs. Config
  > The snakemake rules codify how each tool is executed and how files are created/named.  The *congfig* for each rule is stored seperately and managed via snakemake [profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles). 
  - Presently there are two: `local` and `slurm`. 
    - `local` should run on any '5' generation AWS instance w/at least 256G of ram. 
    - `slurm` requires the matching cluster be created as descrubed in the setup docs.
  - Config may be overriden on the command line.
