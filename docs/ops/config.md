# Configuring Daylily (WIP draft)
Daylily snakemake rules rarely contain any specific config vales.  Instead those are abstracted to snakemake [profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) found in 'config/day_profiles/{slurm,local}]'. Profiles allow running da `.smk` rules where the wrapper for each rule if made more flexible by storing config in `config.yaml` and `rule_config.yaml`. The slurm profile has some additional files to interface with the slurm executor.

## Profile `yaml` Files
### `config.yaml`
Holds config which is used by snakemake when executed.

### `rule_config.yaml`
Holds the per-rule config values

## Profiles
  - `local` should run on any '5' generation AWS instance w/at least 256G of ram. 
  - `slurm` requires the matching cluster be created as descrubed in the setup docs.

## Command line
the `--config` flag may be used multiple times, and takes a key=string|json value. Config set here has precedence over config set elsewhere.

