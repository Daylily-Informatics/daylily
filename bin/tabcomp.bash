#!/usr/bin/env bash

# Tab completion script for Daylily CLI

# Helper function for special config commands
_scc() {
    local sfg="SAMPLES.skipUndetermined:::--config keep_undetermined=0
GENERAL.mode:::--config remote_mode=yes|no
GLIMPSE.proceeedWfailedVCFS:::--config use_existing_finvcfs=yes
HELP.ipythonembed:::--config ipython=yes
BCL2FASTQ.run:::--config ru_path=path/to/ru_dir ru_uid=RU123456 bcl2fq_output_path=/path/to/bcl2fq/output/
SAMPLES.manifest:::--config override_analysis_manifest=/path/to/s.csv
ZZZ.config_usage:::--config make_it_red=yes num_genes=100 more=stuff"

    COMPREPLY+=("press-return")
}

# Completion function for day-activate
_dya() {
    if [[ ${#COMP_WORDS[@]} == 2 ]]; then
        local profiles
        profiles=$(ls -d1 config/day_profiles/*/templates/profile.info 2>/dev/null | awk -F'/' '{print $(NF-2)}')
        COMPREPLY+=($(compgen -W "$profiles" -- "${COMP_WORDS[$COMP_CWORD]}"))
    fi
}

# Completion function for day-deactivate
_dyd() {
    if [[ ${#COMP_WORDS[@]} == 2 ]]; then
        COMPREPLY+=($(compgen -W "help" -- "${COMP_WORDS[$COMP_CWORD]}"))
    fi
}

# Completion function for day-build
_dyb() {
    if [[ ${#COMP_WORDS[@]} == 2 ]]; then
        COMPREPLY+=($(compgen -W "BUILD help" -- "${COMP_WORDS[$COMP_CWORD]}"))
    fi
}


# Completion function for day-run
_dyr() {
    if [[ "${COMP_WORDS[$COMP_CWORD]}" == "" || "${COMP_WORDS[$COMP_CWORD]:0:1}" != "-" ]]; then
        local targets
        # Extract targets, remove trailing colons, and ensure unique values
        targets=$(grep -s -E '^ *rule.*TARGET' workflow/{Snakefile,rules/*.smk} 2>/dev/null | awk '{print $2}' | sed 's/:$//' | sort -u)
        COMPREPLY+=($(compgen -W "$targets" -- "${COMP_WORDS[$COMP_CWORD]}"))
    else


        local options
        options=$(snakemake --help | grep -E '^  --' | awk '{print $1}' | tr -d ',')
        options+=" --keep-temp"  # replace with --notemp bc this is not a snakemake native flag
        COMPREPLY+=($(compgen -W "$options" -- "${COMP_WORDS[$COMP_CWORD]}"))
        
        local options
        options=$(snakemake --help | grep -E '^  --' | awk '{print $1}' | tr -d ',')
        COMPREPLY+=($(compgen -W "$options" -- "${COMP_WORDS[$COMP_CWORD]}"))
    fi
}

# Completion for set genome build
_dyg() {
    local options="hg38 b37"
    COMPREPLY=($(compgen -W "$options" -- "${COMP_WORDS[$COMP_CWORD]}"))
}



# Register completion functions
complete -F _dya day-activate dy-a day_activate
complete -F _dyb day-build day_build dy-b
complete -F _dyd day-deactivate day_deactivate dy-d
complete -F _dyr day-run dy-r
complete -F _scc day-script-cmds
complete -F _dyg gay-set-genome-build dy-g
