#/usr/bin/env bash

_scc() {
    sfg=$(echo "SAMPLES.skipUndetermined(if 0 will not process samples named Undetermined. if 1 will process . default 0):::   '  --config keep_undetermined=0    ',GENERAL.mode(must be set if remote_mode rule is called first. default no):::   '  --config remote_mode=yes|no    ',GLIMPSE.proceeedWfailedVCFS(if re-running a failed glimpse instance will *ONLY* proceed with samples who have produced a VCF file already. default no):::   '  --config use_existing_finvcfs=yes   ',HELP.ipythonembed(open global scope ipy embed. default no):::  ' --config ipython=yes  ',BCL2FASTQ.run(To run bcl2fq w/out manually editing config. these must all be set *AND* a bcl2fastq samplesheet is needed see code. or easier to hack rules_config.yaml really)::: '--config ru_path=path/to/ru_dir  ru_uid=RU123456 bcl2fq_output_path=/path/to/bcl2fq/output/',SAMPLES.manifest(specify an analysis manifest to copy to config/analysis_manifest.csv ....WILL FAIL IF A FILE IS PRESENT ALREADY... an analysis manifest instead of the reccomended use of config/analysis_manifest.csv. default is to stage your sheet to config/analysis_manifest.csv and have it be auto-detected ):::   ' --config  override_analysis_manifest=/path/to/s.csv   ',ZZZ.config_usage(you onlspecify --config 1x, and list all kv pairs after the first, ie:)::: --config make_it_red=yes num_genes=100 more=stuff   '  ")
    python -c "import os, sys; print('\nSpecial case config\n', file=sys.stderr);print(\"\n\n\".join(sorted(\"$sfg\".split(',')))+\"\n\n\", file=sys.stderr)" 1>&2
    echo "
##################################
# You May Hit Enter Safely All This Script Does Is Dump Some Rarely Needed Commmand Reminders
##################################

" 1>&1
    COMPREPLY+=($(compgen -W "press-return" "${COMP_WORDS[$COMP_CWORD]}"))

}

_dya() {
    if [[ ${#COMP_WORDS[@]} == "2" ]]; then
        COMPREPLY+=($(compgen -W "$(ls -d1 config/day_profiles/*/templates/profile.info | cut -d '/' -f 3 | perl -pe 's/\n/ /g;')" "${COMP_WORDS[$COMP_CWORD]}"))
    fi
}

_dyd() {
    if [[ ${#COMP_WORDS[@]} == "2" ]]; then
        COMPREPLY+=($(compgen -W "help" "${COMP_WORDS[$COMP_CWORD]}"))
    fi
}
_dyb() {
    if [[ ${#COMP_WORDS[@]} == "2" ]]; then
        COMPREPLY+=($(compgen -W "BUILD help" "${COMP_WORDS[$COMP_CWORD]}"))
    fi
}


_d() {
    if [[ ${#COMP_WORDS[@]} == "2" ]]; then
        COMPREPLY+=($(compgen -W "help" "${COMP_WORDS[$COMP_CWORD]}"))
    fi
}
_dyr() {
    if [[ ${COMP_WORDS[$COMP_CWORD]} == "" || ${COMP_WORDS[$COMP_CWORD]::1} != "-" ]]; then
        COMPREPLY+=($(compgen -W "$(grep -s -E '^ *rule.*TARGET' workflow{/Snakefile,/rules/*smk} | grep -s TARGET | perl -pe 's/(^.*rule +)(.+)(\: +#.*$)/$2/g;')" "${COMP_WORDS[$COMP_CWORD]}"))
    elif [[ ${COMP_WORDS[$COMP_CWORD]::2} == "--" || ${COMP_WORDS[$COMP_CWORD]::2} == "-" ]]; then
        COMPREPLY+=($(compgen -W "$(snakemake --help | grep -s -E '^  \-\-' | perl -pe 's/(^  *)([a-z|A-Z|0-9|\-]+)( *.*$)/$2/g' | sort | sed 's/\n/ /g')" \-\- "${COMP_WORDS[$COMP_CWORD]}"))
        #elif [[ "${COMP_WORDS[$COMP_CWORD]::1}" == "-" ]]; then
        #    COMPREPLY+=($(compgen -W "$( snakemake --help | grep -s -E '^  \-' | perl -pe 's/(^  *)([a-z|A-Z|0-9|\-]+)( *.*$)/$2/g' | sort | sed 's/\n/ /g'  )" \\- "${COMP_WORDS[$COMP_CWORD]}" ))
    else

        python -c "import os, sys; print('\nSpecial case config\n', file=sys.stderr);print(\"\n\n\".join(sorted(\"$sfg\".split(',')))+\"\n\n\", file=sys.stderr)"
    fi
}

complete -F _d day
complete -F _dya day-activate
complete -F _dya dy-a
complete -F _dya day_activate
complete -F _dyb day-build
complete -F _dyb day_build
complete -F _dyb dy-b
complete -F _dyd day-deactivate
complete -F _dyd day_deactivate
complete -F _dyd dy-d
complete -F _dyr day-run
complete -F _dyr dy-r
complete -F _scc day-script-cmds
