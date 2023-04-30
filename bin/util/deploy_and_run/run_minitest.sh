#!/bin/bash

# Script to pull the repo, build the releavant conda env, then attempt to run a simple pipeline

if [[ "$@" == "" || "${*: 1}" == "-h" || "${*: 1}" == "--help" ]]; then
   echo "
   positional arguments are:
      1)  mod branch to clone for analysis
      2)  the analysis targets -if more than one, must be 'TGT1 TGT2 TGT3' -- IN QUOTES WITH SPACES BETWEEN
      3)  /path/to/directory/to/clone/repo/and/work/
            !! This directory will be create if it does not exist.
            !! Then the tagged version of mod will be cloned here
            !! the multiqc_aln_cov and seqqc modules will be run
      4) if set to 'yes' the base mod environment will be build (and avail for toyr future use.  If anything ele, ir's assumed tou have conda on your path and can activate an environment called DAY as deinfed in the repo yaml def.
   "
   exit 23
fi

dt=$(date +%Y%m%d.%H%M%S)
work_dir=$3/$dt/
echo "making the directory: $work_dir"
sleep 2mkdir -p $work_dir
cd $work_dir
tag=$1
echo "cloning mod $tag"

git clone --branch $tag  git@github.com:iamh2o/mod-pipes.git $tag
cd $tag

if [[ "$4" == "yes" ]]; then
    (bash environment/installconda_setupmod.sh DAY || bash environment/haveconda_setupmod.sh DAY) || echo "t\his should not error...and probably worked... so , lets see"


fi

source mginit
bin/day_activate lcwgs-drmaa
mod-run --jobs 5 --cores 60 -p -k -s workflow/Snakefile $2
mod-deactivate

echo DONE!

exit 0
