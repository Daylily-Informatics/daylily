#!/bin/bash

#Move to environment dir if called from elsewhere
ABSOLUTE_PATH_SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
SCRIPT_DIR=$(python -c "import os; print(os.path.abspath(os.path.dirname(\"$ABSOLUTE_PATH_SCRIPT\")));")
AMBA_DISABLE_LOCKFILE=TRUE
echo "path to environment working dir is $SCRIPT_DIR"
#cd $SCRIPT_DIR

mkdir -p ~$USER/.parallel

################################################################################
# Sets up miniconda and other environments for imputation.
################################################################################
DY_ENVNAME="DAY"


if [[ $1 != "$DY_ENVNAME" ]]; then

    echo """
    Hello! This is the  __ $DY_ENVNAME ___ installation script.

    The DAY env installs the s/w needed to trigger snakemake and run the day (abbreviated to dy-) CLI.  The tools DAY can run have individual environments which are managed by snakemake. You should not need to worry about them(tm).  This code continues to run on both Navops and local

An important detail:

***one*** If you have an existing DAY install, the DAY install will be skipped ==> it must be rebuilt from scratch to save many headaches<==.  To remove the DAY env, activate conda, but not DAY, and run:
   `conda env remove -n DAY`

When that is complete, you may run this script.  If DAY is present and you run this scirpt-- nothing very meaningful should happen.

To run and start the install, rather than see this message, provide just one argument 'DAY'
"""
    return 2
fi

if [[ "$SHELL" != "/bin/bash" ]]; then
    echo "bash is the only supported shell, and you must have a ~/.bashrc file as well. If your bash path is something other than /bin/bash, this script will proceed, but no promises if the shell is not bash."
    sleep 6
fi

export CONDA_DIR="~/miniconda3/"
conda_detected=$(which conda)
if [[ "$conda_detected" != "" ]]; then
    export CONDA_DIR=$(dirname $conda_detected)

    echo "
     conda was detected in your path PATH ::  ( $conda_detected ).  


    "



else
    export CONDA_DIR="$HOME/miniconda3/"
    echo "|||   No conda environment detected.
 Installing to $CONDA_DIR :
    >    wget https://repo.anaconda.com/miniconda/Miniconda3-py39_23.1.0-1-Linux-x86_64.sh
    >    bash Miniconda3-latest-Linux-x86_64.sh -b -p  $CONDA_DIR
    >    rm Miniconda3-latest-Linux-x86_64.sh

    > source $CONDA_DIR/etc/profile.d/conda.sh
    > ~/conda/bin/conda init bash
    > source ~/.bashrc
    > conda activate
    "

    wget wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p  $CONDA_DIR
    rm Miniconda3-latest-Linux-x86_64.sh

    source $CONDA_DIR/etc/profile.d/conda.sh
    ~/conda/bin/conda init bash
    source ~/.bashrc

    conda activate

    echo "CONDA install complete."

fi

$CONDA_DIR/bin/conda init bash  # prob not necessary                                                                                                       


# Update Conda Config
conda config --add channels conda-forge
conda config --add channels bioconda


(conda config --set channel_priority strict && echo '' ) || echo 'failed to set conda pri to strict';
(conda config --set repodata_threads 10 || echo "" ) || echo "repodata_threads not set";
(conda config --set verify_threads 4 && echo '' )   || echo "verify threads not set";
(conda config --set execute_threads 4 && echo '' )   || echo "default threads not set";
(conda config --set always_yes yes && echo 'conda always yes set' ) || echo 'failed to set conda yes';
(conda config --set default_threads 10 && echo "DefThreads=44")      || echo 'failed to set conda default threads';

#  Install DAY If Not Present.

mgcnt=$(ls -d1 "$(dirname $(which conda))"/../envs/DAY | wc -l)
if [[ "$mgcnt" == "0" ]]; then
    echo "Installing DAY \(which will be sourced for you when you \'source dyinit  --project <PROJECT> \' \(see docs\)  you should not need to source MG directly.  Beginning install... "
    echo "..."
    sleep 1.4

    (export PIP_NO_INPUT=1 && export PIP_INDEX_URL="https://pypi.org/simple" && export PIP_EXTRA_INDEX_URL="https://pypi.org/simple"  &&  conda env create   -n DAY -f $SCRIPT_DIR/DAY.yaml && echo "DAY env created") || echo "DAY env not created"

    echo "Install exited with > $? < (if not zero, not the best sign)."
    echo "try the following:

       source dyinit  --project <PROJECT> ; dy-a local; dy-r help;

       #initialized the day CLI, activate the local env settings, runs day to get basic help info on targets
       "
    colr "


     .... And that is that" "red" "ivory" "b"

else
    dn="$(dirname $(which conda))"/../envs/
    echo "
       It appears you have a DAY env (located $dn ), or part of one at least.  You may need to manually remove the conda env dir for DAY, and hack a conf file (google it to be sure).  Then try this again."

fi;


echo "==="
echo "___"
echo 'you should log out of this shell and log back into a fresh shell, run source dyinit  --project <PROJECT> , dy-a *****, etc.'

return 0
