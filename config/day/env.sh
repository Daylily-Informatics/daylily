conda activate || echo NOCONDA
conda_bin=$(which conda)
conda_dir=$(dirname $conda_bin)
if [[ $conda_dir == "" ]]; then
    echo" >>> ERROR: conda has not been found in your PATH.  If you do not have conda installed, you may run the mod installer+build command: source dyinit  --project <PROJECT> && mod-activate BUILD_all"
fi


source $conda_bin/../etc/profile.d/conda.sh 2>/dev/null || source $CONDA_PREFIX_1/etc/profile.d/conda.sh 2>/dev/null || source ~/conda/etc/profile.d/conda.sh 2>/dev/null

conda activate
conda activate DAY
if [[ "$?" != "0" ]]; then
    echo "ERROR----"
    echo "   error-"
    echo "ERROR----"
    echo "The command >conda activate DAY failed.  Have you built the mod env?  >dy-b <."
else
    which colr > /dev/null || sleep 1
    if [[ "$?" != "0" ]]; then
	echo "please hold, pizzazz being added...."
	pip install  --no-input  --extra-index-url https://pypi.org/simple docopt colr > /dev/null 2>&1
	colr "Allright, all set~~" "oldlace" "hotpink" "f"
	sleep .9
    fi;
fi;

. /opt/slurm/etc/slurm.sh

export DAY_ROOT=$PWD
export PATH=$PATH:$CONDA_PREFIX_1/bin/:$PWD/bin/

export DRMAA_LIBRARY_PATH=resources/lib/libdrmaa.so.1 



return 0
