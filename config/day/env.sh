#conda activate || echo NOCONDA
#conda_bin=$(which conda)
#conda_dir=$(dirname $conda_bin)
#if [[ $conda_dir == "" ]]; then
#    echo" >>> ERROR: conda has not been found in your PATH.  If you do not have conda installed, you may run the mod installer+build command: source dyinit  --project <PROJECT> && mod-activate BUILD_all"
#fi




#. /opt/slurm/etc/slurm.sh

#export DAY_ROOT=$PWD
#export PATH=$PATH:$CONDA_PREFIX_1/bin/:$PWD/bin/
#export DRMAA_LIBRARY_PATH=resources/lib/libdrmaa.so.1 



return 0
