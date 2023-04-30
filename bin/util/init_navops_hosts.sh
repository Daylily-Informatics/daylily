# For use in staging the NavOps Cluster hosts

if [[ "$1" == "95" ]]; then
    for i in {1..25}; do ii=$(qsub -pe smp 95 -N nFa -V -cwd ./bin/util/bash_sleep.sh ); echo $ii; done
elif [[ "$1" == "1" ]]; then
    for i in {1..25}; do id=$(qsub -terse -pe smp 1 -N nFb -V -cwd ./bin/util/bash_sleep.sh); echo $id ; done
elif [[ "$1" == "-95" ]]; then
    qstat | grep nFa | cut -d " " -f 6 | parallel 'qdel {}'
elif [[ "$1" == "-1" ]]; then
    qstat | grep nFb | cut -d " " -f 6 | parallel 'qdel {}'
elif [[ "$1" == "qmod" ]]; then
    qstat | cut -d " " -f 6 | parallel 'qmod -c {}'
elif [[ "$1" == "qw" ]]; then
    qstat | grep qw | cut -d " " -f 6 | parallel 'qdel {}'
elif [[ "$1" == "del" ]]; then
    qstat | cut -d " " -f 6 | parallel ' qdel {}'
fi
