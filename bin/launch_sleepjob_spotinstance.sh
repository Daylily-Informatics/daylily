#!/usr/bin/env bash

usage() {
    echo "Launches /opt/slurm/bin/sleep_test.sh on a spot with 1 vCPU held."
    echo "Usage: $0 [-c comment] [-p partition] [-t threads]"
    echo "Defaults:"
    echo "  - comment: RnD"
    echo "  - partition: i8 (options: i8, i128, i192, i192mem, i128mem)"
    echo "  - threads: nproc"
    exit 1
}

# Set default values
COMMENT="RnD"
PARTITION="i8"
THREADS=$(nproc)

while getopts ":c:p:t:h" opt; do
    case ${opt} in
        c ) COMMENT=$OPTARG ;;
        p ) PARTITION=$OPTARG ;;
        t ) THREADS=$OPTARG ;;
        h ) usage ;;
        * ) usage ;;
    esac
done

echo "Submitting job:"
echo "  Comment    : $COMMENT"
echo "  Partition  : $PARTITION"
echo "  Threads    : $THREADS"

sbcmd= """sbatch --cpus-per-task "$THREADS" \
       --comment "$COMMENT" \
       --partition "$PARTITION" \
       /opt/slurm/bin/sleep_test.sh
"""

echo "running : $sbcmd"

$sbcmd


sinfocmd="sinfo"
echo "running sinfo command:"
echo $sinfo
$sinfo

sleep 0.5

sqcmd="squeue"

echo "running queue monitor command:"
echo $sqcmd
$sqcmd


