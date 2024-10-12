#!/bin/bash
# File: manage_xmrig.sh

# Define minimum idle threshold for CPUs to run XMRig
MIN_IDLE_THRESHOLD=1  # Minimum number of idle CPUs needed to start XMRig

# Get total number of CPUs on the system
NCPUS=$(nproc)

# Query Slurm to get the total number of CPUs requested by all jobs running on the node
# This will sum the CPUs for all jobs in the "RUNNING" state on the current node
USED_CPUS=$(squeue --noheader --format="%C" --state=R --nodelist=$(hostname) | awk '{sum += $1} END {print sum}')

# If there are no jobs running, set USED_CPUS to 0
if [ -z "$USED_CPUS" ]; then
  USED_CPUS=0
fi

# Calculate available CPUs for XMRig
MINE_CPUS=$(($NCPUS - $USED_CPUS))

# Ensure that available CPUs are not less than the minimum threshold
if [ "$MINE_CPUS" -lt "$MIN_IDLE_THRESHOLD" ]; then
  MINE_CPUS=0
fi

# Check if XMRig is running
XMRIG_PID=$(pgrep xmrig)

# If there are no available CPUs, stop XMRig
if [ "$MINE_CPUS" -le 0 ]; then
    if [ -n "$XMRIG_PID" ]; then
        kill -9 $XMRIG_PID
        echo "XMRig killed due to no available CPUs (used: $USED_CPUS, total: $NCPUS)"
    fi
else
    # If there are available CPUs and XMRig is not running, start XMRig with available CPUs
    if [ -z "$XMRIG_PID" ]; then
        nohup /fsx/miner/bin/miner_cmd_$(hostname).sh $MINE_CPUS > /tmp/xmrig_$(hostname).log 2>&1 &
        echo "XMRig started with $MINE_CPUS CPUs (used: $USED_CPUS, total: $NCPUS)"
    fi
fi
