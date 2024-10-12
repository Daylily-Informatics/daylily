#!/bin/bash
# File: manage_xmrig.sh

# Define thresholds
MIN_IDLE_THRESHOLD=30  # Minimum idle CPU percentage to start XMRig

# Get the CPU idle percentage (using /proc/stat)
IDLE=$(awk '/^cpu / {idle=$5; total=$2+$3+$4+$5+$6+$7+$8+$9+$10} END {print (idle/total)*100}' /proc/stat)
NCPUS=$(nproc)

# Calculate available CPUs based on idle percentage
MINE_CPUS=$(echo "$NCPUS * ($IDLE / 100)" | bc -l)
MINE_CPUS=${MINE_CPUS%.*}  # Convert to an integer

# Ensure at least one CPU is used when restarting XMRig
if [ "$MINE_CPUS" -lt 1 ]; then
  MINE_CPUS=1
fi

# Check if XMRig is running
XMRIG_PID=$(pgrep xmrig)

if (( $(echo "$IDLE < $MIN_IDLE_THRESHOLD" | bc -l) )); then
    # If idle CPU percentage is too low and XMRig is running, stop it
    if [ -n "$XMRIG_PID" ]; then
        kill -9 $XMRIG_PID
        echo "XMRig killed due to low idle CPU percentage ($IDLE%)"
    fi
elif (( $(echo "$IDLE > $MIN_IDLE_THRESHOLD" | bc -l) )); then
    # If idle CPU percentage is high and XMRig is not running, start it with available CPUs
    if [ -z "$XMRIG_PID" ]; then
        nohup /fsx/miner/bin/miner_cmd_$(hostname).sh $MINE_CPUS > /tmp/xmrig_$(hostname).log 2>&1 &
        echo "XMRig started with $MINE_CPUS CPUs due to high idle CPU percentage ($IDLE%)"
    fi
fi
