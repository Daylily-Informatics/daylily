#!/bin/bash

# Calculate 80% of total memory
TOTAL_MEM=$(grep MemTotal /proc/meminfo | awk '{print $2}')  # in KB
SHM_SIZE_KB=$((TOTAL_MEM * 80 / 100))  # 80% of total memory

# Convert KB to MB
SHM_SIZE_MB=$((SHM_SIZE_KB / 1024))

# Remount /dev/shm with the new size
mount -o remount,size=${SHM_SIZE_MB}M /dev/shm

# Verify new size
echo "/dev/shm resized to 80% of total memory:"
df -h /dev/shm
