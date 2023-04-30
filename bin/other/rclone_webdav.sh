#!/bin/bash
# Hacked script to connect to the webdav server 


source activate DAY
export mnt=$1
export mnt_cache=".day_cache_"$HOSTNAME
echo $mnt
echo $mnt_cache
mkdir $mnt_cache || echo "$mnt_cache exists"

export logging=$2 || export logging=" -v "
config=$3

rclone \
	--config $config \
	--log-file .mod.log \
	--fast-list \
	--no-unicode-normalization \
	--no-check-certificate \
	--async-read=true \
	--use-mmap \
	--multi-thread-streams 1000 \
	--multi-thread-cutoff 40M \
	--tpslimit-burst 4 \
	--use-mmap \
	--no-traverse \
	--max-backlog 100000 \
	--checkers 500 \
	--transfers 500 \
	--buffer-size 768M \
	mount $logging phytoRO:/data/ $mnt \
	--no-modtime \
	--cache-dir $mnt_cache \
	--dir-cache-time 1ms \
	--poll-interval 0h \
	--vfs-cache-max-size 0M \
	--vfs-cache-mode writes \
	--vfs-read-chunk-size 128M \
	--vfs-read-chunk-size-limit 768M \
	--vfs-read-ahead 1024M \
	--read-only &
echo "$!" >$HOSTNAME\_mnt_pid.log
