#!/bin/bash

echo start

sleep_time=600
if [[ "$1" != "" ]]; then
    export sleep_time=$1
fi

sleep $sleep_time

echo end
