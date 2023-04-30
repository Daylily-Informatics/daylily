#!/usr/bin/env bash


sleepwait=$1

if [[ "$1" == '-h' || "$1" == "--help" ]]; then
       echo "This script will display the newest and oldest qstat jobs for $USER, it will take an int as an argument for time between refreshes"
       echo " you need to have colr installed.... 'pip install docopt colr' ."
       echo ""
       exit 0
elif [[ "$1" == "" ]]; then
       sleepwait=5

fi

for i in {1..10000}; do
    setterm -linewrap off && nlines=$(( ($LINES - 8) / 2 )) && topstr="___NEWER___________"$(printf %"$(( $COLUMNS - 19 ))"s |tr " " "^") && sepp=$(printf %"$(( $COLUMNS - 0 ))"s |tr " " "=") && btmstr="__________OLDER___"$(printf %"$(( $COLUMNS - 18 ))"s |tr " " "v") && qstat | head -n 1 &&  qstat  | sort -k 8 | head -n $nlines && echo "" && colr "$topstr" "seagreen2" "maroon3" "b" &&  colr "$sepp" "coral" "chartreuse2" "b" && colr  "$btmstr"  "palegreen" "deeppink" "b" && echo "" && qstat | head -n 1 && qstat  | sort -k 8 | tail -n $nlines && setterm -linewrap on
    sleep $sleepwait
done;

echo ""
echo ""
echo ""
echo "Exiting after "$(( $cycles * $sleepwait ))" seconds"
echo ""

exit 0
