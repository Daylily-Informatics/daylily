#!/usr/bin/env bash

odir2=$1
fnamef=$2
macro_cfg=$3
micro_cfg=$4
name_cfg=$5
gtag=$6
ghash=$7
gbranch=$8
jid=$9
ruu=${10}
exx=${11}
cluster_sample=${12}
ld_pre=${13}
url_root=${14}
search_root=${15}
dkr=${16}
ref=${17}
logb=${18}
url_root2=${19}
mdir=${20}

echo '<html><head><body><h1>A tally of all the data files and plots buried below.</h1>' > $odir2/all_the_links.html; bin/fd -E *multiqc* -E *Multiqc* -E *MULTIQC*  -t f  . results/mod/ | grep -E '(\.png|\.html|\.svg|\.pdf|\.png|\.tsv|\.csv|\.json)' | sort | perl -pe 's/^(.*)$/\<li\>\<a href\=https\:\/\/www\.google\.com\/URL_ROOT_HERE$1\> $1 \<\/a\>\<br\>/g;' >> $odir2/all_the_links.html; echo '</body></html>' >> $odir2/all_the_links.html; ccmd=$(echo "perl -pi -e 's/URL_ROOT_HERE/"$url_root"\//g;' $odir2/all_the_links.html ");  eval $ccmd;


multiqc --interactive   -p -x '*pyc' -x '*.fastq.gz' --sample-filters $name_cfg -n $fnamef -o $odir2 --profile-runtime  -v -c $macro_cfg -c $micro_cfg  $odir2/..

exit 0;
