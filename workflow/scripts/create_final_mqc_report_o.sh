#!/usr/bin/env bash

odir2=$1   # MDIRreportsd,
fnamef=$2  #f"{MDIRreportsd}DAY_{RU[0]}_{EX[0]}_final_multiqc_hcwgs.html",
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


multiqc --interactive   -p -x '*pyc' -x '*.fastq.gz' -f -i "<h3>(ANALYSIS REPORT FOR <a href=https://www.google.com target=ruu >$ruu</a>, occuring in <a href=http://www.google.com >$exx</a>)</h3><br><h4>( vs reference $ref )</h4>" -b "<small><b>$ruu _ $exx  Analysis Set Details____</b><uL>((( gitBranch:$gbranch |||gitTag:$gtag |||gitCommit:$ghash |||DockerImage Describing Envs: $dkr )))<br>Data Links___<ul><li><a href=https://www.google.com >Top Level Analysis Directory For This Run </a><li><a href=https://www.google.com >Raw Summary Data Files For a Subset Of Reports.</a><li><a href=https://wwww.google.com target=aar >A list of every html/png/svg/pdf/jpg found within the inards of this batch- nothing more than a list, not necesarily exhaustive....data in other formats will not be gethers, but ma ybe present.</a></ul>_<ul>[[[<a href=https://google.com >(TBD JOB ID LINKOUT?)</a>||<a href=https://www.google.com >(search  RUN )</a>||<a href=https://www.google.com >(search RUN )</a>]]]</ul>OTHER____<ul><li><a href=https://github.com/Daylily-Informatics/daylily/blob/7e388cfe8698bd615baafa2c09b045d260fc1336/docs/markdown/digital_organization.md>DAY  digital organization </a><li><a href=https://www.google.com >OTHER SEARCH </a></ul>$url_root2 .... $url_root </small></ul>" --sample-filters $name_cfg -n $fnamef -o $odir2 --profile-runtime  -v -c $macro_cfg -c $micro_cfg  $odir2/..

exit 0;
