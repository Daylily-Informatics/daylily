#!/bin/bash

outf=$1
echo '''# DAY PHYTO__RESOURCES Supporting Data

<img src="https://via.placeholder.com/1000x5/5707D9/000000?text=+" valign="bottom" >

''' > $outf;
psize=`du -hc /lll/data/external_data/research_experiments/PHYTO__RESOURCES/ | tail -n 1`;
fandd=`find /lll/data/external_data/research_experiments/PHYTO__RESOURCES/ | wc -l`;
echo "There are:
* $fandd files and top level dirs.
* PHYTO__RESOURCES is $psize

## Under The Hood
" >> $outf
echo "<pre>" >> $outf;
tree  -n --sort=name --dirsfirst -sh    /lll/data/external_data/research_experiments/PHYTO__RESOURCES/  >> $outf ;
echo "</pre>" >> $outf;
