#!/usr/bin/env bash

fofn=$1
ctrld=$2
cd=$3
cvcf=$4
sdf=$5
threads=$6
if [[ "$NUMEXPR_MAX_THREADS" == "" ]]; then
    echo "probably safe to allow threads to go unchanged"
elif [[  "$NUMEXPR_MAX_THREADS" != "" ]]; then
    if [ $threads -ge $NUMEXPR_MAX_THREADS ]; then
	threads=$(( $NUMEXPR_MAX_THREADS - 1 ))
    fi
fi
sample=$7
rtg=$8
ctr_max=$(($9))
if [[ "$ctr_max" == "" || "$ctr_max" == "0" ]]; then
    export ctr_max=$((1))
fi

echo 'sleep 1;' > $cd/concordance.fofn;
ctr = 0
rm -rf $(dirname $fofn) || echo nothingToDel;
mkdir -p $( dirname $fofn ) || sleep .1 ;
mkdir -p $( dirname $fofn )/logs;
export alt_name=$(dirname $ctrld/. | perl -pe 's/^.*\///g;' );
for ctld_i in $(ls -d $ctrld/*/); do
         echo "ZZZZZZZZ--$ctld_i --  ";
         export vcf="$ctld_i/$alt_name.vcf.gz";
         export bed="$ctld_i/$alt_name.bed";
	 echo "$vcf /// $bed";
         export subd=$(echo $ctld_i | perl -pe "s/(.*\/)(.*)(\/$)/\$2/g;");
	 echo "$subd";
	 rm -rf $cd/_$subd || echo nodel;
         export cmd="$rtg vcfeval --vcf-score-field=INFO.RFAQ_ALL --ref-overlap -e $bed -b $vcf -c $cvcf -o $cd/_$subd -t $sdf --threads  $threads";
         export  fin_cmd="workflow/scripts/parse-vcfeval-summary.py $cd/_$subd/summary.txt $sample $bed $subd  $alt_name $cd/_$subd/summary.txt";
         ccmd="($cmd >> $cd/$subd._a.err; $fin_cmd >> $cd/$subd._b.err) || echo FAILED___ &";
         echo "$ccmd" >> $cd/concordance.fofn;
	 if [[ "$ctr" == "$ctr_max" ]]; then
	     echo 'wait;' >> $cd/concordance.fofn;
	     echo 'sleep 10 &' >> $cd/concordance.fofn;
	     ctr=$((0))
	 fi
	 ctr=$(($ctr + 1))
done;
echo 'wait;' >> $cd/concordance.fofn;
