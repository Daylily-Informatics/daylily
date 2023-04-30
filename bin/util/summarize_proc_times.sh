
file_n=$1

ls -1 $PWD/results/day/*/benchmarks/*| parallel -j1 -k """echo -n "{}XyyyX">> $file_n; tail -n 1 {} >> $file_n; """ ; perl -pi -e 's/XyyyX/\t/g' $file_n; perl -pi -e 's/(^.*\/results\/day\/)(.*)(\/benchmarks\/.*_[0-9]\.)(.*)(\.bench\.tsv)(.*$)/$2\t$4\t$6/g;' $file_n ; sort $file_n > $file_n.sort;
