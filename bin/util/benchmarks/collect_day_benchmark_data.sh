#!/usr/bin/env bash

build=$DAY_GENOME_BUILD
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo " run from the DAY_ROOT dir, with no arguments, until day is build aware, then it will accept the build name as the first arg"
    exit 0
fi

if [[ "$build" == "" ]]; then
    echo "No build specified, please set DAY_GENOME_BUILD with dy-g [hg38|b37] and rerun"
    exit 1
fi

(mv results/day/$build/reports/benchmarks_summary.tsv results/day/$build/reports/benchmarks_summary.tsv.old) || echo "no file to move";
echo -n "sample  rule  " > results/day/$build/reports/benchmarks_summary.tsv;
cat $(ls -1 $PWD/results/day/$build/*/benchmarks/* | head -n 1)| head -n 1 >>  results/day/$build/reports/benchmarks_summary.tsv;
fd -p -L -t f -e .bench.tsv . results/day/$build  | parallel -j1 -k "echo -n '{}XyyyX' >> results/day/$build/reports/benchmarks_summary.tsv; tail -n 1 {} >> results/day/$build/reports/benchmarks_summary.tsv; " ;
perl -pi -e 's/XyyyX/\t/g' results/day/$build/reports/benchmarks_summary.tsv ;
perl -pi -e 's/(^.*\/results\/day\/)(.*)(\/benchmarks\/.*_[0-9]\.)(.*)(\.bench\.tsv)(.*$)/$2\t$4\t$6/g;' results/day/$build/reports/benchmarks_summary.tsv ;
perl -pi -e 's/ +/\t/g' results/day/$build/reports/benchmarks_summary.tsv ;
perl -pi -e 's/(^.*\/benchmarks\/)(.+?\.)(.*)(\.bench\.tsv)/$2\t$3/g;' results/day/$build/reports/benchmarks_summary.tsv ;
perl -pi -e 's/_DBC0_0\.//g;'  results/day/$build/reports/benchmarks_summary.tsv ;
echo "

Benchmark Collection Complete, see results/day/$build/reports/benchmarks_summary.tsv"
