ll=$(ls -1 .snakemake/log/*log | head -n 1) && grep 'steps' $ll | grep 'done' | perl -pe 's/\d+ of \d+ steps \((\d+\.*\d*)\%\) done/"\{ \"progress\": ".($1*.01)." \}"/eg;' || sleep 1
