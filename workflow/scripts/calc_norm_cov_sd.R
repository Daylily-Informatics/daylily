Script started on 2024-11-20 23:09:49+00:00 [TERM="screen" TTY="/dev/pts/1" COLUMNS="255" LINES="27"]
[?2004h(base) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ exitsudo su ubuntuexit[K[KR[3mscript workflow/scripts/calc_norm_cov_sd.R results/day/hg38/RIH0_ANA0-HG004-19_DBC0_0/align/bwa2a/alignqc/norm_cov_eveness/RIH0_ANA0-HG004-19_DBC0_0.bwa2a.md.chr$i.regions.[23m[3mb[23m[3med.gz  RIH0_ANA0-HG004-19_DBC0_0_$alnr_chr$i chr$i [23mM[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[Cscript workflow/scripts/calc_norm_cov_sd.R results/day/hg38/RIH0_ANA0-HG004-19_DBC0_0/align/bwa2a/alignqc/norm_cov_eveness/RIH0_ANA0-HG004-19_DBC0_0.bwa2a.md.chr$i.regions.bed.gz  RIH0_ANA0-HG004-19_DBC0_0_$alnr_chr$i chr$i 
[?2004l[?2004h(base) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ Rscript
[?2004lUsage: /path/to/Rscript [--options] [-e expr [-e expr2 ...] | file] [args]

--options accepted are
  --help              Print usage and exit
  --version           Print version and exit
  --verbose           Print information on progress
  --default-packages=list
                      Where 'list' is a comma-separated set
                        of package names, or 'NULL'
or options to R, in addition to --no-echo --no-restore, such as
  --save              Do save workspace at the end of the session
  --no-environ        Don't read the site and user environment files
  --no-site-file      Don't read the site-wide Rprofile
  --no-init-file      Don't read the user R profile
  --restore           Do restore previously saved objects at startup
  --vanilla           Combine --no-save, --no-restore, --no-site-file
                        --no-init-file and --no-environ

'file' may contain spaces but not shell metacharacters
Expressions (one or more '-e <expr>') may be used *instead* of 'file'
See also  ?Rscript  from within R
[?2004h(base) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ Rscript workflow/scripts/calc_norm_cov_sd.R results/day/hg38/RIH0_ANA0-HG004-19_DBC0_0/align/bwa2a/alignqc/norm_cov_eveness/RIH0_ANA0-HG004-19_DBC0_0.bwa2a.md.chr$i.regions.bed.gz  RIH0_ANA0-HG004-19_DBC0_0_$alnr_chr$i chr$i 
[?2004l[?2004h(base) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ Rscript workflow/scripts/calc_norm_cov_sd.R results/day/hg38/RIH0_ANA0-HG004-19_DBC0_0/align/bwa2a/alignqc/norm_cov_eveness/RIH0_ANA0-HG004-19_DBC0_0.bwa2a.md.chr$i.regions.bed.gz  RIH0_ANA0-HG004-19_DBC0_0_$alnr_chr$i chr$i [K[K[K[K[K[K[K[K[K[K[K[K[K[K
[?2004l[?2004h(base) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ Rscript workflow/scripts/calc_norm_cov_sd.R results/day/hg38/RIH0_ANA0-HG004-19_DBC0_0/align/bwa2a/alignqc/norm_cov_eveness/RIH0_ANA0-HG004-19_DBC0_0.bwa2a.md.chr$i.regions.bed.gz  RIH0_ANA0-HG004-19_DBC0_0_$aln[K[K[K[K[K[K[K[K[K[K[K[K
[?2004l[?2004h(base) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ [K(base) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ more [3mresults/day/hg38/RIH0_ANA0-HG004-19_DBC0_0/align/bwa2a/alignqc/norm_cov_eveness/RIH0_ANA0-HG004-19_DBC0_0.bwa2a.md.chr$i.regions.bed.gz [23m[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[Cresults/day/hg38/RIH0_ANA0-HG004-19_DBC0_0/align/bwa2a/alignqc/norm_cov_eveness/RIH0_ANA0-HG004-19_DBC0_0.bwa2a.md.chr$i.regions.bed.gz 
[?2004lmore: cannot open results/day/hg38/RIH0_ANA0-HG004-19_DBC0_0/align/bwa2a/alignqc/norm_cov_eveness/RIH0_ANA0-HG004-19_DBC0_0.bwa2a.md.chr.regions.bed.gz: No such file or directory
[?2004h(base) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ find results/ | grep giabhcr
[?2004lresults/day/MQ/SEQQC_multiqc_data/multiqc_[01;31m[Kgiabhcr[m[K_concordance.txt
results/day/MQ/DAY_final_multiqc_data/multiqc_[01;31m[Kgiabhcr[m[K_concordance.txt
results/day/MQ/SEQQC_multiqc_plots/svg/[01;31m[Kgiabhcr[m[K_concordance.svg
results/day/MQ/SEQQC_multiqc_plots/png/[01;31m[Kgiabhcr[m[K_concordance.png
results/day/MQ/SEQQC_multiqc_plots/pdf/[01;31m[Kgiabhcr[m[K_concordance.pdf
results/day/MQ/DAY_final_multiqc_plots/svg/[01;31m[Kgiabhcr[m[K_concordance.svg
results/day/MQ/DAY_final_multiqc_plots/png/[01;31m[Kgiabhcr[m[K_concordance.png
results/day/MQ/DAY_final_multiqc_plots/pdf/[01;31m[Kgiabhcr[m[K_concordance.pdf
results/day/hg38/other_reports/[01;31m[Kgiabhcr[m[K_concordance_mqc.tsv
results/day/hg38/reports/SEQQC_multiqc_plots/svg/[01;31m[Kgiabhcr[m[K_concordance.svg
results/day/hg38/reports/SEQQC_multiqc_plots/png/[01;31m[Kgiabhcr[m[K_concordance.png
results/day/hg38/reports/SEQQC_multiqc_plots/pdf/[01;31m[Kgiabhcr[m[K_concordance.pdf
results/day/hg38/reports/SEQQC_multiqc_data/multiqc_[01;31m[Kgiabhcr[m[K_concordance.txt
results/day/hg38/reports/DAY_final_multiqc_data/multiqc_[01;31m[Kgiabhcr[m[K_concordance.txt
[?2004h(base) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ find results/ | grep giabhcr[K[K [K[K | grep _mqc
[?2004lresults/day/hg38/other_reports/giab_concordance[01;31m[K_mqc[m[K.tsv
results/day/hg38/other_reports/giabhcr_concordance[01;31m[K_mqc[m[K.tsv
[?2004h(base) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ find results/ | grep giab | grep _mqc | parallel 'ls -lth {}'
[?2004l-rw-rw-r-- 1 ubuntu ubuntu 50K Nov 20 22:37 results/day/hg38/other_reports/giab_concordance_mqc.tsv
-rw-rw-r-- 1 ubuntu ubuntu 3.7K Nov 19 07:52 results/day/hg38/other_reports/giabhcr_concordance_mqc.tsv
[?2004h(base) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ rm [3mresults/day/hg38/other_reports/giabhcr_concordance_mqc.tsv[23mresults/day/hg38/other_reports/giabhcr_concordance_mqc.tsv
[?2004l[?2004h(base) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ rm results/day/hg38/other_reports/giabhcr_concordance_mqc.tsvfind results/ | grep giab | grep _mqc | parallel 'ls -lth {}'[K[C[1P[1P[1P[1P[1@n[1@o[1@r[1@m[1@_[1@C
[?2004l[?2004h(base) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ find results/ | grep norm_C | grep _mqc[1P[1@c
[?2004lresults/day/hg38/other_reports/norm_cov_evenness_combo[01;31m[K_mqc[m[K.tsv
results/day/hg38/other_reports/norm_cov_evenness[01;31m[K_mqc[m[K.tsv
[?2004h(base) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ find results/ | grep norm_c | grep _mqc | parallel 'ls -lth {}'
[?2004l-rw-rw-r-- 1 ubuntu ubuntu 27K Nov 20 21:32 results/day/hg38/other_reports/norm_cov_evenness_combo_mqc.tsv
-rw-rw-r-- 1 ubuntu ubuntu 4.5K Nov 19 08:32 results/day/hg38/other_reports/norm_cov_evenness_mqc.tsv
[?2004h(base) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ find results/ | grep norm_c | grep _mqc | parallel 'ls -lth {}'[C[Krm [3mresults/day/hg38/other_reports/norm_cov_evenness_mqc.tsv[23mresults/day/hg38/other_reports/norm_cov_evenness_mqc.tsv
[?2004l[?2004h(base) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ grep multiqc logs/slurm/multiqc_ff[Kinal_wgs/*
[?2004l[35m[Klogs/slurm/multiqc_final_wgs/multiqc_final_wgs.RIH0_ANA0-HG004-19.1.err[m[K[36m[K:[m[Krule [01;31m[Kmultiqc[m[K_final_wgs:
[35m[Klogs/slurm/multiqc_final_wgs/multiqc_final_wgs.RIH0_ANA0-HG004-19.1.err[m[K[36m[K:[m[K    output: results/day/hg38/reports/DAY_final_[01;31m[Kmultiqc[m[K.html
[35m[Klogs/slurm/multiqc_final_wgs/multiqc_final_wgs.RIH0_ANA0-HG004-19.1.err[m[K[36m[K:[m[K    benchmark: results/day/hg38/benchmarks/DAY_all.final_[01;31m[Kmultiqc[m[K.bench.tsv
[35m[Klogs/slurm/multiqc_final_wgs/multiqc_final_wgs.RIH0_ANA0-HG004-19.1.err[m[K[36m[K:[m[K    reason: Missing output files: results/day/hg38/reports/DAY_final_[01;31m[Kmultiqc[m[K.html
[35m[Klogs/slurm/multiqc_final_wgs/multiqc_final_wgs.RIH0_ANA0-HG004-19.1.err[m[K[36m[K:[m[K        [01;31m[Kmultiqc[m[K -f          --config config/external_tools/[01;31m[Kmultiqc[m[K_header.yaml          --config  config/external_tools/[01;31m[Kmultiqc[m[K_config.yaml          --filename results/day/hg38/reports/DAY_final_[01;31m[Kmultiqc[m[K.html         -i 'Final Multiqc Report'         -b 'Git Info branch:* main tag:yo hash:e61377b'         $(dirname results/day/hg38/logs/report_components_aggregated.done results/day/hg38/other_reports/rules_benchmark_data_mqc.tsv )/../ > results/day/hg38/reports/logs/all__mqc_fin_a.log 2>&1;
[35m[Klogs/slurm/multiqc_final_wgs/multiqc_final_wgs.RIH0_ANA0-HG004-19.1.err[m[K[36m[K:[m[K        ls -lt results/day/hg38/reports/DAY_final_[01;31m[Kmultiqc[m[K.html;
[35m[Klogs/slurm/multiqc_final_wgs/multiqc_final_wgs.RIH0_ANA0-HG004-19.1.out[m[K[36m[K:[m[K-rw-rw-r-- 1 ubuntu ubuntu 11441320 Nov 20 22:37 results/day/hg38/reports/DAY_final_[01;31m[Kmultiqc[m[K.html
[35m[Klogs/slurm/multiqc_final_wgs/multiqc_final_wgs.RIH0_ANA0-HG005-19.1.err[m[K[36m[K:[m[Krule [01;31m[Kmultiqc[m[K_final_wgs:
[35m[Klogs/slurm/multiqc_final_wgs/multiqc_final_wgs.RIH0_ANA0-HG005-19.1.err[m[K[36m[K:[m[K    output: results/day/hg38/reports/DAY_final_[01;31m[Kmultiqc[m[K.html
[35m[Klogs/slurm/multiqc_final_wgs/multiqc_final_wgs.RIH0_ANA0-HG005-19.1.err[m[K[36m[K:[m[K    benchmark: results/day/hg38/benchmarks/DAY_all.final_[01;31m[Kmultiqc[m[K.bench.tsv
[35m[Klogs/slurm/multiqc_final_wgs/multiqc_final_wgs.RIH0_ANA0-HG005-19.1.err[m[K[36m[K:[m[K    reason: Missing output files: results/day/hg38/reports/DAY_final_[01;31m[Kmultiqc[m[K.html, results/day/hg38/benchmarks/DAY_all.final_[01;31m[Kmultiqc[m[K.bench.tsv
[35m[Klogs/slurm/multiqc_final_wgs/multiqc_final_wgs.RIH0_ANA0-HG005-19.1.err[m[K[36m[K:[m[K        ( env bash workflow/scripts/create_final_mqc_report.sh 'results/day/hg38/reports/' 'DAY_final_[01;31m[Kmultiqc[m[K.html'         '/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily/config/external_tools/[01;31m[Kmultiqc[m[K.yaml' '/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily/config/external_tools/[01;31m[Kmultiqc[m[K_final.yaml'         '/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily/config/external_tools/[01;31m[Kmultiqc[m[K_blank.yaml' 'yo' '677105a'         '* main' '' 'RIH0'         'ANA0-HG005-19' 'RIH0_ANA0-HG005-19'         ' ' "$url_root"         'results/day/hg38/'          "none" "hg38" 'results/day/hg38/reports/logs/all_mqc_fin_b.log' "$url_root2"  >> results/day/hg38/reports/logs/all__mqc_fin_a.log  2>&1) || echo "MQC PROBLEM";
[35m[Klogs/slurm/multiqc_final_wgs/multiqc_final_wgs.RIH0_ANA0-HG005-19.1.err[m[K[36m[K:[m[K        ls -lt results/day/hg38/reports/DAY_final_[01;31m[Kmultiqc[m[K.html;
[35m[Klogs/slurm/multiqc_final_wgs/multiqc_final_wgs.RIH0_ANA0-HG005-19.1.out[m[K[36m[K:[m[K-rw-rw-r-- 1 ubuntu ubuntu 6325806 Nov 19 09:02 results/day/hg38/reports/DAY_final_[01;31m[Kmultiqc[m[K.html
[?2004h(base) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ conda actib[Kvate mqc
[?2004l[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ [3mmultiqc -f          --config config/external_tools/multiqc_header.yaml          --config  config/external_tools/multiqc_config.yaml          --filename results/day/hg38/reports/DAY_final_multiq[23m[3mc[23m[3m.html         -i 'Final Multiqc Report'         -b 'Git Info branch:* main tag:yo hash:e61377b'         $(dirname results/day/hg38/logs/report_components_aggregated.done results/day/hg38/other_reports/rules_benchmark_data_mqc.tsv )/../ > results/day/hg38/reports/logs/all__[23m[3mm[23m[3mqc_fin_a.log 2>&1;[23mMM[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[Cmultiqc -f          --config config/external_tools/multiqc_header.yaml          --config  config/external_tools/multiqc_config.yaml          --filename results/day/hg38/reports/DAY_final_multiqc.html         -i 'Final Multiqc Report'         -b 'Git Info branch:* main tag:yo hash:e61377b'         $(dirname results/day/hg38/logs/report_components_aggregated.done results/day/hg38/other_reports/rules_benchmark_data_mqc.tsv )/../ > results/day/hg38/reports/logs/all__mqc_fin_a.log 2>&1;
[?2004l[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ 
[?2004l[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ 
[?2004l[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ multiqc --help | more
[?2004l                                                                                                                                                                                                                                                                                  
 /// MultiQC 🔍 v1.25.1                                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                  
 Usage: multiqc [OPTIONS] [ANALYSIS DIRECTORY]                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                  
 MultiQC aggregates results from bioinformatics analyses across many samples into a single report.                                                                                                                                                                                
 It searches a given directory for analysis logs and compiles an HTML report. It's a general use tool, perfect for summarising the output from numerous bioinformatics tools.                                                                                                     
 To run, supply with one or more directory to scan for analysis results. For example, to run in the current working directory, use 'multiqc .'                                                                                                                                    
                                                                                                                                                                                                                                                                                  
╭─ Main options ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --force            -f  Overwrite any existing reports                                                                                                                                                                                                                          │
│ --config           -c  Specific config file to load, after those in MultiQC dir / home dir / working dir. (PATH)                                                                                                                                                               │
│ --cl-config            Specify MultiQC config YAML on the command line (TEXT)                                                                                                                                                                                                  │
│ --filename         -n  Report filename. Use 'stdout' to print to standard out. (TEXT)                                                                                                                                                                                          │
│ --outdir           -o  Create report in the specified output directory. (TEXT)                                                                                                                                                                                                 │
│ --ignore           -x  Ignore analysis files (GLOB EXPRESSION)                                                                                                                                                                                                                 │
│ --ignore-samples       Ignore sample names (GLOB EXPRESSION)                                                                                                                                                                                                                   │
│ --ignore-symlinks      Ignore symlinked directories and files                                                                                                                                                                                                                  │
│ --file-list        -l  Supply a file containing a list of file paths to be searched, one per row                                                                                                                                                                               │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Choosing modules to run ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --module   -m  Use only this module. Can specify multiple times. (MODULE NAME)                                                                                                                                                                                                 │
│ --exclude  -e  Do not use this module. Can specify multiple times. (MODULE NAME)                                                                                                                                                                                               │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Sample handling ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --dirs           -d   Prepend directory to sample names                                                                                                                                                                                                                        │
[3m--More--[23m│ --dirs-depth     -dd  Prepend n directories to sample names. Negative number to take from start of path. (INTEGER)                                                                                                                                                             │
│ --fullnames      -s   Do not clean the sample names (leave as full file name)                                                                                                                                                                                                  │
│ --fn_as_s_name        Use the log filename as the sample name                                                                                                                                                                                                                  │
│ --replace-names       TSV file to rename sample names during report generation (PATH)                                                                                                                                                                                          │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Report customisation ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --title            -i  Report title. Printed as page header, used for filename if not otherwise specified. (TEXT)                                                                                                                                                              │
│ --comment          -b  Custom comment, will be printed at the top of the report. (TEXT)                                                                                                                                                                                        │
│ --template         -t  Report template to use. (default|gathered|geo|sections|simple)                                                                                                                                                                                          │
│ --sample-names         TSV file containing alternative sample names for renaming buttons in the report (PATH)                                                                                                                                                                  │
│ --sample-filters       TSV file containing show/hide patterns for the report (PATH)                                                                                                                                                                                            │
│ --custom-css-file      Custom CSS file to add to the final report (PATH)                                                                                                                                                                                                       │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Output files ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --flat                    -fp  Use only flat plots (static images)                                                                                                                                                                                                             │
│ --interactive             -ip  Use only interactive plots (in-browser Javascript)                                                                                                                                                                                              │
│ --export                  -p   Export plots as static images in addition to the report                                                                                                                                                                                         │
│ --data-dir/--no-data-dir       Force the parsed data directory to be created.                                                                                                                                                                                                  │
│ --data-format             -k   Output parsed data in a different format. (tsv|csv|json|yaml)                                                                                                                                                                                   │
│ --zip-data-dir            -z   Compress the data directory.                                                                                                                                                                                                                    │
│ --no-report                    Do not generate a report, only export data and plots                                                                                                                                                                                            │
│ --pdf                          Creates PDF report with the 'simple' template. Requires Pandoc to be installed.                                                                                                                                                                 │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ MultiQC behaviour ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --verbose            -v  Increase output verbosity. (INTEGER RANGE)                                                                                                                                                                                                            │
│ --quiet              -q  Only show log warnings                                                                                                                                                                                                                                │
[3m--More--[23m│ --strict                 Don't catch exceptions, run additional code checks to help development.                                                                                                                                                                               │
│ --development,--dev      Development mode. Do not compress and minimise JS, export uncompressed plot data                                                                                                                                                                      │
│ --require-logs           Require all explicitly requested modules to have log files. If not, MultiQC will exit with an error.                                                                                                                                                  │
│ --profile-runtime        Add analysis of how long MultiQC takes to run to the report                                                                                                                                                                                           │
│ --profile-memory         Add analysis of how much memory each module uses. Note that tracking memory will increase the runtime, so the runtime metrics could scale up a few times                                                                                              │
│ --no-megaqc-upload       Don't upload generated report to MegaQC, even if MegaQC options are found                                                                                                                                                                             │
│ --no-ansi                Disable coloured log output                                                                                                                                                                                                                           │
│ --version                Show the version and exit.                                                                                                                                                                                                                            │
│ --help               -h  Show this message and exit.                                                                                                                                                                                                                           │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Options ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --clean-up/--no-clean-up    Remove the temporary directory and log file after finishing                                                                                                                                                                                        │
│ --no-version-check          Disable checking the latest MultiQC version on the server                                                                                                                                                                                          │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
                                                                                                                                                                                                                                                                                  
 See http://multiqc.info for more details.                                                                                                                                                                                                                                        

[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ git commit -a -m "X"
[?2004lgiu[main 6e79e7a] X
 Committer: AWS ParallelCluster user <ubuntu@ip-10-0-0-196.us-west-2.compute.internal>
Your name and email address were configured automatically based
on your username and hostname. Please check that they are accurate.
You can suppress this message by setting them explicitly. Run the
following command and follow the instructions in your editor to edit
your configuration file:

    git config --global --edit

After doing this, you may fix the identity used for this commit with:

    git commit --amend --reset-author

 2 files changed, 131 insertions(+), 27 deletions(-)
 rewrite workflow/scripts/calc_norm_cov_sd.R (99%)
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ giut [K[K[Kt push
[?2004lTo github.com:Daylily-Informatics/daylily.git
 [31m! [rejected]       [m main -> main (fetch first)
[31merror: failed to push some refs to 'github.com:Daylily-Informatics/daylily.git'
[m[33mhint: Updates were rejected because the remote contains work that you do[m
[33mhint: not have locally. This is usually caused by another repository pushing[m
[33mhint: to the same ref. You may want to first integrate the remote changes[m
[33mhint: (e.g., 'git pull ...') before pushing again.[m
[33mhint: See the 'Note about fast-forwards' in 'git push --help' for details.[m
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ git pull
[?2004lremote: Enumerating objects: 12, done.[K
remote: Counting objects:   8% (1/12)[Kremote: Counting objects:  16% (2/12)[Kremote: Counting objects:  25% (3/12)[Kremote: Counting objects:  33% (4/12)[Kremote: Counting objects:  41% (5/12)[Kremote: Counting objects:  50% (6/12)[Kremote: Counting objects:  58% (7/12)[Kremote: Counting objects:  66% (8/12)[Kremote: Counting objects:  75% (9/12)[Kremote: Counting objects:  83% (10/12)[Kremote: Counting objects:  91% (11/12)[Kremote: Counting objects: 100% (12/12)[Kremote: Counting objects: 100% (12/12), done.[K
remote: Compressing objects:  20% (1/5)[Kremote: Compressing objects:  40% (2/5)[Kremote: Compressing objects:  60% (3/5)[Kremote: Compressing objects:  80% (4/5)[Kremote: Compressing objects: 100% (5/5)[Kremote: Compressing objects: 100% (5/5), done.[K
remote: Total 12 (delta 7), reused 12 (delta 7), pack-reused 0 (from 0)[K
Unpacking objects:   8% (1/12)Unpacking objects:  16% (2/12)Unpacking objects:  25% (3/12)Unpacking objects:  33% (4/12)Unpacking objects:  41% (5/12)Unpacking objects:  50% (6/12)Unpacking objects:  58% (7/12)Unpacking objects:  66% (8/12)Unpacking objects:  75% (9/12)Unpacking objects:  83% (10/12)Unpacking objects:  91% (11/12)Unpacking objects: 100% (12/12)Unpacking objects: 100% (12/12), 1.09 KiB | 50.00 KiB/s, done.
From github.com:Daylily-Informatics/daylily
   8b8852f..3226b85  main       -> origin/main
[33mhint: You have divergent branches and need to specify how to reconcile them.[m
[33mhint: You can do so by running one of the following commands sometime before[m
[33mhint: your next pull:[m
[33mhint: [m
[33mhint:   git config pull.rebase false  # merge (the default strategy)[m
[33mhint:   git config pull.rebase true   # rebase[m
[33mhint:   git config pull.ff only       # fast-forward only[m
[33mhint: [m
[33mhint: You can replace "git config" with "git config --global" to set a default[m
[33mhint: preference for all repositories. You can also pass --rebase, --no-rebase,[m
[33mhint: or --ff-only on the command line to override the configured default per[m
[33mhint: invocation.[m
fatal: Need to specify how to reconcile divergent branches.
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ [3mgit config pull.rebase false[23mgit config pull.rebase false
[?2004l[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ git config pull.rebase falsepull[K
[?2004lhint: Waiting for your editor to close the file... [?2004h(B)0[?1049h[1;27r[m[4l[39;49m[?1h=[?1h=[?25l[39;49m[m[H[J[25;131H[0;7m[ Reading... ][m[25;130H[0;7m[ Read 6 lines ][m[H[0;7m  GNU nano 6.2                                                                                            /fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily/.git/MERGE_MSG                                                                                                     [1;273H[m[26d[0;7m^G[m Help[26;19H[0;7m^O[m Write Out[37G[0;7m^W[m Where Is[55G[0;7m^K[m Cut[26;73H[0;7m^T[m Execute[26;91H[0;7m^C[m Location[109G[0;7mM-U[m Undo[26;127H[0;7mM-A[m Set Mark[145G[0;7mM-][m To Bracket    [0;7mM-Q[m Previous[181G[0;7m^B[m Back[26;199H[0;7m^◂[m Prev Word[217G[0;7m^A[m Home[26;235H[0;7m^P[m Prev Line[253G[0;7mM-▴[m Scroll Up[27d[0;7m^X[m Exit[27;19H[0;7m^R[m Read File[37G[0;7m^\[m Replace[27;55H[0;7m^U[m Paste[27;73H[0;7m^J[m Justify[27;91H[0;7m^/[m Go To Line     [0;7mM-E[m Redo[27;127H[0;7mM-6[m Copy[27;145H[0;7m^Q[m Where Was[163G[0;7mM-W[m Next[27;181H[0;7m^F[m Forward[27;199H[0;7m^▸[m Next Word[217G[0;7m^E[m End[27;235H[0;7m^N[m Next Line[253G[0;7mM-▾[m Scroll Down[2dMerge branch 'main' of github.com:Daylily-Informatics/daylily[3d[36m# Please enter a commit message to explain why this merge is necessary,[4d# especially if it merges an updated upstream into a topic branch.[5d#[6d# Lines starting with '#' will be ignored, and an empty message aborts[7d# the commit.[2d[39m[m[34h[?25h[?25l[1;175H[0;7m*[273G[m[34h[?25h[2dxMerge branch 'main' of github.com:Daylily-Informatics/daylilyx[?25l7[2;24r8M[1;27r[2;1Hx[3d[1P[34h[?25h[?25l[25d[0;7mSave modified buffer?                                                                                                                                                                                                                                                             [26;1H Y[m Yes[K[27d[0;7m N[m No  [27;17H[0;7m^C[m Cancel[KMM[34h[?25h[?25l[34h[?25h[?25l[26d[0;7m^G[m Help[26;69H[0;7mM-D[m DOS Format[26;137H[0;7mM-A[m Append[26;205H[0;7mM-B[m Backup File[27d[0;7m^C[m Cancel[17G         [27;69H[0;7mM-M[m Mac Format[27;137H[0;7mM-P[m Prepend[27;205H[0;7m^T[m BrowseMM[0;7mFile Name to Write: /fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily/.git/MERGE_MSG[m[34h[?25h[?25l[25;130H[1K [0;7m[ Writing... ][m[K[1;175H[0;7m [273G[m[25;129H[0;7m[ Wrote 7 lines ][m[J[27d[34h[?25h[27;1H[?1049l[?1l>[?2004l[KMerge made by the 'ort' strategy.
 config/day_profiles/slurm/templates/rule_config.yaml | 2 [32m+[m[31m-[m
 config/external_tools/multiqc_config.yaml            | 2 [32m+[m[31m-[m
 etc/sample_sheet_template.tsv                        | 1 [32m+[m
 workflow/rules/calc_coverage_eveness.smk             | 4 [32m++[m[31m--[m
 4 files changed, 5 insertions(+), 4 deletions(-)
 create mode 100644 etc/sample_sheet_template.tsv
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ git pullconfig pull.rebase falsepull[K
[?2004l^[[A^[[A^[[A^[[AAlready up to date.
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ git pullconfig pull.rebase falsepull[Kshcommit -a -m "X"
[?2004l[main 6f17547] X
 Committer: AWS ParallelCluster user <ubuntu@ip-10-0-0-196.us-west-2.compute.internal>
Your name and email address were configured automatically based
on your username and hostname. Please check that they are accurate.
You can suppress this message by setting them explicitly. Run the
following command and follow the instructions in your editor to edit
your configuration file:

    git config --global --edit

After doing this, you may fix the identity used for this commit with:

    git commit --amend --reset-author

 1 file changed, 41 insertions(+), 1 deletion(-)
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ git commit -a -m "X"pull[Kconfig pull.rebase falsepull[Ksh
[?2004lEnumerating objects: 24, done.
Counting objects:   4% (1/22)Counting objects:   9% (2/22)Counting objects:  13% (3/22)Counting objects:  18% (4/22)Counting objects:  22% (5/22)Counting objects:  27% (6/22)Counting objects:  31% (7/22)Counting objects:  36% (8/22)Counting objects:  40% (9/22)Counting objects:  45% (10/22)Counting objects:  50% (11/22)Counting objects:  54% (12/22)Counting objects:  59% (13/22)Counting objects:  63% (14/22)Counting objects:  68% (15/22)Counting objects:  72% (16/22)Counting objects:  77% (17/22)Counting objects:  81% (18/22)Counting objects:  86% (19/22)Counting objects:  90% (20/22)Counting objects:  95% (21/22)Counting objects: 100% (22/22)Counting objects: 100% (22/22), done.
Delta compression using up to 8 threads
Compressing objects:   7% (1/14)Compressing objects:  14% (2/14)Compressing objects:  21% (3/14)Compressing objects:  28% (4/14)Compressing objects:  35% (5/14)Compressing objects:  42% (6/14)Compressing objects:  50% (7/14)Compressing objects:  57% (8/14)Compressing objects:  64% (9/14)Compressing objects:  71% (10/14)Compressing objects:  78% (11/14)Compressing objects:  85% (12/14)Compressing objects:  92% (13/14)Compressing objects: 100% (14/14)Compressing objects: 100% (14/14), done.
Writing objects:  14% (2/14)Writing objects:  21% (3/14)Writing objects:  28% (4/14)Writing objects:  35% (5/14)Writing objects:  42% (6/14)Writing objects:  50% (7/14)Writing objects:  57% (8/14)Writing objects:  64% (9/14)Writing objects:  71% (10/14)Writing objects:  78% (11/14)Writing objects:  85% (12/14)Writing objects:  92% (13/14)Writing objects: 100% (14/14)Writing objects: 100% (14/14), 5.80 KiB | 2.90 MiB/s, done.
Total 14 (delta 12), reused 0 (delta 0), pack-reused 0
remote: Resolving deltas:   0% (0/12)[Kremote: Resolving deltas:   8% (1/12)[Kremote: Resolving deltas:  16% (2/12)[Kremote: Resolving deltas:  25% (3/12)[Kremote: Resolving deltas:  33% (4/12)[Kremote: Resolving deltas:  41% (5/12)[Kremote: Resolving deltas:  50% (6/12)[Kremote: Resolving deltas:  58% (7/12)[Kremote: Resolving deltas:  66% (8/12)[Kremote: Resolving deltas:  75% (9/12)[Kremote: Resolving deltas:  83% (10/12)[Kremote: Resolving deltas:  91% (11/12)[Kremote: Resolving deltas: 100% (12/12)[Kremote: Resolving deltas: 100% (12/12), completed with 5 local objects.[K
remote: Bypassed rule violations for refs/heads/main:[K
remote: 
remote: - Changes must be made through a pull request.[K
remote: 
remote: - Cannot update this protected ref.[K
remote: 
To github.com:Daylily-Informatics/daylily.git
   3226b85..6f17547  main -> main
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ more .[Kdyinit 
[?2004l# Ensure the script is sourced                                                                                                                          
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "Error: This script must be sourced, not executed directly. Use 'source $0' to run."
    return 3
fi


. "/etc/parallelcluster/cfnconfig"
region=$cfn_region

# Ensure the script is sourced, not executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  echo "Error: This script must be sourced, not executed directly."
  echo "Usage: source $0 --project <project> [--skip-project-check]"
  return 3
fi


# Check if $cfn_region is unset or null
if [[ -z "$region" ]]; then
  echo "Error: cfn_region is not set or is empty."
  return 3
fi

# Function to display usage and return
usage() {
[3m--More--(9%)[23m  echo "This script initializes the Day CLI."
  echo "Region is autodetected as $region"
  echo "Usage: $0 [--project <project_name>] [--skip-project-check]"
  echo " ... if --project is not specified, the defauly budget for the region is used: daylily-omics-analysis-${region}."
  echo "Valid projects: $(grep "^${USER}=" /opt/slurm/etc/projects_list.conf | cut -d'=' -f2)"
  echo "Valid regions: Use a valid AWS region (e.g., us-east-1, us-west-2)."
  return 3
}



# Check if a command exists
command_exists() {
  command -v "$1" &> /dev/null
}

# Ensure AWS CLI is installed
if ! command_exists aws; then
  echo "Error: AWS CLI is not installed. Please install it first."
  return 3
fi

# Parse input arguments

# Parse input arguments
SKIP_PROJECT_CHECK=false
[3m--More--(18%)[23mwhile [[ "$#" -gt 0 ]]; do
  case "$1" in
    --project) PROJECT="$2"; shift 2;;  # Move past the flag and its argument
    --skip-project-check) SKIP_PROJECT_CHECK=true; shift 1;;  # Move past the flag
    *) echo "Unknown parameter passed: $1"; usage; return 3;;  # Handle unknown flags
  esac
done


# Set default project name if not provided
if [[ -z "$PROJECT" ]]; then
    if [[ -z "$region" ]]; then
        echo "Error: Region is not set. Please set the 'region' variable."
        return 3
    fi
    PROJECT="daylily-omics-analysis-${region}"
    echo "Notice: --project not set. Using default project name: $PROJECT"
fi

# Example output to confirm the parsed values
echo "Project: $PROJECT"
echo "Skip Project Check: $SKIP_PROJECT_CHECK"

# Use environment variables if flags are not set

export DAY_PROJECT=$PROJECT
[3m--More--(28%)[23m[K
export DAY_AWS_REGION=$region

# Validate region input
validate_region() {
  if [[ "$1" == "us-west-2" ]]; then
    echo "Region '$1' confirmed as valid."
    return 0
  else
    echo "Warning: Region '$1' is not 'us-west-2'."
    echo "It is recommended to use 'us-west-2' for this operation."
    return 0
  fi
}

# Call the region validation function
validate_region "$region" || return 3

AWS_ACCOUNT_ID=$(aws sts get-caller-identity --query "Account" --output text)

# Skip project check if the flag is provided
if [[ "$SKIP_PROJECT_CHECK" == false ]]; then
  # Ensure the project is valid for the current user
  USER_PROJECTS=$(grep "^${USER}=" /opt/slurm/etc/projects_list.conf | cut -d'=' -f2 | tr -d ' ')
  if [[ ! ",${USER_PROJECTS}," =~ ",${PROJECT}," ]]; then
    echo "Error: Project '$PROJECT' is not valid for user '$USER'."
[3m--More--(38%)[23m[K[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ more dyinit [1P[1P[1P[1P[1@e[1@m[1@a[1@c[1@s
[?2004l[?1049h[34l[?1h=[H[J[27d[K[?1l>[34h[?25h[?1049l[39;49m[?1049h[34l[?1h=[H[J[26d[?25l[3m-UUU:----F1  [0m[39;49m[23m[3m[1m*scratch*   [0m[39;49m[23m[3m   All L1     (Fundamental) -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------[0m[39;49m[23m
M[2d[34h[?25h[34l[27d[?25lLoading /etc/emacs/site-start.d/00debian.el (source)...done[K[H
[34h[?25h[34l[27;33H[?25l50autoconf.el (source)...[K[H
[34h[?25h[34l[27;58H[?25ldone[H
[34h[?25h[34l[27;35H[?25ldictionaries-common.el (source)...[H
[34h[?25h[34l[27;9H[?25ldebian-ispell...[K[H
[34h[?25h[34l[27;9H[?25l/var/cache/dictionaries-common/emacsen-ispell-default.el (source)...[H
[34h[?25h[34l[27;77H[?25ldone[H
[34h[?25h[34l[27;9H[?25ldebian-ispell...done[K[H
[34h[?25h[34l[27;9H[?25l/var/cache/dictionaries-common/emacsen-ispell-dicts.el (source)...[H
[34h[?25h[34l[27;75H[?25ldone[H
[34h[?25h[34l[27;10H[?25letc/emacs/site-start.d/50dictionaries-common[6P[H
[34h[?25h[34l[27;35H[?25ltcsh.el (source)...[K[H
[34h[?25h[34l[27;54H[?25ldone[H
[34h[?25h[34l[>4;1m[?2004h[?1004h[27d[?25lFor information about GNU Emacs and the GNU system, type C-h C-a.[K[H
[34h[?25h[34l[27d[?25lFor information about GNU Emacs and the GNU system, type C-h C-a.[K[H[3mFile Edit Options Buffers Tools Help                                                                                                                                                                                                                                              [0m[39;49m[23m
M
# Ensure the script is sourced[K
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then[K
    echo "Error: This script must be sourced, not executed directly. Use 'source $0' to run."[K
    return 3[K
fi[K
[K
[K
. "/etc/parallelcluster/cfnconfig"[K
region=$cfn_region[K
[K
# Ensure the script is sourced, not executed directly[K
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then[K
  echo "Error: This script must be sourced, not executed directly."[K
  echo "Usage: source $0 --project <project> [--skip-project-check]"[K
  return 3[K
fi[K
[K
[K
# Check if $cfn_region is unset or null[K
if [[ -z "$region" ]]; then[K
  echo "Error: cfn_region is not set or is empty."[K
  return 3[K
fi[K
[K
[3m-UU-:----F1  [0m[39;49m[23m[3m[1mdyinit      [0m[39;49m[23m[3m   Top L1    [0m[39;49m[23m[3mGit-main[0m[39;49m[23m[3m  (Fundamental) --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------[0m[39;49m[23m
M[2d[34h[?25h[34l[27d[?25l[36mI-search: [39;49m[K[H
[34h[?25h[34lM				[?25l[3mIsearch Help[0m[39;49m[23m[26;61H[3m Isearch) [0m[39;49m[23m[H
[34h[?25h[34l[27;11H[?25lA[C[34h[?25h[34l[H
[3;11H[?25l[36m[45mA[39;49m[26;34H[3m2[0m[39;49m[23m[3;12H[34h[?25h[34l[27d[?25lW[C[34h[?25h[34l[3;12HM[?25lif [[ -z "$region" ]]; then[K
  echo "Error: cfn_region is not set or is empty."[4;3Hreturn 3[K
fi[K
[K
# Function to display usage and return
usage() {
  echo "This script initializes the Day CLI."
  echo "Region is autodetected as $region"[11;3Hecho "Usage: $0 [--project <project_name>] [--skip-project-check]"
  echo " ... if --project is not specified, the defauly budget for the region is used: daylily-omics-analysis-${region}."
  echo "Valid projects: $(grep "^${USER}=" /opt/slurm/etc/projects_list.conf | cut -d'=' -f2)"
	Valid regions: Use a valid [36m[45mAW[39;49mS region (e.g., us-east-1, us-west-2)."[15;3Hreturn 3[K
}[K
[K[20;12Ha command exists[K
command_exists() {[K[22;3Hcommand -v "$1" &> /dev/null[K
}[K
[K
# Ensure AWS CLI is installed
[3m 8% L32[0m[39;49m[23m[14;38H[34h[?25h[34l[25;10H[?25l[46mAW[39;49m[14;38H[34h[?25h[34l[27;13H[?25lS[C[34h[?25h[34l[14;38H[?25l[36m[45mS[39;49m[34h[?25h[34l[27;14H[?25l_	[34h[?25h[34l[14;39H[H
[?25l    PROJECT="daylily-omics-analysis-${region}"[3;3H  echo "Notice: --project not set. Using default project name: $PROJECT"
fi[K
[K
# Example output to confirm the parsed values
echo "Project: $PROJECT"[K
echo "Skip Project Check: $SKIP_PROJECT_CHECK"
[K
# Use environment variables if flags are not set
[K
export DAY_PROJECT=$PROJECT[K
[K
export DAY_[36m[45mAWS_[39;49mREGION=$region[K
[K
# Validate region input
validate_region() {[18;3Hif [[ "$1" == "us-west-2" ]]; then[19;5Hecho "Region '$1' confirmed as valid."
    return 0[K
  else[K
  echo "Warning: Region '$1' is not 'us-west-2'."
    echo "It is recommended to use 'us-west-2' for this operation."[24;5Hreturn 0
  fi[K[26;29H[3m25% L80[0m[39;49m[23m[14;16H[34h[?25h[34l[27;15H[?25lA	[34h[?25h[34l[14;16H[2;5H[?25lecho "Region '$1' confirmed as valid."[K[3;5Hreturn 0[K
  else
echo "Warning: Region '$1' is not 'us-west-2'."
    echo "It is recommended to use 'us-west-2' for this operation."
    return 0[K
  fi[K
}
[K
# Call the region validation function
validate_region "$region" || return 3

[36m[45mAWS_A[39;49mCCOUNT_ID=$(aws sts get-caller-identity --query "Account" --output text)[16;3HSkip project check if the flag is provided
if [[ "$SKIP_PROJECT_CHECK" == false ]]; then[18;3H# Ensure the project is valid for the current user[19;3HUSER_PROJECTS=$(grep "^${USER}=" /opt/slurm/etc/projects_list.conf | cut -d'=' -f2 | tr -d ' ')[20;3Hif [[ ! ",${USER_PROJECTS}," =~ ",${PROJECT}," ]]; then[21;3H  echo "Error: Project '$PROJECT' is not valid for user '$USER'."[22;11HValid projects for $USER: $USER_PROJECTS[23;5Hreturn 3[K[24;3Hfi[K
else[26;29H[3m30% L97[0m[39;49m[23m[14;6H[34h[?25h[34l[27;1H[?25lMark saved where search started[14;6H[34h[?25h[34l[1;33H[?25l[3mHelp        [0m[39;49m[23m[14;1HAWS_A[26;35H[3m6   [0m[39;49m[23m[3mGit-main[0m[39;49m[23m[3m  (Fundamental) --------[0m[39;49m[23m[13;1H[34h[?25h[34l[27d[K[26;35H[?25l[3m7[0m[39;49m[23m[14;1H[34h[?25h[34l[?25l[1@e[26;6H[3m**[0m[39;49m[23m[14;2H[34h[?25h[34l[?25l[1@x[34h[?25h[34l[?25l[1@p[34h[?25h[34l[?25l[1@o[34h[?25h[34l[?25l[1@r[34h[?25h[34l[?25l[1@t[34h[?25h[34l[?25l[1@ [34h[?25h[34l[27;1H[?25lSaving file /fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily/dyinit...[14;8H[34h[?25h[34l[27;1H[?25lWrote /fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily/dyinit[K[14;8H[34h[?25h[34l[26;6H[?25l[3m----F1  [0m[39;49m[23m[3m[1mdyinit      [0m[39;49m[23m[3m   30% L97   [0m[39;49m[23m[3mGit:[0m[39;49m[23m[14;8H[34h[?25h[34l[27;1H[?25l[36mI-search: [39;49m[K[14;8H[34h[?25h[34l[1;33H[?25l[3mIsearch Help[0m[39;49m[23m[26;61H[3m Isearch) [0m[39;49m[23m[14;8H[34h[?25h[34l[27;11H[?25lD[C[34h[?25h[34l[14;8H[14;21H[?25l[36m[45mD[39;49m[34h[?25h[34l[27;12H[?25lA[C[34h[?25h[34l[14;22H[H
[?25lecho "AWS Budget for project '$PROJECT' in region '$region':"
echo "  Total: $TOTAL USD"
echo "  Used: $USED USD"
echo "  Percent Used: $PERCENT_USED%"[K
echo "________________________________________________________"[K
sleep 1.3[K
[K
# Check if the budget is exhausted
if (( $(echo "$PERCENT_USED >= 100" | bc -l) )); then
  echo "Warning: Budget for project '$PROJECT' is exhausted!"
fi[K

## export [36m[45mDA[39;49mY_PROJECT="daylily-dev"[K
export APPTAINER_HOME=/fsx/resources/environments/containers/$USER/$(hostname)/
[K
# needed by daycli[K
export DAY_BIOME=AWSPC[K
[K
export DAY_ROOT=$PWD[K
[K
ulimit -n 16384[K
[K
export SENTIEON_TMPDIR='/fsx/scratch/'
[Cxport SENTIEON_LICENSE='/fsx/data/cached_envs/Rare_Cancer_Research_Foundation_eval.lic'[26;29H[3m63% L173[0m[39;49m[23m[14;13H[34h[?25h[34l[18;8H[?25l[46mDA[39;49m

[46mDA[39;49m[14;13H[34h[?25h[34l[27d[?25lY[C[34h[?25h[34l[14;13H[?25l[36m[45mY[39;49m[34h[?25h[34l[18;10H[?25l[46mY[39;49m

[46mY[39;49m[14;14H[34h[?25h[34l[27d[?25l_	[34h[?25h[34l[14;14H[?25l[36m[45m_[39;49m[34h[?25h[34l[27d[?25lP	[34h[?25h[34l[14;15H[?25l[36m[45mP[39;49m[34h[?25h[34l[18;8H[?25lDAY

DAY[14;16H[34h[?25h[34l[27;1H[?25l[36mFailing case-sensitive I-search: [39;49m[41mDAY_P[39;49m[14;16H[34h[?25h[34l[27;1H[?25l[36mWrapped I-search: [39;49mDAY_P[K[14;16H[34h[?25h[34l[H
[?25l        return 3[K
    fi[K
    PROJECT="daylily-omics-analysis-${region}"
    echo "Notice: --project not set. Using default project name: $PROJECT"
fi[K
[K
# Example output to confirm the parsed values
echo "Project: $PROJECT"[K
echo "Skip Project Check: $SKIP_PROJECT_CHECK"[K
[K
# Use environment variables if flags are not set

export [36m[45mDAY_P[39;49mROJECT=$PROJECT[K
[K
export DAY_AWS_REGION=$region
[K
# Validate region input
validate_region() {
  if [[ "$1" == "us-west-2" ]]; then[21;5Hecho "Region '$1' confirmed as valid."
    return 0[K[23;3Helse
    echo "Warning: Region '$1' is not 'us-west-2'."
    echo "It is recommended to use 'us-west-2' for this operation."[K[26;29H[3m25% L78 [0m[39;49m[23m[14;13H[34h[?25h[34l[27;1H[K[14;13H[27;1H[?25lMark saved where search started[14;13H[34h[?25h[34l[27;1H[?25l(No changes need to be saved)[K[14;13H[34h[?25h[34l[1;33H[?25l[3mHelp        [0m[39;49m[23m[14;8HDAY_P[26;61H[3m) --------[0m[39;49m[23m[14;13H[34h[?25h[34l[27;1H[K[26;35H[?25l[3m9[0m[39;49m[23m[15;1H[34h[?25h[34l[26;34H[?25l[3m80[0m[39;49m[23m[16;13H[34h[?25h[34l[26;35H[?25l[3m1[0m[39;49m[23m[17;1H[34h[?25h[34l[27d[?25l[36mI-search: [39;49m[17;1H[34h[?25h[34l[1;33H[?25l[3mIsearch Help[0m[39;49m[23m[26;61H[3m Isearch) [0m[39;49m[23m[17;1H[34h[?25h[34l[27;11H[?25la[C[34h[?25h[34l[17;1H[18;4H[?25l[36m[45ma[39;49m[26;35H[3m2[0m[39;49m[23m[18;5H[34h[?25h[34l[27;12H[?25lc[C[34h[?25h[34l[18;5H[2d[?25lecho "Region '$1' confirmed as valid."[3;5Hreturn 0[4;3Helse[K[5;11HWarning: Region '$1' is not 'us-west-2'."[K
    echo "It is recommended to use 'us-west-2' for this operation."[7;5Hreturn 0
  fi[K
}[K
[K
# Call the region validation function
validate_region "$region" || return 3[K

	AWS_[36m[45mAC[39;49mCOUNT_ID=$(aws sts get-caller-identity --query "Account" --output text)

# Skip project check if the flag is provided
if [[ "$SKIP_PROJECT_CHECK" == false ]]; then
  # Ensure the project is valid for the current user
  USER_PROJECTS=$(grep "^${USER}=" /opt/slurm/etc/projects_list.conf | cut -d'=' -f2 | tr -d ' ')
	! ",${USER_PROJECT[21@S}," =~ ",${PROJECT},[21;11HError: Project '$PROJECT' is not valid for user '$USER'."[22;5Hecho "Valid projects for $USER: $USER_PROJECTS"[23;3H  return 3[24;3Hfi[K
else[K[26;29H[3m30% L97[0m[39;49m[23m[14;14H[34h[?25h[34l[27;13H[?25lc[C[34h[?25h[34l[14;14H[?25l[36m[45mC[39;49m[34h[?25h[34l[14;62H[?25l[46mAcc[39;49m		[34h[?25h[34l[27;14H[?25lo	[34h[?25h[34l[14;15H[?25l[36m[45mO[39;49m[34h[?25h[34l							[?25l[46mo[39;49m		[34h[?25h[34l[27;1H[?25lMark saved where search started[14;16H[34h[?25h[34l[1;33H[?25l[3mHelp        [0m[39;49m[23m[14;12HACCOUNT_ID=$(aws sts get-caller-identity --query "Acco[26;61H[3m) --------[0m[39;49m[23m[14;16H[34h[?25h[34l[27;1H[K[26;35H[?25l[3m8[0m[39;49m[23m[15;1H[34h[?25h[34l[26;35H[?25l[3m9[0m[39;49m[23m[16;1H[34h[?25h[34l[26;34H[?25l[3m100[0m[39;49m[23m[17;1H[34h[?25h[34l[26;36H[?25l[3m1[0m[39;49m[23m[18;1H[34h[?25h[34l[26;36H[?25l[3m2[0m[39;49m[23m[19;1H[34h[?25h[34l[26;36H[?25l[3m3[0m[39;49m[23m[20;1H[34h[?25h[34l[26;36H[?25l[3m4[0m[39;49m[23m[21;1H[34h[?25h[34l[26;36H[?25l[3m5[0m[39;49m[23m[22;1H[34h[?25h[34l[26;36H[?25l[3m6[0m[39;49m[23m[23;1H[34h[?25h[34l[26;36H[?25l[3m7[0m[39;49m[23mMM[34h[?25h[34l[26;36H[?25l[3m8[0m[39;49m[23mM[34h[?25h[34l[1;25r[2;1H[12M[1;27r[14;1H[?25l  echo "Skipping project validation as --skip-project-check was passed."[K
fi[K
[K
# Function to create a new AWS budget[K
create_budget() {[K
[K
  echo "To create budget, please exit and run the following command:"[K
  echo "bash bin/create_budget.sh --amount <AMOUNT> --project <budgetAWSname> --region <REGION>"[K
  echo "...you will then need to add the new project name to /opt/slurm/etc/projects_list.conf"[K
  return 3[K
}[K
[K[26;30H[3m3% L109[0m[39;49m[23m[14;1H[34h[?25h[34l[26;35H[?25l[3m10[0m[39;49m[23m[15;1H[34h[?25h[34l[26;36H[?25l[3m1[0m[39;49m[23m[16;1H[34h[?25h[34l[26;36H[?25l[3m2[0m[39;49m[23m[17;1H[34h[?25h[34l[26;36H[?25l[3m3[0m[39;49m[23m[18;1H[34h[?25h[34l[26;36H[?25l[3m4[0m[39;49m[23m[19;1H[34h[?25h[34l[26;36H[?25l[3m5[0m[39;49m[23m[20;1H[34h[?25h[34l[26;36H[?25l[3m6[0m[39;49m[23m[21;1H[34h[?25h[34l[26;36H[?25l[3m7[0m[39;49m[23m[22;1H[34h[?25h[34l[26;36H[?25l[3m8[0m[39;49m[23m[23;1H[34h[?25h[34l[26;36H[?25l[3m9[0m[39;49m[23mMM[34h[?25h[34l[26;35H[?25l[3m20[0m[39;49m[23mM[34h[?25h[34l[2d[?25l  echo "Skipping project validation as --skip-project-check was passed."[K
fi
[K
# Function to create a new AWS budget[K
create_budget() {[K
[K[8;3Hecho "To create budget, please exit and run the following command:"[9;3Hecho "bash bin/create_budget.sh --amount <AMOUNT> --project <budgetAWSname> --region <REGION>"[10;3Hecho "...you will then need to add the new project name to /opt/slurm/etc/projects_list.conf"[11;3H[2P
}[K
[K
# Query AWS budgets[K
BUDGETS=$(aws budgets describe-budgets --account-id $AWS_ACCOUNT_ID --region "$region" 2>/dev/null)

if [[ -z "$BUDGETS" && "$SKIP_PROJECT_CHECK" != "true" ]]; then
  echo "Error: Unable to retrieve any budgets from AWS. Please check your AWS permissions or configuration."[19;3Hreturn 3
fi[K
[K
# Check if the specified project budget exists[K
MATCHING_BUDGET=$(echo "$BUDGETS" | jq -r ".Budgets[] | select(.BudgetName==\"$PROJECT\")")
[K
if [[ -z "$MATCHING_BUDGET" && "$SKIP_PROJECT_CHECK" != "true" ]]; then[26;29H[3m40% L121[0m[39;49m[23m[14;1H[34h[?25h[34l[26;36H[?25l[3m2[0m[39;49m[23m[15;1H[34h[?25h[34l[26;36H[?25l[3m3[0m[39;49m[23m[16;1H[34h[?25h[34l[26;36H[?25l[3m4[0m[39;49m[23m[17;1H[34h[?25h[34l[26;36H[?25l[3m5[0m[39;49m[23m[18;1H[34h[?25h[34l[26;36H[?25l[3m6[0m[39;49m[23m[19;1H[34h[?25h[34l[26;36H[?25l[3m7[0m[39;49m[23m[20;1H[34h[?25h[34l[26;36H[?25l[3m8[0m[39;49m[23m[21;1H[34h[?25h[34l[26;36H[?25l[3m9[0m[39;49m[23m[22;1H[34h[?25h[34l[26;35H[?25l[3m30[0m[39;49m[23m[23;1H[34h[?25h[34l[26;36H[?25l[3m1[0m[39;49m[23mMM[34h[?25h[34l[26;36H[?25l[3m2[0m[39;49m[23mM[34h[?25h[34l[1;25r[2;1H[12M[1;27r[14;1H[?25l  echo "No matching AWS budget found for project '$PROJECT' in region '$region'."[K
  echo "Available AWS Budgets (which may be specified with the --project flag):"[K
  echo "$BUDGETS" | jq -r '.Budgets[].BudgetName'[K
[K
  read -p "Would you like to create a new budget? (y/n): " RESPONSE[K
  if [[ "$RESPONSE" =~ ^[Yy]$ ]]; then[K
    create_budget || return 3[K
  else[K
    echo "Exiting without creating a budget."[K
    return 3[K
  fi[K
fi[K[26;30H[3m5% L133[0m[39;49m[23m[14;1H[34h[?25h[34l[26;36H[?25l[3m4[0m[39;49m[23m[15;1H[34h[?25h[34l[26;36H[?25l[3m5[0m[39;49m[23m[16;1H[34h[?25h[34l[26;36H[?25l[3m6[0m[39;49m[23m[17;1H[34h[?25h[34l[26;36H[?25l[3m7[0m[39;49m[23m[18;1H[34h[?25h[34l[26;36H[?25l[3m8[0m[39;49m[23m[19;1H[34h[?25h[34l[26;36H[?25l[3m9[0m[39;49m[23m[20;1H[34h[?25h[34l[26;35H[?25l[3m40[0m[39;49m[23m[21;1H[34h[?25h[34l[26;36H[?25l[3m1[0m[39;49m[23m[22;1H[34h[?25h[34l[26;36H[?25l[3m2[0m[39;49m[23m[23;1H[34h[?25h[34l[26;36H[?25l[3m3[0m[39;49m[23mMM[34h[?25h[34l[26;36H[?25l[3m4[0m[39;49m[23mM[34h[?25h[34l[2d[?25l  echo "No matching AWS budget found for project '$PROJECT' in region '$region'."
  echo "Available AWS Budgets (which may be specified with the --project flag):"[K[4;3Hecho "$BUDGETS" | jq -r '.Budgets[].BudgetName'
[K[6;3Hread -p "Would you like to create a new budget? (y/n): " RESPONSE[K[7;3Hif [[ "$RESPONSE" =~ ^[Yy]$ ]]; then
    create_budget || return 3[9;3Helse
    echo "Exiting without creating a budget."[K
    return 3[K[12;3Hfi
fi[K
[K
# Extract budget details for the matching budget[K
TOTAL=$(echo "$MATCHING_BUDGET" | jq -r ".BudgetLimit.Amount")
USED=$(echo "$MATCHING_BUDGET" | jq -r ".CalculatedSpend.ActualSpend.Amount")
[K
[2P	(-z "$TOTAL" || -z "$[41@USED") && "$SKIP_PROJECT_CHECK" != "true"[20;3Hecho "Error: Unable to calculate budget details for project '$PROJECT'.  $SKIP_PROJECT_CHECK"[21;3Hreturn 3
fi[K
[K
# Calculate usage percentage
PERCENT_USED=$(awk "BEGIN {print ($USED / $TOTAL) * 100}")[26;29H[3m51% L145[0m[39;49m[23m[14;1H[34h[?25h[34l[26;36H[?25l[3m6[0m[39;49m[23m[15;1H[34h[?25h[34l[26;36H[?25l[3m7[0m[39;49m[23m[16;1H[34h[?25h[34l[?25l[1@e[26;6H[3m**[0m[39;49m[23m[16;2H[34h[?25h[34l[?25l[1@x[34h[?25h[34l[?25l[1@p[34h[?25h[34l[?25l[1@o[34h[?25h[34l[?25l[1@r[34h[?25h[34l[?25l[1@t[34h[?25h[34l[?25l[1@ [34h[?25h[34l	[C[C[C[C[?25l[1@_[34h[?25h[34l[?25l[1@B[34h[?25h[34l[?25l[1@U[34h[?25h[34l[?25l[1@D[34h[?25h[34l[?25l[1@G[34h[?25h[34l[?25l[1@E[34h[?25h[34l[?25l[1@T[34h[?25h[34l[26;36H[?25l[3m8[0m[39;49m[23m[17;20H[34h[?25h[34l[26;36H[?25l[3m9[0m[39;49m[23m[18;1H[34h[?25h[34l[26;35H[?25l[3m50[0m[39;49m[23m[19;1H[34h[?25h[34l[26;35H[?25l[3m49[0m[39;49m[23m[18;1H[34h[?25h[34l[26;36H[?25l[3m8[0m[39;49m[23m[17;1H[34h[?25h[34l[?25l[1@e[34h[?25h[34l[?25l[1@x[34h[?25h[34l[?25l[1@p[34h[?25h[34l[?25l[1@o[34h[?25h[34l[?25l[1@r[34h[?25h[34l[?25l[1@t[34h[?25h[34l[?25l[1@ [34h[?25h[34l	[C[C[C[?25l[1@_[34h[?25h[34l[?25l[1@B[34h[?25h[34l[?25l[1@U[34h[?25h[34l[?25l[1@D[34h[?25h[34l[?25l[1@G[34h[?25h[34l[?25l[1@E[34h[?25h[34l[?25l[1@T[34h[?25h[34l[26;36H[?25l[3m9[0m[39;49m[23m[18;1H[34h[?25h[34l[26;35H[?25l[3m50[0m[39;49m[23m[19;19H[34h[?25h[34l[C		[C[27;1H[?25lAuto-saving...[19;18H[34h[?25h[34l[27;15H[?25ldone[19;18H[34h[?25h[34l[27;1H[K[19;18H[?25l[1@_[34h[?25h[34l[?25l[1@B[34h[?25h[34l[?25l[1@U[34h[?25h[34l[?25l[1@D[34h[?25h[34l[?25l[1@G[34h[?25h[34l[?25l[1@E[34h[?25h[34l[?25l[1@T[34h[?25h[34l[26;35H[?25l[3m49[0m[39;49m[23m[18;1H[34h[?25h[34l[26;35H[?25l[3m50[0m[39;49m[23m[19;1H[34h[?25h[34l[26;36H[?25l[3m1[0m[39;49m[23m[20;1H[34h[?25h[34l[26;36H[?25l[3m0[0m[39;49m[23m[19;85H[34h[?25h[34l[?25l[1@_[34h[?25h[34l[?25l[1@B[34h[?25h[34l[?25l[1@U[34h[?25h[34l[?25l[1@D[34h[?25h[34l[?25l[1@G[34h[?25h[34l[?25l[1@E[34h[?25h[34l[?25l[1@T[26;30H[3m0[0m[39;49m[23m[19;46H[34h[?25h[34l[27;1H[?25lSaving file /fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily/dyinit...[19;46H[34h[?25h[34l[27;1H[?25lWrote /fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily/dyinit[K[19;46H[34h[?25h[34l[26;6H[?25l[3m--[0m[39;49m[23m[19;46H[34h[?25h[34l[27;1H[?25l[36mI-search: [39;49m[K[19;46H[34h[?25h[34l[1;33H[?25l[3mIsearch Help[0m[39;49m[23m[26;61H[3m Isearch) [0m[39;49m[23m[19;46H[34h[?25h[34l[27;11H[?25lt[C[34h[?25h[34l[19;46H			[?25l[36m[45mT[39;49m[34h[?25h[34l[27;12H[?25lo[C[34h[?25h[34l[19;66H[?25lT[20;23H[36m[45mto[39;49m[26;36H[3m1[0m[39;49m[23m[20;25H[34h[?25h[34l[27;13H[?25lt[C[34h[?25h[34l[20;25H[?25lto[25;44H[36m[45mTOT[39;49m[26;36H[3m6[0m[39;49m[23m[25;47H[34h[?25h[34l[27;14H[?25la	[34h[?25h[34l[25;47H[?25l[36m[45mA[39;49m[34h[?25h[34l[27;15H[?25ll	[34h[?25h[34l[25;48H[?25l[36m[45mL[39;49m[34h[?25h[34l[16;8H[?25l[46mTOTAL[39;49m


[46mTOTAL[39;49m[25;49H[34h[?25h[34l

[?25lMark saved where search startedMM			[34h[?25h[34l[1;33H[?25l[3mHelp        [0m[39;49m[23m[16;8HTOTAL


TOTAL[25;44HTOTAL[26;61H[3m) --------[0m[39;49m[23m[25;49H[34h[?25h[34l

[K[25;49H[?25l[1@_[26;6H[3m**[0m[39;49m[23m[25;50H[34h[?25h[34l[?25l[1@B[34h[?25h[34l[?25l[1@U[34h[?25h[34l[?25l[1@D[34h[?25h[34l[?25l[1@G[34h[?25h[34l[?25l[1@E[34h[?25h[34l[?25l[1@T[34h[?25h[34l[26;36H[?25l[3m5[0m[39;49m[23m[24;29H[34h[?25h[34l[26;36H[?25l[3m6[0m[39;49m[23m[25;28H[34h[?25h[34l[C[C[C		[C[C[C[C[C[C		[?25l[1@_[34h[?25h[34l[?25l[1@B[34h[?25h[34l[?25l[1@U[34h[?25h[34l[?25l[1@D[34h[?25h[34l[?25l[1@G[34h[?25h[34l[?25l[1@E[34h[?25h[34l[?25l[1@T[34h[?25h[34l

[?25lSaving file /fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily/dyinit...[25;47H[34h[?25h[34l

[?25lWrote /fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily/dyinit[K[25;47H[34h[?25h[34l[26;6H[?25l[3m--[0m[39;49m[23m[25;47H[34h[?25h[34l

[?25l[36mI-search: [39;49m[K[25;47H[34h[?25h[34l[1;33H[?25l[3mIsearch Help[0m[39;49m[23m[26;61H[3m Isearch) [0m[39;49m[23m[25;47H[34h[?25h[34l[27;11H[?25l$[C[34h[?25h[34l[25;47H	[C[?25l[36m[45m$[39;49m[34h[?25h[34l[2d[?25l[46m$[39;49mPROJECT' in region '[46m$[39;49m

	[46m$[39;49m


[46m$[39;49mRESPONSE" =~ ^[Yy][46m$[39;49m[16;21H[46m$[39;49m(echo "[46m$[39;49m[17;20H[46m$[39;49m(echo "[46m$[39;49m[19;12H[46m$[39;49mTOTAL_BUDGET" || -z "[46m$[39;49mUSED_BUDGET") && "[46m$[39;49m
		[46m$[39;49mPROJECT'.  [46m$[39;49m[25;14H[46m$[39;49m(awk "BEGIN {print ([46m$[39;49m[25;51H[34h[?25h[34l[27;12H[?25lt[C[34h[?25h[34l[25;51H[?25l[36m[45mT[39;49m[34h[?25h[34l[27;13H[?25lo[C[34h[?25h[34l[25;52H[?25l[36m[45mO[39;49m[34h[?25h[34l[27;14H[?25lt	[34h[?25h[34l[25;53H[?25l[36m[45mT[39;49m[34h[?25h[34l[2;51H[?25l$PROJECT' in region '$

	$


$RESPONSE" =~ ^[Yy]$[16;21H$(echo "$[17;20H$(echo "$[19;13H[46mTOT[39;49mAL_BUDGET" || -z "$USED_BUDGET") && "$
		$PROJECT'.  $[25;14H$(awk "BEGIN {print ($[25;54H[34h[?25h[34l[27;15H[?25la	[34h[?25h[34l[25;54H[?25l[36m[45mA[39;49m[34h[?25h[34l

		[?25ll[C[34h[?25h[34l[25;55H[?25l[36m[45mL[39;49m[34h[?25h[34l[19;16H[?25l[46mAL[39;49m[25;56H[34h[?25h[34l[H
[?25lif [[ (-z "[46m$TOTAL[39;49m_BUDGET" || -z "$USED_BUDGET") && "$SKIP_PROJECT_CHECK" != "true" ]]; then
	Error: Unable to calculate budget details for project '$PROJECT'.  $SKIP_PROJECT_CHECK"[4;3Hreturn 3[K
fi
[K
# Calculate usage percentage[K
PERCENT_USED=$(awk "BEGIN {print ($USED_BUDGET / [46m$TOTAL[39;49m_BUDGET) * 100}")
[K
# Display budget information[K
echo ""[K
echo "________________________________________________________"
echo "AWS Budget for project '$PROJECT' in region '$region':"
echo "  Total: [36m[45m$TOTAL[39;49m USD"
echo "  Used: $USED USD"[K
[Ccho "  Percent Used: $PERCENT_USED%"[K
[Ccho "________________________________________________________"[K
sleep 1.3
[K
# Check if the budget is exhausted[K
if (( $(echo "$PERCENT_USED >= 100" | bc -l) )); then
  echo "Warning: Budget for project '$PROJECT' is exhausted!"
fi
[K
## export DAY_PROJECT="daylily-dev"[K
[3m8% L162[0m[39;49m[23m[14;22H[34h[?25h[34l[27;1H[?25lMark saved where search started[14;22H[34h[?25h[34l[1;33H[?25l[3mHelp        [0m[39;49m[23m[2;12H$TOTAL[8;50H$TOTAL[14;16H$TOTAL[26;61H[3m) --------[0m[39;49m[23m[14;21H[34h[?25h[34l[27;1H[K[14;22H[?25l_ USD"[26;6H[3m**[0m[39;49m[23m[14;23H[34h[?25h[34l[?25lB USD"[34h[?25h[34l[?25lU USD"			[34h[?25h[34l[?25lD USD"[34h[?25h[34l[?25lG USD"[34h[?25h[34l[?25lE USD"[34h[?25h[34l[?25lT USD"[34h[?25h[34l[26;36H[?25l[3m3[0m[39;49m[23m[15;25H[34h[?25h[34l[?25l_ USD"[34h[?25h[34l[?25lB USD"[34h[?25h[34l[?25lU USD"[34h[?25h[34l[?25lD USD"[34h[?25h[34l[?25lG USD"			[34h[?25h[34l[?25lE USD"[34h[?25h[34l[?25lT USD"[34h[?25h[34l[27;1H[?25l[36mI-search: [39;49m[15;27H[34h[?25h[34l[1;33H[?25l[3mIsearch Help[0m[39;49m[23m[26;61H[3m Isearch) [0m[39;49m[23m[15;27H[34h[?25h[34l[27;11H[?25l$[C[34h[?25h[34l[15;27H
[?25l[36m[45m$[39;49m[26;36H[3m4[0m[39;49m[23m[16;24H[34h[?25h[34l[2;12H[?25l[46m$[39;49mTOTAL_BUDGET" || -z "[46m$[39;49mUSED_BUDGET") && "[46m$[39;49m
		[46m$[39;49mPROJECT'.  [46m$[39;49m[8;14H[46m$[39;49m(awk "BEGIN {print ([46m$[39;49mUSED_BUDGET / [46m$[39;49m[13;31H[46m$[39;49mPROJECT' in region '[46m$[39;49m
		[46m$[39;49m
[46m$[39;49m[21;7H[46m$[39;49m(echo "[46m$[39;49m[22;38H[46m$[39;49m[16;24H[34h[?25h[34l[27;12H[?25lu[C[34h[?25h[34l[16;24H[H
[?25lecho "  Total: [46m$[39;49mTOTAL_BUDGET USD"[K
echo "  Used: [46m$[39;49mUSED_BUDGET USD"[K
echo "  Percent Used: [46m$[39;49mPERCENT_USED%"
echo "________________________________________________________"
sleep 1.3
[K
# Check if the budget is exhausted[K
if (( [46m$[39;49m(echo "[46m$[39;49mPERCENT_USED >= 100" | bc -l) )); then
  echo "Warning: Budget for project '[46m$[39;49mPROJECT' is exhausted!"
fi[K
[K
## export DAY_PROJECT="daylily-dev"[K
[Cxport APPTAINER_HOME=/fsx/resources/environments/containers/[36m[45m$U[39;49mSER/$(hostname)/
[K
# needed by daycli[K
[Cxport DAY_BIOME=AWSPC[K
[K
export DAY_ROOT=$PWD
[K
ulimit -n 16384[K
[K
export SENTIEON_TMPDIR='/fsx/scratch/'
export SENTIEON_LICENSE='/fsx/data/cached_envs/Rare_Cancer_Research_Foundation_eval.lic'
export SENTIEON_INSTALL_DIR='/fsx/data/cached_envs/sentieon-genomics-202308.03'[26;29H[3m63% L17[0m[39;49m[23m[14;64H[34h[?25h[34l[27;13H[?25ls[C[34h[?25h[34l[14;64H[?25l[36m[45mS[39;49m[34h[?25h[34l[27;14H[?25le	[34h[?25h[34l[14;65H[?25l[36m[45mE[39;49m[34h[?25h[34l[27;15H[?25lr	[34h[?25h[34l[14;66H[?25l[36m[45mR[39;49m[34h[?25h[34l[27;16H[?25lt[C[34h[?25h[34l[?25l[36mFailing I-search: [39;49m$user[41mt[39;49m[14;67H[34h[?25h[34l[2;16H[?25l$
$[4;23H$[9;7H$(echo "$[10;38H$[14;67H[34h[?25h[34l[27;1H[?25l[36mI-search: [39;49m$user[K[14;67H[34h[?25h[34l[27;15H[K[14;67H[?25lR[34h[?25h[34l[3;15H[?25l[46m$USE[39;49m[14;66H[34h[?25h[34l[27;15H[?25ld	[34h[?25h[34l[?25l[36mFailing I-search: [39;49m$use[41md[39;49m[14;66H[34h[?25h[34l[3;19H[?25l[46mD[39;49m[14;66H[34h[?25h[34l[27;1H[?25l[36mWrapped I-search: [39;49m$used[14;66H[34h[?25h[34l[H
[?25l  if [[ "$RESPONSE" =~ ^[Yy]$ ]]; then
    create_budget || return 3[K
  else[K
    echo "Exiting without creating a budget."[K
    return 3[7;3Hfi
fi[K
[K
# Extract budget details for the matching budget[K
export TOTAL_BUDGET=$(echo "$MATCHING_BUDGET" | jq -r ".BudgetLimit.Amount")
export USED_BUDGET=$(echo "$MATCHING_BUDGET" | jq -r ".CalculatedSpend.ActualSpend.Amount")
[K
if [[ (-z "$TOTAL_BUDGET" || -z "[36m[45m$USED[39;49m_BUDGET") && "$SKIP_PROJECT_CHECK" != "true" ]]; then[15;3Hecho "Error: Unable to calculate budget details for project '$PROJECT'.  $SKIP_PROJECT_CHECK"
  return 3[K
fi[K

# Calculate usage percentage
PERCENT_USED=$(awk "BEGIN {print ($USED_BUDGET / $TOTAL_BUDGET) * 100}")
[K
# Display budget information
[Ccho ""[K
[Ccho "________________________________________________________"[K
[Ccho "AWS Budget for project '$PROJECT' in region '$region':"[K[26;29H[3m54% L150[0m[39;49m[23m[14;39H[34h[?25h[34l[20;35H[?25l[46m$USED[39;49m[14;39H[34h[?25h[34l[?25l[46m$USED[39;49m[20;35H[36m[45m$USED[39;49m[26;36H[3m6[0m[39;49m[23m[20;40H[34h[?25h[34l[2;3H[?25lecho "Error: Unable to calculate budget details for project '$PROJECT'.  $SKIP_PROJECT_CHECK"[3;3Hreturn 3[K
fi[K
[K
# Calculate usage percentage
PERCENT_USED=$(awk "BEGIN {print ([46m$USED[39;49m_BUDGET / $TOTAL_BUDGET) * 100}")
[K
# Display budget information
echo ""[K
[Ccho "________________________________________________________"[K
[Ccho "AWS Budget for project '$PROJECT' in region '$region':"[K
echo "  Total: $TOTAL_BUDGET USD"
echo "  Used: [36m[45m$USED[39;49m_BUDGET USD"[K
echo "  Percent Used: $PERCENT_USED%"[K
echo "________________________________________________________"
sleep 1.3[19;4Hheck if the budget is exhausted
if (( $(echo "$PERCENT_USED >= 100" | bc -l) )); then[K[21;3Hecho "Warning: Budget for project '$PROJECT' is exhausted!"
fi[K
[K
## export DAY_PROJECT="daylily-dev"[K
[Cxport APPTAINER_HOME=/fsx/resources/environments/containers/$USER/$(hostname)/[26;30H[3m9% L163[0m[39;49m[23m[14;20H[34h[?25h[34l[27;1H[?25l[36mFailing wrapped I-search: [39;49m[41m$used[39;49m[14;20H[34h[?25h[34l[27;1H[?25l[36mWrapped I-search: [39;49m$used[K[14;20H[34h[?25h[34l[2;3H[?25lif [[ "$RESPONSE" =~ ^[Yy]$ ]]; then[K[3;3H  create_budget || return 3
  else
echo "Exiting without creating a budget."
    return 3[K
  fi[K
fi
[K
# Extract budget details for the matching budget
[Cxport TOTAL_BUDGET=$(echo "$MATCHING_BUDGET" | jq -r ".BudgetLimit.Amount")
[Cxport USED_BUDGET=$(echo "$MATCHING_BUDGET" | jq -r ".CalculatedSpend.ActualSpend.Amount")
[K
if [[ (-z "$TOTAL_BUDGET" || -z "[36m[45m$USED[39;49m_BUDGET") && "$SKIP_PROJECT_CHECK" != "true" ]]; then
  echo "Error: Unable to calculate budget details for project '$PROJECT'.  $SKIP_PROJECT_CHECK"
  return 3[K
fi[K

[Calculate usage percentage[K
PERCENT_USED=$(awk "BEGIN {print ([46m$USED[39;49m_BUDGET / $TOTAL_BUDGET) * 100}")
[K
# Display budget information
echo ""
echo "________________________________________________________"
[Ccho "AWS Budget for project '$PROJECT' in region '$region':"[K[26;30H[3m4% L150[0m[39;49m[23m[14;39H[34h[?25h[34l[27;1H[?25lAuto-saving...[K[14;39H[34h[?25h[34l[27;1H[?25l[36mWrapped I-search: [39;49m$used[14;39H[34h[?25h[34l[27;1H[?25lMark saved where search started[14;39H[34h[?25h[34l[1;33H[?25l[3mHelp        [0m[39;49m[23m[14;34H$USED[20;35H$USED[26;36H[3m1  [0m[39;49m[23m[3mGit:main[0m[39;49m[23m[3m  (Fundamental) --------[0m[39;49m[23m[15;39H[34h[?25h[34l[27;1H[K[15;40H[26;36H[?25l[3m2[0m[39;49m[23m[16;11H[34h[?25h[34l[26;36H[?25l[3m3[0m[39;49m[23m[17;1H[34h[?25h[34l[26;36H[?25l[3m4[0m[39;49m[23m[18;1H[34h[?25h[34l[26;36H[?25l[3m5[0m[39;49m[23m[19;1H[34h[?25h[34l[26;36H[?25l[3m6[0m[39;49m[23m[20;1H[34h[?25h[34l[?25l[1@e[34h[?25h[34l[?25l[1@x[34h[?25h[34l[?25l[1@p[34h[?25h[34l[?25l[1@o[34h[?25h[34l[?25l[1@r[34h[?25h[34l[?25l[1@t[34h[?25h[34l[?25l[1@ [34h[?25h[34l[27;1H[?25lAuto-saving...[20;8H[34h[?25h[34l[27;15H[?25ldone[20;8H[34h[?25h[34l[1;25r[22;1H[1L[1;27r[27;1H[K[20;7H[K
[?25lPERCENT_USED=$(awk "BEGIN {print ($USED_BUDGET / $TOTAL_BUDGET) * 100}")
[K[26;36H[3m7[0m[39;49m[23m[21;1H[34h[?25h[34l[1;25r[21;1H[1L[1;27r[21;1H[K[26;36H[?25l[3m8[0m[39;49m[23m[22;1H[34h[?25h[34l[1;25r[22;1H[1L[1;27r[22;1H[K[26;36H[?25l[3m9[0m[39;49m[23m[23;1H[34h[?25h[34l[1;25r[22;1H[1M[1;27r[25;1H[?25lecho ""[K[26;36H[3m8[0m[39;49m[23m[22;1H[34h[?25h[34l[1;25r[21;1H[1M[1;27r[25;1H[?25lecho "________________________________________________________"[K[26;36H[3m7[0m[39;49m[23m[21;1H[34h[?25h[34l[1;25r[22;1H[1M[1;27r[20;7H[?25lPERCENT_USED=$(awk "BEGIN {print ($USED_BUDGET / $TOTAL_BUDGET) * 100}")
[K



echo "AWS Budget for project '$PROJECT' in region '$region':"[K[26;36H[3m6[0m[39;49m[23m[20;7H[34h[?25h[34l[?25l[1@ [34h[?25h[34l[27;1H[?25lSaving file /fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily/dyinit...[20;8H[34h[?25h[34l[27;1H[?25lWrote /fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily/dyinit[K[20;8H[34h[?25h[34l[26;6H[?25l[3m--[0m[39;49m[23m[20d[34h[?25h[34l[27;1H[K[20;8H[27;1H[K[?1004l[?2004l[>4m[?1l>[34h[?25h[?1049l[39;49m[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ git commit -a -m "X"
[?2004lg[main 2aa9531] X
 Committer: AWS ParallelCluster user <ubuntu@ip-10-0-0-196.us-west-2.compute.internal>
Your name and email address were configured automatically based
on your username and hostname. Please check that they are accurate.
You can suppress this message by setting them explicitly. Run the
following command and follow the instructions in your editor to edit
your configuration file:

    git config --global --edit

After doing this, you may fix the identity used for this commit with:

    git commit --amend --reset-author

 2 files changed, 484 insertions(+), 8 deletions(-)
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ giut [K[K[Kt push
[?2004lEnumerating objects: 11, done.
Counting objects:   9% (1/11)Counting objects:  18% (2/11)Counting objects:  27% (3/11)Counting objects:  36% (4/11)Counting objects:  45% (5/11)Counting objects:  54% (6/11)Counting objects:  63% (7/11)Counting objects:  72% (8/11)Counting objects:  81% (9/11)Counting objects:  90% (10/11)Counting objects: 100% (11/11)Counting objects: 100% (11/11), done.
Delta compression using up to 8 threads
Compressing objects:  16% (1/6)Compressing objects:  33% (2/6)Compressing objects:  50% (3/6)Compressing objects:  66% (4/6)Compressing objects:  83% (5/6)Compressing objects: 100% (6/6)Compressing objects: 100% (6/6), done.
Writing objects:  16% (1/6)Writing objects:  33% (2/6)Writing objects:  50% (3/6)Writing objects:  66% (4/6)Writing objects:  83% (5/6)Writing objects: 100% (6/6)Writing objects: 100% (6/6), 13.09 KiB | 4.36 MiB/s, done.
Total 6 (delta 4), reused 0 (delta 0), pack-reused 0
remote: Resolving deltas:   0% (0/4)[Kremote: Resolving deltas:  25% (1/4)[Kremote: Resolving deltas:  50% (2/4)[Kremote: Resolving deltas:  75% (3/4)[Kremote: Resolving deltas: 100% (4/4)[Kremote: Resolving deltas: 100% (4/4), completed with 4 local objects.[K
remote: Bypassed rule violations for refs/heads/main:[K
remote: 
remote: - Changes must be made through a pull request.[K
remote: 
remote: - Cannot update this protected ref.[K
remote: 
To github.com:Daylily-Informatics/daylily.git
   6f17547..2aa9531  main -> main
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ git pull
[?2004lAlready up to date.
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ du -hs results/[K
[?2004l239G	results
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ du -hs results[1@s[1@i[1@z[1@e[1@=[1@$[1@([C[C[C[C[C[C[C[C[C[C[C[C[C[C)
[?2004lecho [?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ echo $size
[?2004l239G results
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ emacs dyinit 
[?2004l[?1049h[34l[?1h=[H[J[27d[K[?1l>[34h[?25h[?1049l[39;49m[?1049h[34l[?1h=[H[J[26d[?25l[3m-UUU:----F1  [0m[39;49m[23m[3m[1m*scratch*   [0m[39;49m[23m[3m   All L1     (Fundamental) -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------[0m[39;49m[23m
M[2d[34h[?25h[34l[27d[?25lLoading /etc/emacs/site-start.d/00debian.el (source)...done[K[H
[34h[?25h[34l[27;33H[?25l50autoconf.el (source)...[K[H
[34h[?25h[34l[27;58H[?25ldone[H
[34h[?25h[34l[27;35H[?25ldictionaries-common.el (source)...[H
[34h[?25h[34l[27;9H[?25ldebian-ispell...[K[H
[34h[?25h[34l[27;9H[?25l/var/cache/dictionaries-common/emacsen-ispell-default.el (source)...[H
[34h[?25h[34l[27;77H[?25ldone[H
[34h[?25h[34l[27;9H[?25ldebian-ispell...done[K[H
[34h[?25h[34l[27;9H[?25l/var/cache/dictionaries-common/emacsen-ispell-dicts.el (source)...[H
[34h[?25h[34l[27;75H[?25ldone[H
[34h[?25h[34l[27;10H[?25letc/emacs/site-start.d/50dictionaries-common[6P[H
[34h[?25h[34l[27;35H[?25ltcsh.el (source)...[K[H
[34h[?25h[34l[27;54H[?25ldone[H
[34h[?25h[34l[>4;1m[?2004h[?1004h[27d[?25lFor information about GNU Emacs and the GNU system, type C-h C-a.[K[H
[34h[?25h[34l[27d[?25lFor information about GNU Emacs and the GNU system, type C-h C-a.[K[H[3mFile Edit Options Buffers Tools Help                                                                                                                                                                                                                                              [0m[39;49m[23m
M
# Ensure the script is sourced[K
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then[K
    echo "Error: This script must be sourced, not executed directly. Use 'source $0' to run."[K
    return 3[K
fi[K
[K
[K
. "/etc/parallelcluster/cfnconfig"[K
region=$cfn_region[K
[K
# Ensure the script is sourced, not executed directly[K
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then[K
  echo "Error: This script must be sourced, not executed directly."[K
  echo "Usage: source $0 --project <project> [--skip-project-check]"[K
  return 3[K
fi[K
[K
[K
# Check if $cfn_region is unset or null[K
if [[ -z "$region" ]]; then[K
  echo "Error: cfn_region is not set or is empty."[K
  return 3[K
fi[K
[K
[3m-UU-:----F1  [0m[39;49m[23m[3m[1mdyinit      [0m[39;49m[23m[3m   Top L1    [0m[39;49m[23m[3mGit-main[0m[39;49m[23m[3m  (Fundamental) --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------[0m[39;49m[23m
M[2d[34h[?25h[34l[27d[K[26;34H[?25l[3m2[0m[39;49m[23m[H

[34h[?25h[34l[26;34H[?25l[3m3[0m[39;49m[23m[4;1H[34h[?25h[34l[26;34H[?25l[3m4[0m[39;49m[23m[5;1H[34h[?25h[34l[26;34H[?25l[3m5[0m[39;49m[23m[6;1H[34h[?25h[34l[26;34H[?25l[3m6[0m[39;49m[23m[7;1H[34h[?25h[34l[26;34H[?25l[3m7[0m[39;49m[23m[8;1H[34h[?25h[34l[1;25r[9;1H[1L[1;27r[9;1H[K[26;6H[?25l[3m**--F1  [0m[39;49m[23m[3m[1mdyinit      [0m[39;49m[23m[3m   Top L8[0m[39;49m[23m[9;1H[34h[?25h[34l[1;25r[10;1H[1L[1;27r[10;1H[K[26;34H[?25l[3m9[0m[39;49m[23m[10;1H[34h[?25h[34l[26;34H[?25l[3m8[0m[39;49m[23m[9;1H[34h[?25h[34l[?25le[34h[?25h[34l[?25lx[34h[?25h[34l[?25lp[34h[?25h[34l[?25lo[34h[?25h[34l[?25lr[34h[?25h[34l[?25lt[34h[?25h[34l	[?25lD[34h[?25h[34l[?25lA[34h[?25h[34l[?25lY[34h[?25h[34l[?25l_[34h[?25h[34l[?25lC[34h[?25h[34l[?25lO[34h[?25h[34l[?25lN[34h[?25h[34l[?25lT[34h[?25h[34l[?25lA[34h[?25h[34l[?25lC[34h[?25h[34l[?25lT[34h[?25h[34l[?25l_[34h[?25h[34l[?25lE[34h[?25h[34l[?25lM[34h[?25h[34l[?25lA[34h[?25h[34l[?25lI[34h[?25h[34l[?25lL[34h[?25h[34l[?25l=[34h[?25h[34l[?25lj[34h[?25h[34l[?25lo[34h[?25h[34l[?25lh[34h[?25h[34l[?25ln[34h[?25h[34l[?25l@[34h[?25h[34l[?25ld[34h[?25h[34l[?25la[34h[?25h[34l[?25ly[34h[?25h[34l[?25ll[34h[?25h[34l[?25li[34h[?25h[34l[?25ll[34h[?25h[34l[?25ly[34h[?25h[34l[?25li[34h[?25h[34l[?25ln[34h[?25h[34l[?25lf[34h[?25h[34l[?25lo[34h[?25h[34l[?25lr[34h[?25h[34l[?25lm[34h[?25h[34l[?25la[34h[?25h[34l[?25lt[34h[?25h[34l[?25li[34h[?25h[34l[?25lc[34h[?25h[34l[?25ls[34h[?25h[34l[?25l.[34h[?25h[34l[?25lc[34h[?25h[34l[?25lo[34h[?25h[34l[?25lm[34h[?25h[34l[27;1H[?25lSaving file /fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily/dyinit...[9;53H[34h[?25h[34l[27;1H[?25lWrote /fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily/dyinit[K[9;53H[34h[?25h[34l[26;6H[?25l[3m----F1  [0m[39;49m[23m[3m[1mdyinit      [0m[39;49m[23m[3m   Top L8    [0m[39;49m[23m[3mGit:[0m[39;49m[23m[9;53H[34h[?25h[34l[27;1H[K[9;53H[27;1H[?25l(No changes need to be saved)[9;53H[34h[?25h[34l[27;1H[K[9;53H[27;1H[K[?1004l[?2004l[>4m[?1l>[34h[?25h[?1049l[39;49m[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ ls /fsx/scratch/
[?2004l[0m[01;34mclair3_tmp_20241120090003[0m  [01;34mclair3_tmp_20241120103222[0m  [01;34mclair3_tmp_20241120104201[0m       i192-dy-all-12_spot_price.log  i192-dy-all-2_spot_price.log  i192-dy-all-6_spot_price.log  [01;34mocto_tmp_20241120103747[0m  [01;34mocto_tmp_20241120104157[0m  [01;34mocto_tmp_20241120104208[0m
[01;34mclair3_tmp_20241120102912[0m  [01;34mclair3_tmp_20241120103455[0m  [01;34mdeepvariant_tmp_20241120104201[0m  i192-dy-all-13_spot_price.log  i192-dy-all-3_spot_price.log  i192-dy-all-7_spot_price.log  [01;34mocto_tmp_20241120103757[0m  [01;34mocto_tmp_20241120104200[0m
[01;34mclair3_tmp_20241120102951[0m  [01;34mclair3_tmp_20241120104156[0m  i192-dy-all-10_spot_price.log   i192-dy-all-14_spot_price.log  i192-dy-all-4_spot_price.log  i192-dy-all-8_spot_price.log  [01;34mocto_tmp_20241120103804[0m  [01;34mocto_tmp_20241120104204[0m
[01;34mclair3_tmp_20241120103133[0m  [01;34mclair3_tmp_20241120104200[0m  i192-dy-all-11_spot_price.log   i192-dy-all-1_spot_price.log   i192-dy-all-5_spot_price.log  i192-dy-all-9_spot_price.log  [01;34mocto_tmp_20241120104155[0m  [01;34mocto_tmp_20241120104205[0m
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ more /s[Kfsx/scratch/*price.log
[?2004l::::::::::::::
/fsx/scratch/i192-dy-all-10_spot_price.log
::::::::::::::
2024-11-20 11:00:21 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.464100 USD/hour
2024-11-20 12:20:08 - Region: us-west-2, AZ: us-west-2d, Instance type: c7i.48xlarge, Spot price: 1.929000 USD/hour
2024-11-20 20:57:29 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.449500 USD/hour
[3m--More--(Next file: /fsx/scratch/i192-dy-all-11_spot_price.log)[23m[K[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ more /fsx/scratch/*price.log[C[C[C[C[1P /fsx/scratch/*price.log[1P /fsx/scratch/*price.log[1P /fsx/scratch/*price.log[1P /fsx/scratch/*price.logc /fsx/scratch/*price.loga /fsx/scratch/*price.logt /fsx/scratch/*price.log
[?2004l2024-11-20 11:00:21 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.464100 USD/hour
2024-11-20 12:20:08 - Region: us-west-2, AZ: us-west-2d, Instance type: c7i.48xlarge, Spot price: 1.929000 USD/hour
2024-11-20 20:57:29 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.449500 USD/hour
2024-11-20 12:19:57 - Region: us-west-2, AZ: us-west-2d, Instance type: c7i.48xlarge, Spot price: 1.929000 USD/hour
2024-11-20 20:57:29 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.449500 USD/hour
2024-11-20 08:53:21 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.464100 USD/hour
2024-11-20 12:19:57 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.461900 USD/hour
2024-11-20 21:01:04 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.48xlarge, Spot price: 3.742700 USD/hour
2024-11-20 08:53:29 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.464100 USD/hour
2024-11-20 12:20:05 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.461900 USD/hour
2024-11-20 21:05:09 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.48xlarge, Spot price: 3.742700 USD/hour
2024-11-20 08:53:56 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.464100 USD/hour
2024-11-20 10:03:39 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.464100 USD/hour
2024-11-20 11:47:52 - Region: us-west-2, AZ: us-west-2d, Instance type: r7i.metal-48xl, Spot price: 1.270100 USD/hour
2024-11-20 20:25:22 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.449500 USD/hour
2024-11-20 10:03:40 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.464100 USD/hour
2024-11-20 10:44:39 - Region: us-west-2, AZ: us-west-2d, Instance type: c7i.48xlarge, Spot price: 1.930800 USD/hour
2024-11-20 11:47:50 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.464100 USD/hour
2024-11-20 20:25:14 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.449500 USD/hour
2024-11-20 10:03:39 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.464100 USD/hour
2024-11-20 11:47:58 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.464100 USD/hour
2024-11-20 20:25:14 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.449500 USD/hour
2024-11-20 10:03:49 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.464100 USD/hour
2024-11-20 10:21:54 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.464100 USD/hour
2024-11-20 10:44:39 - Region: us-west-2, AZ: us-west-2d, Instance type: c7i.48xlarge, Spot price: 1.930800 USD/hour
2024-11-20 11:21:58 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.464100 USD/hour
2024-11-20 11:47:50 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.464100 USD/hour
2024-11-20 20:25:24 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.449500 USD/hour
2024-11-20 10:18:52 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.464100 USD/hour
2024-11-20 11:22:07 - Region: us-west-2, AZ: us-west-2d, Instance type: r7i.metal-48xl, Spot price: 1.270100 USD/hour
2024-11-20 12:20:00 - Region: us-west-2, AZ: us-west-2d, Instance type: r7i.metal-48xl, Spot price: 1.270100 USD/hour
2024-11-20 20:25:14 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.449500 USD/hour
2024-11-20 10:44:39 - Region: us-west-2, AZ: us-west-2d, Instance type: r7i.metal-48xl, Spot price: 1.270100 USD/hour
2024-11-20 12:19:57 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.461900 USD/hour
2024-11-20 20:25:14 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.449500 USD/hour
2024-11-20 10:44:39 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.464100 USD/hour
2024-11-20 11:21:57 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.464100 USD/hour
2024-11-20 20:25:14 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.449500 USD/hour
2024-11-20 11:00:21 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.464100 USD/hour
2024-11-20 12:20:10 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.461900 USD/hour
2024-11-20 20:57:37 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.449500 USD/hour
2024-11-20 11:00:21 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.464100 USD/hour
2024-11-20 11:22:09 - Region: us-west-2, AZ: us-west-2d, Instance type: r7i.metal-48xl, Spot price: 1.270100 USD/hour
2024-11-20 12:19:57 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.461900 USD/hour
2024-11-20 20:57:29 - Region: us-west-2, AZ: us-west-2d, Instance type: m7i.metal-48xl, Spot price: 1.449500 USD/hour
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ cat /fsx/scratch/*price.log[Kgitpu[K[K pull
[?2004lAlready up to date.
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ sogit pull
[?2004lAlready up to date.
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ git pull
[?2004lremote: Enumerating objects: 14, done.[K
remote: Counting objects:   7% (1/14)[Kremote: Counting objects:  14% (2/14)[Kremote: Counting objects:  21% (3/14)[Kremote: Counting objects:  28% (4/14)[Kremote: Counting objects:  35% (5/14)[Kremote: Counting objects:  42% (6/14)[Kremote: Counting objects:  50% (7/14)[Kremote: Counting objects:  57% (8/14)[Kremote: Counting objects:  64% (9/14)[Kremote: Counting objects:  71% (10/14)[Kremote: Counting objects:  78% (11/14)[Kremote: Counting objects:  85% (12/14)[Kremote: Counting objects:  92% (13/14)[Kremote: Counting objects: 100% (14/14)[Kremote: Counting objects: 100% (14/14), done.[K
remote: Compressing objects:  20% (1/5)[Kremote: Compressing objects:  40% (2/5)[Kremote: Compressing objects:  60% (3/5)[Kremote: Compressing objects:  80% (4/5)[Kremote: Compressing objects: 100% (5/5)[Kremote: Compressing objects: 100% (5/5), done.[K
remote: Total 14 (delta 9), reused 14 (delta 9), pack-reused 0 (from 0)[K
Unpacking objects:   7% (1/14)Unpacking objects:  14% (2/14)Unpacking objects:  21% (3/14)Unpacking objects:  28% (4/14)Unpacking objects:  35% (5/14)Unpacking objects:  42% (6/14)Unpacking objects:  50% (7/14)Unpacking objects:  57% (8/14)Unpacking objects:  64% (9/14)Unpacking objects:  71% (10/14)Unpacking objects:  78% (11/14)Unpacking objects:  85% (12/14)Unpacking objects:  92% (13/14)Unpacking objects: 100% (14/14)Unpacking objects: 100% (14/14), 2.75 KiB | 80.00 KiB/s, done.
From github.com:Daylily-Informatics/daylily
   2aa9531..1438cb7  main       -> origin/main
Updating 2aa9531..1438cb7
Fast-forward
 bin/proc_spot_price_logs.sh               | 39 [32m+++++++++++++++++++++++++++++++++++++++[m
 config/external_tools/multiqc_config.yaml | 47 [32m+++++++++++++++++++++++++++++++++++++++++++++++[m
 config/external_tools/multiqc_header.yaml | 10 [32m++++++[m[31m----[m
 workflow/rules/multiqc_final_wgs.smk      | 15 [32m+++++++++++++[m[31m--[m
 4 files changed, 105 insertions(+), 6 deletions(-)
 create mode 100644 bin/proc_spot_price_logs.sh
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ source biu[Kn/
all_zones.csv                              create_budget.sh                           day_deactivate                             daylily-cfg-headnode.bash                  init_daycli                                ssh_into_daylily
calcuate_spotprice_for_cluster_yaml.py     create_daylily_omics_analysis_s3.sh        day_run                                    dy-r                                       other/                                     tabcomp.bash
cat_input_files.sh                         day-script-commands                        daylily-cfg-ephemeral-cluster              gen_bcl2fq_ss.py                           proc_spot_price_logs.sh                    test_sleep.sh
check_current_spot_market_by_zones.py      day_activate                               daylily-cfg-ephemeral-cluster.bash         helpers/                                   qc_fastq_files_for_read_order_matching.sh  util/
check_instance_az_presence.py              day_build                                  daylily-cfg-headnode                       init_cloudstackformation.sh                sleep_test.sh                              
(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ source bin/c[Kproc_spot_price_logs.sh 
[?2004lSOURCE ME


Unique Instance Types:
Instance type: c7i.48xlarge
Instance type: m7i.48xlarge
Instance type: m7i.metal-48xl
Instance type: r7i.metal-48xl

Average Spot Price (USD/hour): 1.58141
Median Spot Price (USD/hour): 1.4641
vCPU Cost (USD/min): 0.000137275
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ source bin/proc_spot_price_logs.sh git pull[Kcat /fsx/scratch/*price.loggit pull[K
[?2004l^[[Aremote: Enumerating objects: 4, done.[K
remote: Counting objects:  25% (1/4)[Kremote: Counting objects:  50% (2/4)[Kremote: Counting objects:  75% (3/4)[Kremote: Counting objects: 100% (4/4)[Kremote: Counting objects: 100% (4/4), done.[K
remote: Total 4 (delta 3), reused 4 (delta 3), pack-reused 0 (from 0)[K
Unpacking objects:  25% (1/4)Unpacking objects:  50% (2/4)Unpacking objects:  75% (3/4)Unpacking objects: 100% (4/4)Unpacking objects: 100% (4/4), 471 bytes | 58.00 KiB/s, done.
From github.com:Daylily-Informatics/daylily
   1438cb7..d129b49  main       -> origin/main
^[[AUpdating 1438cb7..d129b49
Fast-forward
 bin/proc_spot_price_logs.sh | 5 [32m+++++[m
 1 file changed, 5 insertions(+)
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ git pullsource bin/proc_spot_price_logs.sh 
[?2004lSOURCE ME


Unique Instance Types:
Instance type: c7i.48xlarge
Instance type: m7i.48xlarge
Instance type: m7i.metal-48xl
Instance type: r7i.metal-48xl

Average Spot Price (USD/hour): 1.58141
Median Spot Price (USD/hour): 1.4641
vCPU Cost (USD/min): 0.000137275
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ source bin/proc_spot_price_logs.sh git pull[K
[?2004l^[[A^[[Aremote: Enumerating objects: 4, done.[K
remote: Counting objects:  25% (1/4)[Kremote: Counting objects:  50% (2/4)[Kremote: Counting objects:  75% (3/4)[Kremote: Counting objects: 100% (4/4)[Kremote: Counting objects: 100% (4/4), done.[K
remote: Total 4 (delta 3), reused 4 (delta 3), pack-reused 0 (from 0)[K
Unpacking objects:  25% (1/4)Unpacking objects:  50% (2/4)Unpacking objects:  75% (3/4)Unpacking objects: 100% (4/4)Unpacking objects: 100% (4/4), 341 bytes | 42.00 KiB/s, done.
From github.com:Daylily-Informatics/daylily
   d129b49..81491ac  main       -> origin/main
Updating d129b49..81491ac
Fast-forward
 bin/proc_spot_price_logs.sh | 1 [32m+[m
 1 file changed, 1 insertion(+)
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ git pullsource bin/proc_spot_price_logs.sh 
[?2004lSOURCE ME


Unique Instance Types:
Instance type: c7i.48xlarge Instance type: m7i.48xlarge Instance type: m7i.metal-48xl Instance type: r7i.metal-48xl
Instance type: c7i.48xlarge
Instance type: m7i.48xlarge
Instance type: m7i.metal-48xl
Instance type: r7i.metal-48xl

Average Spot Price (USD/hour): 1.58141
Median Spot Price (USD/hour): 1.4641
vCPU Cost (USD/min): 0.000137275
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ echo $INSTANCE_TYPES)_[K[K_LINE 
[?2004lInstance type: c7i.48xlarge Instance type: m7i.48xlarge Instance type: m7i.metal-48xl Instance type: r7i.metal-48xl
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ echo $INSTANCE_TYPES_LINE [9@source bin/proc_spot_price_logs.sh[C[Kcat results/day/hg38/other_reports/ry[Kules_benchmark_data_mqc.tsv 
[?2004lcombined_rule	sample	rule	s	h:m:s	max_rss	max_vms	max_uss	max_pss	io_in	io_out	mean_load	cpu_time	rule_prefix	rule_suffix
alignstats_summary-all.	all.	alignstats_summary	0.2405	0:00:00	3.38	7.71	0.12	0.35	0.0	0.0	0.0	0.0	alignstats_summary	NA
final_multiqc-DAY_all.	DAY_all.	final_multiqc	39.5161	0:00:39	544.01	6596.95	444.81	472.93	180.65	17.46	37.96	15.48	final_multiqc	NA
alignstats_smmary_compile-all.	all.	alignstats_smmary_compile	0.7071	0:00:00	11.13	27.76	4.04	4.44	4.48	0.04	0.0	0.0	alignstats_smmary_compile	NA
raw_fastqc-SEQQC-multiqc_.	SEQQC-multiqc_.	raw_fastqc	48.3452	0:00:48	486.65	7625.06	384.82	412.07	2042.91	47.05	85.03	41.23	raw_fastqc	NA
strobe.oct.24-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.oct.24	883.8152	0:14:43	8466.06	13566.77	8495.96	8497.56	140.0	0.1	473.24	4182.7	strobe	oct.24
sent.lfq2.concordance-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.lfq2.concordance	33.2076	0:00:33	13497.05	716739.74	13518.56	13519.51	770.61	159.3	168.12	56.46	sent	lfq2.concordance
sent.sentd.concat.fofn-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.sentd.concat.fofn	0.2802	0:00:00	3.25	7.72	0.12	0.34	0.0	0.0	0.0	0.0	sent	sentd.concat.fofn
bwa2a.tiddit.sv.vcf-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.tiddit.sv.vcf	913.7027	0:15:13	7827.48	49095.2	7701.54	7724.6	5905.44	7456.35	523.06	1027.58	bwa2a	tiddit.sv.vcf
sent.clair3.bcfstat-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.clair3.bcfstat	0.7995	0:00:00	16.5	2561.94	10.21	11.53	3.08	0.0	0.0	1.9	sent	clair3.bcfstat
strobe.tiddit.sv.vcf.sort-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.tiddit.sv.vcf.sort	1.6874	0:00:01	15.0	39.48	10.55	11.07	25.15	57.05	30.81	0.02	strobe	tiddit.sv.vcf.sort
sent.clair3.rtgvcfstats-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.clair3.rtgvcfstats	2.243	0:00:02	874.93	714301.71	872.2	872.54	60.04	0.09	120.8	2.73	sent	clair3.rtgvcfstats
strobe.deep.12-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.12	756.5578	0:12:36	33821.24	204332.97	23247.14	23525.94	1952.54	3946.96	2398.04	68.8	strobe	deep.12
strobe.clair3.22-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.clair3.22	798.6406	0:13:18	10642.89	19498.31	9372.93	9545.34	679.06	522.53	261.35	636.43	strobe	clair3.22
strobe.clair3.peddy-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.clair3.peddy	7.843	0:00:07	143.98	2553.95	141.18	141.64	0.0	8.17	39.39	0.0	strobe	clair3.peddy
strobe.manta-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.manta	440.3811	0:07:20	6260.05	13592.9	6300.4	6303.33	30059.24	144.36	4386.65	4.5	strobe	manta
strobe.deep.bcfstat-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.bcfstat	1.138	0:00:01	13.5	2561.87	12.49	12.67	7.14	0.0	50.09	2.51	strobe	deep.bcfstat
bwa2a.sentd.24-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.sentd.24	95.6385	0:01:35	12200.41	18036.46	9124.57	9131.48	719.9	28.24	1556.96	1140.98	bwa2a	sentd.24
strobe.norm_cov_eveness-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.norm_cov_eveness	628.1167	0:10:28	1938.0	34651.26	1939.17	1939.71	566.03	6960.41	14.19	7.28	strobe	norm_cov_eveness
strobe.lfq2.rtgvcfstats-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.lfq2.rtgvcfstats	0.9907	0:00:00	70.5	712872.53	55.76	56.14	0.0	0.09	0.0	0.21	strobe	lfq2.rtgvcfstats
strobe.clair3.merge-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.clair3.merge	2.3882	0:00:02	19.5	1108.33	18.53	18.68	1.94	20.76	0.0	0.0	strobe	clair3.merge
bwa2a.deep.24-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.deep.24	95.5541	0:01:35	33993.46	204777.78	23243.73	23521.75	1431.62	2209.85	2138.04	1645.94	bwa2a	deep.24
bwa2a.deep.concordance-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.deep.concordance	234.5685	0:03:54	9661.6	705744.48	9656.34	9664.96	1.11	157.75	23.69	56.62	bwa2a	deep.concordance
strobe.clair3.23-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.clair3.23	689.2243	0:11:29	53375.95	91669.16	49094.01	49270.47	1706.38	1159.96	1459.11	379.36	strobe	clair3.23
strobe.deep.13-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.13	558.6623	0:09:18	34124.88	204336.98	23428.7	23707.25	122.52	3405.74	2311.85	43.44	strobe	deep.13
strobe.oct.concordance-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.oct.concordance	245.0714	0:04:05	10750.8	716741.51	10744.72	10753.37	692.83	144.83	25.71	63.12	strobe	oct.concordance
bwa2a.clair3.peddy-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.clair3.peddy	7.7472	0:00:07	140.91	2554.09	113.39	127.89	2.89	8.13	34.69	0.0	bwa2a	clair3.peddy
sent.deep.24-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.deep.24	140.8003	0:02:20	5109.12	107811.88	2944.42	3161.89	1037.91	2319.83	1048.11	36.87	sent	deep.24
strobe.deep.10-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.10	522.8685	0:08:42	34018.38	205422.16	23411.57	23689.83	149.9	4025.61	3732.34	71.51	strobe	deep.10
sent.alignstats-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.alignstats	602.4617	0:10:02	981.0	2481.5	983.96	984.12	0.48	0.09	138.48	834.54	sent	alignstats
strobe.lfq2.concat.fofn-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.lfq2.concat.fofn	0.3482	0:00:00	15.0	39.48	10.56	11.35	0.0	0.0	0.0	0.03	strobe	lfq2.concat.fofn
sent.sentd.rtgvcfstats-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.sentd.rtgvcfstats	3.1798	0:00:03	1388.31	714236.71	1395.08	1395.49	73.83	0.09	110.59	3.68	sent	sentd.rtgvcfstats
bwa2a.oct.24-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.oct.24	698.497	0:11:38	9698.23	15280.52	9711.45	9719.83	0.0	0.12	437.78	3058.65	bwa2a	oct.24
bwa2a.clair3.merge-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.clair3.merge	2.4978	0:00:02	21.0	1108.34	18.07	18.2	1.48	20.75	0.0	0.0	bwa2a	clair3.merge
strobe.deep.11-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.11	513.7636	0:08:33	33796.18	204052.17	23245.89	23526.01	168.19	3822.32	3362.62	71.22	strobe	deep.11
sent.clair3.concat.fofn-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.clair3.concat.fofn	0.2733	0:00:00	3.25	7.72	0.12	0.34	0.0	0.0	0.0	0.0	sent	clair3.concat.fofn
bwa2a.deep.bcfstat-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.deep.bcfstat	1.2742	0:00:01	15.0	2497.87	12.47	12.7	12.41	0.0	43.98	2.49	bwa2a	deep.bcfstat
strobe.deep.16-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.16	672.5973	0:11:12	33951.41	204139.84	23296.39	23574.1	1456.0	3311.49	1435.81	51.46	strobe	deep.16
sent.deep.22-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.deep.22	205.8653	0:03:25	34030.54	205744.16	23305.81	23586.49	0.0	2898.27	2408.29	49.39	sent	deep.22
bwa2a.oct.23-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.oct.23	248.6955	0:04:08	15185.73	19718.62	15049.17	15057.64	1053.65	17.47	1399.02	3479.87	bwa2a	oct.23
sent.lfq2.bcfstat-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.lfq2.bcfstat	0.4255	0:00:00	15.0	39.31	7.02	8.73	0.0	0.0	0.0	0.02	sent	lfq2.bcfstat
bwa2a.mrkdup-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.mrkdup	591.4633	0:09:51	402843.68	438788.68	402967.87	402968.54	12926.22	39356.92	7574.14	136.72	bwa2a	mrkdup
sent.dysgu.sv.vcf-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.dysgu.sv.vcf	1198.6815	0:19:58	2955.44	15423.45	2956.12	2956.8	865.75	9214.18	91.8	1100.58	sent	dysgu.sv.vcf
strobe.deep.17-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.17	615.5508	0:10:15	34112.35	203932.91	23380.13	23658.04	318.01	3481.04	1712.25	53.28	strobe	deep.17
sent.deep.23-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.deep.23	402.1128	0:06:42	33869.57	205115.77	23262.65	23540.77	1290.25	2937.89	3200.04	55.35	sent	deep.23
strobe.sentd.rtgvcfstats-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.sentd.rtgvcfstats	2.5047	0:00:02	1348.08	714171.71	1349.93	1350.44	16.41	0.09	136.38	3.64	strobe	sentd.rtgvcfstats
bwa2a.alignstats-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.alignstats	608.0853	0:10:08	978.0	2289.37	983.84	984.02	1.58	0.06	139.08	845.99	bwa2a	alignstats
bwa2a.oct.22-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.oct.22	601.9971	0:10:01	9882.7	14691.99	9885.76	9892.05	647.45	15.55	1153.29	6943.17	bwa2a	oct.22
bwa2a.sentd.23-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.sentd.23	45.8261	0:00:45	5327.02	14772.46	5449.19	5456.09	1237.8	39.87	1487.92	683.61	bwa2a	sentd.23
strobe.sentd.concat.fofn-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.sentd.concat.fofn	0.357	0:00:00	3.25	7.72	0.13	0.34	0.0	0.0	0.0	0.0	strobe	sentd.concat.fofn
strobe.mrkdup-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.mrkdup	507.4446	0:08:27	341728.19	380372.26	341855.23	341855.9	4372.56	34561.77	8235.46	176.09	strobe	mrkdup
bwa2a.deep.23-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.deep.23	506.436	0:08:26	33970.43	204147.75	23223.61	23503.29	1680.77	2908.81	2449.98	63.11	bwa2a	deep.23
strobe.clair3.24-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.clair3.24	372.4333	0:06:12	9979.55	18931.93	8780.8	8952.42	329.24	331.91	305.91	244.01	strobe	clair3.24
strobe.deep.14-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.14	612.802	0:10:12	33763.32	205909.89	23250.77	23528.69	1562.69	3201.49	1787.54	50.31	strobe	deep.14
sent.qmap-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.qmap	2025.588	0:33:45	25744.41	158766.87	25751.63	25752.21	275.28	1.34	321.47	6511.75	sent	qmap
strobe.oct.22-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.oct.22	756.2383	0:12:36	9631.72	14286.39	9635.09	9640.14	0.01	11.34	912.1	6898.15	strobe	oct.22
bwa2a.sentd.22-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.sentd.22	75.2135	0:01:15	6724.82	15860.46	6829.57	6836.28	620.41	36.95	1058.41	796.02	bwa2a	sentd.22
bwa2a.mrkdup.sort.picard-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.mrkdup.sort.picard	101.5563	0:01:41	1409.15	33939.93	1408.02	1408.91	845.73	1.07	76.22	77.76	bwa2a	mrkdup.sort.picard
strobe.deep.15-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.15	507.5892	0:08:27	34128.83	204440.46	23312.52	23591.01	0.0	3160.39	1650.85	35.6	strobe	deep.15
bwa2a.deep.22-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.deep.22	246.4416	0:04:06	34103.88	205997.57	23279.96	23560.09	957.8	2875.59	1816.53	36.15	bwa2a	deep.22
strobe.oct.23-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.oct.23	267.1789	0:04:27	15540.2	21177.99	15546.87	15553.11	1176.17	25.38	1431.64	3825.48	strobe	oct.23
bwa2a.samt-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.samt	1405.7103	0:23:25	12.0	666.55	7.06	8.45	253.27	0.12	124.38	286.66	bwa2a	samt
sent.deep.merge-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.deep.merge	3.4326	0:00:03	19.5	1044.46	18.23	18.47	0.09	46.89	31.26	1.33	sent	deep.merge
sent.lfq2.24-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.lfq2.24	1301.2736	0:21:41	94.8	268.66	88.39	88.67	65.43	0.9	98.09	1276.4	sent	lfq2.24
strobe.mrkdup.sort.picard-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.mrkdup.sort.picard	100.7299	0:01:40	1409.35	33744.93	1411.66	1412.53	0.0	1.17	78.04	79.72	strobe	mrkdup.sort.picard
sent.deep.peddy-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.deep.peddy	7.9794	0:00:07	2159.0	43346.02	267.94	381.14	0.07	15.78	36.97	0.0	sent	deep.peddy
strobe.sentd.concordance-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.sentd.concordance	251.0148	0:04:11	11600.59	716739.74	11609.88	11610.95	529.46	140.45	22.78	57.28	strobe	sentd.concordance
bwa2a.clair3.24-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.clair3.24	425.9889	0:07:05	9965.56	19238.72	8734.0	8906.74	655.83	445.57	245.09	297.66	bwa2a	clair3.24
strobe.deep.19-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.19	331.8227	0:05:31	33987.37	204299.32	23246.42	23524.73	933.86	3195.36	2882.75	54.01	strobe	deep.19
sent.clair3.22-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.clair3.22	858.5353	0:14:18	10683.89	19374.02	9507.61	9679.21	1197.52	613.94	241.46	671.24	sent	clair3.22
strobe.sentd.24-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.sentd.24	102.8278	0:01:42	11917.11	18036.4	12018.2	12025.11	473.42	30.28	1360.37	0.7	strobe	sentd.24
strobe.deep.18-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.18	521.5303	0:08:41	33842.7	206368.7	23274.56	23552.89	110.15	3027.32	1747.04	42.27	strobe	deep.18
bwa2a.lfq2.24-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.lfq2.24	1350.2754	0:22:30	98.93	273.65	92.35	92.67	130.14	0.93	98.99	1336.62	bwa2a	lfq2.24
sent.mosdepth-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.mosdepth	275.372	0:04:35	1987.5	34783.26	1992.62	1993.26	0.0	4980.78	180.2	497.11	sent	mosdepth
sent.clair3.23-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.clair3.23	702.5175	0:11:42	52908.93	92905.92	48322.57	48499.82	1157.74	1298.36	1214.9	16.58	sent	clair3.23
strobe.samt-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.samt	1371.2919	0:22:51	13.5	670.17	7.04	8.39	1.41	0.12	123.42	280.3	strobe	samt
sent.norm_cov_eveness-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.norm_cov_eveness	632.1117	0:10:32	1941.0	34781.26	1943.85	1944.41	830.36	6853.99	18.17	10.71	sent	norm_cov_eveness
strobe.sentd.22-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.sentd.22	71.7983	0:01:11	5434.71	15371.52	5515.49	5522.38	541.89	35.04	1034.72	743.19	strobe	sentd.22
sent.lfq2.rtgvcfstats-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.lfq2.rtgvcfstats	1.0183	0:00:01	69.0	712872.47	41.16	48.46	0.0	0.09	0.0	0.18	sent	lfq2.rtgvcfstats
sent.clair3.24-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.clair3.24	385.107	0:06:25	10567.1	19326.27	9366.61	9539.74	856.21	442.62	275.76	264.56	sent	clair3.24
bwa2a.lfq2.23-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.lfq2.23	3000.799	0:50:00	170.5	342.65	164.44	164.62	950.57	5.2	99.54	2986.95	bwa2a	lfq2.23
dirsetup-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	dirsetup	0.2252	0:00:00	14.1	31.32	7.04	7.47	0.0	0.0	0.0	0.0	dirsetup	NA
sent.clair3.concordance-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.clair3.concordance	235.347	0:03:55	15402.2	716739.74	15416.57	15417.64	652.84	142.83	23.36	55.07	sent	clair3.concordance
strobe.oct.concat.fofn-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.oct.concat.fofn	0.2648	0:00:00	3.25	7.72	0.12	0.34	0.0	0.0	0.0	0.0	strobe	oct.concat.fofn
bwa2a.alNsort-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.alNsort	2316.4667	0:38:36	293007.5	758153.7	290117.92	290119.71	17030.53	76170.23	10345.24	3567.85	bwa2a	alNsort
strobe.sentd.23-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.sentd.23	46.6646	0:00:46	4914.43	14388.4	5027.65	5034.57	1120.99	40.42	1567.94	731.97	strobe	sentd.23
strobe.lfq2.concordance-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.lfq2.concordance	31.3539	0:00:31	9195.07	716739.74	9210.16	9210.94	90.49	171.2	183.14	58.17	strobe	lfq2.concordance
bwa2a.lfq2.22-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.lfq2.22	2071.5402	0:34:31	81.63	255.2	74.91	75.45	534.05	4.19	99.29	2056.87	bwa2a	lfq2.22
bwa2a.deep.rtgvcfstats-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.deep.rtgvcfstats	1.8641	0:00:01	1277.35	703306.53	1272.05	1272.41	7.05	0.09	150.63	2.99	bwa2a	deep.rtgvcfstats
bwa2a.deep.peddy-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.deep.peddy	7.7611	0:00:07	145.31	2553.96	141.18	141.73	6.05	8.04	38.77	0.0	bwa2a	deep.peddy
sent.alNsort-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.alNsort	1232.0186	0:20:32	86367.91	122352.46	86599.26	86609.33	42839.34	112762.64	11645.63	187.23	sent	alNsort
bwa2a.dysgu.sv.vcf.sort-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.dysgu.sv.vcf.sort	6.5978	0:00:06	82.5	676.93	80.77	81.17	61.88	71.75	79.71	0.99	bwa2a	dysgu.sv.vcf.sort
strobe.oct.rtgvcfstats-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.oct.rtgvcfstats	2.1844	0:00:02	1269.04	714236.71	1234.27	1242.51	8.62	0.11	121.82	2.87	strobe	oct.rtgvcfstats
bwa2a.clair3.23-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.clair3.23	709.8739	0:11:49	58314.98	96778.02	54005.62	54182.73	1377.53	1282.06	1581.74	400.03	bwa2a	clair3.23
bwa2a.deep.concat.fofn-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.deep.concat.fofn	0.3176	0:00:00	15.0	40.47	7.64	8.77	0.0	0.0	0.0	0.03	bwa2a	deep.concat.fofn
sent.lfq2.22-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.lfq2.22	1443.6878	0:24:03	84.63	258.86	78.54	79.36	319.16	4.18	98.74	1425.49	sent	lfq2.22
strobe.lfq2.peddy-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.lfq2.peddy	7.5799	0:00:07	133.5	2544.98	132.96	133.69	0.0	8.02	34.44	0.0	strobe	lfq2.peddy
sent.sentd.concordance-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.sentd.concordance	263.6091	0:04:23	14125.37	716739.74	14122.2	14130.96	1294.54	353.59	41.54	109.94	sent	sentd.concordance
bwa2a.tiddit.sv.vcf.sort-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.tiddit.sv.vcf.sort	0.5269	0:00:00	13.5	39.26	6.84	9.02	0.0	0.0	0.0	0.02	bwa2a	tiddit.sv.vcf.sort
sent.lfq2.concat.fofn-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.lfq2.concat.fofn	0.3345	0:00:00	15.0	39.48	10.53	11.44	0.0	0.0	0.0	0.04	sent	lfq2.concat.fofn
bwa2a.deep.merge-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.deep.merge	3.0601	0:00:03	19.5	980.36	17.25	17.43	2.43	46.31	17.74	0.82	bwa2a	deep.merge
bwa2a.sentd.bcfstat-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.sentd.bcfstat	1.3512	0:00:01	22.5	2504.16	18.95	19.33	5.25	0.0	41.5	2.48	bwa2a	sentd.bcfstat
bwa2a.clair3.22-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.clair3.22	822.2263	0:13:42	9434.64	18404.29	8239.5	8412.84	993.35	603.81	257.17	642.0	bwa2a	clair3.22
fastqc-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	fastqc	2410.4122	0:40:10	5283.01	11140.43	5282.17	5282.43	266.82	1.59	197.96	4771.46	fastqc	NA
sent.goleft-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.goleft	1.6671	0:00:01	48.0	125.09	52.82	53.0	28.44	6.75	74.88	1.28	sent	goleft
strobe.lfq2.merge-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.lfq2.merge	1.0427	0:00:01	13.5	596.76	10.07	10.83	1.16	8.97	0.0	0.39	strobe	lfq2.merge
sent.lfq2.23-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.lfq2.23	2835.9178	0:47:15	169.01	342.65	162.72	163.55	387.0	5.2	99.01	2807.7	sent	lfq2.23
bwa2a.lfq2.bcfstat-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.lfq2.bcfstat	0.4245	0:00:00	13.5	39.25	6.78	9.01	0.0	0.0	0.0	0.02	bwa2a	lfq2.bcfstat
sent.deep.rtgvcfstats-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.deep.rtgvcfstats	1.9933	0:00:01	1185.65	714236.71	1150.88	1159.18	7.14	0.09	137.44	2.92	sent	deep.rtgvcfstats
strobe.deep.3-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.3	820.5713	0:13:40	34115.82	204469.95	23294.8	23573.03	0.01	4542.41	3117.65	81.57	strobe	deep.3
bwa2a.dysgu.sv.vcf-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.dysgu.sv.vcf	1061.0823	0:17:41	2795.56	11843.24	2788.34	2788.59	768.47	8519.56	97.96	1039.35	bwa2a	dysgu.sv.vcf
sent.oct.peddy-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.oct.peddy	8.1185	0:00:08	2117.27	43338.69	253.34	366.07	0.0	15.56	33.54	0.0	sent	oct.peddy
bwa2a.sentd.rtgvcfstats-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.sentd.rtgvcfstats	2.545	0:00:02	1351.5	703176.53	1358.74	1359.32	0.0	0.09	135.7	3.65	bwa2a	sentd.rtgvcfstats
strobe.deep.2-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.2	1006.5007	0:16:46	33996.61	204210.31	23312.0	23590.13	697.56	5053.85	2802.44	91.39	strobe	deep.2
strobe.tiddit.sv.vcf-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.tiddit.sv.vcf	1459.2361	0:24:19	10997.15	52148.21	10826.0	10852.28	6187.47	8227.52	368.24	1951.12	strobe	tiddit.sv.vcf
bwa2a.oct.bcfstat-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.oct.bcfstat	0.9597	0:00:00	13.5	2563.48	13.57	13.75	8.98	0.0	0.0	1.97	bwa2a	oct.bcfstat
strobe.deep.concordance-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.concordance	251.2069	0:04:11	8524.25	716739.74	8520.25	8528.99	1.11	146.14	21.75	55.24	strobe	deep.concordance
bwa2a.norm_cov_eveness-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.norm_cov_eveness	623.7049	0:10:23	1942.5	34651.26	1944.6	1944.8	603.82	6946.96	15.44	7.14	bwa2a	norm_cov_eveness
bwa2a.lfq2.rtgvcfstats-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.lfq2.rtgvcfstats	1.4739	0:00:01	76.5	712936.66	76.12	76.7	56.64	0.09	0.0	0.33	bwa2a	lfq2.rtgvcfstats
sent.oct.merge-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.oct.merge	3.2944	0:00:03	19.5	1748.05	17.55	17.71	8.83	62.52	48.37	1.23	sent	oct.merge
strobe.deep.1-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.1	1367.0145	0:22:47	34154.28	206330.48	23646.23	23926.59	2785.29	5392.32	2311.2	104.14	strobe	deep.1
strobe.sentd.peddy-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.sentd.peddy	7.6012	0:00:07	135.0	2545.0	131.34	132.01	0.12	8.13	34.73	0.0	strobe	sentd.peddy
sent.mrkdup-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.mrkdup	618.7007	0:10:18	416378.5	446481.27	415331.75	415332.42	286.75	40193.29	8535.11	318.03	sent	mrkdup
bwa2a.sentd.concat.fofn-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.sentd.concat.fofn	0.2762	0:00:00	3.38	7.72	0.13	0.35	0.0	0.0	0.0	0.0	bwa2a	sentd.concat.fofn
sent.oct.22-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.oct.22	548.5647	0:09:08	10441.16	15420.54	10444.67	10449.35	797.92	15.52	1238.37	6793.87	sent	oct.22
bwa2a.vb2-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.vb2	301.3024	0:05:01	494.2	514.68	487.65	487.96	13.36	4.21	180.63	544.4	bwa2a	vb2
bwa2a.lfq2.concat.fofn-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.lfq2.concat.fofn	0.3214	0:00:00	13.5	40.47	11.06	11.96	0.0	0.0	0.0	0.03	bwa2a	lfq2.concat.fofn
strobe.deep.24-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.24	130.812	0:02:10	33920.82	205308.52	23357.61	23635.94	582.55	2301.51	1625.72	27.11	strobe	deep.24
strobe.lfq2.bcfstat-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.lfq2.bcfstat	0.3839	0:00:00	3.0	7.71	0.12	0.36	0.0	0.0	0.0	0.0	strobe	lfq2.bcfstat
sent.sentd.bcfstat-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.sentd.bcfstat	1.4589	0:00:01	21.0	2760.21	16.71	17.99	19.4	0.0	43.86	2.55	sent	sentd.bcfstat
strobe.dysgu.sv.vcf-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.dysgu.sv.vcf	1272.1867	0:21:12	3640.16	13359.29	3631.98	3634.07	1241.29	8917.14	96.72	1237.79	strobe	dysgu.sv.vcf
strobe.sentd.merge-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.sentd.merge	4.9265	0:00:04	21.0	1044.59	18.58	18.77	178.49	107.83	68.79	1.4	strobe	sentd.merge
sent.deep.concat.fofn-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.deep.concat.fofn	0.3672	0:00:00	18.0	41.55	11.41	12.02	0.0	0.0	0.0	0.04	sent	deep.concat.fofn
bwa2a.oct.concordance-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.oct.concordance	242.3587	0:04:02	10125.74	705744.48	10138.34	10138.86	1.11	156.76	24.6	60.2	bwa2a	oct.concordance
strobe.alignstats-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.alignstats	541.0421	0:09:01	969.0	1895.95	973.86	974.11	264.58	0.07	135.99	736.04	strobe	alignstats
sent.oct.23-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.oct.23	248.0172	0:04:08	15018.37	19637.96	14814.59	14820.7	0.0	17.51	1395.98	3463.18	sent	oct.23
strobe.deep.22-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.22	203.8466	0:03:23	34251.32	204840.69	23299.6	23578.16	649.75	2772.25	2247.7	44.58	strobe	deep.22
strobe.clair3.rtgvcfstats-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.clair3.rtgvcfstats	1.5206	0:00:01	384.48	714301.71	371.46	379.47	0.0	0.09	145.96	2.37	strobe	clair3.rtgvcfstats
sent.oct.24-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.oct.24	602.0754	0:10:02	9555.14	14838.62	9565.05	9570.33	638.2	0.09	573.34	3452.4	sent	oct.24
strobe.oct.bcfstat-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.oct.bcfstat	0.9379	0:00:00	18.0	2500.17	14.26	14.47	7.65	0.0	0.0	1.99	strobe	oct.bcfstat
sent.tiddit.sv.vcf.sort-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.tiddit.sv.vcf.sort	0.7509	0:00:00	19.5	39.11	17.18	17.48	12.43	12.43	0.0	0.03	sent	tiddit.sv.vcf.sort
bwa2a.oct.merge-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.oct.merge	3.2223	0:00:03	21.0	1620.04	17.83	17.98	4.83	62.34	57.18	1.43	bwa2a	oct.merge
sent.manta-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.manta	466.8884	0:07:46	5819.18	13381.14	5857.78	5861.03	32569.3	79.21	3722.07	19101.78	sent	manta
strobe.deep.7-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.7	764.7734	0:12:44	37064.86	202636.58	23250.9	23530.83	1999.25	4257.79	2565.03	57.56	strobe	deep.7
strobe.vb2-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.vb2	404.2421	0:06:44	530.16	551.25	524.51	524.94	22851.67	4.41	140.3	567.14	strobe	vb2
strobe.deep.23-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.23	355.0479	0:05:55	33810.36	205667.98	23233.89	23511.75	523.8	2851.5	3123.41	41.93	strobe	deep.23
bwa2a.oct.peddy-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.oct.peddy	7.7242	0:00:07	123.0	2534.22	127.27	127.83	0.0	7.93	29.17	0.0	bwa2a	oct.peddy
sent.oct.concordance-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.oct.concordance	235.7916	0:03:55	9576.69	716804.74	9588.04	9589.24	94.51	179.23	25.13	59.83	sent	oct.concordance
strobe.deep.6-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.6	1037.0181	0:17:17	33963.44	205359.76	23298.71	23577.21	1083.24	4308.37	1879.95	73.65	strobe	deep.6
bwa2a.clair3.bcfstat-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.clair3.bcfstat	0.9144	0:00:00	16.5	3009.94	12.57	12.73	3.04	0.0	0.0	1.98	bwa2a	clair3.bcfstat
sent.vb2-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.vb2	385.7103	0:06:25	490.62	513.02	485.88	486.27	9883.49	4.23	153.19	591.09	sent	vb2
sent.clair3.merge-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.clair3.merge	2.5168	0:00:02	18.0	1108.34	17.98	18.13	0.36	17.63	0.0	0.85	sent	clair3.merge
strobe.deep.20-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.20	420.6757	0:07:00	33830.73	204182.68	23232.41	23509.82	109.22	3174.46	2105.56	54.99	strobe	deep.20
bwa2a.mosdepth-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.mosdepth	295.8601	0:04:55	1990.5	34781.26	1995.98	1996.5	0.0	5092.79	181.12	536.29	bwa2a	mosdepth
sent.deep.bcfstat-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.deep.bcfstat	1.0779	0:00:01	16.5	2753.87	12.55	12.77	7.25	0.0	0.0	1.97	sent	deep.bcfstat
bwa2a.sentd.merge-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.sentd.merge	4.8232	0:00:04	19.5	1108.6	18.3	18.61	10.07	107.68	57.91	0.93	bwa2a	sentd.merge
strobe.deep.5-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.5	1155.0313	0:19:15	34227.65	203857.15	23346.97	23626.53	1161.43	4091.93	1657.99	69.51	strobe	deep.5
sent.clair3.peddy-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.clair3.peddy	10.4086	0:00:10	148.02	2559.44	146.8	147.74	131.0	8.08	31.93	0.0	sent	clair3.peddy
strobe.qmap-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.qmap	1969.7384	0:32:49	25738.78	159310.87	25742.9	25745.02	955.96	1.51	324.23	6480.04	strobe	qmap
strobe.deep.21-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.21	481.4569	0:08:01	34753.55	202673.63	23005.68	23283.57	1496.66	2610.64	916.05	37.76	strobe	deep.21
strobe.clair3.concat.fofn-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.clair3.concat.fofn	0.267	0:00:00	3.38	7.72	0.13	0.34	0.0	0.0	0.0	0.0	strobe	clair3.concat.fofn
bwa2a.sentd.peddy-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.sentd.peddy	13.0915	0:00:13	1245.62	38730.8	220.95	283.05	142.66	8.08	26.06	0.0	bwa2a	sentd.peddy
bwa2a.clair3.concordance-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.clair3.concordance	233.4378	0:03:53	9418.8	716739.74	9432.76	9433.28	761.43	166.29	24.7	58.11	bwa2a	clair3.concordance
strobe.deep.4-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.4	852.8652	0:14:12	34215.5	204248.98	23478.81	23756.6	2652.38	4500.17	2778.03	79.82	strobe	deep.4
strobe.clair3.concordance-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.clair3.concordance	235.8064	0:03:55	10431.48	716741.58	10441.59	10442.37	91.39	174.35	26.42	63.11	strobe	clair3.concordance
sent.oct.bcfstat-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.oct.bcfstat	1.084	0:00:01	13.5	2691.57	13.68	13.81	2.93	0.0	0.0	1.94	sent	oct.bcfstat
strobe.deep.peddy-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.peddy	10.7111	0:00:10	1220.66	38681.39	181.32	243.15	135.29	8.09	28.53	0.0	strobe	deep.peddy
strobe.lfq2.24-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.lfq2.24	636.5225	0:10:36	81.77	255.51	75.57	75.75	271.46	0.79	96.68	615.55	strobe	lfq2.24
bwa2a.lfq2.peddy-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.lfq2.peddy	10.1541	0:00:10	1227.97	38730.04	219.44	280.24	128.43	7.98	31.0	0.0	bwa2a	lfq2.peddy
bwa2a.goleft-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.goleft	1.3516	0:00:01	37.5	125.78	41.05	41.29	19.88	5.45	50.75	0.94	bwa2a	goleft
bwa2a.clair3.concat.fofn-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.clair3.concat.fofn	0.2645	0:00:00	3.25	7.72	0.13	0.34	0.0	0.0	0.0	0.0	bwa2a	clair3.concat.fofn
strobe.deep.merge-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.merge	3.0861	0:00:03	19.5	1108.36	17.93	18.16	0.02	46.03	38.52	0.91	strobe	deep.merge
sent.oct.rtgvcfstats-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.oct.rtgvcfstats	2.1665	0:00:02	1324.79	714236.71	1311.84	1319.88	6.13	0.09	125.92	3.01	sent	oct.rtgvcfstats
bwa2a.lfq2.merge-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.lfq2.merge	1.1076	0:00:01	15.0	596.76	10.37	11.29	0.0	0.0	0.0	0.0	bwa2a	lfq2.merge
sent.oct.concat.fofn-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.oct.concat.fofn	0.2722	0:00:00	3.38	7.72	0.12	0.34	0.0	0.0	0.0	0.0	sent	oct.concat.fofn
strobe.dysgu.sv.vcf.sort-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.dysgu.sv.vcf.sort	1.4198	0:00:01	115.5	676.93	105.89	106.49	98.71	113.61	0.0	1.54	strobe	dysgu.sv.vcf.sort
sent.samt-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.samt	1423.0082	0:23:43	13.5	670.17	7.66	8.98	369.39	0.12	124.27	199.67	sent	samt
strobe.deep.8-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.8	563.1193	0:09:23	33935.97	204484.95	23308.86	23587.0	57.12	3792.97	3072.98	51.25	strobe	deep.8
strobe.sentd.1-24-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.sentd.1-24	516.9326	0:08:36	23082.17	27045.78	23084.49	23091.41	36582.17	6684.55	14881.68	38004.98	strobe	sentd.1-24
sent.sentd.22-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.sentd.22	75.694	0:01:15	6626.99	15668.46	6415.23	6422.14	675.53	37.12	1066.0	806.91	sent	sentd.22
kat-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	kat	1117.5975	0:18:37	34053.0	38656.26	34108.16	34108.37	58.29	0.09	753.55	7290.64	kat	NA
strobe.oct.peddy-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.oct.peddy	7.7012	0:00:07	129.0	2535.44	128.03	128.54	0.14	7.98	31.38	0.0	strobe	oct.peddy
bwa2a.qmap-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.qmap	2181.169	0:36:21	26949.88	158309.84	26956.2	26957.56	478.49	1.43	308.5	6727.38	bwa2a	qmap
strobe.goleft-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.goleft	1.2743	0:00:01	34.5	126.53	41.06	41.24	0.0	5.85	54.17	1.01	strobe	goleft
strobe.deep.9-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.9	550.0829	0:09:10	33740.69	204386.27	23277.55	23555.53	403.01	3723.62	2551.05	64.94	strobe	deep.9
bwa2a.clair3.rtgvcfstats-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.clair3.rtgvcfstats	1.4381	0:00:01	459.23	714301.71	441.88	449.86	0.0	0.11	140.71	2.24	bwa2a	clair3.rtgvcfstats
strobe.oct.merge-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.oct.merge	3.329	0:00:03	19.5	1620.04	17.84	17.97	4.82	61.19	48.04	1.51	strobe	oct.merge
sent.sentd.23-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.sentd.23	50.453	0:00:50	8993.22	14580.46	9128.28	9135.19	1223.04	125.26	957.88	1889.38	sent	sentd.23
bwa2a.oct.concat.fofn-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.oct.concat.fofn	0.2939	0:00:00	3.38	7.72	0.12	0.34	0.0	0.0	0.0	0.0	bwa2a	oct.concat.fofn
sent.dysgu.sv.vcf.sort-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.dysgu.sv.vcf.sort	0.9777	0:00:00	121.5	145.16	116.01	116.96	5.39	45.42	0.0	0.25	sent	dysgu.sv.vcf.sort
sent.deep.concordance-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.deep.concordance	233.6439	0:03:53	10153.64	716739.74	10166.12	10166.91	90.51	169.26	25.71	60.82	sent	deep.concordance
sent.mrkdup.sort.picard-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.mrkdup.sort.picard	99.476	0:01:39	1400.32	33679.93	1400.36	1401.36	0.0	1.06	79.24	79.76	sent	mrkdup.sort.picard
sent.sentd.24-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.sentd.24	97.6823	0:01:37	11947.8	17716.46	7853.54	7860.45	560.11	29.25	1498.36	825.42	sent	sentd.24
strobe.mosdepth-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.mosdepth	279.4456	0:04:39	3468.0	34651.26	3473.56	3474.12	116.21	4892.72	161.02	448.53	strobe	mosdepth
bwa2a.sentd.concordance-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.sentd.concordance	242.8063	0:04:02	14188.26	705744.48	14210.16	14210.97	531.4	142.31	21.99	54.67	bwa2a	sentd.concordance
strobe.deep.rtgvcfstats-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.rtgvcfstats	1.9063	0:00:01	1236.93	714171.71	1218.16	1218.87	7.02	0.09	140.24	2.89	strobe	deep.rtgvcfstats
bwa2a.lfq2.concordance-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.lfq2.concordance	31.8483	0:00:31	10348.9	716739.74	10338.88	10339.81	714.58	160.29	177.57	56.96	bwa2a	lfq2.concordance
sent.tiddit.sv.vcf-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.tiddit.sv.vcf	897.2881	0:14:57	8560.78	50201.09	8389.77	8417.66	6163.09	7668.08	513.76	1025.2	sent	tiddit.sv.vcf
bwa2a.manta-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.manta	378.999	0:06:18	5870.23	13382.31	5907.94	5911.11	34352.6	151.01	4605.64	0.09	bwa2a	manta
strobe.alNsort-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.alNsort	1060.1808	0:17:40	300695.95	3203411.86	300864.64	300867.35	11267.18	36477.88	5140.6	4421.54	strobe	alNsort
sent.sentd.peddy-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.sentd.peddy	10.629	0:00:10	2164.04	43339.29	263.22	376.55	128.73	15.87	26.39	0.0	sent	sentd.peddy
sent.lfq2.merge-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.lfq2.merge	1.1204	0:00:01	15.0	596.76	10.37	11.27	6.73	12.0	0.0	0.0	sent	lfq2.merge
strobe.sentd.bcfstat-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.sentd.bcfstat	1.3637	0:00:01	21.0	2503.91	18.63	18.86	21.78	0.0	46.32	2.53	strobe	sentd.bcfstat
strobe.deep.concat.fofn-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.deep.concat.fofn	2.0811	0:00:02	12.0	39.11	6.12	8.29	1.6	0.0	0.0	0.0	strobe	deep.concat.fofn
strobe.lfq2.22-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.lfq2.22	992.7576	0:16:32	68.95	242.3	64.27	64.94	274.94	3.62	98.29	975.91	strobe	lfq2.22
bwa2a.oct.rtgvcfstats-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.oct.rtgvcfstats	2.2115	0:00:02	1273.49	703178.27	1252.04	1252.39	0.0	0.09	122.38	2.82	bwa2a	oct.rtgvcfstats
sent.sentd.merge-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.sentd.merge	4.9823	0:00:04	18.0	1108.48	18.27	18.45	0.12	109.44	54.27	1.28	sent	sentd.merge
sent.lfq2.peddy-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.lfq2.peddy	10.0865	0:00:10	145.5	2554.55	145.2	145.87	128.07	7.93	30.54	0.0	sent	lfq2.peddy
strobe.clair3.bcfstat-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.clair3.bcfstat	0.9557	0:00:00	16.5	2561.94	12.45	12.59	0.0	0.0	0.0	1.94	strobe	clair3.bcfstat
strobe.lfq2.23-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.lfq2.23	2035.6979	0:33:55	168.95	342.51	162.46	163.18	403.18	5.08	99.46	2024.91	strobe	lfq2.23
sent.manta-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.manta	434.808	0:07:14	5884.28	13403.05	5914.86	5918.28	34528.44	80.17	3965.18	17622.75	sent	manta
bwa2a.clair3.24-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.clair3.24	251.7865	0:04:11	6387.99	13365.47	5998.9	6156.98	734.98	154.88	135.36	197.8	bwa2a	clair3.24
sent.deep.peddy-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.deep.peddy	8.0981	0:00:08	1594.9	33101.06	140.11	228.6	7.04	15.78	32.09	0.0	sent	deep.peddy
bwa2a.oct.concat.fofn-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.oct.concat.fofn	0.2649	0:00:00	3.38	7.72	0.12	0.34	0.0	0.0	0.0	0.0	bwa2a	oct.concat.fofn
sent.clair3.concordance-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.clair3.concordance	246.0538	0:04:06	9535.64	705744.48	9522.8	9531.41	74.37	141.86	24.2	60.35	sent	clair3.concordance
strobe.vb2-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.vb2	270.4267	0:04:30	493.23	515.68	488.62	488.92	280.37	4.22	167.87	454.11	strobe	vb2
strobe.oct.bcfstat-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.oct.bcfstat	1.13	0:00:01	18.0	2499.92	13.94	14.13	17.39	0.0	52.22	2.5	strobe	oct.bcfstat
sent.deep.concordance-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.deep.concordance	239.0035	0:03:59	19140.28	716739.74	19151.45	19151.79	1.11	154.77	24.3	58.58	sent	deep.concordance
strobe.sentd.rtgvcfstats-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.sentd.rtgvcfstats	3.3421	0:00:03	1354.05	714171.71	1359.85	1360.36	76.11	0.1	123.52	4.16	strobe	sentd.rtgvcfstats
sent.deep.merge-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.deep.merge	2.9688	0:00:02	19.5	1108.48	18.16	18.33	2.96	44.64	37.64	0.93	sent	deep.merge
strobe.lfq2.concordance-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.lfq2.concordance	33.0632	0:00:33	19054.52	716743.44	19068.08	19069.0	0.09	166.88	178.17	60.31	strobe	lfq2.concordance
bwa2a.dysgu.sv.vcf.sort-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.dysgu.sv.vcf.sort	7.3793	0:00:07	88.5	112.88	85.83	86.63	68.68	72.5	84.21	0.12	bwa2a	dysgu.sv.vcf.sort
strobe.sentd.concat.fofn-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.sentd.concat.fofn	0.2673	0:00:00	3.38	7.72	0.13	0.34	0.0	0.0	0.0	0.0	strobe	sentd.concat.fofn
strobe.deep.merge-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.deep.merge	3.1045	0:00:03	19.5	1108.49	17.95	18.18	2.37	44.75	35.73	0.94	strobe	deep.merge
bwa2a.goleft-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.goleft	1.5184	0:00:01	48.0	125.59	28.51	29.09	92.4	3.82	42.69	0.86	bwa2a	goleft
sent.clair3.23-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.clair3.23	1125.7148	0:18:45	20456.11	53737.38	18281.13	18458.33	2384.27	1897.74	471.08	987.5	sent	clair3.23
strobe.clair3.bcfstat-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.clair3.bcfstat	1.0476	0:00:01	15.0	3329.86	12.6	12.95	13.43	0.0	0.0	1.87	strobe	clair3.bcfstat
strobe.sentd.24-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.sentd.24	95.5958	0:01:35	11878.09	18709.89	11977.79	11984.69	268.99	23.83	698.08	1058.19	strobe	sentd.24
bwa2a.mrkdup.sort.picard-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.mrkdup.sort.picard	105.7973	0:01:45	1402.63	33744.93	1403.53	1404.59	850.01	1.07	72.98	77.76	bwa2a	mrkdup.sort.picard
strobe.deep.peddy-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.deep.peddy	10.7067	0:00:10	1225.73	38726.98	223.84	285.57	128.54	8.08	29.02	0.0	strobe	deep.peddy
sent.clair3.22-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.clair3.22	857.7927	0:14:17	11649.54	21261.57	10260.88	10435.01	879.6	601.62	234.33	681.82	sent	clair3.22
bwa2a.manta-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.manta	385.2358	0:06:25	5807.0	13405.67	5858.6	5861.48	34628.61	124.28	4064.31	4.71	bwa2a	manta
bwa2a.oct.rtgvcfstats-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.oct.rtgvcfstats	2.3472	0:00:02	1350.62	703241.53	1341.2	1346.74	11.84	0.11	145.09	3.61	bwa2a	oct.rtgvcfstats
bwa2a.deep.concat.fofn-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.deep.concat.fofn	0.2984	0:00:00	13.5	39.26	10.28	10.82	0.0	0.0	0.0	0.02	bwa2a	deep.concat.fofn
bwa2a.dysgu.sv.vcf-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.dysgu.sv.vcf	1004.717	0:16:44	2858.63	11916.11	2853.61	2854.31	4027.76	8465.6	97.64	980.99	bwa2a	dysgu.sv.vcf
bwa2a.lfq2.bcfstat-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.lfq2.bcfstat	0.4484	0:00:00	15.0	39.25	6.68	8.94	0.0	0.0	0.0	0.02	bwa2a	lfq2.bcfstat
strobe.oct.rtgvcfstats-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.oct.rtgvcfstats	2.295	0:00:02	1340.34	714236.71	1324.86	1328.29	0.0	0.09	137.56	3.38	strobe	oct.rtgvcfstats
strobe.sentd.23-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.sentd.23	44.1119	0:00:44	5094.85	14324.32	5222.85	5229.75	1818.77	53.05	2091.83	923.19	strobe	sentd.23
sent.clair3.24-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.clair3.24	243.0873	0:04:03	3824.4	13476.41	2855.94	3017.19	361.84	150.48	114.04	193.25	sent	clair3.24
strobe.sentd.22-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.sentd.22	82.5451	0:01:22	8997.62	15051.37	9128.76	9135.69	568.5	103.4	811.8	1364.21	strobe	sentd.22
sent.sentd.concordance-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.sentd.concordance	244.8569	0:04:04	16721.64	716804.74	16731.89	16732.75	771.14	174.22	25.27	62.43	sent	sentd.concordance
strobe.clair3.concordance-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.clair3.concordance	240.44	0:04:00	16393.73	716739.74	16405.89	16406.89	714.42	174.23	25.48	61.75	strobe	clair3.concordance
sent.sentd.bcfstat-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.sentd.bcfstat	1.5329	0:00:01	19.5	2566.67	17.49	17.68	19.08	0.0	41.06	2.66	sent	sentd.bcfstat
strobe.sentd.peddy-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.sentd.peddy	11.9181	0:00:11	1226.33	38728.38	221.56	283.61	128.66	8.13	25.51	0.0	strobe	sentd.peddy
bwa2a.clair3.22-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.clair3.22	849.5138	0:14:09	12380.39	22116.22	11096.45	11265.54	320.64	592.52	231.56	685.79	bwa2a	clair3.22
sent.qmap-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.qmap	2057.7732	0:34:17	24940.88	159177.85	24950.56	24951.31	764.82	1.55	318.46	6553.42	sent	qmap
bwa2a.deep.rtgvcfstats-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.deep.rtgvcfstats	1.8634	0:00:01	1258.45	703241.53	1236.0	1241.55	0.0	0.09	144.29	2.88	bwa2a	deep.rtgvcfstats
strobe.oct.concat.fofn-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.oct.concat.fofn	0.2972	0:00:00	3.25	7.72	0.12	0.34	0.0	0.0	0.0	0.0	strobe	oct.concat.fofn
bwa2a.clair3.23-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.clair3.23	1110.8055	0:18:30	20465.97	54219.57	18185.55	18356.85	2533.05	1907.79	523.68	17.01	bwa2a	clair3.23
strobe.sentd.merge-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.sentd.merge	5.5521	0:00:05	19.5	1108.63	18.59	18.82	110.42	119.83	68.28	1.54	strobe	sentd.merge
sent.sentd.peddy-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.sentd.peddy	8.2615	0:00:08	151.63	2559.34	146.42	147.0	16.46	15.87	39.56	0.0	sent	sentd.peddy
strobe.alignstats-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.alignstats	541.6631	0:09:01	970.5	2023.91	973.94	974.14	1.58	0.06	134.16	726.92	strobe	alignstats
bwa2a.sentd.24-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.sentd.24	72.1826	0:01:12	12551.19	18932.35	12631.32	12638.22	206.97	20.05	1085.18	783.35	bwa2a	sentd.24
strobe.oct.concordance-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.oct.concordance	253.3781	0:04:13	14427.57	716739.74	14421.29	14422.49	70.38	135.8	20.34	52.34	strobe	oct.concordance
strobe.lfq2.bcfstat-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.lfq2.bcfstat	0.4361	0:00:00	13.5	39.25	7.78	9.48	0.0	0.0	0.0	0.02	strobe	lfq2.bcfstat
sent.tiddit.sv.vcf-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.tiddit.sv.vcf	912.0512	0:15:12	8030.34	49315.37	7882.76	7912.63	7972.71	7569.63	537.84	989.17	sent	tiddit.sv.vcf
strobe.dysgu.sv.vcf-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.dysgu.sv.vcf	1199.345	0:19:59	3554.96	16030.58	3531.79	3546.53	1827.81	8869.92	97.65	1171.22	strobe	dysgu.sv.vcf
sent.sentd.merge-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.sentd.merge	5.6387	0:00:05	19.5	1108.64	18.22	18.35	119.41	121.14	58.72	1.37	sent	sentd.merge
strobe.qmap-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.qmap	1904.1384	0:31:44	25552.92	159050.86	25563.36	25563.89	237.53	1.51	335.01	6377.67	strobe	qmap
sent.vb2-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.vb2	347.7907	0:05:47	489.26	511.47	484.32	484.68	22522.93	4.25	148.96	518.27	sent	vb2
sent.sentd.rtgvcfstats-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.sentd.rtgvcfstats	2.5107	0:00:02	1357.44	714171.71	1346.54	1352.21	0.0	0.09	135.01	3.57	sent	sentd.rtgvcfstats
strobe.clair3.rtgvcfstats-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.clair3.rtgvcfstats	1.7426	0:00:01	1331.57	714301.71	1334.84	1335.54	5.21	0.09	162.11	3.04	strobe	clair3.rtgvcfstats
sent.oct.24-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.oct.24	406.2766	0:06:46	4700.55	9747.58	4713.32	4721.9	467.45	0.03	373.34	1517.17	sent	oct.24
bwa2a.oct.bcfstat-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.oct.bcfstat	1.1443	0:00:01	16.5	2628.13	14.3	14.52	5.44	0.0	51.45	2.53	bwa2a	oct.bcfstat
sent.sentd.22-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.sentd.22	83.0118	0:01:23	8995.18	14836.36	9128.89	9135.83	562.2	109.8	856.3	1367.81	sent	sentd.22
strobe.clair3.concat.fofn-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.clair3.concat.fofn	0.2704	0:00:00	3.25	7.72	0.13	0.34	0.0	0.0	0.0	0.0	strobe	clair3.concat.fofn
sent.sentd.concat.fofn-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.sentd.concat.fofn	0.2674	0:00:00	3.38	7.72	0.12	0.34	0.0	0.0	0.0	0.0	sent	sentd.concat.fofn
strobe.goleft-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.goleft	1.2975	0:00:01	36.0	126.03	40.66	40.94	19.88	5.98	51.69	0.96	strobe	goleft
bwa2a.deep.concordance-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.deep.concordance	246.1816	0:04:06	9078.65	705744.48	9086.3	9086.84	468.38	145.86	23.94	59.07	bwa2a	deep.concordance
sent.alignstats-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.alignstats	565.1889	0:09:25	979.5	2481.38	983.71	984.35	514.81	0.18	140.43	793.95	sent	alignstats
sent.sentd.23-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.sentd.23	45.4339	0:00:45	4294.47	14324.36	4420.65	4427.58	1876.47	53.38	2045.21	929.86	sent	sentd.23
bwa2a.tiddit.sv.vcf.sort-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.tiddit.sv.vcf.sort	0.7141	0:00:00	28.5	43.0	25.51	25.65	13.04	10.41	0.0	0.04	bwa2a	tiddit.sv.vcf.sort
bwa2a.sentd.peddy-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.sentd.peddy	7.8468	0:00:07	150.04	2559.57	146.56	147.13	0.0	8.08	33.38	0.0	bwa2a	sentd.peddy
sent.mrkdup-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.mrkdup	582.0518	0:09:42	356213.04	396080.9	353419.05	353419.71	15212.51	39859.95	8616.62	262.09	sent	mrkdup
sent.oct.22-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.oct.22	543.3463	0:09:03	10114.04	14911.59	9699.85	9701.43	803.56	15.44	1244.42	6762.05	sent	oct.22
sent.clair3.rtgvcfstats-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.clair3.rtgvcfstats	1.6523	0:00:01	361.93	703306.53	369.02	369.33	3.31	0.09	128.48	2.25	sent	clair3.rtgvcfstats
strobe.oct.peddy-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.oct.peddy	7.6554	0:00:07	133.5	2539.77	130.68	131.25	0.0	7.97	34.1	0.0	strobe	oct.peddy
sent.deep.rtgvcfstats-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.deep.rtgvcfstats	1.8503	0:00:01	1281.31	703306.53	1276.85	1277.37	0.0	0.09	148.82	2.93	sent	deep.rtgvcfstats
sent.sentd.24-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.sentd.24	77.1197	0:01:17	12457.12	19071.71	12551.39	12558.31	411.8	43.12	1041.71	1214.24	sent	sentd.24
sent.clair3.peddy-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.clair3.peddy	8.3326	0:00:08	1230.07	38728.55	221.52	283.51	0.0	8.08	41.16	0.0	sent	clair3.peddy
strobe.sentd.concordance-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.sentd.concordance	246.8468	0:04:06	23506.79	716739.74	23529.35	23529.79	526.39	144.86	23.28	57.73	strobe	sentd.concordance
sent.deep.bcfstat-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.deep.bcfstat	1.1258	0:00:01	16.5	2689.86	12.43	12.67	7.45	0.0	48.85	2.51	sent	deep.bcfstat
bwa2a.sentd.merge-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.sentd.merge	5.4868	0:00:05	19.5	1108.64	18.31	18.46	110.28	119.41	68.48	1.59	bwa2a	sentd.merge
strobe.oct.merge-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.oct.merge	4.1149	0:00:04	21.0	1940.02	19.01	19.2	1.62	76.41	52.24	1.47	strobe	oct.merge
sent.oct.23-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.oct.23	294.3104	0:04:54	16500.68	20656.81	16510.07	16511.62	73.17	44.34	2303.16	6779.07	sent	oct.23
sent.clair3.merge-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.clair3.merge	2.9511	0:00:02	19.5	1108.31	17.9	18.04	1.46	25.22	20.68	0.25	sent	clair3.merge
strobe.lfq2.rtgvcfstats-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.lfq2.rtgvcfstats	1.0307	0:00:01	73.5	712872.53	56.36	58.74	0.0	0.11	0.0	0.23	strobe	lfq2.rtgvcfstats
strobe.norm_cov_eveness-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.norm_cov_eveness	623.717	0:10:23	1941.0	34781.26	1943.25	1943.76	837.77	6922.37	19.08	6.74	strobe	norm_cov_eveness
sent.oct.bcfstat-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.oct.bcfstat	1.3608	0:00:01	15.0	2563.72	13.87	14.04	12.07	0.0	43.6	2.46	sent	oct.bcfstat
bwa2a.sentd.22-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.sentd.22	81.8454	0:01:21	9000.62	14900.35	9129.14	9136.06	550.01	110.68	858.63	1377.86	bwa2a	sentd.22
strobe.lfq2.concat.fofn-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.lfq2.concat.fofn	0.2956	0:00:00	13.5	39.48	10.54	11.43	0.0	0.0	0.0	0.03	strobe	lfq2.concat.fofn
bwa2a.lfq2.peddy-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.lfq2.peddy	7.6952	0:00:07	130.5	2533.94	123.32	124.0	0.21	7.97	33.39	0.0	bwa2a	lfq2.peddy
bwa2a.sentd.23-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.sentd.23	44.3188	0:00:44	4696.21	14388.35	4822.68	4829.58	1829.5	52.93	2107.17	935.7	bwa2a	sentd.23
bwa2a.qmap-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.qmap	2130.8777	0:35:30	25347.17	159374.87	25353.0	25353.99	310.47	1.57	310.57	6619.73	bwa2a	qmap
sent.deep.concat.fofn-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.deep.concat.fofn	0.3467	0:00:00	13.5	39.26	10.26	10.74	0.0	0.0	0.0	0.02	sent	deep.concat.fofn
bwa2a.lfq2.merge-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.lfq2.merge	1.2505	0:00:01	15.0	596.76	10.47	11.41	1.2	14.7	0.0	0.37	bwa2a	lfq2.merge
sent.clair3.concat.fofn-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.clair3.concat.fofn	0.2983	0:00:00	3.25	7.72	0.12	0.34	0.0	0.0	0.0	0.0	sent	clair3.concat.fofn
bwa2a.oct.concordance-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.oct.concordance	249.5625	0:04:09	9032.91	705874.48	9032.21	9033.08	0.11	146.81	23.7	59.99	bwa2a	oct.concordance
strobe.clair3.24-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.clair3.24	306.8337	0:05:06	3864.99	13177.38	3280.05	3445.83	441.28	87.54	114.66	248.96	strobe	clair3.24
sent.mrkdup.sort.picard-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.mrkdup.sort.picard	103.3632	0:01:43	1407.45	33744.93	1406.58	1407.91	816.0	1.07	75.4	78.33	sent	mrkdup.sort.picard
sent.norm_cov_eveness-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.norm_cov_eveness	651.2996	0:10:51	1554.0	2727.5	1557.6	1557.86	588.26	6915.38	14.2	10.27	sent	norm_cov_eveness
strobe.oct.24-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.oct.24	494.1646	0:08:14	5316.73	10644.23	5344.83	5346.5	452.5	0.08	685.13	3385.94	strobe	oct.24
sent.lfq2.rtgvcfstats-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.lfq2.rtgvcfstats	1.0635	0:00:01	72.0	712872.53	55.34	57.67	0.0	0.11	0.0	0.22	sent	lfq2.rtgvcfstats
bwa2a.deep.24-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.deep.24	108.2757	0:01:48	23274.84	192923.06	12322.61	12576.32	1120.23	2087.14	702.26	484.8	bwa2a	deep.24
strobe.lfq2.24-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.lfq2.24	88.4304	0:01:28	81.16	254.27	73.7	74.4	0.0	0.64	85.25	75.61	strobe	lfq2.24
bwa2a.sentd.concordance-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.sentd.concordance	244.8413	0:04:04	17090.56	716739.74	17097.3	17098.21	529.36	153.95	22.92	56.69	bwa2a	sentd.concordance
strobe.deep.rtgvcfstats-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.deep.rtgvcfstats	1.9805	0:00:01	1120.44	714301.71	1104.16	1112.39	7.36	0.09	135.72	2.91	strobe	deep.rtgvcfstats
sent.goleft-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.goleft	1.5366	0:00:01	46.5	125.34	51.87	52.11	92.4	3.86	42.92	0.87	sent	goleft
sent.lfq2.24-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.lfq2.24	173.6333	0:02:53	90.27	263.61	83.26	84.01	0.0	0.7	95.41	165.76	sent	lfq2.24
strobe.deep.concat.fofn-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.deep.concat.fofn	0.3003	0:00:00	16.5	40.47	11.07	11.69	0.0	0.0	0.0	0.03	strobe	deep.concat.fofn
strobe.tiddit.sv.vcf-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.tiddit.sv.vcf	1329.248	0:22:09	9881.42	51170.37	9734.18	9762.11	1593.0	8098.22	416.94	1769.85	strobe	tiddit.sv.vcf
bwa2a.vb2-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.vb2	344.9786	0:05:44	533.18	557.57	530.34	531.04	24171.97	4.38	118.19	407.88	bwa2a	vb2
bwa2a.alignstats-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.alignstats	570.7356	0:09:30	982.5	2353.47	983.88	984.04	1.58	0.07	139.11	794.16	bwa2a	alignstats
bwa2a.samt-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.samt	1360.1502	0:22:40	15.0	670.17	10.48	11.57	517.54	0.12	124.43	196.57	bwa2a	samt
sent.lfq2.concat.fofn-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.lfq2.concat.fofn	0.2932	0:00:00	13.5	40.47	11.05	11.92	0.0	0.0	0.0	0.03	sent	lfq2.concat.fofn
bwa2a.clair3.concordance-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.clair3.concordance	237.2277	0:03:57	15524.97	716739.74	15542.44	15543.22	87.3	177.21	25.86	62.55	bwa2a	clair3.concordance
strobe.lfq2.merge-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.lfq2.merge	1.1807	0:00:01	13.5	596.76	10.37	11.24	1.06	13.84	0.0	0.45	strobe	lfq2.merge
sent.alNsort-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.alNsort	1208.56	0:20:08	84127.73	118220.46	84358.49	84368.62	44952.82	112170.0	11233.2	367.54	sent	alNsort
sent.lfq2.22-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.lfq2.22	2032.0927	0:33:52	67.42	241.8	63.64	63.85	461.77	4.15	99.74	2026.87	sent	lfq2.22
bwa2a.sentd.bcfstat-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.sentd.bcfstat	1.3567	0:00:01	16.5	2694.79	17.65	17.86	18.7	0.0	46.38	2.32	bwa2a	sentd.bcfstat
strobe.mosdepth-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.mosdepth	274.6434	0:04:34	1987.5	34781.26	1993.41	1993.94	0.0	5000.75	164.09	451.04	strobe	mosdepth
strobe.lfq2.peddy-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.lfq2.peddy	7.5989	0:00:07	133.5	2541.0	131.97	132.84	0.0	8.02	34.35	0.0	strobe	lfq2.peddy
strobe.mrkdup.sort.picard-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.mrkdup.sort.picard	105.2796	0:01:45	1397.47	33744.93	1396.87	1397.77	525.92	1.11	71.17	75.46	strobe	mrkdup.sort.picard
sent.lfq2.23-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.lfq2.23	4184.2604	1:09:44	167.93	342.53	164.43	164.7	1432.91	7.27	99.17	4150.19	sent	lfq2.23
strobe.samt-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.samt	1348.9558	0:22:28	13.5	670.17	6.9	7.96	4.63	0.12	123.45	230.24	strobe	samt
sent.lfq2.peddy-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.lfq2.peddy	7.5934	0:00:07	136.5	2540.98	133.19	133.86	0.09	7.93	34.51	0.0	sent	lfq2.peddy
strobe.oct.22-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.oct.22	576.7008	0:09:36	9952.8	14589.52	9956.19	9962.48	1038.8	11.23	1202.29	6934.0	strobe	oct.22
bwa2a.clair3.peddy-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.clair3.peddy	10.0368	0:00:10	1227.54	38726.45	222.38	284.42	128.43	8.12	32.13	0.0	bwa2a	clair3.peddy
strobe.mrkdup-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.mrkdup	504.4296	0:08:24	388082.8	421854.35	388136.53	388137.19	4452.12	34305.81	8489.68	98.28	strobe	mrkdup
strobe.clair3.22-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.clair3.22	774.2664	0:12:54	10301.96	19191.76	9113.95	9287.82	1137.95	515.27	238.66	624.08	strobe	clair3.22
bwa2a.deep.23-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.deep.23	351.174	0:05:51	34529.52	203998.46	23252.89	23526.02	2842.98	3323.34	3865.66	39.39	bwa2a	deep.23
sent.dysgu.sv.vcf.sort-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.dysgu.sv.vcf.sort	6.8275	0:00:06	118.5	143.53	114.57	115.03	67.62	10.19	84.6	6.0	sent	dysgu.sv.vcf.sort
strobe.lfq2.23-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.lfq2.23	4292.458	1:11:32	168.89	342.5	162.56	163.44	1082.53	7.12	99.53	4272.34	strobe	lfq2.23
strobe.alNsort-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.alNsort	1058.6153	0:17:38	293395.98	3203613.86	293426.17	293428.86	11266.63	36845.22	4846.91	4260.92	strobe	alNsort
bwa2a.deep.bcfstat-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.deep.bcfstat	1.1712	0:00:01	15.0	2689.86	12.44	12.6	7.38	0.0	48.62	2.51	bwa2a	deep.bcfstat
sent.lfq2.merge-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.lfq2.merge	1.23	0:00:01	15.0	596.76	10.36	11.24	1.21	14.77	0.0	0.4	sent	lfq2.merge
bwa2a.lfq2.concordance-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.lfq2.concordance	33.6332	0:00:33	22827.07	716739.74	22840.11	22841.02	0.0	167.89	178.02	60.94	bwa2a	lfq2.concordance
strobe.oct.23-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.oct.23	393.3559	0:06:33	17761.94	23110.15	17767.41	17773.6	2041.9	40.01	2397.0	9428.98	strobe	oct.23
bwa2a.clair3.bcfstat-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.clair3.bcfstat	1.0927	0:00:01	16.5	2497.87	12.38	12.57	2.65	0.02	47.62	2.45	bwa2a	clair3.bcfstat
strobe.clair3.23-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.clair3.23	1070.7691	0:17:50	19525.62	52093.12	17207.05	17377.95	2181.96	1683.05	491.89	985.34	strobe	clair3.23
sent.oct.concordance-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.oct.concordance	244.099	0:04:04	21606.17	716739.74	21621.62	21622.5	763.91	157.97	23.64	58.25	sent	oct.concordance
bwa2a.clair3.merge-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.clair3.merge	2.693	0:00:02	21.0	1108.47	17.94	18.42	256.94	20.82	0.0	0.28	bwa2a	clair3.merge
strobe.lfq2.22-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.lfq2.22	1119.0826	0:18:39	67.39	241.77	61.7	62.51	0.06	3.52	97.96	1096.5	strobe	lfq2.22
bwa2a.deep.22-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.deep.22	211.5813	0:03:31	34084.64	204018.22	23289.62	23569.32	932.87	2841.65	3002.19	35.54	bwa2a	deep.22
sent.tiddit.sv.vcf.sort-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.tiddit.sv.vcf.sort	0.7711	0:00:00	16.5	39.48	14.86	15.17	14.86	12.22	0.0	0.02	sent	tiddit.sv.vcf.sort
bwa2a.lfq2.concat.fofn-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.lfq2.concat.fofn	0.322	0:00:00	13.5	40.47	10.68	11.57	0.0	0.0	0.0	0.03	bwa2a	lfq2.concat.fofn
bwa2a.deep.peddy-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.deep.peddy	7.8782	0:00:07	149.73	2559.46	118.89	133.38	0.0	8.05	39.47	0.0	bwa2a	deep.peddy
sent.oct.concat.fofn-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.oct.concat.fofn	0.2767	0:00:00	3.25	7.72	0.12	0.33	0.0	0.0	0.0	0.0	sent	oct.concat.fofn
sent.deep.24-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.deep.24	76.7707	0:01:16	34046.56	205034.4	23391.36	23672.08	382.41	2044.98	921.52	637.29	sent	deep.24
sent.oct.merge-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.oct.merge	3.9295	0:00:03	19.5	1940.01	18.88	19.35	13.59	78.36	64.11	1.35	sent	oct.merge
bwa2a.oct.24-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.oct.24	267.6974	0:04:27	4620.69	9001.02	4633.36	4639.49	367.06	0.04	462.64	1238.98	bwa2a	oct.24
bwa2a.deep.merge-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.deep.merge	2.87	0:00:02	19.5	1108.48	17.9	18.06	0.57	44.14	36.6	0.99	bwa2a	deep.merge
sent.oct.peddy-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.oct.peddy	11.2832	0:00:11	144.86	2559.35	131.39	132.19	139.94	15.56	29.57	0.0	sent	oct.peddy
strobe.clair3.merge-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.clair3.merge	3.1766	0:00:03	19.5	1033.38	15.23	16.51	2.43	31.39	0.0	0.88	strobe	clair3.merge
strobe.deep.bcfstat-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.deep.bcfstat	1.2227	0:00:01	16.5	3457.86	12.65	13.03	52.17	0.0	45.87	2.48	strobe	deep.bcfstat
strobe.clair3.peddy-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.clair3.peddy	7.6221	0:00:07	132.0	2539.68	105.25	119.27	0.0	8.18	34.01	0.0	strobe	clair3.peddy
bwa2a.lfq2.rtgvcfstats-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.lfq2.rtgvcfstats	1.0644	0:00:01	72.0	712872.53	36.6	40.65	0.0	0.11	0.0	0.22	bwa2a	lfq2.rtgvcfstats
bwa2a.norm_cov_eveness-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.norm_cov_eveness	623.8618	0:10:23	1939.5	34716.26	1944.86	1945.49	991.79	6908.28	10.36	7.01	bwa2a	norm_cov_eveness
sent.mosdepth-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.mosdepth	278.8555	0:04:38	1990.5	34716.26	1993.78	1993.94	52.14	4940.71	176.9	493.66	sent	mosdepth
fastqc-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	fastqc	2382.1171	0:39:42	5278.64	11140.43	5276.79	5277.78	272.89	1.62	197.81	4712.17	fastqc	NA
sent.clair3.bcfstat-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.clair3.bcfstat	0.8587	0:00:00	16.5	2561.87	10.14	11.45	0.48	0.0	0.0	1.94	sent	clair3.bcfstat
strobe.manta-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.manta	436.0292	0:07:16	6177.93	13548.91	6213.0	6216.24	29875.52	90.13	4428.54	19561.4	strobe	manta
bwa2a.lfq2.24-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.lfq2.24	131.8473	0:02:11	95.71	269.98	88.66	89.17	4.86	0.5	79.89	105.51	bwa2a	lfq2.24
sent.oct.rtgvcfstats-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.oct.rtgvcfstats	3.1407	0:00:03	1352.13	703241.53	1356.45	1356.75	57.01	0.09	66.92	3.41	sent	oct.rtgvcfstats
strobe.deep.24-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.deep.24	239.7487	0:03:59	9395.36	139629.47	4079.54	4316.77	620.77	2147.29	549.16	25.22	strobe	deep.24
bwa2a.tiddit.sv.vcf-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.tiddit.sv.vcf	863.8708	0:14:23	8081.67	49520.44	7911.04	7940.9	6150.77	7369.96	464.32	997.15	bwa2a	tiddit.sv.vcf
sent.samt-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.samt	1396.307	0:23:16	15.0	670.17	10.55	11.17	803.8	0.12	112.24	164.79	sent	samt
bwa2a.lfq2.23-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.lfq2.23	4011.0135	1:06:51	170.39	342.52	162.55	163.21	1158.81	7.24	99.07	3974.55	bwa2a	lfq2.23
strobe.deep.23-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.deep.23	384.2031	0:06:24	33944.67	205357.48	23302.29	23576.06	2771.36	3248.59	3408.0	42.91	strobe	deep.23
strobe.sentd.bcfstat-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.sentd.bcfstat	1.441	0:00:01	21.0	2502.05	14.51	15.78	2.25	0.0	41.65	2.48	strobe	sentd.bcfstat
bwa2a.clair3.concat.fofn-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.clair3.concat.fofn	0.3018	0:00:00	3.38	7.72	0.12	0.34	0.0	0.0	0.0	0.0	bwa2a	clair3.concat.fofn
sent.lfq2.concordance-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.lfq2.concordance	33.4002	0:00:33	22876.03	716739.74	22892.45	22893.38	84.04	168.88	179.03	61.03	sent	lfq2.concordance
bwa2a.oct.peddy-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.oct.peddy	7.7314	0:00:07	141.61	2553.94	141.04	141.63	5.75	7.93	34.02	0.0	bwa2a	oct.peddy
strobe.deep.22-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.deep.22	190.4539	0:03:10	33988.04	204904.59	23306.3	23584.69	916.38	2756.73	1848.56	32.88	strobe	deep.22
bwa2a.lfq2.22-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.lfq2.22	1687.1047	0:28:07	67.39	241.8	61.64	62.35	469.32	4.09	98.62	1663.62	bwa2a	lfq2.22
strobe.tiddit.sv.vcf.sort-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.tiddit.sv.vcf.sort	1.5174	0:00:01	13.5	39.11	9.45	10.35	26.08	17.18	33.63	0.75	strobe	tiddit.sv.vcf.sort
bwa2a.sentd.rtgvcfstats-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.sentd.rtgvcfstats	2.7127	0:00:02	1357.8	703241.53	1361.45	1361.81	0.0	0.1	125.5	0.0	bwa2a	sentd.rtgvcfstats
bwa2a.oct.merge-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.oct.merge	3.937	0:00:03	19.5	1940.01	18.88	19.1	11.74	77.99	56.54	1.37	bwa2a	oct.merge
bwa2a.alNsort-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.alNsort	2252.4194	0:37:32	290005.31	758564.59	292177.99	292179.76	17031.86	75901.47	10212.69	3559.58	bwa2a	alNsort
strobe.deep.concordance-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.deep.concordance	244.6609	0:04:04	10033.66	716739.74	10035.46	10044.14	79.15	155.91	24.15	60.1	strobe	deep.concordance
strobe.dysgu.sv.vcf.sort-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	strobe.dysgu.sv.vcf.sort	1.2406	0:00:01	177.0	205.24	173.8	174.16	96.32	113.18	0.0	0.08	strobe	dysgu.sv.vcf.sort
bwa2a.sentd.concat.fofn-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.sentd.concat.fofn	0.2625	0:00:00	3.38	7.72	0.13	0.35	0.0	0.0	0.0	0.0	bwa2a	sentd.concat.fofn
sent.deep.22-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.deep.22	210.8791	0:03:30	36526.01	202729.91	23313.83	23593.8	815.39	2845.04	2812.0	35.76	sent	deep.22
bwa2a.oct.23-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.oct.23	314.9373	0:05:14	16450.68	21004.8	16446.77	16453.02	1908.99	32.2	2262.83	7126.91	bwa2a	oct.23
sent.lfq2.bcfstat-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.lfq2.bcfstat	0.4477	0:00:00	13.5	39.31	7.02	9.26	0.0	0.0	0.0	0.02	sent	lfq2.bcfstat
sent.dysgu.sv.vcf-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.dysgu.sv.vcf	1110.1057	0:18:30	2934.15	15401.91	2930.09	2930.5	3836.18	9153.46	98.09	1089.02	sent	dysgu.sv.vcf
bwa2a.mosdepth-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.mosdepth	305.9275	0:05:05	2919.0	34716.26	2927.87	2928.08	474.53	5012.72	173.17	530.17	bwa2a	mosdepth
bwa2a.mrkdup-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.mrkdup	695.1713	0:11:35	385543.84	426724.16	384971.5	384972.17	12832.43	39035.5	12390.8	204.2	bwa2a	mrkdup
sent.deep.23-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	sent.deep.23	374.8576	0:06:14	34028.09	203999.69	23281.43	23554.97	2852.8	3349.42	3510.91	43.29	sent	deep.23
bwa2a.oct.22-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.oct.22	522.8109	0:08:42	9473.68	14363.78	9295.46	9301.56	788.39	7.48	1290.44	6747.17	bwa2a	oct.22
bwa2a.clair3.rtgvcfstats-RIH0_ANA0-HG004-19	RIH0_ANA0-HG004-19	bwa2a.clair3.rtgvcfstats	1.5753	0:00:01	431.89	714301.71	414.06	422.35	0.0	0.09	138.78	2.35	bwa2a	clair3.rtgvcfstats
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ cat results/day/hg38/other_reports/rules_benchmark_data_mqc.tsv  | head
[?2004lcombined_rule	sample	rule	s	h:m:s	max_rss	max_vms	max_uss	max_pss	io_in	io_out	mean_load	cpu_time	rule_prefix	rule_suffix
alignstats_summary-all.	all.	alignstats_summary	0.2405	0:00:00	3.38	7.71	0.12	0.35	0.0	0.0	0.0	0.0	alignstats_summary	NA
final_multiqc-DAY_all.	DAY_all.	final_multiqc	39.5161	0:00:39	544.01	6596.95	444.81	472.93	180.65	17.46	37.96	15.48	final_multiqc	NA
alignstats_smmary_compile-all.	all.	alignstats_smmary_compile	0.7071	0:00:00	11.13	27.76	4.04	4.44	4.48	0.04	0.0	0.0	alignstats_smmary_compile	NA
raw_fastqc-SEQQC-multiqc_.	SEQQC-multiqc_.	raw_fastqc	48.3452	0:00:48	486.65	7625.06	384.82	412.07	2042.91	47.05	85.03	41.23	raw_fastqc	NA
strobe.oct.24-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	strobe.oct.24	883.8152	0:14:43	8466.06	13566.77	8495.96	8497.56	140.0	0.1	473.24	4182.7	strobe	oct.24
sent.lfq2.concordance-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.lfq2.concordance	33.2076	0:00:33	13497.05	716739.74	13518.56	13519.51	770.61	159.3	168.12	56.46	sent	lfq2.concordance
sent.sentd.concat.fofn-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.sentd.concat.fofn	0.2802	0:00:00	3.25	7.72	0.12	0.34	0.0	0.0	0.0	0.0	sent	sentd.concat.fofn
bwa2a.tiddit.sv.vcf-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	bwa2a.tiddit.sv.vcf	913.7027	0:15:13	7827.48	49095.2	7701.54	7724.6	5905.44	7456.35	523.06	1027.58	bwa2a	tiddit.sv.vcf
sent.clair3.bcfstat-RIH0_ANA0-HG005-19	RIH0_ANA0-HG005-19	sent.clair3.bcfstat	0.7995	0:00:00	16.5	2561.94	10.21	11.53	3.08	0.0	0.0	1.9	sent	clair3.bcfstat
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ cat results/day/hg38/other_reports/rules_benchmark_data_mqc.tsv  | head[K[38Pecho $INSTANCE_TYPES_LINE[C[9@source bin/proc_spot_price_logs.sh[Cgit pull[K
[?2004lremote: Enumerating objects: 10, done.[K
remote: Counting objects:  10% (1/10)[Kremote: Counting objects:  20% (2/10)[Kremote: Counting objects:  30% (3/10)[Kremote: Counting objects:  40% (4/10)[Kremote: Counting objects:  50% (5/10)[Kremote: Counting objects:  60% (6/10)[Kremote: Counting objects:  70% (7/10)[Kremote: Counting objects:  80% (8/10)[Kremote: Counting objects:  90% (9/10)[Kremote: Counting objects: 100% (10/10)[Kremote: Counting objects: 100% (10/10), done.[K
remote: Compressing objects:  50% (1/2)[Kremote: Compressing objects: 100% (2/2)[Kremote: Compressing objects: 100% (2/2), done.[K
remote: Total 10 (delta 8), reused 10 (delta 8), pack-reused 0 (from 0)[K
Unpacking objects:  10% (1/10)Unpacking objects:  20% (2/10)Unpacking objects:  30% (3/10)Unpacking objects:  40% (4/10)Unpacking objects:  50% (5/10)Unpacking objects:  60% (6/10)Unpacking objects:  70% (7/10)Unpacking objects:  80% (8/10)Unpacking objects:  90% (9/10)Unpacking objects: 100% (10/10)Unpacking objects: 100% (10/10), 1.25 KiB | 36.00 KiB/s, done.
From github.com:Daylily-Informatics/daylily
   81491ac..34d84cf  main       -> origin/main
Updating 81491ac..34d84cf
Fast-forward
 bin/proc_benchmark_runtime.sh             | 16 [32m++++++++++++++++[m
 config/external_tools/multiqc_header.yaml |  6 [32m++++[m[31m--[m
 workflow/rules/multiqc_final_wgs.smk      |  4 [32m++++[m
 3 files changed, 24 insertions(+), 2 deletions(-)
 create mode 100644 bin/proc_benchmark_runtime.sh
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ source bin/proc_benchmark_runtime.sh 
[?2004lSOURCE ME

Total Runtime (min): 2473.02
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ more bin/proc_benchmark_runtime.sh 
[?2004lecho "SOURCE ME"
echo ""

# Path to the input file
input_file="results/day/hg38/other_reports/rules_benchmark_data_mqc.tsv"

# Extract the 's' column, sum its values, and divide by 60 to get the total runtime in minutes
total_runtime=$(awk -F '\t' 'NR>1 {sum += $4} END {print sum / 60}' "$input_file")
export TOTAL_RUNTIME_MIN=$total_runtime

# Calculate total vCPU minutes (runtime in minutes * 192)
total_vcpu_min=$(awk "BEGIN {print $total_runtime * 192}")
export TOTAL_VCPU_MIN=$total_vcpu_min

# Output the results
echo "Total Runtime (min): $TOTAL_RUNTIME_MIN"
[?2004h(mqc) ubuntu@ip-10-0-0-196:/fsx/analysis_results/ubuntu/hg38/hg005hg007/daylily$ more bin/proc_benchmark_runtime.sh [2@sourc[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[Cgit pull[Ksource bin/proc_benchmark_runtime.sh 