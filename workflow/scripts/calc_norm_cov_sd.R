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
export SENTIEON_INSTALL_DIR='/fsx/data/cached_envs/sentieon-genomics-202308.03'[26;29H[3m63% L17[0m[39;49m[23m[14;64H[34h[?25h[34l[27;13H[?25ls[C[34h[?25h[34l[14;64H[?25l[36m[45mS[39;49m[34h[?25h[34l[27;14H[?25le	[34h[?25h[34l[14;65H[?25l[36m[45mE[39;49m[34h[?25h[34l[27;15H[?25lr	[34h[?25h[34l[14;66H[?25l[36m[45mR[39;49m[34h[?25h[34l[27;16H[