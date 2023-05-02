# Analysis Manifest

- This is the 'sample sheet' required to run daylily.  An example can be found that will run the included test data [here](../../.test_data/data/0.01xwgs_HG002.samplesheet.csv).
- Another example which runs the 7 GIAB google-brain 30x dataset can be found [here](giab_30x_b37_analysis_manifest.csv). NOTE- the full giab dataset should only be run via the slurm profile.
- Copy this to a file named `analysis_manifest.csv` in the `config` top level directory, and it will be picked up when running `dy-r`.
- Optional column headers allow for sub sampling the input fastqs, as well as setting the `-k` bwa flag.
- If specified, truth VCF and BED files will be auto-detected and create concordance analysis reports.
