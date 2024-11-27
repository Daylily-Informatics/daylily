# Running Analysis on an Ephemeral Cluster

The following assumes you have completed the cluster creation steps, and ideally you have also created a PCUI dashboard to monitor and interact with your clusters.

## Step 1: 
### Option A: ssh into cluster headnode
    
```bash
export AWS_PROFILE=your_profile

source bin/ssh_into_daylily --profile your_profile

# you will be prompted to select the region your cluster is running in, to choose the PEM file used to create the cluster (which you must have a copt of in ~/.ssh/), and finally the cluster you want to connect to, and you'll be connected to the headnode of that cluster.

# always work in a tmux or scren session in case your connection is lost
tmux new -s first_analysis

```

### Option B: Use the PCUI to ssh into the cluster

- Visit the PCUI dashboard you created for your cluster, select the region your cluster is running in, and then select the cluster.  Wait a few seconds, and you'll be able to click the 'shell' button to open a ssh terminal in your browser.
- From the new terminal session:
```bash

# become the ubuntu user and activate the shell
sudo su ubuntu
bash

# always work in a tmux or scren session in case your connection is lost
tmux new -s first_analysis

```

## Step 2: Prepare The Analysis Directory
- It is assumed that the headnode was successfully configured when the cluster was created. If you do not have conda available as ubuntu and there is no `~/projects/daylily`, you might need to run the `bin/daylily-cfg-headnode` script from your remote machine first.
- Assuming the headnode config was successfuk... from your tmux session as the ubuntu user on the headnode:

```bash 
cd /fsx/analysis_results/ubuntu
mkdir first_analysis
cd first_analysis

# NOTE!  __yes__ you are generally advised to clone the github repo to run any analysis w/in.  This might seem like overkill, but is a good practice (and has little downside beyond bent senses of propriety). For why, please read the snakemake rolling paper on reproducible and sustainable analysis: https://doi.org/10.12688/f1000research.16740.1 . 

git clone git@github.com:Daylily-Informatics/daylily.git
cd daylily

cp .test/

# initialize the DAY conda env + other helpful tools
. dyinit  # you can speciy --project AWSBUDGETNAME to change from the default region specific budget. New budgets created in AWS must also be added to the allowed ubuntu projects in /opt/slurm/etc/projects_list.conf to be recognized by the cluster. jobs will not launch if the budget is not recognized.

# move an analysis manifest into place, which haas a tiny toy dataset.
cp .test_data/data/0.01xwgs_HG002_hg38.samplesheet.csv confiig/analysis_manifest.csv # YOU WILL USE THE ./bin/create_analysis_manifest script to create your own manifest for your own data when the time comes.

# Activate the local execution profile
dy-a local

## NOTE: TO CHANGE THE REFERENCE BUILD, you must edit the config/day_profiles/local/rule_config.yaml file to point to the desired reference build.  This is done on line 4 of the yaml.  The default is hg38. If using the slurm profile, this will need to be changed as well.

# run a toy local analysis to make sure everything is working
dy-r seqqc -j 1  -p -k -n  # -n runs in dryrun mode, -p prints the commands, -k will keep going where possible if there ar errors  -j 1 limits one job at a time

# Run the analysis for real
dy-r seqqc -j 1  -p -k

# All results will be located here following any analysis run
find results/ 

# Next, run the same analysis on the cluster, activate the slurm profile
dy-a slurm

# clean up the local results
rm -rf results .snakemake logs

y-r seqqc -j 3  -p -k -n  # -n runs in dryrun mode, -p prints the commands, -k will keep going where possible if there ar errors  -j 1 limits one job at a time

# Run the analysis for real
dy-r seqqc -j 3  -p -k

# Jobs will limited to 3 at a time.
# No compute nodes should be running, so it may take 5-15+min to get spot instnaes running and available for your jobs.  This is normal, and varies depending on the spot market.

## Your tmux terminal will be blocked as things log to stdout.  

# USE PCUI to monitor the cluster jobs -- or -- move to another tmux pane or terminal on the headnode

sinfo # reports on the overall compute fleet status

squeue # reports on jobs running and pending

watch $sq_cmd # fancy formatted squeue with all the infoz presented

# ... ultimately things will succeed or fail in the original window.

# Logs can be found in ./logs/slurm/*

# The seqqc results are a multiqc report, which can be found

ls results/day/{hg38,b37}/reports/*.html

```

## Step 2.5: (optional) Run a test analysis with GIAB data
- Each GIAB sample has 30x  2x150 ILMN google brain data available for use be default.  There are pre-configured analysis manifests for these data sets, which also include the info to automatically create concordance reports with the GIAB truth sets (both for `b37` and `hg38`).  To run `bwa mem2 + doppelmark + deepvariant` on the GIAB samples vs `hg38`, you can do the following.
- From a tmux session on the headnode, in the github repo directory prepared above.

### hg38
Copy the pre-configured giab analysis manifest for hg38, and run the analysis. This manifest points to the fastq.gz files for each, as well as the truth data set needed to produce the concordance reports.
```bash

. dyinit 
dy-a slurm

vi config/day_profiles/slurm/rule_config.yaml # the default build is hg38, so no changes are needed here 

cp .test_data/data/giab_30x_hg38_analysis_manifest.csv config/analysis_manifest.csv

dy-r produce_snv_concordances -p  -k -j 9 -n

dy-r produce_snv_concordances  -p -k -j 9

# The concordance reports will be found in results/day/hg38/other_reports/
ls results/day/hg38/other_reports/*giab*

```

### b37
You will need to modify the rule_config.yaml to use b37, and copy the pre-configured giab analysis manifest for b37.
```bash

. dyinit 
dy-a slurm

# you will need to chang ethe reference build to b37 from hg38
vi config/day_profiles/slurm/rule_config.yaml 

cp .test_data/data/giab_30x_b37_analysis_manifest.csv config/analysis_manifest.csv

dy-r produce_snv_concordances -p  -k -j 9 -n

dy-r produce_snv_concordances  -p -k -j 9

# The concordance reports will be found in results/day/b37/other_reports/
ls results/day/b37/other_reports/*giab*

```



## Step 3: Prepare Your Own Data // Create an Analysis Manifest
_the analysis_manifest.csv_ is being fully re-worked ASAP, so please bear with the current process for the moment.


### FIRST: stage your input fastqs 
You must have your fastq.gz files stored in the S3 bucket which you used to attach to the /fsx filesystem mounted to this headnode (and will be mounted to all compute nodes).
- FASTQ files must be located in this dir, or a sub-dir `s3://<YOURPREFIX>-omics-analysis-<REGION>/data/tmp_input_sample_data/`, and the files should be visible from the headnode/compute nodes with `ls /fsx/data/tmp_input_sample_data/*`.
  - You will be referencing the input data via the `/fsx/data/tmp_input_sample_data/` path.
    - Fastqs should be deleted from this location once analysis is complete.  I advise keeping the BAM/CRAM as your fastq source and ditching the fastqs entirely. If you wish to also keep the fastq, move them someplace else once you are done with the analysis-- this is not meant to be fastq storage.

### SECOND: create a fastq manifest
Create a `fastq_manifest.tsv` with the following columns (all are required):
```text
sample_id   fastq1_fsx_path  fastq2_fsx_path  expected_copies_x  expected_copies_y
HG001-10x   /fsx/data/tmp_input_sample_data/HG001-10x_R1.fastq.gz /fsx/data/tmp_input_sample_data/HG001-10x_R2.fastq.gz  2 0
HG002-10x   /fsx/data/tmp_input_sample_data/HG002-10x_R1.fastq.gz /fsx/data/tmp_input_sample_data/HG002-10x_R2.fastq.gz  1 1
HG001-30x   /fsx/data/tmp_input_sample_data/HG001-30x_R1.fastq.gz /fsx/data/tmp_input_sample_data/HG001-30x_R2.fastq.gz  2 0
```
- sample_id should only contain alphanumeric characters and dashes. No underscores, spaces or other special characters. This sample_id should be succinct, as it will be included in the output file names (file names are also in line to be re-worked along with the analysis_manifest).
- expected_copies_x and expected_copies_y are used to assess if expected X and Y ploidy is observerd in qc reports (if requested).

### THIRD: create an analysis manifest from your fastq manifest
_among the first changes coming to this part of the process, moving from a `csv` to `tsv`_
```bash
bash bin/create_analysis_manifest -f fastq_manifest.tsv -o config/analysis_manifest.csv
```

- This will create a correctly formatted `config/analysis_manifest.csv` file for your analysis to use. The `analysis_manifest.csv` file has a lot going on, and much of it is deprecated.  I do not advise tinkering with it unless very eager to do so.
- Every time you run `dy-r COMMAND ARGS`, this manifest will be used to determine which samples to run. You can re-run after making changes to this manifest, or other parts of the code, but please take the time to read the snakemake docs and understand which flags will be needed to achieve the desired behavior (it is very easy to blow away work if you are not careful... get in the habit of always using `-n`).

## Step 4: Run Your Analysis (bwa mem2, doppelmark, deepvariant)!
- From the github repo directory on the headnode, in a tmux terminal, with your desired `config/analysis_manifest.csv` in place:
```bash

. dyinit
dy-a slurm

# If you wish to enable an aligner other than bwa mem2 or a variant caller other than deep variant, you will need to edit the config/day_profiles/slurm/rule_config.yaml before proceeding.  

# See what the execution plan looks like
dy-r produce_snv_concordances produce_manta produce_tiddit-p -k -j 9 -n  # NOTE::: your spot quotas will probably not allow running more than 2-3 192vcpu instances at a time. You can send more jobs than your quotas allow, but watch your quotas and the cloudwatch dashboard if you are not seeing as many instances spin up as you expect.


# Run the analysis for real
dy-r produce_snv_concordances produce_manta produce_tiddit-p -k -j 9

# Monitor the jobs with squeue, watch, or the PCUI dashboard

# results can be found in ./results/ and logs in ./logs/slurm/*/*

find results/ | grep -E 't.vcf.gz|bam.bai'

```


## Step 5: (bonus) Build The MultiQC Report
- From the github repo directory on the headnode, in a tmux terminal, with your desired `config/analysis_manifest.csv` in place:
```bash

. dyinit
dy-a slurm

dy-r produce_multiqc_final_wgs -p -k -j 9 -n

dy-r produce_multiqc_final_wgs -p -k -j 9

# When complete, generate the benchmark metrics file

bash bin/util/benchmarks/collect_benchmarks.sh hg38 # or b37 if you are using b37, the script will indicate where to find the output file

# You'll find the multiqc report here

ls results/day/{hg38,b37}/reports/*.html

```
