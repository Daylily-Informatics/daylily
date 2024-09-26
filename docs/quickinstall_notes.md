# Overview

Capturing the steps I took to get the daylily framework up and running in a `classic` AWS environment.

# Before Beginning
- `daylily` was tuned to run in the `us-west-2c` AZ, for cost and availabiliyu of some specific compute resources not found in all AZs. The install should all build and proceed using the `us-west-2c` AZ.
- A cloudstack formation template will run which will create some important resources for the epheral cluster, namely: the public VPC, public subnet, private subnet, and a policy to allow tagging and budget tracking of resources by pcluster.
- The install steps should work for both `bash` and `zsh` shells, but the `bash` shell is the default for the `pcluster` environment.

# Steps


## Prepare Your Local Shell
The shell (probably your laptop) you will be using to ssh into the clusters you'll be creating.

### Clone The `daylily` Repository
```bash
git clone git@github.com:Daylily-Informatics/daylily.git
cd daylily
```

### Install Conda (if not already installed)
_bash and zsh are known to work_
- [Conda, miniconda specifically, is required to be setup in your shell to proceed.](docs/install/Install.md#run-daylily-init)
 

## Only Necessary To Follow Once Per AWS Acount
### AWS 
#### 1. Install / Configure The AWS CLI
_Create And Save AWS CLI Credentials In `.aws`_
- For your terminal/shell account these files `~/.aws/{credentials,config}` should exist.

**config** should look like:
```
[default]
region = us-west-2
output = json
```

**credentials** should look like:
```
[default]
aws_access_key_id = <ACCESS_KEY>
aws_secret_access_key = <SECRET_ACCESS_KEY>
```

#### 2. Create A New SSH Key Pair (type: ed25519) & Downloaded `.pem` File
Move, copy and chmod: 
```bash
mv ~/Downloads/<yourkey>.pem  ~/.ssh/<yourkey>.pem 
chmod 400 ~/.ssh/<yourkey>.pem`
```


#### 3. Create `YOURPREFIX-omics-analysis` s3 Bucket
- Your new bucket name should end in `-omics-analysis` and be unique to your account. This allows easier auto-detection latter in the `daycli` setup.
- And it should be created in `us-west-2`.

##### Copy The `daylily` Public References/Resources To Your New Bucket (choose the most recent version)
From your terminal/shell (and to be safe, in a screen or tmux session as it may run for a while).

Move to the root of the cloned repo, please run the following code (which will prompt you to select the version of the reference data you would like to copy to your bucket).

```bash

echo -n "enter your bucket prefix (ie: myorg): "
read your_bucket_prefix

s3_reference_bucket_url="s3://rcrf-omics-analysis/daylily"

# Prompt user to select an S3 reference data version
echo "Select the S3 reference data version from the list below:"

# Check if the version file exists
version_file="config/day_cluster/s3_ref_data.versions"
if [[ ! -f "$version_file" ]]; then
    echo "Error: Version file '$version_file' does not exist. Exiting."
    return 1
fi


# Read versions from the file into an array using a while loop for compatibility
versions=()
while IFS= read -r line; do
    versions+=("$line")
done < "$version_file"


# Check if any versions were read
if [[ ${#versions[@]} -eq 0 ]]; then
    echo "Error: No versions found in '$version_file'. Exiting."
    return 1
fi

# Present the list of versions for selection
select version_choice in "${versions[@]}"; do
    if [[ -n "$version_choice" ]]; then
        echo "You selected version: $version_choice"
        s3_reference_data_version=$version_choice
        break
    else
        echo "Invalid selection, please try again."
    fi
done



cmda="aws s3 cp s3://daylily-references-public/data/$s3_reference_data_version s3://${your_bucket_prefix}-omics-analysis/data --recursive --request-payer requester "
cmdb="aws s3 cp s3://daylily-references-public/cluster_boot_config/$s3_reference_data_version s3://${your_bucket_prefix}-omics-analysis/cluster_boot_config --recursive  --request-payer requester "

echo " Now, run the following commands in your shell:\n\n\t"
echo "  $cmda & $cmdb & echo 'waiting for both to complete' & wait \n\n"

```

The above is composing the commands to copy the current versions of the reference and config files to your local s3 bucket. This is the gist of it (using `v0.9`):
```bash
aws s3 cp s3://daylily-references-public/data/v0.9 s3://<YOURPREFIX>-omics-analysis/data --recursive --request-payer requester 
aws s3 cp s3://daylily-references-public/cluster_boot_config/v0.9 s3://<YOURPREFIX>-omics-analysis/cluster_boot_config 
```



#### 4. Create A New `pcluster` Stack
From the `daylily` repository root:
```bash
bin/init_cloudstackformation.sh ./config/day_cluster/pcluster_env.yml <short-resource-prefix-to-use>
```
- This will run and create a VPC, private and public subnet as well as an IAM policy which allows pcluster to tag resources for budget tracking and allow budget tools to see these tags.
- The script should block the terminal while running, and report back on success.  If failure, go to the cloudformation console and look at the logs for the stack to see what went wrong & delete the stack before attempting again.
- **please record/take note** of the `Public Subnet ID`, `Private Subnet ID`, and `Policy ARN` will be reported here and needed in the `daycli` setup.


## Local DAYCLI Setup & Ephemeral Cluster Initialization (only need to do this once per new cluster)
Your local/laptop shell.
Assumes all of the setup steps above have been completed.

### Local `pcluster` DAYCLI Setup & Ephemeral Cluster Initialization

#### 1. Install The `DAYCLI` Conda Environment (if not present) & Initialize The Ephemeral Cluster
- From the `daylily` repository root:
  - Your AWS cli credentials are automatically detected via the aws cli, and you will be informed of success/fail for this step. You can simply hit enter to proceed through this.
  - You will be prompted to select the correct PEM file.
  - 
  - **full path to your .ssh/<mykey>.pm**,  _so on a Mac `/Users/<you>/.ssh/<mykey>.pem`_

    ```bash
    source bin/daylily-cfg-ephemeral-cluster
    ``` 
    - You will be prompted to:
      - select the full path to your $HOME/.ssh/<mykey>.pem (from detected .pem files)
      - select the `s3` bucket you created and seeded, options presented will be any with names ending in `-omics-analysis`. Or you may select `1` and manually enter a s3 url.
      - select the `Public Subnet ID` created when the cloudstack formation script was run earlier.
      - select the `Private Subnet ID` created when the cloudstack formation script was run earlier.from the cloudformation stack output.
      - select the `Policy ARN` created when the cloudstack formation script was run earlier.
      - enter a name to asisgn your new ephemeral cluster (ie: `<myorg>-omics-analysis`)
      - enter the path to a pcluser config file (ie: `config/day_cluster/pcluster_config.yaml`), enter selects the default.
  - The configuration script will then create  `_cluster.yaml` and `cluster_init_vals.txt` files in `~/.config/daylily/` and these will be used to create the new ephemeral cluster.
  - First, a dryrun `pcluster create-clsuter` command will be run to ensure the configuration is correct. You will be informed of success/fail for this step, and asked if you wish to proceed in actual creation of the cluster, enter `y`.
  - The ephemeral cluster creation will begin and a monitoring script will watch for its completion.  This may take 10-20min. You may exit the process at this time, and can check on the cluster via the pcluster cli, or the cloudformation console.

##### Head Node Configuration (only needed to be done once per new cluster)
_**todo**: bake the headnode as an AMI and use this AMI in the cluster congifuration, making this step unnecessary_
- Headnode configuration will proceed automatically following the above successful cluster creation  (you may manually trigger running `source bin/daylily-init-headnode $pem_file $region` if not started automatically).  
- This process will:
  -  prompt you to select the name of the cluster just created.
  -  will be presented an `rsa` key you will be asked to save as a github ssh/gpg key. When you acknowledge you have saved this key, the `daylily` repo will be cloned to the headnode into `~/projects/daylily`.
  -  The `daylily` cli will be installed on the headnode, which includes creation of a `DAY` conda environment.
  -  Upon completion, you will see something like:
  
```text
You can now SSH into the head node with the following command:
ssh -i /Users/daylily/.ssh/omics-analysis-b.pem ubuntu@52.24.138.65
Once logged in, as the 'ubuntu' user, run the following commands:
  cd ~/projects/daylily
  source dyinit
  source bin/day_activate local
  dy-r help
 
Setup complete. You can now start working with Daylily on the head node.
```
  - You are ready to roll.

### Cluster Initialization Complete ( head node confguration not necessarily complete )
You can confirm the cluster creation was successful with the following command.
```bash
pcluster list-clusters --region us-west-2
```

#### COST OF IDLE CLUSTER
> The cluster, once created, will *MINIMALLY* cost whatever the hourly on-demand cost of the head node you chose (the default instance type should run ~$1/hr). And then, the cost of the fsx filesystem. 

> When the cluster is doing actual work, spot instances are created as needed, and released if idle for longer than (configurable, but default 5min). Spot instance pools are configured to greatly increase finding very cheap spot instances, but this is not guaranteed. You pay per min for each spot instance running. MONITOR YOUR COSTS AND USE.

> The intended use case is not to leave a cluster idle for too long, but spin the thing up and down as needed. Leaving the cluster idle for two+ days or so is probably starting to be too long idle.

> 

### Working With The Ephemeral Clusters (via `pcluster`)
_**note**: AWS quotas limit many of the resoureces used by pcluster, please be aware of these quotas and monitor your usage. You may need to request increases to these quotas to run larger clusters._
#### Activate The DAYCLI Conda Environment
```bash
conda activate DAYCLI
```

#### `pcluster` CLI Usage
```text
pcluster -h  
usage: pcluster [-h]
                {list-clusters,create-cluster,delete-cluster,describe-cluster,update-cluster,describe-compute-fleet,update-compute-fleet,delete-cluster-instances,describe-cluster-instances,list-cluster-log-streams,get-cluster-log-events,get-cluster-stack-events,list-images,build-image,delete-image,describe-image,list-image-log-streams,get-image-log-events,get-image-stack-events,list-official-images,configure,dcv-connect,export-cluster-logs,export-image-logs,ssh,version}
                ...

pcluster is the AWS ParallelCluster CLI and permits launching and management of HPC clusters in the AWS cloud.

options:
  -h, --help            show this help message and exit

COMMANDS:
  {list-clusters,create-cluster,delete-cluster,describe-cluster,update-cluster,describe-compute-fleet,update-compute-fleet,delete-cluster-instances,describe-cluster-instances,list-cluster-log-streams,get-cluster-log-events,get-cluster-stack-events,list-images,build-image,delete-image,describe-image,list-image-log-streams,get-image-log-events,get-image-stack-events,list-official-images,configure,dcv-connect,export-cluster-logs,export-image-logs,ssh,version}
    list-clusters       Retrieve the list of existing clusters.
    create-cluster      Create a managed cluster in a given region.
    delete-cluster      Initiate the deletion of a cluster.
    describe-cluster    Get detailed information about an existing cluster.
    update-cluster      Update a cluster managed in a given region.
    describe-compute-fleet
                        Describe the status of the compute fleet.
    update-compute-fleet
                        Update the status of the cluster compute fleet.
    delete-cluster-instances
                        Initiate the forced termination of all cluster compute nodes. Does not work with AWS Batch clusters.
    describe-cluster-instances
                        Describe the instances belonging to a given cluster.
    list-cluster-log-streams
                        Retrieve the list of log streams associated with a cluster.
    get-cluster-log-events
                        Retrieve the events associated with a log stream.
    get-cluster-stack-events
                        Retrieve the events associated with the stack for a given cluster.
    list-images         Retrieve the list of existing custom images.
    build-image         Create a custom ParallelCluster image in a given region.
    delete-image        Initiate the deletion of the custom ParallelCluster image.
    describe-image      Get detailed information about an existing image.
    list-image-log-streams
                        Retrieve the list of log streams associated with an image.
    get-image-log-events
                        Retrieve the events associated with an image build.
    get-image-stack-events
                        Retrieve the events associated with the stack for a given image build.
    list-official-images
                        List Official ParallelCluster AMIs.
    configure           Start the AWS ParallelCluster configuration.
    dcv-connect         Permits to connect to the head node through an interactive session by using NICE DCV.
    export-cluster-logs
                        Export the logs of the cluster to a local tar.gz archive by passing through an Amazon S3 Bucket.
    export-image-logs   Export the logs of the image builder stack to a local tar.gz archive by passing through an Amazon S3 Bucket.
    ssh                 Connects to the head node instance using SSH.
    version             Displays the version of AWS ParallelCluster.

For command specific flags, please run: "pcluster [command] --help"
```

##### List Clusters
```bash
pcluster list-clusters --region us-west-2
```

##### Describe Cluster

```bash
pcluster describe-cluster -n $cluster_name --region us-west-2
```

ie: to get the public IP of the head node.
```bash
pcluster describe-cluster -n $cluster_name --region us-west-2 | grep 'publicIpAddress' | cut -d '"' -f 4
```

##### SSH Into Cluster Headnode
From your local shell, you can ssh into the head node of the cluster using the following command.
```bash
ssh -i $pem_file ubuntu@$cluster_ip_address 
```

###### First Time Logging Into Head Node
**Confirm `daylily` CLI Is Working**
```bash
cd ~/projects/daylily
. dyinit # inisitalizes the daylily cli
dy-a local # activates the local config
dy-r help # should show the help menu
```
This should produce a magenta `WORKFLOW SUCCESS` message and `RETURN CODE: 0` at the end of the output.  If so, you are set.

**Confirm `/fsx/` directories are present**
```bash
ls -lth /fsx/

total 130K
drwxrwxrwx 3 root root 33K Sep 26 09:22 environments
drwxr-xr-x 5 root root 33K Sep 26 08:58 data
drwxrwxrwx 5 root root 33K Sep 26 08:35 analysis_results
drwxrwxrwx 3 root root 33K Sep 26 08:35 resources
```

**run a local test**
```bash
. dyinit
dy-a local
dy-r produce_deduplicated_bams -p -j 2 -n # dry run
dy-r produce_deduplicated_bams -p -j 2 
```
This will produce a job plan, and then begin executing. The sample manifest can be found in `.test_data/data/0.01x_3_wgs_HG002.samplesheet.csv` (i am aware this is not a `.tsv` :-P ). Runtime on the default small test data runnin locally on the default headnode instance type should be ~5min.
```text

NOTE! NOTE !! NOTE !!! ---- The Undetermined Sample Is Excluded. Set --config keep_undetermined=1 to process it.
Building DAG of jobs...
Creating conda environment workflow/envs/vanilla_v0.1.yaml...
Downloading and installing remote packages.
Environment for /home/ubuntu/projects/daylily/workflow/rules/../envs/vanilla_v0.1.yaml created (location: ../../../../fsx/resources/environments/conda/ubuntu/ip-10-0-0-37/f7b02dfcffb9942845fe3a995dd77dca_)
Creating conda environment workflow/envs/strobe_aligner.yaml...
Downloading and installing remote packages.
Environment for /home/ubuntu/projects/daylily/workflow/rules/../envs/strobe_aligner.yaml created (location: ../../../../fsx/resources/environments/conda/ubuntu/ip-10-0-0-37/a759d60f3b4e735d629d60f903591630_)
Using shell: /usr/bin/bash
Provided cores: 16
Rules claiming more threads will be scaled down.
Provided resources: vcpu=16
Job stats:
job                          count    min threads    max threads
-------------------------  -------  -------------  -------------
doppelmark_dups                  1             16             16
pre_prep_raw_fq                  1              1              1
prep_results_dirs                1              1              1
produce_deduplicated_bams        1              1              1
stage_supporting_data            1              1              1
strobe_align_sort                1             16             16
workflow_staging                 1              1              1
total                            7              1             16
```
This should exit with a magenta success message and `RETURN CODE: 0`. Results can be found in `results/day/{hg38,b37}`.

**run a small slurm test**
The following will submit jobs to the slurm scheduler on the headnode, and spot instances will be spun up to run the jobs (modulo limits imposed by config and quotas).

First, create a working directory on the `/fsx/` filesystem.
```bash
# create a working analysis directory
mkdir -p /fsx/analysis_results/ubuntu/init_test
cd /fsx/analysis_results/ubuntu/init_test
git clone git@github.com:Daylily-Informatics/daylily.git
cd daylily

#  prepare to run the test
tmux new -s slurm_test
. dyinit
dy-a slurm

# create a test manifest for one giab sample only
head -n 2 .test_data/data/giab_30x_hg38_analysis_manifest.csv  > config/analysis_manifest.csv

# run the test
dy-r produce_snv_concordances -p -k -n
```

Which will produce a plan that looks like
```text

Job stats:
job                           count    min threads    max threads
--------------------------  -------  -------------  -------------
deep_concat_fofn                  1              2              2
deep_concat_index_chunks          1              4              4
deepvariant                      24             64             64
doppelmark_dups                   1            192            192
dv_sort_index_chunk_vcf          24              4              4
pre_prep_raw_fq                   1              1              1
prep_deep_chunkdirs               1              1              1
prep_for_concordance_check        1             32             32
prep_results_dirs                 1              1              1
produce_snv_concordances          1              1              1
stage_supporting_data             1              1              1
strobe_align_sort                 1            192            192
workflow_staging                  1              1              1
total                            59              1            192
```

Run the test with
```bash
dy-r produce_snv_concordances -p -k 
```
note: the first time you run a pipeline, if the docker images are not cached, there can be a delay in starting jobs as the docker images are cached. They are only pulled 1x per cluster lifetime, so subsequent runs will be faster.


Once jobs begin to be submitted, you can monitor from another shell on the headnode with:
```bash
squeue

# or 
sinfo

# I find helpful
watch squeue

# also for the headnode
glances
```


##### SSH Into Compute Nodes
You can not access compute nodes directly, but can access them via the head node. From the head node, you can determine if there are running compute nodes with `squeue`, and use the node names to ssh into them.
```bash
ssh <nodename>
```



##### Delete Cluster
**warning**: this will delete all resources created for the ephemeral cluster, importantly, including the fsx filesystem. You must export any analysis results created in `/fsx/analysis_results` from the `fsx` filesystem  back to `s3` before deleting the cluster. 

###### Export `fsx` Analysis Results Back To S3
- Go to the 'fsx' AWS console and select the filesystem for your cluster.
- Under the `Data Repositories` tab, select the `fsx` filesystem and click `Export to S3`. Export can only currently be carried out back to the same s3 which was mounted to the fsx filesystem. 
- Specify the export path as `analysis_results` (or be more specific to an `analysis_results/subdir`), the path you enter is named relative to the mountpoint of the fsx filesystem on the cluster head and compute nodes, which is `/fsx/`. Start the export. This can take 10+ min.  When complete, confirm the data is now visible in the s3 bucket which was exported to. Once you confirm the export was successful, you can delete the cluster (which will delete the fsx filesystem).

###### Delete The Cluster, For Real
_note: this will not modify/delete the s3 bucket mounted to the fsx filesystem, nor will it delete the policyARN, or private/public subnets used to config the ephemeral cluster._
```bash
pcluster delete-cluster-instances -n <cluster-name> --region us-west-2
pcluster delete-cluster -n <cluster-name> --region us-west-2
``` 
- You can monitor the status of the cluster deletion using `pcluster list-clusters --region us-west-2` and/or `pcluster describe-cluster -n <cluster-name> --region us-west-2`. Deletion can take ~10min depending on the complexity of resources created and fsx filesystem size.
