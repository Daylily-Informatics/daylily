# Daylily AWS Ephemeral Cluster Setup (0.7.131d)


**beta release**

Daylily is a framework for setting up ephemeral AWS clusters optimized for genomics data analysis. It leverages AWS ParallelCluster and provides automated scripts for cluster creation, management, and teardown.


<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/0000002.png" valign="bottom" ></a></p>


## Table of Contents
- [Before Beginning](#before-beginning)
- [Strongly Suggested: AWS Parallel Cluster UI](#strongly-suggested--pcui)
- [Installation Steps](#installation-steps)
  - [Prepare Your Local Shell](#prepare-your-local-shell)
  - [Clone The `daylily` Repository](#clone-the-daylily-repository)
  - [Install Conda](#install-conda-if-not-already-installed)
- [AWS Setup](#aws-setup)
- [Cluster Initialization](#cluster-initialization)
- [Working with the Ephemeral Cluster](#working-with-the-ephemeral-cluster)
- [Monitoring Tools](#monitoring-tools)
- [Known Issues](#known-issues)
- [Future Development](#future-development)

---

# What's It All About?

## BFAIR: Bioinformatics FAIR Principles 
_this is a title rough idea! - not likely final, but is the gist of it_

## Comprehensive Cost Transparency & Predictability

> 30x FASTQ->VCF, in 1hr, for <img src="docs/images/000000.png" valign="bottom" width=2 >~$3.00<img src="docs/images/000000.png" valign="bottom" width=2 > @ F-score 0.99xx

> Be up and running in a few hours from reading this sentence. 
> - Have your first GIAB 30x Fastq->VCF (both snv and sv) ~60min later. 
> - The (_onetime_) cost of staging data is ~$20, analysis will be ` ~$3.00 to $5.00 ` (pricing is established dynamically at cluster creation, and you can inspect the max bound on spot prices which are possible, this sets your upper bound... as does complexity of pipeline, but more on that latter).

**Time to result, Cost of analysis, Accuracy && Integrated Concordance/Comparison**: These are key elements required in making solid analysis decisions, and in setting the stage for analysis decisions which can improve in tandem as the field advances. 

**Cost Optimization & Predictability**:  With benchmarked, reproducible analysis on stable and reproducible computing platforms, it is possible to optimze for the most beneficial compute locale && to predict expected (and bound highest) per-sample cost to analyze.

**Cost Transparency & Management**:  IRT views into what is being spent, and where.  Not just compute, but data transfer, sroage, and other _ALL_ other costs. No specialized hardware investment needed, no contracts, pay for what you use. 

**Pro Open Source**: All out of the box functionality here is open source, and does not require any investment in software liscneces, etc. This is important for both future proof reproducibility and ongoing cost management. This said, daylily is not hostile to s/w that requires liscences (if selection of closed s/w is made understanding tradeoffs, if any, in long term reproducibility), but you will need to purchase those separately.

- https://github.com/aws/aws-parallelcluster makes the cluster creation/management possible.
- snakemake
- all the tools!

<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff00ff><img src="docs/images/000000.png" valign="bottom" ></a></p>

<p valign="middle"><img src="docs/images/000000.png" valign="bottom" ></p>

# Installation -- PREREQUISITES

## AWS 

### Quotas
There are a handful of quotas which will greatly limit (or block) your ability to create and run an ephemeral cluster.  These quotas are set by AWS and you must request increases. The `daylily-cfg=ephemeral-cluster` script will check these quotas for you, and warn if it appears they are too low,  but you should be aware of them and [request increases proactively // these requests have no cost](https://console.aws.amazon.com/servicequotas/home).

**dedicated instances**
_pre region quotas_
- `Running Dedicated r7i Hosts` >= 1 **!!(AWS default is 0) !!**
- `Running On-Demand Standard (A, C, D, H, I, M, R, T, Z) instances` must be >= 9 **!!(AWS default is 5) !!** just to run the headnode, and will need to be increased further for ANY other dedicated instances you (presently)/(will) run.

**spot instances**
_per region quotas_
- `All Standard (A, C, D, H, I, M, R, T, Z) Spot Instance Requests` must be >= 310 (and preferable >=2958) **!!(AWS default is 5) !!**

**fsx lustre**
_per region quotas_
- should minimally allow creation of a FSX Lustre filesystem with >= 4.8 TB storage, which should be the default.

**other quotas**
May limit you as well, keep an eye on `VPC` & `networking` specifically.


### Activate Cost Allocation Tags (optional, but strongly suggested)
The cost management built into daylily requires use of budgets and cost allocation tags. Someone with permissions to do so will need to activate these tags in the billing console. *note: if no clusters have yet been created, these tags may not exist to be activeted until the first cluster is running. Activating these tags can happen at any time, it will not block progress on the installation of daylily if this is skipped for now*. See [AWS cost allocation tags](https://us-east-1.console.aws.amazon.com/billing/home#/tags)

The tags to activate are:
  ```text
  aws-parallelcluster-jobid
  aws-parallelcluster-username
  aws-parallelcluster-project
  aws-parallelcluster-clustername
  aws-parallelcluster-enforce-budget
  ```


### CLI User Account
#### CLI Credentials
- Whoever will be creating and managing the daylily ephemeral clusters will need to have aws cli credentials saved in their ~/.aws directory. Please ensure these are saved before proceeding.

#### Inline Policies
AWS Parellel Cluster requires a number of permissions and policies be attached to the user who will be creating and managing the clusters. 
[Their documentation is here](https://docs.aws.amazon.com/parallelcluster/v2/ug/iam.html). 
[I compiled the inline policies and permissions necessary to proceed here](config/aws/daylily-omics-analysis.json), this daylily policy `json` includes a few other policies needed for advanced cost management and to build a PCUI.

- As an admin, you can apply these to the user with the following:

```bash
# Define variables
AWS_USERNAME="<your_username>" # Replace with your actual AWS IAM username
AWS_ACCOUNT_ID="<your_account_id>" # Replace with your actual AWS account number

# Copy the template JSON file to a temporary file
cp config/aws/daylily-omics-analysis.json ./daylily-omics-analysis.tmp.json

# Replace placeholders in the temporary file
perl -pi -e "s/AWS_USERNAME/$AWS_USERNAME/g" ./daylily-omics-analysis.tmp.json
perl -pi -e "s/AWS_ACCOUNT_ID/$AWS_ACCOUNT_ID/g" ./daylily-omics-analysis.tmp.json

# Apply the IAM policy
aws iam put-user-policy --user-name "$AWS_USERNAME" --policy-name daylily-omics-analysis --policy-document file://./daylily-omics-analysis.tmp.json

# Clean up the temporary file
rm -f ./daylily-omics-analysis.tmp.json
```
  - Or, copy the policy from `./daylily-omics-analysis.tmp.json` and apply it to the user via the AWS IAM console.

##### A Note On Budgets
- The access necesary to *view* budgets is beyond the scope of this config, please work with your admin to set that up. If you are able to create clusters and whatnot, then the budgeting infrastructure should be working.

#### SSH Key Pair(s)
You will need one keypair per region you intend to run in.
- via the EC2 Networking Dashboard
- The key name must include `-omics-analysis` in the name.
- The key must be of type: **ed25519**)
-  Downloaded `.pem` locally & be sure to have a copy on the local machine you intend to manage things from, you can't download it again!

##### Place .pem File & Set Permissions

```bash
chmod 700 ~/.ssh

mv ~/Downloads/<yourkey>.pem  ~/.ssh/<yourkey>.pem 
chmod 400 ~/.ssh/<yourkey>.pem`
```

---


## Prerequisites (On Your Local Machine) 
Local machine development has been exclusively on a mac, using the `zsh` shell (I am moving everything to bash slowly, so there is a mix of bash and zsh at the moment). _suggestion: run things in tmux or screen_

Very good odds this will work on any mac and most Linux distros (ubuntu 22.04 are what the cluster nodes run). Windows, I would not hazard a guess.

### System Packages
Install with `brew` or `apt-get` or ... :
- `python3`
- `git`
- `jq`
- `wget`
- `tmux` (optional, but suggested)
- `emacs` (optional, I guess, but I'm not sure how you live without it)



### AWS CLI
#### Install
Use `brew install awscli` or `sudo apt-get install awscli` or see [Install the AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/install-cliv2.html).

#### AWS CLI Configure - Opt 1
Run `aws configure` and enter your credentials and region.

#### AWS CLI Configure - Opt 2
Create the files and directories manually:

```bash
mkdir ~/.aws
chmod 700 ~/.aws

touch ~/.aws/credentials
chmod 600 ~/.aws/credentials

touch ~/.aws/config
chmod 600 ~/.aws/config
```

Edit `~/.aws/config`, which should look like:

```yaml
[default]
region = us-west-2
output = json
```

Edit `~/.aws/credentials`, and add your deets, which should look like:
```yaml
[default]
aws_access_key_id = <ACCESS_KEY>
aws_secret_access_key = <SECRET_ACCESS_KEY>
region = <REGION>
```

### `daylily` Repository
Clone the `daylily` repository to your local machine.

```bash
git clone https://github.com/Daylily-Informatics/daylily.git  # or, if you have set ssh keys with github and intend to make changes:  git clone git@github.com:Daylily-Informatics/daylily.git
cd daylily
```


### Miniconda
Install with:
```bash
source bin/install_miniconda.sh
```

### DAYCLI Environment
```bash
source bin/init_daycli 

# DAYCLI should now be active... did it work?
colr  'did it work?' 0,100,255 255,100,0

```


## [daylily-references-public](#daylily-references-public-bucket-contents) Reference Bucket(s)
- The `daylily-references-public` bucket is preconfigured with all the necessary reference data to run the various pipelines, as well as including GIAB reads for automated concordance. 
- This bucket will need to be cloned to a new bucket with the name `<YOURPREFIX>-omics-analysis-<REGION>`, one for each region you intend to run in.
- These S3 buckets are tightly coupled to the `Fsx lustre` filesystems (which allows 1000s of concurrnet spot instances to read/write to the shared filesystem, making reference and data management both easier and highly performant). 
- [Continue for more on this topic,,,](#s3-reference-bucket--fsx-filesystem).
- This will cost you ~$23 to clone w/in `us-west-2`, up to $110 across regions. _(one time, per region, cost)_ 
- The bust will cost ~$14.50/mo to keep hot in `us-west-2`. It is not advised, but you may opt to remove unused reference data to reduce the monthly cost footprint by up to 65%. _(monthly ongoing cost)_

### Clone `daylily-references-public` to YOURPREFIX-omics-analysis-REGION
- `YOURPREFIX` will be used as the bucket name prefix. Please keep it short. The new bucket name will be `YOURPREFIX-omics-analysis-REGION` and created in the region you specify. You may name the buckets in other ways, but this will block you from using the `daylily-cfg-ephemeral-cluster` script, which is largely why you're here.
- Cloning it will take 1 to many hours.
  
**Use the following script**
_from a tmux or screen session, as the copy will take 1 to many hours_

```bash
conda activate DAYCLI
# help
bash ./bin/create_daylily_omics_analysis_s3.sh -h

AWS_PROFILE=<your_profile>
BUCKET_PREFIX=<your_prefix>

# dryrun
bash ./bin/create_daylily_omics_analysis_s3.sh  --disable-warn --region us-west-2 --profile daylily --bucket-prefix $BUCKET_PREFIX

# run for real
bash ./bin/create_daylily_omics_analysis_s3.sh  --disable-warn --region us-west-2 --profile daylily --bucket-prefix $BUCKET_PREFIX --disable-dryrun

```

---

# Installation -- Daylily Ephemeral Cluster Creation 
You'll need to have an S3 bucket and a keypair available in any region you choose. You may run concurrently in multiple AZs w/in the same region, however, you are likely to hit quota restrictions if you are not proactive in requesting increases. For any cluster crashes (after creating one successfully, see the cloud monitoring dashboard in the region for the clusters, or stack formation logs- these often have clues re: quota issues).

## Select The Most Cost Effective AZ (based on current pricing, not guaruntees... though, daylily does bound the max cost possible)
```text
| Zone       | Total Instances Considered | Instances with Pricing | Median Spot Price | Mean Spot Price | Min Spot Price | Max Spot Price | Avg of Lowest 3 Prices | Harmonic Mean Price | Estimated EC2 Cost / Genome (no FSX or S3 costs included here) | Stability (Max-Min Spread) |
|------------|----------------------------|------------------------|-------------------|----------------|---------------|---------------|------------------------|--------------------|-----------------------------|----------------------------|
| us-east-1d | 6                          | 6                      | **2.2443**         | **2.6635**      | **1.2944**    | 6.0426        | **1.5202**              | **2.0483**         | **4.5605**                  | 4.7482                     |
| us-east-1b | 6                          | 6                      | 2.5682             | 3.1314         | 1.8549         | 6.8788        | 1.9998                 | 2.5703              | 5.9994                      | 5.0239                     |
| us-west-2c | 6                          | 6                      | 3.2456             | 3.4618         | 1.3648         | **5.6738**    | 2.2753                 | 2.8179              | 6.8259                      | **4.3090**                 |
| us-west-2d | 6                          | 6                      | 3.0717             | 3.3662         | 1.7559         | 5.8005        | 2.2841                 | 2.8928              | 6.8523                      | 4.0446                     |
| us-east-1c | 6                          | 6                      | 3.9424             | 3.8837         | 2.1061         | 5.4064        | 3.1122                 | 3.5621              | 9.3365                      | **3.3003**                 |
| us-west-2b | 6                          | 6                      | 3.9970             | 4.4728         | 2.3100         | 9.1402        | 2.8915                 | 3.6808              | 8.6745                      | 6.8302                     |
| us-east-1a | 6                          | 6                      | 3.5854             | 4.0252         | 2.9681         | 6.8254        | 3.2490                 | 3.7470              | 9.7470                      | 3.8573                     |
| us-west-2a | 6                          | 6                      | 3.9594             | 4.8392         | 3.4149         | 8.8822        | 3.4376                 | 4.2998              | 10.3128                     | 5.4673                     |
```


## Create An Ephemeral Cluster In The AZ You Select
This script is doing quite a lot. Please check it out to understand it.  It will be confirming that as many pre-requisites that have been anicipated are met, and creating some resources where needed (VPCs, subnets, budgets).

```bash
AWS_PROFILE=<your_profile>
region_az=<region-az>
source bin/daylily-cfg-ephemeral-cluster --region-az $region_az --profile <your_profile>

```    
- Your aws credentials will be auto-detected and used to query appropriate resources to select from to proceed. You will be prompted to:
      - select the full path to your $HOME/.ssh/<mykey>.pem (from detected .pem files)
      - select the `s3` bucket you created and seeded, options presented will be any with names ending in `-omics-analysis`. Or you may select `1` and manually enter a s3 url.
      - select the `Public Subnet ID` created when the cloudstack formation script was run earlier.
      - select the `Private Subnet ID` created when the cloudstack formation script was run earlier.from the cloudformation stack output.
      - select the `Policy ARN` created when the cloudstack formation script was run earlier.
      - enter a name to asisgn your new ephemeral cluster (ie: `<myorg>-omics-analysis`)
      - enter the path to a pcluser config file (ie: `config/day_cluster/pcluster_config.yaml`), enter selects the default.
  
  - A script will be run to check the current spot pricing for all instance types defined in the clusder config yaml which was composed for you, and will put `max-spot-price` values for each group of instances, which limits the max price you will pay for spot instances. This is a hard limit, and if the spot price exceeds this, the instance will be terminated. This is a safety measure to prevent runaway costs.
  - The configuration script will then create  `CLUSTERNAME_cluster.yaml` and `CLUSTERNAME_cluster_init_vals.txt` files in `~/.config/daylily/` and these will be used to create the new ephemeral cluster.
  - First, a dryrun `pcluster create-cluster` command will be run to ensure the configuration is correct. You will be informed of success/fail for this step, and asked if you wish to proceed in actual creation of the cluster, enter `y`.
  - The ephemeral cluster creation will begin and a monitoring script will watch for its completion. **this can take from 20m to an hour to complete**, depending on the region, size of Fsx requested, S3 size, etc.  There is a max timeout set in the cluster config yaml of 1hr, which will cause a failure if the cluster is not up in that time. 

The terminal will block, and if all is well, you will see the following message:

```text
You can now SSH into the head node with the following command:
ssh -i /Users/daylily/.ssh/omics-analysis-b.pem ubuntu@52.24.138.65
Once logged in, as the 'ubuntu' user, run the following commands:
  cd ~/projects/daylily
  source dyinit -h
  source dyinit  --project PROJECT
  dy-a local
  dy-r help
 
Setup complete. You can now start working with Daylily on the head node.
```
  - You are ready to roll.


### Review Clusters
You can confirm the cluster creation was successful with the following command (alternatively, use the PCUI console).

```bash
pcluster list-clusters --region $REGION
```

### Confirm The Headnode Is Configured
[See the instructions here](#first-time-logging-into-head-node) to confirm the headnode is configured and ready to run the daylily pipeline.


---

# Costs
## Monitoring (tags and budgets)
- Every resource created by daylily is tagged to allow in real time monitoring of costs, to whatever level of granularity you desire. This is intended as a tool for not only managing costs, but as a very important metric to track in assessing various tools utility moving ahead (are the costs of a tool worth the value of the data produced by it, and how does this tool compare with others in the same class?)


## Regulating Usage via Budgets
- During setup of each ephemeral cluster, each cluster can be configured to enforce budgets. Meaning, job submission will be blocked if the budget specifiecd has been exceeded.
  

## OF HOT & IDLE CLUSTER ( ~$1.68 / hr )
For default configuration, once running, the hourly cost will be ~ **$1.68** (_note:_ the cluster is intended to be up only when in use, not kept hot and inactive).
The cost drivers are:
1. `r7i.2xlarge` on-demand headnode = `$0.57 / hr`.
2. `fsx` filesystem = `$1.11 / hr` (for 4.8TB, which is the default size for daylily. You do not pay by usage, but by size requested). 
3.  No other EC2 or storage (beyond the s3 storage used for the ref bucket and your sample data storage) costs are incurred.

## OF RUNNING CLUSTER ( >= $2.50 / hr )
There is the idle hourly costs, plus...

### Spot instances ( ~$1.20 / hr per 192vcpu instance )
For v192 spots, the cost is generally $1 to $3 per hour _(if you are discriminating in your AZ selection, the cost should be closer to $1/hr)_.
- You pay for spot instances as they spin up and until they are shut down (which all happens automatically). The max spot price per resource group limits the max costs (as does the max number of instances allowed per group, and your quotas).

### Data transfer, during analysis ( ~$0.00 )
There are no anticipated or observed costs in runnin the default daylily pipeline, as all data is already located in the same region as the cluster. The cost to reflect data from Fsx back to S3 is effectively $0.

### Data transfer, staging and moving off cluster ( ~$0.00 to > $0.00/hr )
Depending on your data management strategy, these costs will be zero or more than zero.
- You can use  Fsx to bounce results back to the mounted S3 bucket, then move the results elsewhre, or move them from the cluster to another bucket (the former I know has no impact on performance, the latter might interfere with network latency?).

### Storage, during analysis ( ~$0.00 )
- You are paying for the fsx filesystem, which are represented in the idle cluster hourly cost. There are no costs beyond this for storage.
- *HOWEVER*, you are responsible for sizing the Fsx filesystem for your workload, and to be on top of moving data off of it as analysis completes. Fsx is not intended as a long term storage location, but is very much a high performance scratch space.


## OF DELETED CLUSTER -- compute and Fsx ( ~$0.00 )
If its not being used, there is no cost incurred.

## OF REFERENCE DATA in S3 ( $14.50 / month )
- The reference bucket will cost ~$14.50/mo to keep available in `us-west-2`, and one will be needed in any AZ you intend to run in. 
- You should not store your sample or analysis data here long term.

## OF SAMPLE / READ DATA in S3 ( $0.00 to $A LOT / month )
- I argue that it is unecessary to store `fastq` files once bams are (properly) created, as the bam can reconstitute the fastq. So, the cost of storing fastqs beyond initial analysis, should be `$0.00`.

## OF RESULTS DATA in S3 ( $Varies, are you storing BAM or CRAM, vcf.gz or gvcf.gz? )
I suggest:
- CRAM (which is ~1/3 the size of BAM, costs correspondingly less in storage and transfer costs).
- gvcf.gz, which are bigger than vcf.gz, but contain more information, and are more useful for future analysis. _note, bcf and vcf.gz sizes are effectively the same and do not justify the overhad of managing the additional format IMO._





---

# PCUI (technically optional, but you will be missing out)
[Install instructions here](https://docs.aws.amazon.com/parallelcluster/latest/ug/install-pcui-v3.html#install-pcui-steps-v3), launch it using the public subnet created in your cluster, and the vpcID this public net belongs to. These go in the `ImageBuilderVpcId` and `ImageBuilderSubnetId` respectively.

- You should be sure to enable SSM which allows remote access to the nodes from the PCUI console. https://docs.aws.amazon.com/systems-manager/latest/userguide/session-manager-getting-started-ssm-user-permissions.html


## PCUI Costs ( ~ $1.00 / month )
- *[< $1/month>](https://docs.aws.amazon.com/parallelcluster/latest/ug/install-pcui-costs-v3.html)*


---

# Working With The Ephemeral Clusters

## PCUI
- Visit your url created when you built a PCUI

## DAYCLI & AWS Parallel Cluster CLI (pcluster)

### Activate The DAYCLI Conda Environment
```bash
conda activate DAYCLI
```

### `pcluster` CLI Usage
**WARNING:**  you are advised to run `aws configure set region <REGION>` to set the region for use with the pcluster CLI, to avoid the errors you will cause when the `--region` flag is omitted.

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

#### List Clusters
```bash
pcluster list-clusters --region us-west-2
```

#### Describe Cluster

```bash
pcluster describe-cluster -n $cluster_name --region us-west-2
```

ie: to get the public IP of the head node.
```bash
pcluster describe-cluster -n $cluster_name --region us-west-2 | grep 'publicIpAddress' | cut -d '"' -f 4
```

#### SSH Into Cluster Headnode
##### Basic
From your local shell, you can ssh into the head node of the cluster using the following command.
```bash
ssh -i $pem_file ubuntu@$cluster_ip_address 
```
##### Facilitated
```bash
AWS_PROFILE=<profile_name>
bin/daylily-ssh-into-headnode 
```

##### Confirm Headnode Configuration Is Complete

**Is `daylily` CLI Available & Working**
```bash
cd ~/projects/daylily
. dyinit # inisitalizes the daylily cli
dy-a local # activates the local config
dy-r help # should show the help menu
```

This should produce a magenta `WORKFLOW SUCCESS` message and `RETURN CODE: 0` at the end of the output.  If so, you are set. If not, see the next section.

###### Headnode Confiugration Incomplete

If there is no `~/projects/daylily` directory, or the `dyinit` command is not found, the headnode configuration is incomplete. 

**Attempt To Complete Headnode Configuration**
From your remote terminal that you created the cluster with, run the following commands to complete the headnode configuration.
```bash
conda activate DAYCLI

source ./bin/daylily-cfg-headnode $PATH_TO_PEM $CLUSTER_AWS_REGION $AWS_PROFILE
```

> If the problem persists, ssh into the headnode, and attempt to run the commands as the ubuntu user which are being attempted by the `daylily-cfg-headnode` script.

##### Confirm Headnode /fsx/ Directory Structure

**Confirm `/fsx/` directories are present**
```bash
ls -lth /fsx/

total 130K
drwxrwxrwx 3 root root 33K Sep 26 09:22 environments
drwxr-xr-x 5 root root 33K Sep 26 08:58 data
drwxrwxrwx 5 root root 33K Sep 26 08:35 analysis_results
drwxrwxrwx 3 root root 33K Sep 26 08:35 resources
```

##### Run A Local Test Workflow

```bash
. dyinit  --project PROJECT
dy-a local

head -n 2 .test_data/data/giab_30x_hg38_analysis_manifest.csv

dy-r produce_deduplicated_bams -p -j 2 -n # dry run
dy-r produce_deduplicated_bams -p -j 2 
```

###### More On The `-j` Flag
The `-j` flag specified in `dy-r` limits the number of jobs submitted to slurm. For out of the box settings, the advised range for `-j` is 1 to 10. You may omit this flag, and allow submitting all potnetial jobs to slurm, which slurm, /fsx, and the instances can handle growing to the 10s or even 100 thousands of instances... however, various quotas will begin causing problems before then.  The `local` defauly is set to `-j 1` and `slurm` is set to `-j 10`, `-j` may be set to any int > 0.

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


##### Run A Slurm Test Workflow
The following will submit jobs to the slurm scheduler on the headnode, and spot instances will be spun up to run the jobs (modulo limits imposed by config and quotas).

First, create a working directory on the `/fsx/` filesystem.
```bash
# create a working analysis directory
mkdir -p /fsx/analysis_results/ubuntu/init_test
cd /fsx/analysis_results/ubuntu/init_test
git clone https://github.com/Daylily-Informatics/daylily.git  # or, if you have set ssh keys with github and intend to make changes:  git clone git@github.com:Daylily-Informatics/daylily.git
cd daylily

#  prepare to run the test
tmux new -s slurm_test
. dyinit 
dy-a slurm

# create a test manifest for one giab sample only, which will run on the 0.01x test dataset
head -n 2 .test_data/data/0.01xwgs_HG002_hg38.samplesheet.csv > config/analysis_manifest.csv

# run the test, which will auto detect the analysis_manifest.csv file & will run this all via slurm
dy-r produce_snv_concordances -p -k -j 3 -n
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
dy-r produce_snv_concordances -p -k -j 6  #  -j 6 will run 6 jobs in parallel max, which is done here b/c the test data runs so quickly we do not need to spin up one spor instance per deepvariant job & since 3 dv jobs can run on a 192 instance, this flag will limit creating only  2 instances at a time.
```

note1: the first time you run a pipeline, if the docker images are not cached, there can be a delay in starting jobs as the docker images are cached. They are only pulled 1x per cluster lifetime, so subsequent runs will be faster.

note2: The first time a cold cluster requests spot instances, can take some time (~10min) to begin winning spot bids and running jobs. Hang tighe, and see below for monitoring tips.

#### (RUN ON A FULL 30x WGS DATA SET)
**ALERT** The `analysis_manifest.csv` is being re-worked to be more user friendly. The following will continue to work, but will be repleaced with a less touchy method soon.

##### Specify A Single Sample Manifest
You may repeat the above, and use the pre-existing analysis_manifest.csv template `.test_data/data/giab_30x_hg38_analysis_manifest.csv`.

```bash
tmux new -s slurm_test_30x_single

mkdir /fsx/analysis_results/ubuntu/slurm_single
cd /fsx/analysis_results/ubuntu/slurm_single

git clone https://github.com/Daylily-Informatics/daylily.git  # or, if you have set ssh keys with github and intend to make changes:  git clone git@github.com:Daylily-Informatics/daylily.git
cd daylily

. dyinit  --project PROJECT 
dy-a slurm

# TO create a single sample manifest
head -n 2 .test_data/data/giab_30x_hg38_analysis_manifest.csv > config/analysis_manifest.csv

dy-r produce_snv_concordances -p -k -j 10 -n  # dry run

dy-r produce_snv_concordances -p -k -j 10 # run jobs, and wait for completion
```


##### Specify A Multi-Sample Manifest (in this case, all 7 GIAB samples)
```bash

tmux new -s slurm_test_30x_multi

mkdir /fsx/analysis_results/ubuntu/slurm_multi_30x_test
cd /fsx/analysis_results/ubuntu/slurm_multi_30x_test

git clone https://github.com/Daylily-Informatics/daylily.git  # or, if you have set ssh keys with github and intend to make changes:  git clone git@github.com:Daylily-Informatics/daylily.git
cd daylily

. dyinit  --project PROJECT 
dy-a slurm

# copy full 30x giab sample template to config/analysis_manifest.csv
cp .test_data/data/giab_30x_hg38_analysis_manifest.csv  config/analysis_manifest.csv

dy-r produce_snv_concordances -p -k -j 10 -n  # dry run

dy-r produce_snv_concordances -p -k -j 10 # run jobs, and wait for completion
```

#### Monitor Slurm Submitted Jobs

Once jobs begin to be submitted, you can monitor from another shell on the headnode(or any compute node) with:
```bash
# The compute fleet, only nodes in state 'up' are running spots. 'idle' are defined pools of potential spots not bid on yet.
sinfo
PARTITION AVAIL  TIMELIMIT  NODES  STATE NODELIST
i8*          up   infinite     12  idle~ i8-dy-gb64-[1-12]
i32          up   infinite     24  idle~ i32-dy-gb64-[1-8],i32-dy-gb128-[1-8],i32-dy-gb256-[1-8]
i64          up   infinite     16  idle~ i64-dy-gb256-[1-8],i64-dy-gb512-[1-8]
i96          up   infinite     16  idle~ i96-dy-gb384-[1-8],i96-dy-gb768-[1-8]
i128         up   infinite     28  idle~ i128-dy-gb256-[1-8],i128-dy-gb512-[1-10],i128-dy-gb1024-[1-10]
i192         up   infinite      1  down# i192-dy-gb384-1
i192         up   infinite     29  idle~ i192-dy-gb384-[2-10],i192-dy-gb768-[1-10],i192-dy-gb1536-[1-10]

# running jobs, usually reflecting all running node/spots as the spot teardown idle time is set to 5min default.
squeue
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
                 1      i192 D-strobe   ubuntu PD       0:00      1 (BeginTime)
# ST = PD is pending
# ST = CF is a spot has been instantiated and is being configured
# PD and CF sometimes toggle as the spot is configured and then begins running jobs.

 squeue
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
                 1      i192 D-strobe   ubuntu  R       5:09      1 i192-dy-gb384-1
# ST = R is running


# Also helpful
watch squeue

# also for the headnode
glances
```


#### SSH Into Compute Nodes
You can not access compute nodes directly, but can access them via the head node. From the head node, you can determine if there are running compute nodes with `squeue`, and use the node names to ssh into them.

```bash
ssh i192-dy-gb384-1
```


#### Delete Cluster
**warning**: this will delete all resources created for the ephemeral cluster, importantly, including the fsx filesystem. You must export any analysis results created in `/fsx/analysis_results` from the `fsx` filesystem  back to `s3` before deleting the cluster. 
- During cluster config, you will choose if Fsx and the EBS volumes auto-delete with cluster deletion. If you disable auto-deletion, these idle volumes can begin to cost a lot, so keep an eye on this if you opt for retaining on deletion.

##### Export `fsx` Analysis Results Back To S3
- Go to the 'fsx' AWS console and select the filesystem for your cluster.
- Under the `Data Repositories` tab, select the `fsx` filesystem and click `Export to S3`. Export can only currently be carried out back to the same s3 which was mounted to the fsx filesystem. 
- Specify the export path as `analysis_results` (or be more specific to an `analysis_results/subdir`), the path you enter is named relative to the mountpoint of the fsx filesystem on the cluster head and compute nodes, which is `/fsx/`. Start the export. This can take 10+ min.  When complete, confirm the data is now visible in the s3 bucket which was exported to. Once you confirm the export was successful, you can delete the cluster (which will delete the fsx filesystem).

##### Delete The Cluster, For Real
_note: this will not modify/delete the s3 bucket mounted to the fsx filesystem, nor will it delete the policyARN, or private/public subnets used to config the ephemeral cluster._
```bash
pcluster delete-cluster-instances -n <cluster-name> --region us-west-2
pcluster delete-cluster -n <cluster-name> --region us-west-2
``` 
- You can monitor the status of the cluster deletion using `pcluster list-clusters --region us-west-2` and/or `pcluster describe-cluster -n <cluster-name> --region us-west-2`. Deletion can take ~10min depending on the complexity of resources created and fsx filesystem size.

# Other Monitoring Tools

## PCUI
... Let me plug it again!

## Quick SSH Into Headnode
(also, can be done via pcui)

`bin/daylily-ssh-into-headnode`

_alias it for your shell:_ `alias goday="source ~/git_repos/daylily/bin/daylily-ssh-into-headnode"`


## AWS Cloudwatch
- The AWS Cloudwatch console can be used to monitor the cluster, and the resources it is using.  This is a good place to monitor the health of the cluster, and in particular the slurm and pcluster logs for the headnode and compute fleet.
- Navigate to your `cloudwatch` console, then select `dashboards` and there will be a dashboard named for the name you used for the cluster. Follow this link (be sure you are in the `us-west-2` region) to see the logs and metrics for the cluster.
- Reports are not automaticaly created for spot instances, but you may extend this base report as you like.  This dashboard is automatically created by `pcluster` for each new cluster you create (and will be deleted when the cluster is deleted).



 
# And There Is More

## S3 Reference Bucket & Fsx Filesystem

### PREFIX-omics-analysis-REGION Reference Bucket
Daylily relies on a variety of pre-built reference data and resources to run. These are stored in the `daylily-references-public` bucket. You will need to clone this bucket to a new bucket in your account, once per region you intend to operate in.  

> This is a design choice based on leveraging the `FSX` filesystem to mount the data to the cluster nodes. Reference data in this S3 bucket are auto-mounted an available to the head and all compute nodes (*Fsx supports 10's of thousands of concurrent connections*), further, as analysis completes on the cluster, you can choose to reflect data back to this bucket (and then stage elsewhere). Having these references pre-arranged aids in reproducibility and allows for the cluster to be spun up and down with negligible time required to move / create refernce data. 

> BONUS: the 7 giab google brain 30x ILMN read sets are included with the bucket to standardize benchmarking and concordance testing.

> You may add / edit (not advised) / remove data (say, if you never need one of the builds, or don't wish to use the GIAB reads) to suit your needs.

#### Reference Bucket Metrics
*Onetime* cost of between ~$27 to ~$108 per region to create bucket.
*monthly S3 standard* cost of ~$14/month to continue hosting it.
- Size: 617.2GB, and contains 599 files.
- Source bucket region: `us-west-2`
- Cost to store S3 (standard: $14.20/month, IA: $7.72/month, Glacier: $2.47 to $0.61/month)
- Data transfer costs to clone source bucket
  - within us-west-2: ~$3.40
  - to other regions: ~$58.00
- Accelerated transfer is used for the largest files, and adds ~$24.00 w/in `us-west-2` and ~$50 across regions.
- Cloning w/in `us-west-2` will take ~2hr, and to other regions ~7hrs.
- Moving data between this bucket and the FSX filesystem and back is not charged by size, but by number of objects, at a cost of `$0.005 per 1,000 PUT`. The cost to move 599 objecsts back and forth once to Fsx is `$0.0025`(you do pay for Fsx _when it is running, which is only when you choose to run analysus_).

### The `YOURPREFIX-omics-analysis-REGION` s3 Bucket
- Your new bucket name needs to end in `-omics-analysis-REGION` and be unique to your account.
- One bucket must be created per `REGION` you intend to run in.
- The reference data version is currently `0.7`, and will be replicated correctly using the script below.
- The total size of the bucket will be 779.1GB, and the cost of standard S3 storage will be ~$30/mo.
- Copying the daylily-references-public bucket will take ~7hrs using the script below.

#### daylily-references-public Bucket Contents
- `hg38` and `b37` reference data files (including supporting tool specific files).
- 7 google-brain ~`30x` Illunina 2x150 `fastq.gz` files for all 7 GIAB samples (`HG001,HG002,HG003,HG004,HG005,HG006,HG007`).
- snv and sv truth sets (`v4.2.1`) for all 7 GIAB samples in both `b37` and `hg38`.
- A handful of pre-built conda environments and docker images (for demonstration purposes, you may choose to add to your own instance of this bucket to save on re-building envs on new eclusters).
- A handful of scripts and config necessary for the ephemeral cluster to run.

_note:_ you can choose to eliminate the data for `b37` or `hg38` to save on storage costs. In addition, you may choose to eliminate the GIAB fastq files if you do not intend to run concordance or benchmarking tests (which is advised against as this framework was developed explicitly to facilitate these types of comparisons in an ongoing way).

##### Top Level Diretories
See the secion on [shared Fsx filesystem](#shared-fsx-filesystem) for more on hos this bucket interacts with these ephemeral cluster region specific S3 buckets.
```text
.
├── cluster_boot_config  # used to configure the head and compute nodes in the ephemeral cluster, is not mounted to cluster nodes
└── data # this directory is mounted to the head and compute nodes under /fsx/data as READ-ONLY. Data added to the S3 bucket will become available to the fsx mount, but can not be written to via FSX
    ├── cached_envs
    │   ├── conda
    │   └── containers
    ├── genomic_data
    │   ├── organism_annotations
    │   │   └── H_sapiens
    │   │       ├── b37
    │   │       └── hg38
    │   ├── organism_reads
    │   │   └── H_sapiens
    │   │       └── giab
    │   └── organism_references
    │       └── H_sapiens
    │           ├── b37
    │           └── hg38
    └── tool_specific_resources
        └── verifybam2
            ├── exome
            └── test
```




# Fsx Filesystem
Are region specific, and may only intereact with `S3` buckets in the same region as the filesystem. There are region specific quotas to be aware of.
- Fsx filesystems are extraordinarily fast, massively scallable (both in IO operations as well as number of connections supported -- you will be hard pressed to stress this thing out until you have 10s of thousands of concurrent connected instances).  It is also a pay-to-play product, and is only cost effective to run while in active use.  
- Daylily uses a `scratch` type instance, which auto-mounts the region specific `s3://PREFIX-omics-analysis-REGION/data` directory to the fsx filesystem as `/fsx/data`.  `/fsx` is available to the head node and all compute nodes.  
- When you delete a cluster, the attached `Fsx Lustre` filesystem will be deleted as well.  
- > **BE SURE YOU REFLECT ANALYSIS REUSLTS BACK TO S3 BEFORE DELETING YOUR EPHEMERAL CLUSTER** ... do this via the Fsx dashboard and create a data export task to the same s3 bucket you used to seed the fsx filesystem ( you will probably wish to define exporting `analysis_results`, which will export back to `s3://PREFIX-omics-analysis-REGION/FSX-export-DATETIME/` everything in `/fsx/analysis_results` to this new FSX-export directory.  **do not export the entire `/fsx` mount, this is not tested and might try to duplicate your reference data as well!** ).  This can take 10+ min to complete, and you can monitor the progress in the fsx dashboard & delete your cluster once the export is complete.
- Fsx can only mount one s3 bucket at a time, the analysis_results data moved back to S3 via the export should be moved again to a final destination (w/in the same region ideally) for longer term storage.  
- All of this handling of data is amendable to being automated, and if someone would like to add a cluster delete check which blocks deletion if there is unexported data still on /fsx, that would be awesome.
- Further, you may write to any path in `/fsx` from any instance it is mounted to, except `/fsx/data` which is read only and will only update if data mounted from the `s3://PREFIX-omics-analysis-REGION/data` is added/deleted/updated (not advised).

## Fsx Directory Structure
The following directories are created and accessible via `/fsx` on the headnode and compute nodes.
```text

/fsx/
├── analysis_results
│   ├── cromwell_executions  ## in development
│   ├── daylily   ## deprecated
│   └── ubuntu  ## <<<< run all analyses here <<<<
├── data  ## mounted to the s3 bucket PREFIX-omics-analysis-REGION/data
│   ├── cached_envs
│   ├── genomic_data
│   └── tool_specific_resources
├── miners  ## experimental & disabled by default
├── resources
│   └── environments  ## location of cached conda envs and docker images. so they are only created/pulled once per cluster lifetime.
├── scratch  ## scratch space for high IO tools
└── tmp ## tmp used by slurm by default
```




# In Progress // Future Development

## Re-enable Sentieon Workflows & Include in Benchmarking
- I have a demo lisc, and old working workflows (but they are ~2yrs out of date at this point).  I will be updating these workflows and including them in the benchmarking results.

## Add Strobe Aligner To Benchmarking
- The aligner is already included, but I have not been running it as my $ resources are v. limited.

## Using Data From Benchmarking Experiments, Complete The Comprehensive Cost Caclulator
- Rough draft script is running already, with best guesses for things like compute time per-x coverage, etc.

## Break Daylily Into 2 Parts: 1) Ephermal Cluster Manager 2) Analysis Pipeline
- The `daylily` repo grew from an analysis pipeline, and has co-mingled the ephmeral cluster infrastructure (which is not tied to any particular pipeline orchestrator). Breaking it into 2 parts will make things more modular and easier to maintain.

## Update Analysis Pipeline To Run With Snakemake v8.*
- A branch has been started for this work, which is reasonably straightforward. Tasks include:
  -  The AWS parallel cluster slurm snakemake executor, [pcluster-slurm](https://github.com/Daylily-Informatics/snakemake-executor-plugin-pcluster-slurm)  is written, but needs some additional features and to be tested at scale.
  -  Migrate from the current `daylily` `analysis_manifest.csv` to the snakemake `v8.*` `config/samples/units` format (which is much cleaner than the current manifest).
  -  The actual workflow files should need very little tweaking.

## Cromwell & WDL's
- Running Cromwell WDL's is in early stages, and preliminary & still lightly documented work can be found [here](config/CROMWELL/immuno/workflow.sh) ( using the https://github.com/wustl-oncology as starting workflows ).



<p valign="middle"><img src="docs/images/000000.png" valign="bottom" ></p>


# General Components Overview

  > Before getting into the cool informatics business going on, there is a boatload of complex ops systems running to manage EC2 spot instances, navigate spot markets, as well as mechanisms to monitor and observe all aspects of this framework. [AWS ParallelCluster](https://docs.aws.amazon.com/parallelcluster/latest/ug/what-is-aws-parallelcluster.html) is the glue holding everything together, and deserves special thanks.
  
![DEC_components_v2](https://user-images.githubusercontent.com/4713659/236144817-d9b26d68-f50b-423b-8e46-410b05911b12.png)

# Managed Genomics Analysis Services
The system is designed to be robust, secure, auditable, and should only take a matter of days to stand up. [Please contact me for further details](https://us21.list-manage.com/contact-form?u=434d42174af0051b1571c6dce&form_id=23d28c274008c0829e07aff8d5ea2e91).


![daylily_managed_service](https://user-images.githubusercontent.com/4713659/236186668-6ea2ec81-9fe4-4549-8ed0-6fcbd4256dd4.png)


<p valign="middle"><img src="docs/images/000000.png" valign="bottom" ></p>



## Some Bioinformatics Bits, Big Picture

### The DAG For 1 Sample Running Through The `BWA-MEM2ert+Doppelmark+Deepvariant+Manta+TIDDIT+Dysgu+Svaba+QCforDays` Pipeline
#### NOTE: *each* node in the below DAG is run as a self-contained job. Each job/node/rule is distributed to a suitable EC2 spot(or on demand if you prefer) instance to run. Each node is a packaged/containerized unit of work. This dag represents jobs running across sometimes thousands of instances at a time. Slurm and Snakemake manage all of the scaling, teardown, scheduling, recovery and general orchestration: cherry on top: killer observability & per project resource cost reporting and budget controls!
   ![](docs/images/assets/ks_rg.png)
   
   - The above is actually a compressed view of the jobs managed for a sample moving through this pipeline. This view is of the dag which properly reflects parallelized jobs.
   
     ![](docs/images/assets/ks_dag.png)




### Daylily Framework, Cont.

#### [Batch QC HTML Summary Report](http://daylilyinformatics.com:8082/reports/DAY_final_multiqc.html)
> The batch is comprised of google-brain Novaseq 30x HG002 fastqs, and again downsampling to: 25,20,15,10,5x.     
[Example report](http://daylilyinformatics.com:8082/reports/DAY_final_multiqc.html).
![](docs/images/assets/day_qc_1.png)

![](docs/images/assets/day_qc_2.png)
    
    
### [Consistent + Easy To Navigate Results Directory & File Structure](/docs/ops/dir_and_file_scheme.md)
   - A visualization of just the directories (minus log dirs) created by daylily
    _b37 shown, hg38 is supported as well_
![](docs/images/assets/tree_structure/tree.md)
     - [with files](docs/ops/tree_full.md   
    
### [Automated Concordance Analysis Table](http://daylilyinformatics.com:8081/components/daylily_qc_reports/other_reports/giabhcr_concordance_mqc.tsv)
  > Reported faceted by: SNPts, SNPtv, INS>0-<51, DEL>0-51, Indel>0-<51.
  > Generated when the correct info is set in the analysis_manifest.


#### [Performance Monitoring Reports]()
  > Picture and  list of tools

#### [Observability w/CloudWatch Dashboard](https://us-east-2.console.aws.amazon.com/cloudwatch/home?region=us-east-2#)
  > ![](docs/images/assets/cloudwatch.png)
  > ![](/docs/images/assets.cloudwatch_2.png)
  > ![](/docs/images/assets.cloudwatch3.png)

#### [Cost Tracking and Budget Enforcement](https://aws.amazon.com/blogs/compute/using-cost-allocation-tags-with-aws-parallelcluster/)
  > ![](https://d2908q01vomqb2.cloudfront.net/1b6453892473a467d07372d45eb05abc2031647a/2020/07/23/Billing-console-projects-grouping.png)
  - ![](docs/images/assets/costs1.png)
  - ![](docs/images/assets/costs2.png)
  
  
<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/000000.png" valign="bottom" ></a></p>



# Metrics Required To Make Informed Decisions About Choosing An Analysis Pipeline
To make informed decisions about choosing an analysis pipeline, there are four key metrics to consider: accuracy(as generally measured via Fscore), user run time, cost of analysis and reproducibility. Further consideration should then be given to the longevity of results (how results are stored, costs associated with storage, and the ease of access to results). All of these can not be optimized simultaneously, and trade-offs must be made. Making these tradeoffs with intention is the primary goal of the daylily framework.

## Accuracy / Precision / Recall / Fscore
- what is the pipelines perofrmance?

## User Run Time
- how long does it take to run the pipeline?

## Cost Of Analysis
### Init Cost
- is
### Compute Cost
- is
### Storage Cost (for computation)
- is
### Other Costs (ie: data transfer)
- is

## Cost of Storage
- is

## Reproducibility
- what is the asserted reproducibility of the pipeline? For 1 month? 1 year? 5 years? 20 years?
- And how is this tested?

## Longevity of Results
- how are results stored? What are the costs and access mechanisms for these results?

# Sentieon Tools & License
To activate sentieon bwa and sentieon DNA scope, edit the `config/day_profiles/{local,slurm}/templates/rule_config.yaml` file to uncomment the following:
```yaml
active_snv_callers:
    - deep
#    - sentd

active_aligners:
    - bwa2a:
        mkdup: dppl
#    - sent:  # uncomment to run aligner (same deduper must be used presently in all aligner outputs, unless bundled into the align rule)
#        mkdup: dppl  # One ddup'r per analysis run.  dppl and sent currently are active
#    - strobe:

```
- This will enable running these two tools.  You will also need a liscence file from sentieon in order to run these tools.  
- Please [contact them](https://www.sentieon.com/company/) to obtain a valid liscense . 
- Once you have a lisence file, edit the `dyinit` file to include the */fsx/* relative path to this file where `export SENTIEON_LICENSE=` is found.  
- Save the liscence file in the each region specific S3 reference bucket, ie: `s3://PREFIX-omics-analysis-REGION/data/cached_envs/`. When this bucket is mounted to the fsx filesystem, the liscence file will be available to all instances at `/fsx/data/cached_envs/`.



# Contributing
- [Contributing Guidelines](CONTRIBUTING.md)

# Versioning

Daylily uses [Semantic Versioning](https://semver.org/). For the versions available, see the [tags on this repository](https://github.com/Daylily-Informatics/daylily/tags).

# Known Issues

## _Fsx Mount Times Out During Headnode Creation & Causes Pcluster `build-cluster` To Fail_
If the `S3` bucket mounted to the FSX filesystem is too large (the default bucket is close to too large), this can cause Fsx to fail to create in time for pcluster, and pcluster time out fails.  The wait time for pcluster is configured to be much longer than default, but this can still be a difficult to identify reason for cluster creation failure. Probability for failure increases with S3 bucket size, and also if the imported directories are being changed during pcluster creation. Try again, try with a longer timeount, and try with a smaller bucket (ie: remove one of the human reference build data sets, or move to a different location in the bucket not imported by Fsx)

## Cloudstack Formation Fails When Creating Clusters In >1 AZ A Region (must be manually sorted ATM)
The command `bin/init_cloudstackformation.sh ./config/day_cluster/pcluster_env.yml "$res_prefix" "$region_az" "$region" $AWS_PROFILE` does not yet gracefully handle being run >1x per region.  The yaml can be edited to create the correct scoped resources for running in >1 AZ in a region (this all works fine when running in 1AZ in >1 regions), or you can manually create the pub/private subnets, etc for running in multiple AZs in a region. The fix is not difficult, but is not yet automated.


# Compliance / Data Security
Is largely in your hands. AWS Parallel Cluster is as secure or insecure as you set it up to be. https://docs.aws.amazon.com/parallelcluster/v2/ug/security-compliance-validation.html

<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/000000.png" valign="bottom" ></a></p>

<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/0000002.png" valign="bottom" ></a></p>

# [DAY](https://en.wikipedia.org/wiki/Margaret_Oakley_Dayhoff)![](https://placehold.co/60x35/ff03f3/fcf2fb?text=LILLY)

_named in honor of Margaret Oakley Dahoff_ 
 
 
 
