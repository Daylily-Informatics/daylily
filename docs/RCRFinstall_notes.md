# Overview

Capturing the steps I took to get the daylily framework up and running in a `classic` AWS environment.

# Before Beginning
- `daylily` was tuned to run in the `us-west-2c` AZ, for cost and availabiliyu of some specific compute resources not found in all AZs. The install should all build and proceed using the `us-west-2c` AZ.
- A cloudstack formation template will run which will create some important resources for the epheral cluster, namely: the public VPC, public subnet, private subnet, and a policy to allow tagging and budget tracking of resources by pcluster.
- The install steps should work for both `bash` and `zsh` shells, but the `bash` shell is the default for the `pcluster` environment.

# Steps


## Clone The `daylily` Repository
```bash
git clone git@github.com:Daylily-Informatics/daylily.git
cd daylily
```

### Install Conda (if not already installed)
_bash and zsh are known to work_
- [Conda, miniconda specifically, is required to be setup in your shell to proceed.](docs/install/Install.md#run-daylily-init)
 

## Only Necessary To Follow Once Per New Installation

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
- Move, copy and chmod: `mv ~/Downloads/<yourkey>.pem  ~/.ssh/<yourkey>.pem && chmod 400 ~/.ssh/<yourkey>.pem`.
- **please record/take note** of the full path to this file, as it will be needed in the `daycli` setup.

#### 3. Create `YOUR-omics-analysis` s3 Bucket
- In `us-west-2`.
- **please record/take note** of the bucket name, as it will be needed in the `daycli` setup.

##### Copy The `daylily` Public References/Resources To Your New Bucket
From your terminal/shell (and to be safe, in a screen or tmux session as it may run for a while):
**dryrun suggested**
```bash
aws s3 cp s3://rcrf-omics-analysis/daylily/ s3://<YOUR>-omics-analysis/daylily/ --recursive --dryrun  --request-payer requester
```

**actual copy** .. _note, requestor pays is set, which should only cost a few $_
```bash
aws s3 cp s3://daylily-references-public/ s3://<YOUR>-omics-analysis/ --recursive  --request-payer requester
```


#### 4. Create A New `pcluster` Stack
From the `daylily` repository root:
```bash
bin/init_cloudstackformation.sh ./config/day_cluster/pcluster_env.yml <short-resource-prefix-to-use>
```
- This will run and create a VPC, private and public subnet as well as an IAM policy which allows pcluster to tag resources for budget tracking and allow budget tools to see these tags.
- The script should block the terminal while running, and report back on success.  If failure, go to the cloudformation console and look at the logs for the stack to see what went wrong & delete the stack before attempting again.
- **please record/take note** of the `public-subnet`, `private-subnet`, and `arn-tag-policy` will be reported here and needed in the `daycli` setup.


### Local `pcluster` DAYCLI Setup

#### 2. Install The `daycli` CLI
From the `daylily` repository root:
```bash
source bin/daylily-init
```
- Your AWS cli credentials are automatically detected via the aws cli, and you will be informed of success/fail for this step. You can simply hit enter to proceed through this.
- Then,you will be prompted to enter the following (which is collected in the prior steps):
- **full path to your .ssh/<mykey>.pm**,  _so on a Mac `/Users/<you>/.ssh/<mykey>.pem`_
- 

## Steps To Follow For Each New Cluster