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

# Test if the directories exist in the S3 bucket
data_dir="$bucket_url/data/$s3_reference_data_version"
cluster_boot_config_dir="$bucket_url/cluster_boot_config/$s3_reference_data_version"

echo "Checking if directories exist in the S3 bucket..."

# Check if data directory exists
if aws s3 ls "$data_dir/" &> /dev/null; then
    echo "Data directory '$data_dir' exists."
else
    echo "Error: Data directory '$data_dir' does not exist. Exiting."
    return 1
fi

# Check if cluster_boot_config directory exists
if aws s3 ls "$cluster_boot_config_dir/" &> /dev/null; then
    echo "Cluster boot config directory '$cluster_boot_config_dir' exists."
else
    echo "Error: Cluster boot config directory '$cluster_boot_config_dir' does not exist. Exiting."
    return 1
fi


cmda="aws s3 cp s3://rcrf-omics-analysis/daylily/data/$s3_reference_data_version s3://${your_bucket_prefix}-omics-analysis/data --recursive --request-payer requester "
cmdb="aws s3 cp s3://rcrf-omics-analysis/daylily/cluster_boot_config/$s3_reference_data_version s3://${your_bucket_prefix}-omics-analysis/cluster_boot_config --recursive  --request-payer requester "

echo "Now, run the following commands in your shell."
echo "  $cmda & $cmdb & echo 'waiting for both to complete' & wait "

```



#### 4. Create A New `pcluster` Stack
From the `daylily` repository root:
```bash
bin/init_cloudstackformation.sh ./config/day_cluster/pcluster_env.yml <short-resource-prefix-to-use>
```
- This will run and create a VPC, private and public subnet as well as an IAM policy which allows pcluster to tag resources for budget tracking and allow budget tools to see these tags.
- The script should block the terminal while running, and report back on success.  If failure, go to the cloudformation console and look at the logs for the stack to see what went wrong & delete the stack before attempting again.
- **please record/take note** of the `Public Subnet ID`, `Private Subnet ID`, and `Policy ARN` will be reported here and needed in the `daycli` setup.


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
    - You will be prompted to select:
      - the full path to your $HOME/.ssh/<mykey>.pem (from detected .pem files)
      - the `s3` bucket you created and seeded from the `daylily-references-public` bucket. If your bucket name contains `omics-analysis`, it will be suggested in a list of buckets, otherwise, select `1` and enter a valid s3 bucket url.

      - the `Public Subnet ID` from the cloudformation stack output.
      - the `Private Subnet ID` from the cloudformation stack output.
      - the `Policy ARN` from the cloudformation stack output.
      - the `YOUR-omics-analysis` s3 bucket name.
      - the short resource prefix you used in the cloudformation stack creation.

## Steps To Follow For Each New Cluster