#!/bin/bash

# Check if the script is running in Zsh and emulate Zsh behavior
if [ -n "$ZSH_VERSION" ]; then
  emulate -L zsh  # Ensure Zsh behaves like Zsh (if required)
fi

# Default configuration
pass_on_warn=false  # Default for warnings causing failures

# Function to display help
usage() {
    echo "Usage: source $0 [--region-az REGION-AZ] [--pass-on-warn]"
    echo "       --region-az REGION-AZ   Specify the AWS region and availability zone (default: us-west-2d)"
    echo "       --pass-on-warn          Allow warnings to pass without failure (default: false)"
}

# Check if the script is sourced
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "Error: This script must be sourced, not executed directly. Use 'source $0' to run."
    return 3
fi

# Check if no arguments are provided and show usage
if [[ $# -eq 0 ]]; then
    usage
    return 0
fi

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --region-az)
            region_az="$2"
            region="${region_az:0:-1}"  # Remove the last character (Bash and Zsh compatible)
            shift 2
            ;;
        --pass-on-warn)
            pass_on_warn=true
            shift
            ;;
        -h|--help)
            usage
            return 0
            ;;
        *)
            echo "Error: Unknown option: $1"
            usage
            return 3
            ;;
    esac
done

# Function to handle warnings
handle_warning() {
    local message="$1"
    echo -e "\nWARNING: $message"
    echo "Consider setting --pass-on-warn to continue without failure."
    if ! $pass_on_warn; then
        return 3
    fi
}

# Validate AWS region
validate_region() {
    echo "RGN: $1"
    echo "NOTE: Run 'aws configure set region $1' to set the region for pcluster CLI usage."
    sleep 2.3

    if [[ "$1" == "us-west-2" ]]; then
        echo "Region '$1' confirmed as valid."
        return 0
    else
        echo "Warning: Region '$1' is not 'us-west-2'. It is recommended to use 'us-west-2'."
    fi
}

# Main logic
echo "Welcome to the Daylily CLI Setup"
echo "Daylily is configured to run in region $region (AZ: $region_az)."
echo "Ensure your ~/.aws/config file matches this region."

# Ensure Conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: Conda is not installed. Please install Miniconda."
    return 3
fi

# Activate or create the 'DAYCLI' Conda environment
if conda env list | grep -q "^DAYCLI "; then
    echo "Activating existing 'DAYCLI' environment."
else
    echo "Creating 'DAYCLI' environment."
    conda create -y -n DAYCLI -c conda-forge python parallel nodejs==18.15.0 aws-parallelcluster==3.10.1 flask==2.2.5
fi

conda activate DAYCLI
pip install colr==0.9.1 pyyaml==6.0.2 tabulate==0.9.0

# Check if GNU Parallel is installed
if ! parallel --version &> /dev/null; then
    echo "Error: GNU Parallel is not installed. Please install it via 'conda install -c conda-forge parallel'."
    return 3
else
    parallel --citation <<< "will cite"
fi

# Check AWS CLI version
aws_version=$(aws --version 2>&1 | cut -d ' ' -f 1 | sed 's/aws-cli\///g')
if [[ "$aws_version" != "1.27.92" ]]; then
    handle_warning "Expected AWS CLI version 1.27.92, but found $aws_version."
fi

# Validate AWS credentials
echo "Verifying AWS credentials..."
if ! aws sts get-caller-identity --region "$region" &> /dev/null; then
    echo "Error: AWS credentials are invalid or do not have access to the AWS account in region $region."
    return 3
else
    echo "AWS credentials verified successfully."
fi

# Call region validation
validate_region "$region" || return 3


ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)

# Get AWS caller identity
caller_identity_arn=$(aws sts get-caller-identity --query 'Arn' --output text --region "$region")
required_policies=("pcluster-omics-analysis-fleet")

# Check policies function
check_policies() {
    local entity_name="$1"
    local entity_type="$2"
    local policies=()

    if [[ "$entity_type" == "user" ]]; then
        policies=$(aws iam list-user-policies --user-name "$entity_name" --query 'PolicyNames' --output text --region "$region")
    else
        policies=$(aws iam list-role-policies --role-name "$entity_name" --query 'PolicyNames' --output text --region "$region")
    fi

    for policy_name in "${required_policies[@]}"; do
        if ! echo "$policies" | grep -wq "$policy_name"; then
            missing_policies+=("$policy_name")
        fi
    done
}

# IAM policy check (continue on warning for root users)
echo "Checking required policies for your AWS user/role..."
check_policies "$caller_identity_arn" "role"
if [[ ${#missing_policies[@]} -gt 0 ]]; then
    handle_warning "Missing required inline policies: ${missing_policies[*]}"
else
    echo "All required policies are attached."
fi

# AWS Quota Check Function
check_quota() {
    local resource_name="$1"
    local quota_code="$2"
    local recommended_count="$3"
    local sccode="$4"
    
    echo "Checking quota for $resource_name..."
    quota=$(aws service-quotas get-service-quota --service-code $sccode --quota-code "$quota_code" --region "$region" 2>/dev/null)

    if [[ $? -ne 0 ]]; then
        handle_warning "Error: Unable to retrieve quota for $resource_name."
        return 3
    fi

    quota_value=$(echo "$quota" | jq -r '.Quota.Value')

    if (( $(echo "$quota_value < $recommended_count" | bc -l) )); then
        handle_warning "Warning: $resource_name quota is below the recommended $recommended_count. Current quota: $quota_value."
    else
        echo "$resource_name quota is sufficient: $quota_value."
    fi
}

# Quota Checks for Key AWS Resources
echo "Performing AWS quota checks..."

# Dedicated Instances
check_quota "Spot vCPU Max" "L-1216C47A" 20 ec2

# Spot vCPU Max Quota 
check_quota "Spot vCPU Max" "L-34B43A08" 310 ec2

# VPC Quota (Default: 5 per region)
check_quota "VPCs" "L-F678F1CE" 4 vpc

# Elastic IP Quota (Default: 5 per region)
check_quota "Elastic IPs" "L-0263D0A3" 4 ec2

# NAT Gateway Quota (Default: 5 per AZ)
check_quota "NAT Gateways" "L-FE5A380F" 4 vpc

# Internet Gateway Quota (Default: 5 per region)
check_quota "Internet Gateways" "L-A4707A72" 4 vpc

# Check for 'pcluster-omics-analysis' policy and create if necessary
echo "Checking for required IAM policy 'pcluster-omics-analysis'..."
policy_name="pcluster-omics-analysis"
policy_arn=$(aws iam list-policies --query "Policies[?PolicyName=='$policy_name'].Arn" --output text --region "$region")
if [[ -z "$policy_arn" ]]; then
    echo "Creating IAM policy '$policy_name'..."
    policy_document='{
      "Version": "2012-10-17",
      "Statement": [
        {
          "Effect": "Allow",
          "Action": "iam:CreateServiceLinkedRole",
          "Resource": "*",
          "Condition": {
            "StringLike": {
              "iam:AWSServiceName": "spot.amazonaws.com"
            }
          }
        }
      ]
    }'
    aws iam create-policy --policy-name "$policy_name" --policy-document "$policy_document" --region "$region"
fi


## PEM File Check
# Fetch all key pairs in the given region
key_data=$(aws ec2 describe-key-pairs --region "$region" --query 'KeyPairs[*].[KeyName, KeyType]' --output text | grep omics)

# Filter for keys with the required suffix
filtered_keys=()
while IFS=$'\t' read -r key_name key_type; do
    if [[ $key_type == "ed25519" ]]; then
        pem_file_f=$HOME/.ssh/"${key_name}.pem"
        pem_status="pem - NOT FOUND"
        [[ -f "$pem_file_f" ]] && pem_status="pem - FOUND"
        filtered_keys+=("$key_name|$key_type|$pem_status|$pem_file_f")
    fi
done <<< "$key_data"

# Present valid keys to the user
echo "Available ED25519 keys with matching PEM files:"
valid_keys=()
index=1
for key_info in "${filtered_keys[@]}"; do
    IFS='|' read -r key_name key_type pem_status pem_file <<< "$key_info"
    if [[ $pem_status == "pem - FOUND" ]]; then
        valid_keys+=("$key_name")
        echo "[$index] Key: $key_name | Type: $key_type | PEM: $pem_status"
        ((index++))
    fi
done

# Check if no valid keys were found
if [[ ${#valid_keys[@]} -eq 0 ]]; then
    echo "No valid ED25519 keys with matching PEM files were found."
    return 3
fi
# Present valid keys to the user using 'select'
echo "Select a valid ED25519 key with a matching PEM file:"
PS3="Enter the corresponding number: "
select key_info in "${filtered_keys[@]}"; do
    if [[ -n "$key_info" ]]; then
        # Extract the key name from the selected item
        selected_key=$(echo "$key_info" | awk -F'|' '{print $1}')
        echo "You selected: $selected_key"
        break
    else
        echo "Invalid selection. Please try again."
    fi
done


# Proceed with further operations using the selected key
# Example: connect to an instance or perform other operations
pem_file=$HOME/.ssh/"${selected_key}.pem"

# Check if the file exists
if [[ -e "$pem_file" ]]; then
    echo "File $pem_file exists. Proceeding with further operations."
else
    echo "File $pem_file does not exist. Exiting."
    return 3  # Use 'exit 1' if outside a function
fi
echo "Using key: $selected_key with PEM file $pem_file"


# Region Check
echo "Checking if the 'daylily-omics-analysis-${region}' budget exists..."
budget_name="daylily-omics-analysis-${region}"

# Query AWS for the budget
budget_exists=$(aws budgets describe-budgets \
    --query "Budgets[?BudgetName=='$budget_name'] | [0].BudgetName" \
    --output text --region "$region" --account-id "$ACCOUNT_ID")

# Check if the output is empty
if [[ -z "$budget_exists" || "$budget_exists" == "None" ]]; then
    echo "Budget '$budget_name' not found."
    echo -n "Enter an email address to receive alerts: "
    read user_email
    echo -n "Enter a budget amount (default: 500): "
    read amount
    amount=${amount:-500}

    echo "Creating the budget..."
    bin/create_budget.sh -p "$budget_name" -r "$region" -t 75 -e "$user_email" -a "$amount"
else
    echo "Default budget '$budget_name' exists. Proceeding"
fi


# Query available S3 buckets and present only those in the correct region
echo "Fetching S3 buckets matching 'omics-analysis' in $region..."
buckets=$(aws s3api list-buckets --query "Buckets[].Name" --output text | tr '\t' '\n' | \
while read -r bucket; do \
  region=$(aws s3api get-bucket-location --bucket "$bucket" --query "LocationConstraint" --output text); \
  echo "$bucket ($region)" | grep "omics-analysis"; \
done)

matching_buckets=()
while IFS= read  bucket; do
    bucket_region=$(echo "$bucket" | awk '{print $NF}')
    if [[ "$bucket_region" == "($region)" ]]; then
        matching_buckets+=("$bucket")
    fi
done <<< "$buckets"

if [[ ${#matching_buckets[@]} -eq 0 ]]; then
    echo "No matching buckets found in region $region."
    return 1
fi

echo "Select a matching S3 bucket:"
select bucket_choice in "${matching_buckets[@]}"; do
    selected_bucket=$(echo "$bucket_choice" | awk '{print $1}')
    echo "You selected: $selected_bucket"
    break
done

bucket=$(echo $bucket_choice | cut -f 1 -d " " )


# Check if the required folders exist in the bucket
check_s3_folder() {
    local bucket=$1
    local folder=$2

    if aws s3 ls "s3://$bucket/$folder/" > /dev/null 2>&1; then
        echo "Folder s3://$bucket/$folder/ exists."
    else
        echo "ERROR: Folder s3://$bucket/$folder/ does not exist."
        return 1
    fi
}

echo "Validating required folders in bucket: $selected_bucket"
check_s3_folder "$bucket" "data" || return 3
check_s3_folder "$bucket" "cluster_boot_config" || return 3 

echo "All required folders exist in $bucket."

## TODO: add a check to confirim the bucket has the expected s3://%{bucket}/data and s3://%{bucket}/cluster_boot_config folders

# Query for public and private subnets and ARN
echo "Checking subnets and ARN in the specified AZ: $region_az..."

public_subnets=$(aws ec2 describe-subnets \
  --query "Subnets[?AvailabilityZone=='$region_az'].{ID: SubnetId, Name: Tags[?Key=='Name']|[0].Value}" \
  --region "$region" --output text | grep "Public Subnet"
  )
  
private_subnets=$(aws ec2 describe-subnets \
  --query "Subnets[?AvailabilityZone=='$region_az'].{ID: SubnetId, Name: Tags[?Key=='Name']|[0].Value}" \
  --region "$region" --output text | grep "Private Subnet"
)

arn_count=$(aws iam list-policies --query 'Policies[?PolicyName==`pclusterTagsAndBudget`].Arn' --output text --region "$region" | wc -l)

if [[ -z "$public_subnets" && -z "$private_subnets" && $arn_count -eq 0 ]]; then
    echo "All required resources are missing. Running the initialization script..."
    bin/init_cloudstackformation.sh ./config/day_cluster/pcluster_env.yml "daylily-cloudstack-$region_az" "$region_az" "$region"
elif [[ -z "$public_subnets" || -z "$private_subnets" || $arn_count -eq 0 ]]; then
    echo "Error: Incomplete setup. Ensure public and private subnets and the required ARN are available."
    return 1
else
    echo "All required resources are available."
fi

echo "Setup complete. Proceeding with the remaining steps..."

# Select public and private subnets
public_subnets_new=$(aws ec2 describe-subnets --query "Subnets[*].[SubnetId, Tags[?Key=='Name'].Value | [0], AvailabilityZone]" --region "$region" --output text | grep "Public Subnet" | grep $region_az)
echo "Select a Public Subnet:"
public_subnet_choices=()
while read -r subnet_id subnet_name; do
    public_subnet_choices+=("$subnet_name ($subnet_id)")
done <<< "$public_subnets_new"

select public_subnet_choice in "${public_subnet_choices[@]}"; do
    public_subnet=$(echo "$public_subnet_choice" | sed -n 's/.*(\(.*\)).*/\1/p')
    echo "You selected Public Subnet: $public_subnet"
    break
done

echo "Select a Private Subnet:"
private_subnets_new=$(aws ec2 describe-subnets --query "Subnets[*].[SubnetId, Tags[?Key=='Name'].Value | [0], AvailabilityZone]" --region "$region" --output text | grep "Public Subnet" | grep $region_az)
private_subnet_choices=()
while read -r subnet_id subnet_name; do
    private_subnet_choices+=("$subnet_name ($subnet_id)")
done <<< "$private_subnets_new"

select private_subnet_choice in "${private_subnet_choices[@]}"; do
    private_subnet=$(echo "$private_subnet_choice" | sed -n 's/.*(\(.*\)).*/\1/p')
    echo "You selected Private Subnet: $private_subnet"
    break
done

# Select an IAM policy ARN for 'pclusterTagsAndBudget'
policy_arns=($(aws iam list-policies --query 'Policies[?PolicyName==`pclusterTagsAndBudget`].Arn' --output text --region "$region"))
echo "Select an IAM policy ARN for 'pclusterTagsAndBudget':"
select arn_policy_id in "${policy_arns[@]}"; do
    echo "You selected: $arn_policy_id"
    break
done

# Cluster name input
echo -n "Enter the name for your cluster (alphanumeric and '-'): "
read cluster_name
if [[ ! "$cluster_name" =~ ^[a-zA-Z0-9\-]+$ ]]; then
    echo "Error: Invalid cluster name."
    return 3
fi

# Cluster config file input
echo -n "Enter the path to the Daylily cluster config YAML file [press Enter to use default 'config/day_cluster/prod_cluster.yaml']: "
read cluster_yaml

if [[ -z "$cluster_yaml" ]]; then
    cluster_yaml="config/day_cluster/prod_cluster.yaml"
    echo "Using default: $cluster_yaml"
fi

if [[ ! -f "$cluster_yaml" ]]; then
    echo "Error: YAML file '$cluster_yaml' does not exist."
    return 3
fi 

xmr_pool_url=na
xmr_wallet=na

# XMR Mining Setup
echo "Enable XMR mining with idle compute? (1) No (2) Yes:"
select enable_xmr in "No" "Yes"; do
    if [[ "$enable_xmr" == "Yes" ]]; then
        echo -n "Specify a mining pool URL (default=pool.supportxmr.com:3333): "
        read xmr_pool_url
        [[ -z "$xmr_pool_url" ]] && xmr_pool_url="pool.supportxmr.com:3333"

        echo -n "Specify a destination wallet (default wallet will donate coins to health research): "
        read xmr_wallet
        [[ -z "$xmr_wallet" ]] && xmr_wallet="42s1tbj3qp6T2QBw9BzYxxRn1e9afiZhcRr8yoe2QYLn5ZG6Fk9XzEh1719moDeUgY6tReaJDF464ZhbEyg4cXMYNRXJve2"

        if ! nc -z $(echo "$xmr_pool_url" | cut -d: -f1) $(echo "$xmr_pool_url" | cut -d: -f2); then
            handle_warning "Cannot reach the mining pool URL: $xmr_pool_url"
        else
            echo "Mining pool URL is reachable."
        fi
    else
        enable_xmr="0"
    fi
    break
done


# Prepare configuration
mkdir -p $HOME/.config/daylily
target_conf=$HOME/.config/daylily/${cluster_name}_cluster.yaml.init
target_conf_fin=$HOME/.config/daylily/${cluster_name}_cluster.yaml
regsub_vals=$HOME/.config/daylily/${cluster_name}_cluster_init_vals.txt

cp "$cluster_yaml" "$target_conf"
pem_name=$(basename "$pem_file" | cut -d '.' -f 1)

bucket_url="s3://$selected_bucket"
bucket_name=$bucket

# Write variables to config
cat <<EOF > $regsub_vals
REGSUB_REGION=$region
REGSUB_PUB_SUBNET=$public_subnet
REGSUB_KEYNAME=$pem_name
REGSUB_S3_BUCKET_INIT=$bucket_url
REGSUB_S3_BUCKET_NAME=$bucket_name
REGSUB_S3_IAM_POLICY=$arn_policy_id
REGSUB_PRIVATE_SUBNET=$private_subnet
REGSUB_S3_BUCKET_REF=$bucket_url
REGSUB_XMR_MINE=$enable_xmr
REGSUB_XMR_POOL_URL=$xmr_pool_url
REGSUB_XMR_WALLET=$xmr_wallet
EOF

bash bin/other/regsub_yaml.sh $regsub_vals $target_conf


echo "Calculating max spot bid prices per partition resource group..."
echo "PYTHON: $(which python)"

python bin/calcuate_spotprice_for_cluster_yaml.py -i $target_conf -o $target_conf_fin --az $region_az --profile $AWS_PROFILE

# Check the exit status of the previous command
if [[ $? -ne 0 ]]; then
    echo "Error: Failed to calculate spot bid prices."
    echo "This is likely due to a policy/permissions issue in AWS."
    echo "Please refer to the quickstart guide for more information on required policies."
    return 3  # Use 'exit 3' if outside a function or script
fi

# Dry run cluster creation
echo "Running a cluster creation dry run..."
json_response=$(pcluster create-cluster -n "$cluster_name" -c "$target_conf_fin" --dryrun true --region "$region")

message=$(echo "$json_response" | jq -r '.message')

# Validate the message
if [[ "$message" == "Request would have succeeded, but DryRun flag is set." ]]; then
    echo "Dry run successful. Proceeding with cluster creation."
elif [[ "$DAY_BREAK" == "1" ]]; then
    echo "DAY_BREAK == '1'.  Exiting without creating the cluster."
    return 0
else
    echo "Error: Dry run failed (  pcluster create-cluster -n $cluster_name -c $target_conf_fin --dryrun true --region $region ). Exiting."
    return 3
fi

# Create the cluster
echo "Creating the cluster '$cluster_name' in region $region"
pcluster create-cluster -n "$cluster_name" -c "$target_conf_fin" --region "$region"

# Monitor cluster creation
echo "Monitoring cluster creation... this may take 10-15 minutes."
python bin/helpers/watch_cluster_status.py "$region"

echo "\nOnce the cluster is complete, configure the head node by running 'source ./bin/daylily-cfg-headnode $pem_file $region'."
echo "\n\t....trying this now!\n"
sleep 2
source ./bin/daylily-cfg-headnode $pem_file $region
echo "\n\t....fin!\n"