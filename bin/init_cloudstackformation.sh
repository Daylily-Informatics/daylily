#!/bin/bash

echo "\n\n\t ALERT 1:::>>> You will need valid aws credentials in your ~/.aws directory in order to proceed.\n\n"
echo "\n\n\tALERT 2:::>>> Please also be sure to have created a key pair in the AWS console which is ed25519 type & save the resulting .pem file (& chmod 400 it).  This will be needed for creating the cluster after this step completes.\n\n... pausing"
sleep 3.3
echo "\n\n\t... continuing\n\n"

# Function to display usage information
function usage() {
  echo "Usage: $0 TEMPLATE_FILE RESOURCE_PREFIX"
  echo
  echo "Arguments:"
  echo "  TEMPLATE_FILE      Full path to the CloudFormation template file (e.g., ./config/day_cluster/pcluster_env.yml)"
  echo "  RESOURCE_PREFIX    A prefix for all resources created by the stack (alphabets and dashes only)"
  echo
  echo "Options:"
  echo "  --help             Display this help message and exit"
}

# Check if help is requested
if [[ "$1" == "--help" ]]; then
  usage
  exit 0
fi

# Check for missing arguments
if [[ -z "$1" || -z "$2" ]]; then
  echo "Error: Missing required arguments."
  usage
  exit 1
fi

TEMPLATE_FILE=$1 # full path to the CloudFormation template file, probably ./config/day_cluster/pcluster_env.yml
RESOURCES_PREFIX=$2 # a prefix for all resources created by the stack, only alphamum and dash alloed, be succinct


# Validate if the template file exists
if [[ ! -f "$TEMPLATE_FILE" ]]; then
  echo "Error: The template file $TEMPLATE_FILE does not exist."
  exit 1
fi

# Validate resource prefix: only allow alphabets and dashes
if [[ ! "$RESOURCES_PREFIX" =~ ^[a-zA-Z-]+$ ]]; then
  echo "Error: The resource prefix can only contain alphabets and dashes."
  exit 1
fi

# Variables
STACK_NAME="pcluster-vpc-stack"
REGION="us-west-2"  # resources are tuned for us-west-2c, it is not advised to change zones
AVAILABILITY_ZONE="us-west-2c"
VPC_CIDR="10.0.0.0/16"
PUBLIC_SUBNET_CIDR="10.0.0.0/24"
PRIVATE_SUBNET_CIDR="10.0.1.0/24"


# Validate that the availability zone starts with the region
if [[ "$AVAILABILITY_ZONE" != "$REGION"* ]]; then
  echo "Error: Availability zone ($AVAILABILITY_ZONE) does not match the region ($REGION)."
  exit 1
fi

# Create the CloudFormation Stack
aws cloudformation create-stack \
  --stack-name $STACK_NAME \
  --template-body file://$TEMPLATE_FILE \
  --region $REGION \
  --parameters ParameterKey=EnvironmentName,ParameterValue=$RESOURCES_PREFIX \
               ParameterKey=VpcCIDR,ParameterValue=$VPC_CIDR \
               ParameterKey=PublicSubnetCIDR,ParameterValue=$PUBLIC_SUBNET_CIDR \
               ParameterKey=PrivateSubnetCIDR,ParameterValue=$PRIVATE_SUBNET_CIDR \
               ParameterKey=AvailabilityZone,ParameterValue=$AVAILABILITY_ZONE \
  --capabilities CAPABILITY_NAMED_IAM

# Wait for the stack to finish creating (success or failure)
echo "Waiting for stack creation to complete..."
aws cloudformation wait stack-create-complete --stack-name $STACK_NAME --region $REGION

# Check if stack creation succeeded or failed
STACK_STATUS=$(aws cloudformation describe-stacks --stack-name $STACK_NAME --region $REGION --query "Stacks[0].StackStatus" --output text)

if [[ "$STACK_STATUS" == "CREATE_COMPLETE" ]]; then
  echo "Stack creation succeeded."

  # Get the Policy ARN, Public Subnet ID, and Private Subnet ID
  POLICY_ARN=$(aws cloudformation describe-stacks --stack-name $STACK_NAME --region $REGION --query "Stacks[0].Outputs[?OutputKey=='PclusterPolicy'].OutputValue" --output text)
  PUBLIC_SUBNET_ID=$(aws cloudformation describe-stacks --stack-name $STACK_NAME --region $REGION --query "Stacks[0].Outputs[?OutputKey=='PublicSubnets'].OutputValue" --output text)
  PRIVATE_SUBNET_ID=$(aws cloudformation describe-stacks --stack-name $STACK_NAME --region $REGION --query "Stacks[0].Outputs[?OutputKey=='PrivateSubnet'].OutputValue" --output text)

  echo "Policy ARN: $POLICY_ARN"
  echo "Public Subnet ID: $PUBLIC_SUBNET_ID"
  echo "Private Subnet ID: $PRIVATE_SUBNET_ID"

else
  echo "Stack creation failed with status: $STACK_STATUS"

  # Get the stack events to show failure reason
  aws cloudformation describe-stack-events --stack-name $STACK_NAME --region $REGION --query "StackEvents[?ResourceStatus=='CREATE_FAILED'].{Resource:LogicalResourceId,Reason:ResourceStatusReason}" --output table
fi
