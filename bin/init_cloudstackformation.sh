#!/bin/bash


TEMPLATE_FILE=$1

# Variables
STACK_NAME="pcluster-vpc-stack"
REGION="us-west-2"
ENVIRONMENT_NAME="MyEnvironmentName"
VPC_CIDR="10.0.0.0/16"
PUBLIC_SUBNET_CIDR="10.0.0.0/24"
PRIVATE_SUBNET_CIDR="10.0.1.0/24"

# Create the CloudFormation Stack
aws cloudformation create-stack \
  --stack-name $STACK_NAME \
  --template-body file://$TEMPLATE_FILE \
  --region $REGION \
  --parameters ParameterKey=EnvironmentName,ParameterValue=$ENVIRONMENT_NAME \
               ParameterKey=VpcCIDR,ParameterValue=$VPC_CIDR \
               ParameterKey=PublicSubnetCIDR,ParameterValue=$PUBLIC_SUBNET_CIDR \
               ParameterKey=PrivateSubnetCIDR,ParameterValue=$PRIVATE_SUBNET_CIDR \
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
