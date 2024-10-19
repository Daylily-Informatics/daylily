#!/usr/bin/env bash

# Function to display usage information
print_help() {
  echo "Usage: $0 -p <project_name> -a <amount> -r <region> [-h]"
  echo ""
  echo "Options:"
  echo "  -p    Project name (required)"
  echo "  -a    Budget amount in USD (required)"
  echo "  -r    AWS region (required)"
  echo "  -h    Display this help message"
  echo ""
  echo "Example:"
  echo "  $0 -p my-project -a 5000 -r us-west-2"
}

# Parse command-line options
while getopts "p:a:r:h" opt; do
  case ${opt} in
    p) PROJECT_NAME=$OPTARG ;;
    a) AMOUNT=$OPTARG ;;
    r) REGION=$OPTARG ;;
    h) 
      print_help
      exit 0
      ;;
    *)
      print_help
      exit 1
      ;;
  esac
done

# Check if all required parameters are provided
if [[ -z "$PROJECT_NAME" || -z "$AMOUNT" || -z "$REGION" ]]; then
  echo "ERROR: Project name, amount, and region are required."
  print_help
  exit 1
fi

# Define the budget JSON template
TEMPLATE='{
    "BudgetLimit": {
        "Amount": "<amount>",
        "Unit": "USD"
    },
    "BudgetName": "<project_name>",
    "BudgetType": "COST",
    "CostFilters": {
        "TagKeyValue": [
            "user:aws-parallelcluster-project$<project_name>"
        ]
    },
    "CostTypes": {
        "IncludeCredit": true,
        "IncludeDiscount": true,
        "IncludeOtherSubscription": true,
        "IncludeRecurring": true,
        "IncludeRefund": true,
        "IncludeSubscription": true,
        "IncludeSupport": true,
        "IncludeTax": true,
        "IncludeUpfront": true,
        "UseBlended": false
    },
    "TimeUnit": "MONTHLY"
}'

# Replace placeholders with input values
BUDGET_JSON=$(echo "$TEMPLATE" | \
    sed "s/<amount>/$AMOUNT/g" | \
    sed "s/<project_name>/$PROJECT_NAME/g")

# Save the JSON to a temporary file
TEMP_JSON=$(mktemp)
echo "$BUDGET_JSON" > "$TEMP_JSON"

ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)

# Create the budget using AWS CLI
echo "Creating budget for project: $PROJECT_NAME with amount: $AMOUNT USD in region: $REGION"
aws budgets create-budget --account-id $ACCOUNT_ID \
    --budget file://"$TEMP_JSON" --region "$REGION"

# Check if the budget creation succeeded
if [[ $? -eq 0 ]]; then
  echo "Successfully created budget for project: $PROJECT_NAME with amount: $AMOUNT USD"
  echo "You will need to add this budget name to the projects file found in the s3:*-omics-analysis/cluster_boot_config/projects_list.conf so it will be included in new clusters. For already running clusters, you will need to edit the file found in /opt/slrum/etc/projects_list.conf as well."
else
  echo "ERROR: Failed to create budget. Please check your AWS credentials, permissions, and region."
fi

# Clean up temporary JSON file
rm "$TEMP_JSON"
