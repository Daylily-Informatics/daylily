#!/usr/bin/env bash

# Function to display usage information
print_help() {
  echo "Usage: $0 -p <project_name> -a <amount> -r <region> -e <email> -t <thresholds> [-h]"
  echo ""
  echo "Options:"
  echo "  -p    Project name (required)"
  echo "  -a    Budget amount in USD (required)"
  echo "  -r    AWS region (required)"
  echo "  -e    Email address for budget alerts (required)"
  echo "  -t    Comma-separated alert thresholds (required, e.g., 50,80,100)"
  echo "  -h    Display this help message"
  echo ""
  echo "Example:"
  echo "  $0 -p my-project -a 5000 -r us-west-2 -e user@example.com -t 50,80,100"
}

# Parse command-line options
while getopts "p:a:r:e:t:h" opt; do
  case ${opt} in
    p) PROJECT_NAME=$OPTARG ;;
    a) AMOUNT=$OPTARG ;;
    r) REGION=$OPTARG ;;
    e) EMAIL=$OPTARG ;;
    t) THRESHOLDS=$OPTARG ;;
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
if [[ -z "$PROJECT_NAME" || -z "$AMOUNT" || -z "$REGION" || -z "$EMAIL" || -z "$THRESHOLDS" ]]; then
  echo "ERROR: Project name, amount, region, email, and thresholds are required."
  print_help
  exit 1
fi

# Define the budget JSON template
BUDGET_TEMPLATE='{
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
BUDGET_JSON=$(echo "$BUDGET_TEMPLATE" | \
    sed "s/<amount>/$AMOUNT/g" | \
    sed "s/<project_name>/$PROJECT_NAME/g")

# Save budget JSON to a temporary file
TEMP_BUDGET_JSON=$(mktemp)
echo "$BUDGET_JSON" > "$TEMP_BUDGET_JSON"

ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)

# Create the budget using AWS CLI
echo "Creating budget for project: $PROJECT_NAME with amount: $AMOUNT USD in region: $REGION"
aws budgets create-budget --account-id "$ACCOUNT_ID" \
    --budget file://"$TEMP_BUDGET_JSON" --region "$REGION"

# Check if the budget creation succeeded
if [[ $? -eq 0 ]]; then
  echo "Successfully created budget for project: $PROJECT_NAME"

  # Create budget alerts using the correct AWS CLI syntax
  echo "Creating budget alerts..."
  IFS=',' read -ra THRESHOLD_ARRAY <<< "$THRESHOLDS"
  for threshold in "${THRESHOLD_ARRAY[@]}"; do
    echo P:$PROJECT_NAME T:$threshold E:$EMAIL R:$REGION A:$AMOUNT AI:$ACCOUNT_ID 
    aws budgets create-notification --account-id "$ACCOUNT_ID" \
      --budget-name "$PROJECT_NAME" \
      --notification '{"ComparisonOperator":"GREATER_THAN","NotificationType":"ACTUAL","Threshold":'"$threshold"',"ThresholdType":"PERCENTAGE"}' \
      --region "$REGION"   --subscribers '[{"Address":"'"$EMAIL"'","SubscriptionType":"EMAIL"}]' 

    if [[ $? -eq 0 ]]; then
      echo "Successfully created alert for threshold: $threshold%"
    else
      echo "ERROR: Failed to create alert for threshold: $threshold%. Please check your AWS credentials and permissions."
    fi
  done

  echo "You will need to add this budget name to the projects file found in the s3:*-omics-analysis/cluster_boot_config/projects_list.conf so it will be included in new clusters. For already running clusters, you will need to edit the file found in /opt/slurm/etc/projects_list.conf as well."
else
  echo "ERROR: Failed to create budget. Please check your AWS credentials, permissions, and region."
fi

# Clean up temporary JSON file
rm "$TEMP_BUDGET_JSON"
