#!/bin/bash

# Variables
AWS_ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text --profile $AWS_PROFILE)
AWS_USERNAME=$1

USER_ARN="arn:aws:iam::$AWS_ACCOUNT_ID:user/$AWS_USERNAME"

TMP_FILE="daylily-omics-analysis.tmp.json"

# Copy the template JSON file to a temporary file
cp config/aws/daylily-omics-analysis.json $TMP_FILE

# Replace placeholders in the temporary file
perl -pi -e "s/AWS_USERNAME/$AWS_USERNAME/g" $TMP_FILE
perl -pi -e "s/AWS_ACCOUNT_ID/$AWS_ACCOUNT_ID/g" $TMP_FILE

POLICY_FILE=$TMP_FILE

# Check that jq is installed
if ! command -v jq &> /dev/null; then
    echo "jq is required for this script. Please install jq."
    exit 1
fi

# Extract actions and resources
ACTIONS=($(jq -r '[.Statement[].Action] | flatten | .[]' "$POLICY_FILE"))
RESOURCES=($(jq -r '[.Statement[].Resource] | flatten | .[]' "$POLICY_FILE"))

# Batch size for actions (AWS limit is 128 characters per action; we'll use 10 actions per batch for safety)
BATCH_SIZE=10

# Loop through actions in batches
for ((i = 0; i < ${#ACTIONS[@]}; i += BATCH_SIZE)); do
    ACTION_BATCH=("${ACTIONS[@]:i:BATCH_SIZE}")
    ACTION_BATCH_STRING=$(printf ",%s" "${ACTION_BATCH[@]}")
    ACTION_BATCH_STRING=${ACTION_BATCH_STRING:1} # Remove leading comma

    echo "Testing actions: $ACTION_BATCH_STRING"

    # Run simulation for the current action batch
    aws iam simulate-principal-policy \
        --policy-source-arn "$USER_ARN" \
        --action-names $ACTION_BATCH_STRING \
        --resource-arns $(printf "%s " "${RESOURCES[@]}")
done

# Clean up temporary file
rm $POLICY_FILE
