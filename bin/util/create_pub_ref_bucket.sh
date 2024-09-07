
bucket_name="daylily-references-public"
region="us-west-2"
bucket_user_pays_policy="user_pays.json"

aws s3api create-bucket --bucket $bucket_name --region $region --acl private

aws s3api put-bucket-policy --bucket daylily-references-public --policy file://user_pays.json

aws s3api put-bucket-request-payment --bucket daylily-references-public --request-payment-configuration Payer=Requester

# aws s3 cp s3://daylily-references-private2/data/  s3://daylily-references-public/



