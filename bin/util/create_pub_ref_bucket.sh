


aws s3api create-bucket --bucket daylily-references-public --region us-west-2 --create-bucket-configuration LocationConstraint=us-west-2 --acl private

aws s3api put-public-access-block --bucket daylily-references-public --public-access-block-configuration BlockPublicAcls=false,IgnorePublicAcls=false,BlockPublicPolicy=false,RestrictPublicBuckets=false


aws s3api put-bucket-policy --bucket daylily-references-public --policy file://bin/util/user_pays.json

aws s3api put-bucket-request-payment --bucket daylily-references-public --request-payment-configuration Payer=Requester

# aws s3 cp s3://daylily-references-private2/data/  s3://daylily-references-public/



