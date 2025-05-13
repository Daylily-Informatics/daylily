# Creating New Cluster

## Setup service user (as admin)
- `rcrf-daylily-service`
- `https://265676937214.signin.aws.amazon.com/console`
- `3c3Fd@y!!`

## Get security credentials and setup ssh key (as service user)
- Credentials (saved to ~/.aws/{credentials|config}) as the profile name `rcrf-daylily-service`
- From `us-west-2` EC2 dashboard, created key named `rcrf-omics-analysis-service-usw2.pem`, and

```bash
mv ~/Downloads/rcrf-omics-analysis-service-usw2.pem ~/.ssh/
chmod 700 ~/.ssh
chmod 400 ~/.ssh/rcrf-omics-analysis-service-usw2.pem
```

## setup an rcrf-daylily-service handle w/rclone

```bash
rclone config
# select new
# enter name rcrf-daylily-service
# choose '4' Amazon S3 for storage type
# choose '1' AWS for storage provider
# choose '1' to enter credential by hand
# enter access key
# enter secret access key
# enter '4' for us-west-2 for region
# leave endpint blank, hit enter
# leave location constraint blank, hit enter
# hit enter to accept default ACL
# hit enter for server side encryption default
# hit enter for sse key default
# hit enter for default service class
# enter 'n' for editing advanced config
# enter 'y' to save remote'
# the remote should be displayed w/any others in the Name col, enter 'q' to exit config

rclone lsd rcrf-daylily-service:  # will list all of the buckets

rclone lsd rcrf-daylily-service:rcrf-omics-analysis-us-west-2 # will list all directories in the bucket

rclone ls rcrf-daylily-service:rcrf-omics-analysis-us-west-2/cluster_boot_config # will recursively list all files i
```

## Create references bucket (on laptop, using rcrf-daylily-service profile)
_if a references bucket already exists in the region, this can be skipped_

```bash

AWS_PROFILE=rcrf-daylily-service
BUCKET_PREFIX=rcrf-ds
REGION=us-west-2

bash ./bin/create_daylily_omics_analysis_s3.sh  --disable-warn --region $REGION --profile $AWS_PROFILE --bucket-prefix $BUCKET_PREFIX

bash ./bin/create_daylily_omics_analysis_s3.sh  --disable-warn --region $REGION --profile $AWS_PROFILE --bucket-prefix $BUCKET_PREFIX --disable-dryrun
```

- Once completed, I copied the sentieon liscence from an existing bucket to the new one.
```bash
 rclone copy rcrf-daylily-service:rcrf-ds-omics-analysis-us-west-2/data/cached_envs/Rare_Cancer_Research_Foundation_eval.lic rcrf-daylily-service:rcrf-omics-analysis-us-west-2/data/cached_envs/
 ```

## Scan regions for spot pricing
_we are going to assume we will run in us-west-2, so lets just look at the AZ pricing there._


Run this command:
```bash
AWS_PROFILE=rcrf-daylily-service
REGION=us-west-2          
OUT_TSV=./init_daylily_cluster.tsv

./bin/check_current_spot_market_by_zones.py -o $OUT_TSV --profile $AWS_PROFILE   

```

- I chose `us-west-2c`, however, pricing in us-west* is higher than I usually see it, with min spot cost ~`$2.80` while `ap-south-1*` pricing is ~`$1.01`.  We do not have reference buckets ready in `ap-south-1`, and moving data there is pretty slow and our quotas are pretty low there.  So, will proceed with us-west-2 for now.
  

## Creating an ephemeral cluster

> For cluster named ** `(rcrf-oa)(-)(usw2)-(jem)`**:
> - $1 = rcrf omics analysis cluster
> - $3 = region-az # I've found helpful to know this at a glance
> - $5 = short descriptor, ie: my initials here

```zsh
AWS_PROFILE=rcrf-daylily-service
REGION_AZ=us-west-2c

source bin/daylily-create-ephemeral-cluster --region-az $REGION_AZ --profile $AWS_PROFILE --pass-on-warn

# SSH KEY Selected
2) rcrf-omics-analysis-service-usw2|ed25519|pem - FOUND|/Users/daylily/.ssh/rcrf-omics-analysis-service-usw2.pem

# S3 Bucket selection
1) rcrf-ds-omics-analysis-us-west-2 (us-west-2)

# Pub Subnet
1) rcrf-omics-analysis Public Subnet    us-west-2c (subnet-050d2932318547278)

# Priv Subnet
1) rcrf-omics-analysis Private Subnet   us-west-2c (subnet-00ae8f786365a423e)

# ARN for 'pclusterTagsAndBudget'
1) arn:aws:iam::265676937214:policy/pclusterTagsAndBudget

# cluster name
rcrf-oa-usw2-jem

# string of allowed user names
hit enter for defaults

# email address for budget alerts
yours@rcrf.org

# Budget alert amount
200

# Enforce budget
1) No

# Path to cluster config.yaml
hit enter for default

# FSX Size
1) 4800

# Enable Detailed Monitoring
1) No

# Delete local vol on term
2) Yes

# Retain/Delete Fsx on cluster delete
2) Delete
```

-The cluster creation script will then go on to calculate spot price caps for the cluster config.yaml which was created from the values you entered above.  This has been saved to `~/.config/daylily/CLUSTERNAME.yaml`.
-If these calculations are successful, the script will attempt a dry-run with the `.yaml`, and if the `dry-run` is successful, move to creating the cluster.  Which will run for ~20min, and block the terminal until completed.  
-Once complete, it will go on to run the headnode config script, and then to run a few tests on the headnode to confirm things are running.


# Parallel Cluster UI

I create a parallel cluster UI and rooted it in `us-west-2`.  You can create these from the [templates found here](https://docs.aws.amazon.com/parallelcluster/latest/ug/install-pcui-v3.html). They should cost a few $ a month to run, so we might just want to leave the one I spun up in place. When creating, you need to enter the admin email (who gets an email with credentials to login the first time), the VPC to run in, and a subnet to run it in.  This subnet should ideally not be a VPC or subnet created by parallel cluster, as they get torn down and re-created often.  I chose to run this in the VPC and public subnet dewey is in presently.  You have each been emailed login credentials, please login and use them to gain access to the PCUI.  This interface is a really nice way to see at a high level what is going on with various clusters.

Our current PCUI is here: [https://c8pfjtp4zf.execute-api.us-west-2.amazonaws.com/pcui/clusters](https://c8pfjtp4zf.execute-api.us-west-2.amazonaws.com/pcui/clusters).

## Open a terminal on the headnode and become `ubuntu`
- Login to the PCUI.
- Select the cluster.
- Click `shell`, which will open a web-terminal in the browser with a weird SSM shell.
- Login to ubuntu user `udo su - ubuntu`.
- As ubuntu user, init bash and open tmux (it's helpful to immediately go to tmux, b/c the shells via ssm have some weirdnesses)
  ```
  bash
  tmux new -s blah
  ```


# Trigger the remote test run
_this command is printed for you at the end of successful cluster creation_

`source bin/daylily-run-ephemeral-cluster-remote-tests /Users/daylily/.ssh/rcrf-omics-analysis-service-usw2.pem us-west-2  rcrf-daylily-service`

- Run this, and enter `y` to run this in a tmux session you'll be able to connect to on the tmux server.

**ssh into the headnode and check the tmux session** _or_ **use PCUI to open a shell on the headnode**
Then monitor the job (which will likely be busy building some cache respources if this is a new cluster)

```bash
tmux ls
tmux a -t cluster_t
```

# Running the WGS test pipeline on HG001 and HG007

```bash
ssh -i /Users/daylily/.ssh/rcrf-omics-analysis-service-usw2.pem ubuntu@50.112.171.179

mkdir -p /fsx/analysis_results/ubuntu/controls
cd cd /fsx/analysis_results/ubuntu/controls

git clone git@github.com:Daylily-Informatics/daylily.git
cd daylily

cp .test_data/data/giab_30x_hg38_analysis_manifest.csv config/analysis_manifest.csv
emacs config/analysis_manifest.csv # remove all but the HG001 and HG007 sample rows, leave header row

. dyinit

dy-a slurm hg38

# dry run
dy-r produce_multiqc_final_wgs produce_snv_concordances -p -k -j 40 --config genome_build=hg38 aligners=['bwa2a','sent','strobe']  dedupers=['dppl']  snv_callers=['oct','clair3','deep','sentd','lfq2'] sv_callers=['manta','tiddit','dysgu'] -n

# real deal
dy-r produce_multiqc_final_wgs produce_snv_concordances -p -k -j 40 --config genome_build=hg38 aligners=['bwa2a','sent','strobe']  dedupers=['dppl']  snv_callers=['oct','clair3','deep','sentd','lfq2'] sv_callers=['manta','tiddit','dysgu']

# monitor, from another tmux
cd /fsx/analysis_results/ubuntu/controls/daylily
. dyinit

watch $sq_cmd
```


# Running the WGS pipeline on our 2 patient samples

## Createa an analysis_manifest.csv from a sample_sheet.tsv

- Add `rcrf-daylily-service` to the headnode `~/.aws/credentials` and `~/.aws/config`.
- make a copy of the `etc/analysis_samples_template.tsv` to `./etc/my_analysis.tsv`
- Hack the new tsv to reflect your data.

Run using your new samplesheet/manifest:

```bash
python bin/daylily-analysis-samples-to-manifest.py local ./etc/my_analysis.tsv rcrf-daylily-service
cp  ./R0_analysis_manifest.csv config/analysis_manifest.csv
dy-a local hg38
dy-r produce_snv_concordances -p -k -j 1 --config aligners=[\"bwa2a\"] dedupers=[\"dppl\"] genome_build=\"hg38\" snv_callers=[\"deep\"] -n
```


# S3 Ref bucket policy modification

{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Principal": {
                "AWS": "arn:aws:iam::265676937214:root"
            },
            "Action": [
                "s3:GetBucketLocation",
                "s3:ListBucket",
                "s3:GetObject"
            ],
            "Resource": [
                "arn:aws:s3:::rcrf-ds-omics-analysis-us-west-2",
                "arn:aws:s3:::rcrf-ds-omics-analysis-us-west-2/*"
            ]
        }
    ]
}
