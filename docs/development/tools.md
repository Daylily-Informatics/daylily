# Pcluster

## Head Node

### Run Test Job
`aws s3 cp s3://daylily-references/cluster_boot_config/sleep_test.sh - | sbatch --comment RandD --partition i4-5`