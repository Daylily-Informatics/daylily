#!/bin/bash

# MIT No Attribution
# Copyright 2020 Amazon.com, Inc. or its affiliates. All Rights Reserved.

# The script configures the Slurm cluster after the deployment.

touch /tmp/$HOSTNAME.postinstallBEGIN

. "/etc/parallelcluster/cfnconfig"

region="$1"
bucket="$2"  # specified in the cluster yaml, bucket-name, no s3:// prefix
miner_pool="$3"  # specified in the cluster yaml, miner_pool
wallet="$4"  # specified in the cluster yaml, wallet

aws configure set region $region


# Function to log spot price
log_spot_price() {

  TOKEN=$(curl -X PUT 'http://169.254.169.254/latest/api/token' -H 'X-aws-ec2-metadata-token-ttl-seconds: 21600')
  instance_type=$(curl -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-type)
  spot_price=$(aws ec2 describe-spot-price-history \
    --instance-types "$instance_type" \
    --region "$region" \
    --product-description "Linux/UNIX" \
    --query 'SpotPriceHistory[0].SpotPrice' \
    --output text)
  echo "Spot price for instance type $instance_type: $spot_price USD/hour" >> /fsx/scratch/$HOSTNAME_spot_price.log
}

mkdir -p /tmp/jobs
chmod -R a+wrx /tmp/jobs

# Configure hugepages and namespaces (common to both head and compute nodes)
echo "vm.nr_hugepages=2048" | sudo tee -a /etc/sysctl.conf
echo "vm.hugetlb_shm_group=27" | sudo tee -a /etc/sysctl.conf
echo "kernel.unprivileged_userns_clone=1" | sudo tee /etc/sysctl.d/00-local-userns.conf
echo "user.max_user_namespaces=15076" | sudo tee -a /etc/sysctl.conf
sudo sysctl -p

# Ensure the user exists
sudo adduser --uid 1002 --disabled-password --gecos "" daylily || echo "daylily user add failed"

# Compute Node Configuration
if [ "${cfn_node_type}" == "ComputeFleet" ]; then
  echo "Configuring compute node..."

  # Set up cron job for monitoring
  echo "* * * * * /opt/slurm/sbin/check_tags.sh" | sudo tee /var/spool/cron/crontabs/root

  log_spot_price

else
  # Head Node Configuration
  echo "Configuring head node..."

  # Log the spot price
  log_spot_price

  # Create and configure the check_tags.sh script for head node
  cat <<'EOF' | sudo tee /opt/slurm/sbin/check_tags.sh
#!/bin/bash
source /etc/profile


TOKEN=$(curl -X PUT 'http://169.254.169.254/latest/api/token' -H 'X-aws-ec2-metadata-token-ttl-seconds: 21600')
itype=$(curl -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-type)

update=0
tag_userid=""
tag_jobid=""
tag_project=""

if [ ! -f /tmp/jobs/jobs_users ] || [ ! -f /tmp/jobs/jobs_ids ]; then
  exit 0
fi

active_users=$(cat /tmp/jobs/jobs_users | sort | uniq)
active_jobs=$(cat /tmp/jobs/jobs_ids | sort)
echo $active_users > /tmp/jobs/tmp_jobs_users
echo $active_jobs > /tmp/jobs/tmp_jobs_ids

if [ -f /tmp/jobs/jobs_projects ]; then
  active_projects=$(cat /tmp/jobs/jobs_projects | sort | uniq)
  echo $active_projects > /tmp/jobs/tmp_jobs_projects
fi

if [ ! -f /tmp/jobs/tag_userid ] || [ ! -f /tmp/jobs/tag_jobid ]; then
  echo $active_users > /tmp/jobs/tag_userid
  echo $active_jobs > /tmp/jobs/tag_jobid
  echo $active_projects > /tmp/jobs/tag_project
  update=1
else
  active_users=$(cat /tmp/jobs/tmp_jobs_users)
  active_jobs=$(cat /tmp/jobs/tmp_jobs_ids)
  if [ -f /tmp/jobs/tmp_jobs_projects ]; then
    active_projects=$(cat /tmp/jobs/tmp_jobs_projects)
  fi

  tag_userid=$(cat /tmp/jobs/tag_userid)
  tag_jobid=$(cat /tmp/jobs/tag_jobid)
  if [ -f /tmp/jobs/tag_project ]; then
    tag_project=$(cat /tmp/jobs/tag_project)
  fi

  if [ "${active_users}" != "${tag_userid}" ]; then
    tag_userid="${active_users}"
    echo ${tag_userid} > /tmp/jobs/tag_userid
    update=1
  fi

  if [ "${active_jobs}" != "${tag_jobid}" ]; then
    tag_jobid="${active_jobs}"
    echo ${tag_jobid} > /tmp/jobs/tag_jobid
    update=1
  fi

  if [ "${active_projects}" != "${tag_project}" ]; then
    tag_project="${active_projects}"
    echo ${tag_project} > /tmp/jobs/tag_project
    update=1
  fi
fi

if [ ${update} -eq 1 ]; then
  TOKEN=$(curl -X PUT "http://169.254.169.254/latest/api/token" -H "X-aws-ec2-metadata-token-ttl-seconds: 21600")
  MyInstID=$(curl -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)
  aws ec2 create-tags --resources ${MyInstID} --tags Key=aws-parallelcluster-username,Value="${tag_userid}" --region ${region}
  aws ec2 create-tags --resources ${MyInstID} --tags Key=aws-parallelcluster-jobid,Value="${tag_jobid}" --region ${region}
  aws ec2 create-tags --resources ${MyInstID} --tags Key=aws-parallelcluster-project,Value="${tag_project}" --region ${region}
fi
EOF  

  sudo chmod +x /opt/slurm/sbin/check_tags.sh

  # Prolog and Epilog Scripts for SLURM (Head Node Only)
  cat <<'EOF' | sudo tee /opt/slurm/sbin/prolog.sh
#!/bin/bash
export SLURM_ROOT=/opt/slurm
echo "${SLURM_JOB_USER}" >> /tmp/jobs/jobs_users
echo "${SLURM_JOBID}" >> /tmp/jobs/jobs_ids


Project=$($SLURM_ROOT/bin/scontrol show job ${SLURM_JOBID} | grep Comment | awk -F"=" '{print $2}')
if [ ! -z "${Project}" ]; then
  echo "${Project}" >> /tmp/jobs/jobs_projects
fi
EOF

  cat <<'EOF' | sudo tee /opt/slurm/sbin/epilog.sh
#!/bin/bash
export SLURM_ROOT=/opt/slurm
sed -i "0,/${SLURM_JOB_USER}/d" /tmp/jobs/jobs_users
sed -i "0,/${SLURM_JOBID}/d" /tmp/jobs/jobs_ids

Project=$($SLURM_ROOT/bin/scontrol show job ${SLURM_JOB_ID} | grep Comment | awk -F"=" '{print $2}')
if [ ! -z "${Project}" ]; then
  sed -i "0,/${Project}/d" /tmp/jobs/jobs_projects
fi
EOF

  sudo chmod +x /opt/slurm/sbin/prolog.sh /opt/slurm/sbin/epilog.sh

  # Configure SLURM with Prolog and Epilog
  echo "PrologFlags=Alloc" | sudo tee -a /opt/slurm/etc/slurm.conf
  echo "Prolog=/opt/slurm/sbin/prolog.sh" | sudo tee -a /opt/slurm/etc/slurm.conf
  echo "Epilog=/opt/slurm/sbin/epilog.sh" | sudo tee -a /opt/slurm/etc/slurm.conf


  aws s3 cp s3://${bucket}/cluster_boot_config/sbatch /opt/slurm/bin/sbatch
  sudo chmod +x /opt/slurm/bin/sbatch

  (sudo mv /opt/slurm/bin/srun /opt/slurm/sbin/srun) || echo "mv failed: mv /opt/slurm/bin/srun /opt/slurm/sbin/srun"
  sudo ln -s /opt/slurm/bin/sbatch /opt/slurm/bin/srun

  aws s3 cp s3://${bucket}/cluster_boot_config/projects_list.conf /opt/slurm/etc/projects_list.conf

  # Restart SLURM Controller
  sudo systemctl restart slurmctld
  touch /tmp/$HOSTNAME.postslurmcfg

fi

# Common Steps for All Nodes

# Update and install necessary packages
export DEBIAN_FRONTEND=noninteractive
sudo apt update -y
sudo apt install -y tmux emacs rclone parallel atop htop glances fd-find docker.io \
                    build-essential libssl-dev uuid-dev libgpgme-dev squashfs-tools \
                    libseccomp-dev pkg-config openjdk-11-jdk wget unzip nasm yasm isal \
                    fuse2fs gocryptfs cpulimit

# Docker setup
sudo groupadd docker
sudo usermod -aG docker ubuntu root
sudo systemctl enable docker && sudo systemctl start docker

# Install Apptainer
export AVERSION=1.3.1
cp /fsx/data/tool_specific_resources/apptainer-${AVERSION}.tar.gz .
tar -xzf apptainer-${AVERSION}.tar.gz && cd apptainer-${AVERSION}
sudo ./mconfig && sudo make -C builddir && sudo make -C builddir install

# Install Cromwell and Go (using cached versions)
cp /fsx/data/tool_specific_resources/cromwell_87.jar /usr/local/bin/cromwell.jar
cp /fsx/data/tool_specific_resources/womtool_87.jar /usr/local/bin/womtool.jar
cp /fsx/data/tool_specific_resources/go1.20.4.linux-amd64.tar.gz .
sudo tar -xzvf go1.20.4.linux-amd64.tar.gz -C /usr/local

sudo ln -s /usr/local/go/bin/{go,gofmt} /usr/bin/

# Create directories and set permissions
for dir in /fsx/{tmp,analysis_results/{daylily,ubuntu,cromwell_executions},resources/{environments/{conda,containers}},scratch,miners/bin,miners/logs}; do
  mkdir -p $dir && chmod -R a+wrx $dir
done

# Mining Setup (if miner_pool is specified)
if [ "$miner_pool" != "na" ]; then
  echo "Starting mining..."
  /fsx/miners/bin/$(hostname)_miner.sh "$miner_pool" "$wallet" &
else
  echo "No miner pool specified, skipping mining."
fi

# Custom Shutdown Script (Head Node Only)
if [ "${cfn_node_type}" != "ComputeFleet" ]; then
  cat <<'EOF' | sudo tee /etc/systemd/system/custom-shutdown.service
[Unit]
Description=Custom Shutdown Script
Before=shutdown.target reboot.target

[Service]
Type=oneshot
ExecStart=/opt/slurm/etc/shutdown_script_ubuntu_head.sh
RemainAfterExit=true

[Install]
WantedBy=halt.target reboot.target
EOF

  sudo systemctl daemon-reload && sudo systemctl enable custom-shutdown.service
fi

# Copy cached data from S3
cp -r /fsx/data/cached_envs/conda/* /fsx/resources/environments/conda/ubuntu/$USER/$(hostname)/
cp -r /fsx/data/cached_envs/containers/* /fsx/resources/environments/containers/$USER/$(hostname)/


# Finalization
touch /tmp/$HOSTNAME.postinstallcomplete
echo "Post-installation complete."
exit 0
