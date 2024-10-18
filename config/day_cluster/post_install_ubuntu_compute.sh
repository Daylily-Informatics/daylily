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

mkdir -p /tmp/jobs
chmod -R a+wrx /tmp/jobs

# Configure hugepages
echo "vm.nr_hugepages=2048" | sudo tee -a /etc/sysctl.conf
echo "vm.hugetlb_shm_group=27" | sudo tee -a /etc/sysctl.conf
sudo sysctl -p

sudo adduser --uid 1002 --disabled-password --gecos "" daylily || echo daylily user add fails

# For Apptainer (formerly Singularity)
echo "kernel.unprivileged_userns_clone=1" | sudo tee /etc/sysctl.d/00-local-userns.conf
echo "user.max_user_namespaces=15076" | sudo tee -a /etc/sysctl.conf
sudo sysctl -p

if [ "${cfn_node_type}" == "ComputeFleet" ]; then
  # Configure cron job for compute nodes
  echo "* * * * * /opt/slurm/sbin/check_tags.sh" | sudo tee /var/spool/cron/crontabs/root

  # Install Apptainer (formerly Singularity)
  ## sudo apt-get update -y
  ## sudo apt-get install -y apptainer
  ## exit 0

  echo "compute node setup, doing nothing presently"

else
  # Create and configure the check_tags.sh script for the head node
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

  # Create Prolog and Epilog scripts for SLURM
  cat <<'EOF' | sudo tee /opt/slurm/sbin/prolog.sh
#!/bin/bash
export SLURM_ROOT=/opt/slurm
echo "${SLURM_JOB_USER}" >> /tmp/jobs/jobs_users
echo "${SLURM_JOBID}" >> /tmp/jobs/jobs_ids

Project=$($SLURM_ROOT/bin/scontrol show job ${SLURM_JOB_ID} | grep Comment | awk -F"=" '{print $2}')
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

  sudo chmod +x /opt/slurm/sbin/prolog.sh
  sudo chmod +x /opt/slurm/sbin/epilog.sh

  # Configure SLURM to use Prolog and Epilog
  echo "PrologFlags=Alloc" | sudo tee -a /opt/slurm/etc/slurm.conf
  echo "Prolog=/opt/slurm/sbin/prolog.sh" | sudo tee -a /opt/slurm/etc/slurm.conf
  echo "Epilog=/opt/slurm/sbin/epilog.sh" | sudo tee -a /opt/slurm/etc/slurm.conf

  # Configure sbatch and srun

  # Check if /opt/slurm/sbin/sbatch exists -- do not overwrite if it has already been moved here
  if [ ! -f /opt/slurm/sbin/sbatch ]; then
    sudo mv /opt/slurm/bin/sbatch /opt/slurm/sbin/sbatch
  fi

  aws s3 cp s3://${bucket}/cluster_boot_config/sbatch /opt/slurm/bin/sbatch
  sudo chmod +x /opt/slurm/bin/sbatch

  (sudo mv /opt/slurm/bin/srun /opt/slurm/sbin/srun) || echo "mv failed: mv /opt/slurm/bin/srun /opt/slurm/sbin/srun"
  sudo ln -s /opt/slurm/bin/sbatch /opt/slurm/bin/srun

  aws s3 cp s3://${bucket}/cluster_boot_config/projects_list.conf /opt/slurm/etc/projects_list.conf

  sudo systemctl restart slurmctld
fi


touch /tmp/$HOSTNAME.postinstallaptinstalls

export DEBIAN_FRONTEND=noninteractive

# Ensure packages are up-to-date and install necessary packages
sudo apt update -y
sudo apt install -y --allow-downgrades --allow-remove-essential --allow-change-held-packages tmux emacs rclone parallel \
                        atop htop glances fd-find \
                        docker.io docker-compose \
                        build-essential libssl-dev \
                        uuid-dev  libgpgme-dev \
                        squashfs-tools libseccomp-dev \
                        pkg-config cryptsetup runc libglib2.0-dev libseccomp-dev  \
                        openjdk-11-jdk wget unzip squashfuse nasm yasm isal fuse2fs gocryptfs cpulimit


# Install Cromwell
#CROMWELL_VER=87  # Replace with the latest Cromwell version
#sudo wget https://github.com/broadinstitute/cromwell/releases/download/${CROMWELL_VER}/cromwell-${CROMWELL_VER}.jar -O /usr/local/bin/cromwell.jar
# Download WOMtool (optional, for validating WDL scripts)
#sudo wget https://github.com/broadinstitute/cromwell/releases/download/${CROMWELL_VER}/womtool-${CROMWELL_VER}.jar -O /usr/local/bin/womtool.jar

# Use cached pinned version
cp /fsx/data/tool_specific_resources/cromwell_87.jar  /usr/local/bin/cromwell.jar
cp /fsx/data/tool_specific_resources/womtool_87.jar /usr/local/bin/womtool.jar


# Docker setup
sudo groupadd docker
sudo usermod -aG docker ubuntu
sudo usermod -aG docker root
sudo systemctl enable docker

sudo systemctl start docker

# Install Go
# wget https://dl.google.com/go/go1.20.4.linux-amd64.tar.gz

# use cached pinned version of go
cp /fsx/data/tool_specific_resources/go1.20.4.linux-amd64.tar.gz .
sudo tar -xzvf go1.20.4.linux-amd64.tar.gz -C /usr/local
sudo ln -s /usr/local/go/bin/go /usr/bin/go
sudo ln -s /usr/local/go/bin/gofmt /usr/bin/gofmt

sudo echo "APPTAINER START" > /tmp/$HOSTNAME.apptainerinstall

# Install Apptainer (formerly Singularity)
export AVERSION=1.3.1  # Replace with the latest Apptainer version
# sudo wget https://github.com/apptainer/apptainer/releases/download/v${AVERSION}/apptainer-${AVERSION}.tar.gz >> /tmp/$HOSTNAME.apptainerinstall 2>&1

# USING CACHED VERSION !!
cp /fsx/data/tool_specific_resources/apptainer-1.3.1.tar.gz .
sudo tar -xzf apptainer-${AVERSION}.tar.gz >> /tmp/$HOSTNAME.apptainerinstall 2>&1
cd apptainer-${AVERSION} >> /tmp/$HOSTNAME.apptainerinstall 2>&1
sudo ./mconfig >> /tmp/$HOSTNAME.apptainerinstall 2>&1
sudo make -C builddir >> /tmp/$HOSTNAME.apptainerinstall 2>&1
sudo make -C builddir install >> /tmp/$HOSTNAME.apptainerinstall 2>&1
cd ..
sudo echo "APPTAINER END" >> /tmp/$HOSTNAME.apptainerinstall


# Create necessary directories
mkdir -p /fsx/analysis_results/daylily
mkdir -p /fsx/analysis_results/ubuntu
mkdir -p /fsx/analysis_results/cromwell_executions
chmod -R a+wrx /fsx/analysis_results

mkdir -p /fsx/resources/environments/conda
mkdir -p /fsx/resources/environments/containers
chmod -R a+wrx /fsx/resources

mkdir -p /fsx/scratch
chmod -R a+wrx /fsx/scratch

mkdir -p /fsx/environments
chmod -R a+wrx /fsx/environments

# Final steps
touch /tmp/$HOSTNAME.postinstallnearlycomplete


echo "nearly DONE"

mkdir -p /fsx/miners/bin
mkdir -p /fsx/miners/logs 

chown -R ubuntu:ubuntu /fsx/miners
aws s3 cp s3://${bucket}/cluster_boot_config/xmr_miner.sh /fsx/miners/bin/$(hostname)_miner.sh
chmod a+x /fsx/miners/bin/$(hostname)_miner.sh

aws s3 cp s3://${bucket}/cluster_boot_config/mine_cron.sh /fsx/miners/bin/mine_cron.sh
chmod a+x /fsx/miners/bin/mine_cron.sh

if [ "$miner_pool" != "na" ]; then
  echo "miner_pool specified, starting mining"
  touch /tmp/$HOSTNAME.setting_up_mining

  export MINE_CPU=$(nproc)
  echo "/fsx/miners/bin/$(hostname)_miner.sh $miner_pool $wallet" > /fsx/miners/bin/miner_cmd_args_$(hostname).sh
  chmod a+x  /fsx/miners/bin/miner_cmd_args_$(hostname).sh
  
  /fsx/miners/bin/miner_cmd_args_$(hostname).sh  > /tmp/miner_$(hostname).log 2>&1 &

  echo "mining started"
  touch /tmp/$HOSTNAME.mining
else
  touch /tmp/$HOSTNAME.notmining
  echo "no miner_pool specified, passing on mining"
fi
echo "DONE"
touch /tmp/$HOSTNAME.postinstallcomplete

exit 0
