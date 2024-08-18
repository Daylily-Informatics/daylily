#!/bin/bash

# MIT No Attribution
# Copyright 2020 Amazon.com, Inc. or its affiliates. All Rights Reserved.

# The script configures the Slurm cluster after the deployment.

. "/etc/parallelcluster/cfnconfig"

bucket="$1"  # specified in the cluster yaml, bucket-name, no s3:// prefix

mkdir -p /tmp/jobs
chmod -R a+wrx /tmp/jobs

# Configure hugepages
echo "vm.nr_hugepages=2048" | sudo tee -a /etc/sysctl.conf
echo "vm.hugetlb_shm_group=27" | sudo tee -a /etc/sysctl.conf
sudo sysctl -p

# For Apptainer (formerly Singularity)
sudo adduser --uid 1002 --disabled-password --gecos "" daylily || echo daylily user add fails
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

region=$(curl -s http://169.254.169.254/latest/meta-data/placement/region)
aws configure set region $region

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
  MyInstID=$(curl -s http://169.254.169.254/latest/meta-data/instance-id)
  aws ec2 create-tags --resources ${MyInstID} --tags Key=aws-parallelcluster-username,Value="${tag_userid}"
  aws ec2 create-tags --resources ${MyInstID} --tags Key=aws-parallelcluster-jobid,Value="${tag_jobid}"
  aws ec2 create-tags --resources ${MyInstID} --tags Key=aws-parallelcluster-project,Value="${tag_project}"
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

export DEBIAN_FRONTEND=noninteractive

# Ensure packages are up-to-date and install necessary packages
sudo apt-get update -y
sudo apt-get install -y --allow-downgrades --allow-remove-essential --allow-change-held-packages tmux emacs rclone parallel \
                        atop htop glances fd-find \
                        docker.io docker-compose \
                        build-essential libssl-dev \
                        uuid-dev  libgpgme-dev \
                        squashfs-tools libseccomp-dev \
                        pkg-config cryptsetup runc libglib2.0-dev libseccomp-dev 


# Docker setup
sudo groupadd docker
sudo usermod -aG docker ubuntu
sudo usermod -aG docker root
sudo systemctl enable docker

sudo systemctl start docker

# Install Go
wget https://dl.google.com/go/go1.17.7.linux-amd64.tar.gz
sudo tar -xzvf go1.17.7.linux-amd64.tar.gz -C /usr/local
sudo ln -s /usr/local/go/bin/go /usr/bin/go
sudo ln -s /usr/local/go/bin/gofmt /usr/bin/gofmt

# Install Singularity
git clone --recurse-submodules https://github.com/sylabs/singularity.git
cd singularity
git checkout --recurse-submodules v3.10.0
./mconfig --prefix=/opt/singularity
cd builddir
make
sudo make install

# Create necessary directories
mkdir -p /fsx/analysis_results/daylily
mkdir -p /fsx/analysis_results/ubuntu
chmod -R a+wrx /fsx/analysis_results

mkdir -p /fsx/resources/environments/conda
mkdir -p /fsx/resources/environments/containers
chmod -R a+wrx /fsx/resources

mkdir -p /fsx/environments
chmod -R a+wrx /fsx/environments

# Final steps
touch /tmp/$HOSTNAME.postinstallcomplete
echo "DONE"
exit 0
