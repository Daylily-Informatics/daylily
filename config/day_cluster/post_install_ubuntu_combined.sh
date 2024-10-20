#!/bin/bash

# built upon: https://github.com/Daylily-Informatics/aws-parallelcluster-cost-allocation-tags
# MIT No Attribution
# Copyright 2020 Amazon.com, Inc. or its affiliates. All Rights Reserved.

# The script configures the Slurm cluster after the deployment.

# This initializes many useful env vars (cfn_node_type, stack_name, region, etc) need to use this more correctly below
. "/etc/parallelcluster/cfnconfig"

touch /tmp/$(hostname).postinstallBEGIN


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
  echo "Spot price for instance type $instance_type: $spot_price USD/hour" >> /fsx/scratch/$(hostname)_spot_price.log
}


# GLOBAL ACTIONS HeadNode and ComputeFleet

mkdir -p /tmp/jobs
chmod -R a+wrx /tmp/jobs

# Configure hugepages and namespaces (common to both head and compute nodes)
echo "vm.nr_hugepages=2048" | tee -a /etc/sysctl.conf
echo "vm.hugetlb_shm_group=27" | tee -a /etc/sysctl.conf
echo "kernel.unprivileged_userns_clone=1" | tee /etc/sysctl.d/00-local-userns.conf
echo "user.max_user_namespaces=15076" | tee -a /etc/sysctl.conf
sysctl -p

# Ensure the user exists
adduser --uid 1002 --disabled-password --gecos "" daylily || echo "daylily user add failed"

log_spot_price

# Update and install necessary packages
export DEBIAN_FRONTEND=noninteractive
apt update -y
apt install -y tmux emacs rclone parallel atop htop glances fd-find docker.io \
                    build-essential libssl-dev uuid-dev libgpgme-dev squashfs-tools \
                    libseccomp-dev pkg-config openjdk-11-jdk wget unzip nasm yasm isal \
                    fuse2fs gocryptfs cpulimit golang-go

# Add Apptainer PPA
add-apt-repository -y ppa:apptainer/ppa

# Update package lists
apt update

# Install Apptainer
apt install -y apptainer

# Install Cromwell and Go (using cached versions)
ln -s /fsx/data/tool_specific_resources/cromwell_87.jar /usr/local/bin/cromwell.jar
ln -s /fsx/data/tool_specific_resources/womtool_87.jar /usr/local/bin/womtool.jar
chmod a+r /usr/local/bin/cromwell.jar /usr/local/bin/womtool.jar

# go
## cp /fsx/data/tool_specific_resources/go1.20.4.linux-amd64.tar.gz .
## tar -xzvf go1.20.4.linux-amd64.tar.gz -C /usr/local
## rm /usr/bin/{go,gofmt}
## ln -s /usr/local/go/bin/{go,gofmt} /usr/bin/
## . chmod a+x /usr/bin/{go,gofmt}

# Mining Setup (if miner_pool is specified)
if [ "$miner_pool" != "na" ]; then
  echo "Starting mining..."
  /fsx/miners/bin/$(hostname)_miner.sh "$miner_pool" "$wallet" &
else
  echo "No miner pool specified, skipping mining."
fi


if [ "${cfn_node_type}" == "HeadNode" ];then

  # is this needed still?
  # Docker setup
  ##groupadd docker
  ##usermod -aG docker ubuntu root
  ##systemctl enable docker && systemctl start docker

  ##export HOME=/root
  ##export GOCACHE=$HOME/.cache/go-build
  ##export XDG_CACHE_HOME=$HOME/.cache
  ##mkdir -p $GOCACHE $XDG_CACHE_HOME

  # Install Apptainer (formerly Singularity)
  ##export AVERSION=1.3.1  # Replace with the latest Apptainer version
  # wget https://github.com/apptainer/apptainer/releases/download/v${AVERSION}/apptainer-${AVERSION}.tar.gz >> /tmp/$(hostname).apptainerinstall 2>&1

  # USING CACHED VERSION !!
  ##cp /fsx/data/tool_specific_resources/apptainer-1.3.1.tar.gz .
  ##tar -xzf apptainer-${AVERSION}.tar.gz >> /tmp/$(hostname).apptainerinstall 2>&1
  
  ##chmod -R a+rx apptainer-${AVERSION}
  ##chmod +x apptainer-${AVERSION}/scripts/*

  ##cd apptainer-${AVERSION} 

  ##./mconfig >> /tmp/$(hostname).apptainerinstall 2>&1
  ##make -C builddir -j1  >> /tmp/$(hostname).apptainerinstall 2>&1
  ##make -C builddir install  >> /tmp/$(hostname).apptainerinstall 2>&1
  ##cd ..
  ##echo "APPTAINER END" >> /tmp/$(hostname).apptainerinstall

  # Create necessary directories
  mkdir -p /fsx/analysis_results/cromwell_executions  
  mkdir -p /fsx/analysis_results/ubuntu  
  mkdir -p /fsx/analysis_results/daylily              
  mkdir -p /fsx/miners/logs  
  mkdir -p /fsx/tmp
  mkdir -p /fsx/miners/bin               
  mkdir -p /fsx/scratch
  mkdir -p /fsx/resources/environments/containers/{ubuntu,daylily}/$(hostname)/
  mkdir -p /fsx/resources/environments/conda/{ubuntu,daylily}/$(hostname)/
  chmod -R a+wrx /fsx/analysis_results
  chmod -R a+wrx /fsx/scratch
  chmod -R a+wrx /fsx/miners
  chmod -R a+wrx /fsx/tmp
  chmod -R a+wrx /fsx/resources
  chown -R ubuntu:ubuntu /fsx/miners


  # Copy cached data from S3

  ln -s /fsx/data/cached_envs/conda/* /fsx/resources/environments/conda/ubuntu/$(hostname)/
  ln -s /fsx/data/cached_envs/containers/* /fsx/resources/environments/containers/ubuntu/$(hostname)/
  ln -s /fsx/data/cached_envs/conda/* /fsx/resources/environments/conda/daylily/$(hostname)/
  ln -s /fsx/data/cached_envs/containers/* /fsx/resources/environments/containers/daylily/$(hostname)/


  mv /opt/slurm/bin/sbatch /opt/slurm/sbin/sbatch
  aws s3 cp s3://${bucket}/cluster_boot_config/sbatch /opt/slurm/bin/sbatch
  chmod +x /opt/slurm/bin/sbatch

  mv /opt/slurm/bin/srun /opt/slurm/sbin/srun
  ln -s /opt/slurm/bin/sbatch /opt/slurm/bin/srun

  aws s3 cp s3://${bucket}/cluster_boot_config/projects_list.conf /opt/slurm/etc/projects_list.conf

  aws s3 cp s3://${bucket}/cluster_boot_config/sleep_test.sh /opt/slurm/bin/sleep_test.sh
  chmod a+x /opt/slurm/bin/sleep_test.sh


  # Restart SLURM Controller
  systemctl restart slurmctld
  touch /tmp/$(hostname).postslurmcfg

fi



# Tagging and Budget Bits

if [ "${cfn_node_type}" == "ComputeFleet" ];then

  # Create the folder used to save jobs information

  mkdir -p /tmp/jobs

  # Configure the script to run every minute
  echo "
* * * * * /opt/slurm/sbin/check_tags.sh
" | crontab -
  exit 0
else
  # HEADNODE ONLY

  # Cron script used to update the instance tags

  cat <<'EOF' > /opt/slurm/sbin/check_tags.sh
#!/bin/bash

source /etc/profile
TOKEN=$(curl -X PUT "http://169.254.169.254/latest/api/token" -H "X-aws-ec2-metadata-token-ttl-seconds: 21600")
region=$(curl -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/placement/region)
aws configure set region $region

update=0
tag_userid=""
tag_jobid=""
tag_project=""

if [ ! -f /tmp/jobs/jobs_users ] || [ ! -f /tmp/jobs/jobs_ids ]; then
  exit 0
fi

active_users=$(cat /tmp/jobs/jobs_users | sort | uniq )
active_jobs=$(cat /tmp/jobs/jobs_ids | sort )
echo $active_users > /tmp/jobs/tmp_jobs_users
echo $active_jobs > /tmp/jobs/tmp_jobs_ids
if [ -f /tmp/jobs/jobs_projects ]; then
  active_projects=$(cat /tmp/jobs/jobs_projects | sort | uniq )
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
  
  # Instance ID

  TOKEN=$(curl -X PUT "http://169.254.169.254/latest/api/token" -H "X-aws-ec2-metadata-token-ttl-seconds: 21600")
  MyInstID=$(curl -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)
  tag_userid=$(cat /tmp/jobs/tag_userid)
  tag_jobid=$(cat /tmp/jobs/tag_jobid)
  tag_project=$(cat /tmp/jobs/tag_project)
  aws ec2 create-tags --resources ${MyInstID} --tags Key=aws-parallelcluster-username,Value="${tag_userid}" --region ${region}
  aws ec2 create-tags --resources ${MyInstID} --tags Key=aws-parallelcluster-jobid,Value="${tag_jobid}" --region ${region}
  aws ec2 create-tags --resources ${MyInstID} --tags Key=aws-parallelcluster-project,Value="${tag_project}" --region ${region}  
fi

EOF

   chmod a+x /opt/slurm/sbin/check_tags.sh
   
   # Create Prolog and Epilog to tag the instances
   cat <<'EOF' > /opt/slurm/sbin/prolog.sh
#!/bin/bash

#slurm directory
export SLURM_ROOT=/opt/slurm
echo "${SLURM_JOB_USER}" >> /tmp/jobs/jobs_users
echo "${SLURM_JOBID}" >> /tmp/jobs/jobs_ids

#load the comment of the job.
Project=$($SLURM_ROOT/bin/scontrol show job ${SLURM_JOB_ID} | grep Comment | awk -F'=' '{print $2}')
Project_Tag=""
if [ ! -z "${Project}" ];then
  echo "${Project}" >> /tmp/jobs/jobs_projects
fi

EOF

   cat <<'EOF' > /opt/slurm/sbin/epilog.sh
#!/bin/bash
#slurm directory
export SLURM_ROOT=/opt/slurm
sed -i "0,/${SLURM_JOB_USER}/d" /tmp/jobs/jobs_users
sed -i "0,/${SLURM_JOBID}/d" /tmp/jobs/jobs_ids

#load the comment of the job.
Project=$($SLURM_ROOT/bin/scontrol show job ${SLURM_JOB_ID} | grep Comment | awk -F'=' '{print $2}')
Project_Tag=""
if [ ! -z "${Project}" ];then
  sed -i "0,/${Project}/d" /tmp/jobs/jobs_projects
fi

EOF

   chmod a+x /opt/slurm/sbin/prolog.sh
   chmod a+x /opt/slurm/sbin/epilog.sh
   
   # Configure slurm to use Prolog and Epilog
   echo "PrologFlags=Alloc" >> /opt/slurm/etc/slurm.conf
   echo "Prolog=/opt/slurm/sbin/prolog.sh" >> /opt/slurm/etc/slurm.conf
   echo "Epilog=/opt/slurm/sbin/epilog.sh" >> /opt/slurm/etc/slurm.conf
   
   systemctl restart slurmctld
fi

# Finalization
touch /tmp/$(hostname).postinstallcomplete
echo "Post-installation complete."
exit 0
