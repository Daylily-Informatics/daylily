

sudo adduser --uid 1002 --disabled-password --gecos "" daylily || echo "daylily user add failed"

sudo apt update -y

sudo apt install -y tmux emacs rclone parallel atop htop glances fd-find docker.io   build-essential libssl-dev uuid-dev libgpgme-dev squashfs-tools   libseccomp-dev pkg-config openjdk-11-jdk wget unzip nasm yasm  fuse2fs gocryptfs  golang-go


sudo add-apt-repository -y ppa:apptainer/ppa

# Update package lists
sudp apt update

# Install Apptainer
sudo apt install -y apptainer




ssh-keygen

head ~/.ssh/id_rsa.pub



#curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
curl "https://awscli.amazonaws.com/awscli-exe-linux-aarch64.zip" -o "awscliv2.zip"

unzip awscliv2.zip
sudo ./aws/install
aws --version



MACHINE="linux_arm"

echo "Autoinstalling Miniconda for Linux ARM..."
wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh -O Miniconda3-Linux-aarch64.sh
chmod +x Miniconda3-Linux-aarch64.sh
./Miniconda3-Linux-aarch64.sh -b -p ~/miniconda3
rm Miniconda3-Linux-aarch64.sh
~/miniconda3/bin/conda init $SHELL

$SHELL

mkdir ~/projects

cd ~/projects

git clone 


