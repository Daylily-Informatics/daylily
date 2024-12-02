

sudo adduser --uid 1002 --disabled-password --gecos "" daylily || echo "daylily user add failed"

sudo apt update -y

sudo apt install -y tmux emacs rclone parallel atop htop glances fd-find docker.io   build-essential libssl-dev uuid-dev libgpgme-dev squashfs-tools   libseccomp-dev pkg-config openjdk-11-jdk wget unzip nasm yasm  fuse2fs gocryptfs  golang-go isal


sudo add-apt-repository -y ppa:apptainer/ppa

# Update package lists
sudp apt update

# Install Apptainer
sudo apt install -y apptainer




ssh-keygen

head ~/.ssh/id_rsa.pub



curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"


unzip awscliv2.zip
sudo ./aws/install
aws --version





echo "Autoinstalling Miniconda for Linux ARM..."
wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda3-Linux-x86_64.sh
chmod +x  Miniconda3-Linux-x86_64.sh
./Miniconda3-Linux-x86_64.sh -b -p ~/miniconda3
rm  Miniconda3-Linux-x86_64.sh
~/miniconda3/bin/conda init $SHELL

$SHELL

mkdir ~/projects

cd ~/projects

git clone git@github.com:Daylily-Informatics/github_markdown_text_colorizer.git

git clone git@github.com:Daylily-Informatics/bloom.git

git clone git@github.com:Daylily-Informatics/daylily

git clone git@github.com:Daylily-Informatics/slim_goodie.git

git clone git@github.com:Daylily-Informatics/zebra_day.git

git clone git@github.com:Daylily-Informatics/regulatory_capture.git

git clone git@github.com:Daylily-Informatics/fedex_tracking_day.git

git clone git@github.com:Daylily-Informatics/daylily-informatics.github.io.git

git clone git@github.com:Daylily-Informatics/daylily-web.git

git clone git@github.com:Daylily-Informatics/daylily-web-mobile.git


sudo apt update
sudo apt install -y  apache2

sudo a2enmod proxy proxy_http
sudo systemctl restart apache2


sudo emacs  /etc/apache2/sites-available/bloom.conf

<VirtualHost *:80>
    ServerName bloom.dyly.bio
    ServerAlias www.bloom.dyly.bio

    ProxyPreserveHost On
    ProxyPass / http://127.0.0.1:8912/
    ProxyPassReverse / http://127.0.0.1:8912/

    ErrorLog ${APACHE_LOG_DIR}/bloom-error.log
    CustomLog ${APACHE_LOG_DIR}/bloom-access.log combined
</VirtualHost>


sudo emacs /etc/apache2/sites-available/gtc.conf

<VirtualHost *:80>
    ServerName gtc.dyly.bio
    ServerAlias www.gtc.dyly.bio

    ProxyPreserveHost On
    ProxyPass / http://127.0.0.1:8911/
    ProxyPassReverse / http://127.0.0.1:8911/

    ErrorLog ${APACHE_LOG_DIR}/gtc-error.log
    CustomLog ${APACHE_LOG_DIR}/gtc-access.log combined
</VirtualHost>


sudo emacs /etc/apache2/sites-available/day.conf

<VirtualHost *:80>
    ServerName day.dyly.bio
    ServerAlias www.gtc.dyly.bio

    ProxyPreserveHost On
    ProxyPass / http://127.0.0.1:8914/
    ProxyPassReverse / http://127.0.0.1:8914/

    ErrorLog ${APACHE_LOG_DIR}/gtc-error.log
    CustomLog ${APACHE_LOG_DIR}/gtc-access.log combined
</VirtualHost>




sudo a2ensite bloom.conf
sudo a2ensite gtc.conf
sudo a2ensite day.conf
sudo systemctl reload apache2


# wait and curl to test