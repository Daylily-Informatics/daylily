

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
# Redirect HTTP to HTTPS
<VirtualHost *:80>
    ServerName bloom.dyly.bio
    ServerAlias www.bloom.dyly.bio

    Redirect permanent / https://bloom.dyly.bio/

    ErrorLog ${APACHE_LOG_DIR}/bloom-http-error.log
    CustomLog ${APACHE_LOG_DIR}/bloom-http-access.log combined
</VirtualHost>

# HTTPS Configuration: Proxy to Port 8912
<VirtualHost *:443>
    ServerName bloom.dyly.bio
    ServerAlias www.bloom.dyly.bio

    # Enable SSL
    SSLEngine On
    SSLCertificateFile /etc/letsencrypt/live/dyly.bio/fullchain.pem
    SSLCertificateKeyFile /etc/letsencrypt/live/dyly.bio/privkey.pem

    # Enable SSL Proxy Engine
    SSLProxyEngine On

    # Proxy Configuration to Backend via HTTPS
    ProxyPreserveHost On
    ProxyPass / https://127.0.0.1:8912/
    ProxyPassReverse / https://127.0.0.1:8912/

    # Enforce Strong SSL Settings
    SSLProtocol all -SSLv3 -TLSv1 -TLSv1.1
    SSLCipherSuite HIGH:!aNULL:!MD5:!3DES
    SSLHonorCipherOrder On

    # HSTS Header
    Header always set Strict-Transport-Security "max-age=31536000; includeSubDomains; preload"

    ErrorLog ${APACHE_LOG_DIR}/bloom-https-error.log
    CustomLog ${APACHE_LOG_DIR}/bloom-https-access.log combined
</VirtualHost>



sudo a2enmod proxy proxy_http headers ssl








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



sudo a2ensite bloom.conf
sudo a2ensite gtc.conf
sudo a2ensite day.conf

sudo a2enmod ssl
sudo a2enmod proxy
sudo a2enmod proxy_http
sudo a2enmod headers

sudo ln -s /etc/apache2/sites-available/bloom.conf /etc/apache2/sites-enabled/
#above might fail, ok

sudo apachectl configtest
 #if ok
sudo systemctl reload apache2
sudo systemctl restart apache2



# wait and curl to test


sudo apt update
sudo apt install certbot
sudo apt install certbot python3-certbot-dns-route53


# ADD .aws credentials to root ~/.aws/credentials and config

sudo certbot certonly --dns-route53 -d "*.dyly.bio" -d "dyly.bio" -v

sudo groupadd ssl-users
sudo usermod -aG ssl-users ubuntu

#logout of shell, relogin

sudo chown root:ssl-users /etc/letsencrypt/live/dyly.bio/privkey.pem
sudo chown root:ssl-users /etc/letsencrypt/live/dyly.bio/fullchain.pem


sudo chmod 640 /etc/letsencrypt/live/dyly.bio/privkey.pem
sudo chmod 640 /etc/letsencrypt/live/dyly.bio/fullchain.pem
sudo chmod 750 /etc/letsencrypt/live/dyly.bio
sudo chmod 750 /etc/letsencrypt/live
sudo chmod 710 /etc/letsencrypt
sudo chmod 750 /etc/letsencrypt/live
sudo chmod 750 /etc/letsencrypt/live/dyly.bio
sudo chown -R root:ssl-users /etc/letsencrypt

sudo chmod 640 /etc/letsencrypt/archive/dyly.bio/cert1.pem
sudo chmod 640 /etc/letsencrypt/archive/dyly.bio/chain1.pem
sudo chmod 640 /etc/letsencrypt/archive/dyly.bio/fullchain1.pem
sudo chmod 640 /etc/letsencrypt/archive/dyly.bio/privkey1.pem
sudo chmod 750 /etc/letsencrypt
sudo chmod 710 /etc/letsencrypt
sudo chmod 750 /etc/letsencrypt/archive
sudo chmod 750 /etc/letsencrypt/live
sudo chmod 750 /etc/letsencrypt/renewal
sudo chmod 750 /etc/letsencrypt/renewal-hooks
sudo chmod 750 /etc/letsencrypt/accounts


sudo chown -h root:ssl-users /etc/letsencrypt/live/dyly.bio/*
sudo chmod 750 /etc/letsencrypt/live/dyly.bio/
sudo chmod 640 /etc/letsencrypt/live/dyly.bio/*
sudo chown root:ssl-users /etc/letsencrypt/archive/dyly.bio/*

sudo chmod 640 /etc/letsencrypt/archive/dyly.bio/*



#    --ssl-certfile  /etc/letsencrypt/live/dyly.bio/fullchain.pem \
#    --ssl-keyfile /etc/letsencrypt/live/dyly.bio/privkey.pem


# confirm certbot renewal
sudo systemctl list-timers | grep certbot


sudo nano /etc/letsencrypt/renewal-hooks/post/fix-permissions.sh
#!/bin/bash

# Fix permissions for archive files
chmod 640 /etc/letsencrypt/archive/dyly.bio/*
chown root:ssl-users /etc/letsencrypt/archive/dyly.bio/*

# Fix permissions for live files
chmod 750 /etc/letsencrypt/live/dyly.bio/
chmod 640 /etc/letsencrypt/live/dyly.bio/*

# Restart the server to apply the new certificate
# assuming a running tmux serssion bloom_server
sudo -u ubuntu tmux send-keys -t bloom_server C-c
sleep 10
sudo -u ubuntu tmux send-keys -t b "source /home/ubuntu/miniconda3/bin/activate BLOOM && /home/ubuntu/projects/bloom/run_bloomui.sh --host 0.0.0.0 --port 8912 --mode dev" Enter


sudo chmod +x /etc/letsencrypt/renewal-hooks/post/fix-permissions.sh
sudo certbot renew --dry-run

