

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
    ServerAlias www.bloom.dyly.bio www.bloom.daylilyinformatics.com bloom.daylilyinformatics.com www.bloom.daylily.bio bloom.daylily.bio www.bloom.daylilyinformatics.bio bloom.daylilyinformatics.bio 

    Redirect permanent / https://bloom.dyly.bio/

    ErrorLog ${APACHE_LOG_DIR}/bloom-http-error.log
    CustomLog ${APACHE_LOG_DIR}/bloom-http-access.log combined
</VirtualHost>

# HTTPS Configuration: Proxy to Port 8912
<VirtualHost *:443>
    ServerName bloom.dyly.bio
    ServerAlias www.bloom.dyly.bio www.bloom.daylilyinformatics.com bloom.daylilyinformatics.com www.bloom.daylily.bio bloom.daylily.bio www.bloom.daylilyinformatics.bio bloom.daylilyinformatics.bio 

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


uvicorn main:app --reload --log-level trace --port 8915 --timeout-keep-alive 303 --host 0.0.0.0 \ --ssl-certfile /etc/letsencrypt/live/dyly.bio/fullchain.pem --ssl-keyfile /etc/letsencrypt/live/dyly.bio/privkey.pem
Usage: uvicorn [OPTIONS] APP
Try 'uvicorn --help' for help.

Error: Got unexpected extra arguments ( --ssl-certfile /etc/letsencrypt/live/dyly.bio/fullchain.pem)
sudo a2enmod proxy proxy_http headers ssl


sudo emacs /etc/apache2/sites-available/gtc.conf
<VirtualHost *:443>
    ServerName gtc.dyly.bio
    ServerAlias www.gtc.dyly.bio gtc.dyly.bio www.gtc.daylily.bio gtc.daylily.bio www.gtc.daylilyinformatics.com gtc.daylilyinformatics.com www.daylilyinformatics.bio gtc.daylilyinformatics.bio

    # Enable SSL
    SSLEngine On
    SSLCertificateFile /etc/letsencrypt/live/dyly.bio/fullchain.pem
    SSLCertificateKeyFile /etc/letsencrypt/live/dyly.bio/privkey.pem

    # Proxy Configuration
    ProxyPreserveHost On
    ProxyPass / https://127.0.0.1:8911/
    ProxyPassReverse / https://127.0.0.1:8911/

    # Enforce Strong SSL Settings
    SSLProtocol all -SSLv3 -TLSv1 -TLSv1.1
    SSLCipherSuite HIGH:!aNULL:!MD5:!3DES
    SSLHonorCipherOrder On

    # HSTS Header
    Header always set Strict-Transport-Security "max-age=31536000; includeSubDomains; preload"

    ErrorLog ${APACHE_LOG_DIR}/gtc-https-error.log
    CustomLog ${APACHE_LOG_DIR}/gtc-https-access.log combined
</VirtualHost>
<VirtualHost *:80>
    ServerName gtc.dyly.bio
    ServerAlias www.gtc.dyly.bio gtc.dyly.bio www.gtc.daylily.bio gtc.daylily.bio www.gtc.daylilyinformatics.com gtc.daylilyinformatics.com www.daylilyinformatics.bio gtc.daylilyinformatics.bio

    Redirect permanent / https://gtc.dyly.bio/


    ErrorLog ${APACHE_LOG_DIR}/gtc-80-error.log
    CustomLog ${APACHE_LOG_DIR}/gtc-80-access.log combined
</VirtualHost>


sudo emacs /etc/apache2/sites-available/day-www.conf
<VirtualHost *:80>
    ServerName dyly.bio
    ServerAlias dyly.bio
    Redirect permanent / https://dyly.bio/
    ErrorLog ${APACHE_LOG_DIR}/dyly-www-http-error.log
    CustomLog ${APACHE_LOG_DIR}/dyly-www-http-access.log combined
</VirtualHost>

<VirtualHost *:443>
    ServerName dyly.bio
    ServerAlias dyly.bio

    # Enable SSL
    SSLEngine On
    SSLCertificateFile /etc/letsencrypt/live/dyly.bio/fullchain.pem
    SSLCertificateKeyFile /etc/letsencrypt/live/dyly.bio/privkey.pem

    # Enable SSL Proxy Engine
    SSLProxyEngine On

    # Proxy Configuration to Backend via HTTPS
    ProxyPreserveHost On
    ProxyPass / https://127.0.0.1:8915/
    ProxyPassReverse / https://127.0.0.1:8915/

    # Enforce Strong SSL Settings
    SSLProtocol all -SSLv3 -TLSv1 -TLSv1.1
    SSLCipherSuite HIGH:!aNULL:!MD5:!3DES
    SSLHonorCipherOrder On

    # HSTS Header
    Header always set Strict-Transport-Security "max-age=31536000; includeSubDomains; preload"

    ErrorLog ${APACHE_LOG_DIR}/dyly-www-https-error.log
    CustomLog ${APACHE_LOG_DIR}/dyly-www-https-access.log combined
</VirtualHost>






sudo emacs /etc/apache2/sites-available/day.conf
<VirtualHost *:443>
    ServerName day.dyly.bio
    ServerAlias www.day.dyly.bio day.dyly.bio www.day.daylilyinformatics.com day.daylilyinformatics.com www.day.daylily.bio day.daylily.bio day.daylilyinformatics.bio www.day.daylillyinformatics.bio

    # Enable SSL
    SSLEngine On
    SSLCertificateFile /etc/letsencrypt/live/dyly.bio/fullchain.pem
    SSLCertificateKeyFile /etc/letsencrypt/live/dyly.bio/privkey.pem

    # Proxy Configuration
    ProxyPreserveHost On
    ProxyPass / https://127.0.0.1:8914/
    ProxyPassReverse / https://127.0.0.1:8914/

    # Enforce Strong SSL Settings
    SSLProtocol all -SSLv3 -TLSv1 -TLSv1.1
    SSLCipherSuite HIGH:!aNULL:!MD5:!3DES
    SSLHonorCipherOrder On

    # HSTS Header
    Header always set Strict-Transport-Security "max-age=31536000; includeSubDomains; preload"

    ErrorLog ${APACHE_LOG_DIR}/day-https-443-error.log
    CustomLog ${APACHE_LOG_DIR}/day-https-443-access.log combined
</VirtualHost>
<VirtualHost *:80>
    ServerName day.dyly.bio
    ServerAlias www.day.dyly.bio day.dyly.bio www.day.daylilyinformatics.com day.daylilyinformatics.com www.day.daylily.bio day.daylily.bio day.daylilyinformatics.bio www.day.daylillyinformatics.bio

    Redirect permanent / https://day.dyly.bio/

    ErrorLog ${APACHE_LOG_DIR}/day-error.log
    CustomLog ${APACHE_LOG_DIR}/day-access.log combined
</VirtualHost>


sudo a2ensite bloom.conf
sudo a2ensite gtc.conf
sudo a2ensite day.conf
sudo a2ensite day-www.conf


sudo a2enmod ssl
sudo a2enmod proxy
sudo a2enmod proxy_http
sudo a2enmod headers

sudo ln -s /etc/apache2/sites-available/bloom.conf /etc/apache2/sites-enabled/
#above might fail, ok

sudo apachectl configtest
 #if ok
sudo systemctl daemon-reload

sudo systemctl reload apache2
sudo systemctl restart apache2

sudo systemctl status apache2
sudo apachectl configtest


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
sudo -u ubuntu tmux send-keys -t dyly_bloom C-c
sudo -u ubuntu tmux send-keys -t dyly_gtmc C-c
sudo -u ubuntu tmux send-keys -t dyly_web C-c
sudo -u ubuntu tmux send-keys -t dyly_dewey C-c
sudo -u ubuntu tmux send-keys -t dyly_day C-c

sleep 1


sudo lsof -i :8911 | grep python | awk '{print $2}' | xargs -r sudo kill -9 > /dev/null 2>&1 || echo "kill failed"
sleep 1
sudo lsof -i :8912 | grep gunicorn | awk '{print $2}' | xargs -r sudo kill -9 > /dev/null 2>&1 || echo "kill failed"
sleep 11
sudo lsof -i :8913 | grep gunicorn | awk '{print $2}' | xargs -r sudo kill -9 > /dev/null 2>&1 || echo "kill failed"
sleep 11
sudo lsof -i :8914 | grep python | awk '{print $2}' | xargs -r sudo kill -9 > /dev/null 2>&1 || echo "kill failed"
sleep 1
sudo lsof -i :8915 | grep gunicorn | awk '{print $2}' | xargs -r sudo kill -9 > /dev/null 2>&1 || echo "kill failed"
sleep 1

 sudo systemctl daemon-reload

sudo systemctl reload apache2
sudo systemctl restart apache2

sudo systemctl status apache2

# GTMC
sudo -u ubuntu tmux send-keys -t dyly_gtmc "source /home/ubuntu/miniconda3/bin/activate GMTC && python /home/ubuntu/projects/github_markdown_text_colorizer/github_markdown_text_colorizer/bin/gitmdtxtclr3.py 0.0.0.0 8911 http://54.190.46.215" Enter
sleep 1

# BLOOM
sudo -u ubuntu tmux send-keys -t dyly_bloom "source /home/ubuntu/miniconda3/bin/activate BLOOM && /home/ubuntu/projects/bloom/run_bloomui.sh --host 0.0.0.0 --port 8912 --mode prod" Enter
sleep 1

# DAY- 8913- daylilyUI
#sudo -u ubuntu tmux send-keys -t dyly_day "source /home/ubuntu/miniconda3/bin/activate DYLYD && /home/ubuntu/projects/dyly-day/run_dyly_day.sh --host 0.0.0.0 --port 8914 --mode prod" Enter
#sleep 1

# dyly-web

sudo -u ubuntu tmux send-keys -t dyly_web "source /home/ubuntu/miniconda3/bin/activate DYLY && ./run_dyly_www.sh --mode dev --ssl-certfile /etc/letsencrypt/\
live/dyly.bio/fullchain.pem --ssl-keyfile /etc/letsencrypt/live/dyly.bio/privkey.pem 2>&1 | tee logs/dyly_www.log" Enter
sleep 1


tmux new -s dyly_bloom
tmux new -s dyly_day
tmux new -s dyly_gtmc
tmux new -s dyly_web
tmux new -s dyly_dewey

sudo chmod +x /etc/letsencrypt/renewal-hooks/post/fix-permissions.sh
sudo certbot renew --dry-run


sudo certbot certonly --dns-route53 \
    -d "*.dyly.bio" -d "dyly.bio" \
    -d "*.daylily.bio" -d "daylily.bio" \
    -d "*.daylilyinformatics.com" -d "daylilyinformatics.com" \
    -d "*.daylilyinformatics.bio" -d "daylilyinformatics.bio" \
    -d "*.daylily.cloud" -d "daylily.cloud" \
    -v



<VirtualHost *:443>
    ServerName day.dyly.bio
    ServerAlias www.day.dyly.bio day.dyly.bio www.day.daylilyinformatics.com day.daylilyinformatics.com www.day.daylily.bio day.daylily.bio day.daylilyinformatics.bio www.day.daylillyinformatic\
s.bio

    # Enable SSL
    SSLEngine On
    SSLCertificateFile /etc/letsencrypt/live/dyly.bio/fullchain.pem
    SSLCertificateKeyFile /etc/letsencrypt/live/dyly.bio/privkey.pem

    # Proxy Configuration
    ProxyPreserveHost On
    ProxyPass / https://127.0.0.1:8914/
    ProxyPassReverse / https://127.0.0.1:8914/

    # Enforce Strong SSL Settings
    SSLProtocol all -SSLv3 -TLSv1 -TLSv1.1
    SSLCipherSuite HIGH:!aNULL:!MD5:!3DES
    SSLHonorCipherOrder On

    # HSTS Header
    Header always set Strict-Transport-Security "max-age=31536000; includeSubDomains; preload"

    ErrorLog ${APACHE_LOG_DIR}/day-https-443-error.log
    CustomLog ${APACHE_LOG_DIR}/day-https-443-access.log combined
</VirtualHost>
<VirtualHost *:80>
    ServerName day.dyly.bio
    ServerAlias www.day.dyly.bio day.dyly.bio www.day.daylilyinformatics.com day.daylilyinformatics.com www.day.daylily.bio day.daylily.bio day.daylilyinformatics.bio www.day.daylillyinformatic\
s.bio

    Redirect permanent / https://day.dyly.bio/

    ErrorLog ${APACHE_LOG_DIR}/day-error.log
    CustomLog ${APACHE_LOG_DIR}/day-access.log combined
    </VirtualHost>


    
sudo systemctl daemon-reload

sudo systemctl reload apache2
sudo systemctl restart apache2

sudo systemctl status apache2
sudo apachectl configtest






# Add to
sudo emacs /etc/apache2/apache2.conf
<Location "/wordpress">
    Require all denied
</Location>

<Files "setup-config.php">
    Require all denied
</Files>


sudo systemctl restart apache2   # Debian/Ubuntu
