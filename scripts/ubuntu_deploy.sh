#!/bin/sh

# ubuntu_deploy.sh
#
# Installs the reference server with the example data on an AWS
# Ubuntu instance

sudo apt-get install python-dev zlib1g-dev python-virtualenv
sudo apt-get install apache2 libapache2-mod-wsgi
sudo a2enmod wsgi
sudo mkdir /var/cache/apache2/python-egg-cache
sudo chown www-data:www-data /var/cache/apache2/python-egg-cache/
sudo mkdir /srv/ga4gh
sudo chown $USER /srv/ga4gh
cd /srv/ga4gh
virtualenv ga4gh-server-env
source ga4gh-server-env/bin/activate
pip install --pre ga4gh
deactivate
wget http://www.well.ox.ac.uk/~jk/ga4gh-example-data.tar
tar -xf ga4gh-example-data.tar
echo 'from ga4gh.frontend import app as application
import ga4gh.frontend as frontend
frontend.configure("/srv/ga4gh/config.py")' > /srv/ga4gh/application.wsgi
echo 'DATA_SOURCE = "/srv/ga4gh/ga4gh-example-data"' > /srv/ga4gh/config.py
sudo echo '<VirtualHost *:80>
        ServerAdmin webmaster@localhost
        DocumentRoot /var/www/html

        ErrorLog ${APACHE_LOG_DIR}/error.log
        CustomLog ${APACHE_LOG_DIR}/access.log combined

        WSGIDaemonProcess ga4gh \
            python-path=/srv/ga4gh/ga4gh-server-env/lib/python2.7/site-packages \
            python-eggs=/var/cache/apache2/python-egg-cache
        WSGIScriptAlias /ga4gh /srv/ga4gh/application.wsgi

        <Directory /srv/ga4gh>
            WSGIProcessGroup ga4gh
            WSGIApplicationGroup %{GLOBAL}
            Require all granted
        </Directory>
</VirtualHost>' > /etc/apache2/sites-enabled/000-default.conf
sudo service apache2 restart
wget -qO- http://localhost/ga4gh &> /dev/null
if [ $? -ne 0 ]; then
    echo "FAILURE!"
else
    echo "SUCCESS!"
fi
