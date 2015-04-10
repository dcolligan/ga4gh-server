#!/bin/sh

#
# Deploys a ngnix/uwsgi-based ga4gh server on an Ubuntu AWS instance
#

# TODO at this point we assume the apache instructions were carried out.
# figure out which of these were essential to the nginx install

# TODO problem: 
# flask is configured by the cli at startup, but under wsgi we don't
# go through this phase...
# probably should do something like using the WSGI file in the apache setup

# install nginx
sudo apt-get install nginx

# install uwsgi
pip install uwsgi

# edit the nginx conf file
vi /etc/nginx/nginx.conf

# configure uwsgi
echo '[uwsgi]
socket = /tmp/uwsgi.sock
module = ga4gh.frontend:app' > ga4gh-uwsgi.ini

# start uwsgi
uwsgi ga4gh-uwsgi.ini

# reload the nginx conf file
nginx -s reload  # (or sudo service nginx restart)
