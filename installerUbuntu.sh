#!/bin/sh
apt-get -y install build-essential python-dev
apt-get -y install pygtk2
apt-get -y install scipy
apt-get -y install python-matplotlib

pip install -r requirements.txt
