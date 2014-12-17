#!/bin/sh
ARCH=$(uname -m | sed 's/x86_//;s/i[3-6]86/32/')

if [ -f /etc/lsb-release ]; then
    . /etc/lsb-release
    OS=$DISTRIB_ID
    VER=$DISTRIB_RELEASE
	echo " no auto installation for detected system available"

elif [ -f /etc/debian_version ]; then
    OS=debian-ubuntu 
    VER=$(cat /etc/debian_version)
    echo " %OS detected, starting installation .."
    
    apt-get -y install build-essential python-dev
	apt-get -y install pygtk2
	apt-get -y install scipy
	apt-get -y install python-matplotlib

elif [ -f /etc/redhat-release ]; then
    # TODO add code for Red Hat and CentOS here
    OS=redhat-fedora
    echo $OS "detected, starting installation .."
    
    yum -y install gcc
	yum -y install gcc-gfortran
	yum -y install gcc-c++
	yum -y install python-devel
	yum -y install pygtk2
	yum -y install scipy
	yum -y install python-matplotlib
    
else
    OS=$(uname -s)
    VER=$(uname -r)
    echo " no auto installation for detected system available"
fi

pip install -r requirements.txt
