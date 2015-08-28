#!/bin/sh
OS=$(lsb_release -si)
ARCH=$(uname -m | sed 's/x86_//;s/i[3-6]86/32/')
VER=$(lsb_release -sr)

if [ $OS == "Ubuntu" ]; then
    echo " %OS detected, starting installation via apt-get.."
    apt-get -y install python-pip
    apt-get -y install build-essential 
    apt-get -y install python-dev
    apt-get -y install python-gtk2
    apt-get -y install python-scipy
    apt-get -y install python-matplotlib
    apt-get -y install graphviz
    apt-get -y install libhdf5-dev
    apt-get -y install libxml2-dev
    apt-get -y install libxslt-dev

elif [ $OS == "Fedora" ]; then
    if [ $VER -lt 21 ]; then
        echo $OS "version 21 or earlier detected , starting installation via yum.."
        yum -y install gcc
        yum -y install gcc-gfortran
        yum -y install gcc-c++
        yum -y install python-devel
        yum -y install pygtk2
        yum -y install scipy
        yum -y install python-matplotlib
        yum -y install graphviz
        yum -y install hdf5-devel
        yum -y install libxml-devel
        yum -y install libxslt-devel

    else
        echo $OS "version 22 or later detected, starting installation via dnf.."
        dnf -y install gcc
        dnf -y install gcc-gfortran
        dnf -y install gcc-c++
        dnf -y install python-devel
        dnf -y install pygtk2
        dnf -y install scipy
        dnf -y install python-matplotlib
        dnf -y install graphviz
        dnf -y install hdf5-devel
        dnf -y install libxml-devel
        dnf -y install libxslt-devel
fi
else
    echo "$OS $VER detected. No auto installation for detected system available. Attempting to install python dependencies only."
fi
echo "Installing python dependencies via pip"
pip install -r requirements.txt
