#!/bin/sh
yum -y install gcc
yum -y install gcc-gfortran
yum -y install gcc-c++
yum -y install python-devel
yum -y install pygtk2
yum -y install scipy
yum -y install python-matplotlib

pip install -r requirements.txt
