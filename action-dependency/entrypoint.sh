#!/bin/sh -l

sh -c "echo Hello world my name is $MY_NAME"

apt-get -y update
sh -c "updated tree"
apt-get -y install libnetcdf-dev
apt-get -y install cmake
