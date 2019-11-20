#!/bin/sh -l

sh -c "echo Hello world my name is $MY_NAME"

apt-get -y update
sh -c "echo updated tree"
apt-get -y install libnetcdf-dev
sh -c "echo netcdf stuff"
sh -c "nc-config --libdir"
sh -c "nc-config --includedir"
