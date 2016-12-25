#!/bin/bash
echo ========================================
echo Setting up Genetic Data analysis Library and compiling utility tools
echo ========================================
cd util
sh INSTALL.sh
cd ..
echo ========================================
echo Downloading and extracting example datasets
echo ========================================
wget http://giant.csail.mit.edu/opal/data.tar.bz2
tar -xjf data.tar.bz2
rm data.tar.bz2
