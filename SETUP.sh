#!/bin/bash
echo ========================================
echo Downloading and extracting example datasets
echo ========================================
wget http://giant.csail.mit.edu/opal/data.tar.bz2
tar -xjf data.tar.bz2
rm data.tar.bz2
