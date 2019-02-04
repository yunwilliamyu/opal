#!/bin/bash
# Installs Opal dependencies and downloads some data files.

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd "$DIR"

#===================
# version comparison function from https://stackoverflow.com/questions/4023830/how-to-compare-two-strings-in-dot-separated-version-format-in-bash
#===================
vercomp () {
    if [[ $1 == $2 ]]
    then
        return 0
    fi
    local IFS=.
    local i ver1=($1) ver2=($2)
    # fill empty fields in ver1 with zeros
    for ((i=${#ver1[@]}; i<${#ver2[@]}; i++))
    do
        ver1[i]=0
    done
    for ((i=0; i<${#ver1[@]}; i++))
    do
        if [[ -z ${ver2[i]} ]]
        then
            # fill empty fields in ver2 with zeros
            ver2[i]=0
        fi
        if ((10#${ver1[i]} > 10#${ver2[i]}))
        then
            return 1
        fi
        if ((10#${ver1[i]} < 10#${ver2[i]}))
        then
            return 2
        fi
    done
    return 0
}

echo ========================================
echo Testing dependency: vowpal-wabbit command \'vw\'
if [ -x "$(command -v vw)" ]; then
    echo "vowpal-wabbit found. Testing version..."
    vercomp $(vw --version) 8.1.1
    code=$?
    if [ "$code" -eq 0 ] || [ "$code" -eq 1 ]; then
        echo "vowpal-wabbit version >=8.1.1 found."
    else
		echo "ERROR: vowpal-wabbit version <8.1.1. Please update."
        echo "e.g. on Ubuntu, run \"sudo apt-get install vowpal-wabbit\"."
        echo "On other OS, see install instructions here: https://github.com/VowpalWabbit/vowpal_wabbit"
        echo "Note: Python pip package vowpalwabbit does not install the vw command, so will not work."
        exit 1
    fi
else
    echo "ERROR: vowpal-wabbit not found. Please install."
    echo "e.g. on Ubuntu, run \"sudo apt-get install vowpal-wabbit\"."
    echo "On other OS, see install instructions here: https://github.com/VowpalWabbit/vowpal_wabbit"
    echo "Note: Python pip package vowpalwabbit does not install the vw command, so will not work."
    exit 1
fi
echo ========================================
echo Testing Python package dependencies:

if ! [ -x "$(command -v python)" ]; then
    echo "Python 2 must be installed. See your package manager for details."
    exit 1
fi

python -c "import pandas, sklearn"
if [ "$?" -eq 1 ]; then
    echo "Python packages pandas and/or sklearn not found. Please install."
    echo "e.g. if using pip installer, use \"pip install pandas sklearn\""
    echo "Note, if you are using Ubuntu and do not have pip installed, you can install using \"sudo apt-get install python-pip\"."
    echo "Otherwise, see instructions here: https://pip.pypa.io/en/stable/installing/"
    exit 1
fi


echo ========================================
echo Downloading and extracting example datasets
echo ========================================
if [ ! -f data/examples.installed ]; then
    wget -c http://giant.csail.mit.edu/opal/data.tar.bz2  || curl -O  http://giant.csail.mit.edu/opal/data.tar.bz2
    tar -xjf data.tar.bz2
    rm data.tar.bz2
    touch data/examples.installed
else
    echo "Examples already installed"
fi

echo ========================================
echo Installation complete!
echo Test out Opal on the exmaple data with:
echo ./opal.py simulate data/A1/test data/A1/train out_dir

totalk=$(awk '/^MemTotal:/{print $2}' /proc/meminfo)
if [ "$totalk" -lt 3355000 ]; then
    echo ========================================
    echo WARNING: System has less than 32 GiB of RAM.
    echo Opal may run out of memory for the model using standard parameters.
    echo ========================================
fi
