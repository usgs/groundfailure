#!/bin/bash

echo "Path:"
echo $PATH

# Name of new environment (must also change this in .yml files)
VENV=gf

# Are the reset/travis flags set?
reset=0
travis=0
while getopts rt FLAG; do
  case $FLAG in
    r)
        reset=1;;
    t)
	travis=1;;
  esac
done


# Is conda installed?
conda=$_CONDA_EXE

# If not, install miniconda
if [ ! "$conda" ] ; then
    echo "No conda detected, installing miniconda"
    curl https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
         -o miniconda.sh;
    echo "Install directory: $HOME/miniconda"
    bash miniconda.sh -f -b -p $HOME/miniconda
    rm -f miniconda.sh
fi

# Source bash startup file
#if [ -f $HOME/.bash_profile ]; then
#    echo 'Sourcing .bash_profile'
#    . $HOME/.bash_profile
#    cat $HOME/.bash_profile
#    echo ""
#fi

#if [ -f $HOME/.bashrc ]; then
#    echo 'Sourcing .bashrc'
#    . $HOME/.bashrc
#    cat $HOME/.bashrc
#    echo ""
#fi

# Need this to get conda into path
. $HOME/miniconda/etc/profile.d/conda.sh

echo "PATH:"
echo $PATH
echo ""

# Choose an environment file based on platform
unamestr=`uname`
if [ "$unamestr" == 'Linux' ]; then
    env_file=environment_linux.yml
elif [ "$unamestr" == 'FreeBSD' ] || [ "$unamestr" == 'Darwin' ]; then
    env_file=environment_osx.yml
fi

# If the user has specified the -r (reset) flag, then create an
# environment based on only the named dependencies, without
# any versions of packages specified.
if [ $reset == 1 ]; then
    echo "Ignoring platform, letting conda sort out dependencies..."
    env_file=environment.yml
fi

# Start in conda base environment
echo "Activate base virtual environment"
conda activate base

# Create a conda virtual environment
echo "Creating the $VENV virtual environment"
conda env create -f $env_file --force

# Activate the new environment
echo "Activating the $VENV virtual environment"
source activate $VENV

# This package
echo "Installing $VENV"
pip install -e .

# Tell the user they have to activate this environment
echo "Type 'conda activate $VENV' to use this new virtual environment."
