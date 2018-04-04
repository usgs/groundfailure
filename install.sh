#!/bin/bash

echo "Path:"
echo $PATH

# Name of new environment (must also change this in .yml files)
VENV=gf

# Is the reset flag set?
reset=0
while getopts r FLAG; do
  case $FLAG in
    r)
        reset=1
        
      ;;
  esac
done

# Source bash startup file, not sure why it doesn't do this automatically
if [ -f ~/.bash_profile ]; then
    . ~/.bash_profile
fi
if [ -f ~/.bashrc ]; then
    . ~/.bashrc
fi

# Is conda installed?
conda=$_CONDA_EXE
if [ ! "$conda" ] ; then
    echo "No conda detected, installing miniconda..."
    curl https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
        -o miniconda.sh;
    bash miniconda.sh -f -b -p $HOME/miniconda
fi


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
conda activate base

# Create a conda virtual environment
echo "Creating the $VENV virtual environment:"
conda env create -f $env_file --force

# Activate the new environment
echo "Activating the $VENV virtual environment"
conda activate $VENV

# This package
echo "Installing $VENV..."
pip install -e .

# Tell the user they have to activate this environment
echo "Type 'conda activate $VENV' to use this new virtual environment."
