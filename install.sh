#!/bin/bash

unamestr=`uname`
if [ "$unamestr" == 'Linux' ]; then
    prof=~/.bashrc
    mini_conda_url=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    matplotlibdir=~/.config/matplotlib
elif [ "$unamestr" == 'FreeBSD' ] || [ "$unamestr" == 'Darwin' ]; then
    prof=~/.bash_profile
    mini_conda_url=https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    matplotlibdir=~/.matplotlib
else
    echo "Unsupported environment. Exiting."
    exit
fi

source $prof

# Name of new environment (must also change this in .yml files)
VENV=gf
# Python version
py_ver=3.8

# Set to 1 if you are a developer and want ipython etc. installed
developer=0

# Is conda installed?
conda --version
if [ $? -ne 0 ]; then
    echo "No conda detected, installing miniconda..."

    curl -L $mini_conda_url -o miniconda.sh;

    # if curl fails, bow out gracefully
    if [ $? -ne 0 ];then
        echo "Failed to create download miniconda installer shell script. Exiting."
        exit 1
    fi
    
    echo "Install directory: $HOME/miniconda"

    bash miniconda.sh -f -b -p $HOME/miniconda

    # if miniconda.sh fails, bow out gracefully
    if [ $? -ne 0 ];then
        echo "Failed to run miniconda installer shell script. Exiting."
        exit 1
    fi
    
    . $HOME/miniconda/etc/profile.d/conda.sh
else
    echo "conda detected, installing $VENV environment..."
fi

# # make defaults higher priority, set that priority to strict
# conda config --add channels 'conda-forge'
# conda config --add channels 'defaults'
# conda config --set channel_priority strict

# echo "PATH:"
# echo $PATH
# echo ""

# add source command to profile file if it isn't already there
grep "/etc/profile.d/conda.sh" $prof
if [ $? -ne 0 ]; then
    echo ". $_CONDA_ROOT/etc/profile.d/conda.sh" >> $prof
fi

#env_file=environment.yml

# Start in conda base environment
echo "Activate base virtual environment"
conda activate base

# Remove existing gf environment if it exists
conda remove -y -n $VENV --all

conda install mamba -y -n base -c conda-forge

dev_list=(
    "ipython"
    "spyder"
    "sphinx"
    "sphinx-argparse"
    "jupyterlab"
)

# Package list:
package_list=(
      "python=$py_ver"
      "configobj"
      "descartes"
      "fiona"
      "folium"
      "gdal=3.0"
      "impactutils"
      "libcomcat"
      "mapio"
      "matplotlib-base"
      "numpy"
      "pytables"
      "pytest"
      "pytest-cov"
      "pytest-faulthandler"
      "rasterio"
      "scipy"
      "simplekml"
)

if [ $developer == 1 ]; then
    package_list=( "${package_list[@]}" "${dev_list[@]}" )
    echo ${package_list[*]}
fi

# Create a conda virtual environment
echo "Creating the $VENV virtual environment"
# conda env create -f $env_file --force
mamba create -y -n $VENV -c defaults -c conda-forge \
      --strict-channel-priority ${package_list[*]}

# Bail out at this point if the conda create command fails.
# Clean up zip files we've downloaded
if [ $? -ne 0 ]; then
    echo "Failed to create conda environment.  Resolve any conflicts, then try again."
    exit
fi


# Activate the new environment
echo "Activating the $VENV virtual environment"
conda activate $VENV

# if conda activate fails, bow out gracefully
if [ $? -ne 0 ];then
    echo "Failed to activate ${VENV} conda environment. Exiting."
    exit 1
fi

# upgrade pip, mostly so pip doesn't complain about not being new...
pip install --upgrade pip

# if pip upgrade fails, complain but try to keep going
if [ $? -ne 0 ];then
    echo "Failed to upgrade pip, trying to continue..."
    exit 1
fi

# This package
echo "Installing $VENV"
pip install -e .

# if pip install fails, bow out gracefully
if [ $? -ne 0 ];then
    echo "Failed to pip install this package. Exiting."
    exit 1
fi

# Tell the user they have to activate this environment
echo "Type 'conda activate $VENV' to use this new virtual environment."
