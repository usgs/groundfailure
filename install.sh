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
      "configobj>=5.0"
      "descartes>=1.1"
      "fiona>=1.8"
      "folium>=0.12"
      "gdal>=3.4"
      "hdf5>=1.10"
      "impactutils>=0.8"
      "libcomcat>=2.0"
      "libgdal>=3.4"
      "mapio>=0.7"
      "matplotlib-base>=3.5"
      "numpy>=1.20"
      "pytables>=3.6"
      "pytest>=6.2"
      "pytest-cov>=3.0"
      "pytest-faulthandler>=2.0"
      "rasterio>=1.1"
      "scipy>=1.8"
      "simplekml>=1.3"
)
# package_list=(
#       "python=$py_ver"
#       "configobj=5.0.6"
#       "descartes=1.1.0"
#       "fiona=1.8.13"
#       "folium=0.12.1"
#       "gdal=3.0.2"
#       "hdf5=1.10.6"
#       "impactutils=0.8.32"
#       "libcomcat=2.0.16"
#       "libgdal=3.0.2"
#       "mapio=0.7.31"
#       "matplotlib-base=3.5.1"
#       "numpy=1.20.3"
#       "pytables=3.6.1"
#       "pytest=6.2.5"
#       "pytest-cov=3.0.0"
#       "pytest-faulthandler=2.0.1"
#       "rasterio=1.1.2"
#       "scipy=1.7.3"
#       "simplekml=1.3.6"
# )

if [ $developer == 1 ]; then
    package_list=( "${package_list[@]}" "${dev_list[@]}" )
    echo ${package_list[*]}
fi

# Create a conda virtual environment
echo "Creating the $VENV virtual environment"
conda config --add channels 'defaults'
conda config --add channels 'conda-forge'
conda config --set channel_priority flexible
# conda env create -f $env_file --force
conda create -y -n $VENV ${package_list[*]}

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
