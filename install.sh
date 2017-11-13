#!/bin/bash
echo $PATH

VENV=gf
PYVER=3.5

# Is conda installed?
conda=$(which conda)
if [ ! "$conda" ] ; then
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
        -O miniconda.sh;
    bash miniconda.sh -f -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
fi

conda update -q -y conda
conda config --prepend channels conda-forge
conda config --append channels digitalglobe # for rasterio v 1.0a9
conda config --append channels ioos # for rasterio v 1.0a2


DEPARRAY=(pytables \
          numpy==1.13 \
          scipy=0.19.1 \
          pip \
    	  matplotlib=2.0.2 \
          rasterio=1.0 \
    	  fiona=1.7.8 \
    	  basemap=1.1.0 \
    	  basemap-data-hires=1.1.0 \
          shapely=1.5.17 \
    	  h5py=2.7.0 \
    	  gdal=2.1.4 \
    	  descartes=1.1.0 \
    	  pytest=3.2.0 \
    	  pytest-cov=2.5.1 \
          pytest-mpl=0.7 \
    	  configobj=5.0.6 \
    	  pandas=0.20.3 \
          scikit-learn=0.18.2 \
    	  scikit-image=0.13.0 \
          branca=0.2.0 \
    	  paramiko=2.1.2)

# Is the Travis flag set?
travis=0
while getopts t FLAG; do
  case $FLAG in
    t)
      travis=1
      ;;
  esac
done

# Append additional deps that are not for Travis CI
if [ $travis == 0 ] ; then
    DEPARRAY+=(ipython=6.1.0 spyder=3.2.1 jupyter=1.0.0 seaborn=0.8.0 \
        sphinx=1.6.3)
fi

# Turn off whatever other virtual environment user might be in
source deactivate

# Remove any previous virtual environments
CWD=`pwd`
cd $HOME;
conda remove --name $VENV --all -y
cd $CWD

# Create a conda virtual environment
echo "Creating the $VENV virtual environment"
echo "with the following dependencies:"
echo ${DEPARRAY[*]}
conda create --name $VENV -y python=$PYVER ${DEPARRAY[*]}

# Activate the new environment
echo "Activating the $VENV virtual environment"
source activate $VENV

# psutil
conda install -y psutil=5.2.1

# MapIO and impact-utils
echo "Installing MapIO..."
pip -q install https://github.com/usgs/MapIO/archive/master.zip
echo "Installing impact-utils..."
pip -q install \
    https://github.com/usgs/earthquake-impact-utils/archive/master.zip

# OpenQuake v2.5.0
echo "Downloading OpenQuake v2.5.0..."
curl --max-time 60 --retry 3 -L \
    https://github.com/gem/oq-engine/archive/v2.5.0.zip -o openquake.zip
pip -q install --no-deps openquake.zip
rm openquake.zip

pip install sphinx_rtd_theme
pip install folium

# This package
echo "Installing groundfailure-prod..."
pip install -e .


# Tell the user they have to activate this environment
echo "Type 'source activate $VENV' to use this new virtual environment."
