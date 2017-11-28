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

DEPARRAY=(basemap=1.1.0 \
          basemap-data-hires=1.1.0 \
          beautifulsoup4=4.6.0\
          branca=0.2.0 \
          configobj=5.0.6 \
          descartes=1.1.0 \
          gdal=2.1.4 \
          h5py=2.7.1 \
          lxml=4.1.1 \
          matplotlib=2.0.2 \
          numpy==1.13 \
          pandas=0.20.3 \
          paramiko=2.3.1 \
          pip \
          pytables \
          pytest=3.2.5 \
          pytest-cov=2.5.1 \
          pytest-mpl=0.7 \
          scipy=0.19.1 \
          scikit-learn=0.18.2 \
          scikit-image=0.13.0 \
)

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


# Force incompatible versions of fiona, shapely, rasterio
# and their dependencies
conda install -y -f fiona=1.7.8 
conda install -y -f shapely=1.5.17 
conda install -y -f rasterio=0.36  
conda install -y -f affine=2.1.0
conda install -y -f click=6.7
conda install -y -f click-plugins=1.0.3

# psutil
conda install -y psutil=5.2.1

# Sphinx extensions
echo "Installing Sphinx extensions..."
pip install sphinxcontrib-napoleon
pip install sphinx-argparse
pip install sphinx_rtd_theme

# MapIO and impact-utils
echo "Installing MapIO..."
curl --max-time 60 --retry 3 -L \
     https://github.com/usgs/MapIO/archive/master.zip -o mapio.zip
pip -q install --no-deps mapio.zip
rm mapio.zip

echo "Installing impact-utils..."
curl --max-time 60 --retry 3 -L \
     https://github.com/usgs/earthquake-impact-utils/archive/master.zip\
     -o impactutils.zip
pip -q install --no-deps impactutils.zip
rm impactutils.zip


echo "Installing folium..."
pip install folium
echo "Installing pelican markdown..."

pip install pelican markdown

# This package
echo "Installing groundfailure..."
pip install -e .


# Tell the user they have to activate this environment
echo "Type 'source activate $VENV' to use this new virtual environment."

