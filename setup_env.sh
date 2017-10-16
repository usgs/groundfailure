#!/bin/sh
echo $PATH

VENV=gf
PYVER=3.5


DEPARRAY=(pytables numpy==1.12.1 scipy=0.19.0 pip matplotlib=1.5.3 jupyter=1.0.0 \
          rasterio=0.36.0 fiona=1.7.6 basemap=1.1.0 basemap-data-hires=1.1.0 \
          shapely=1.5.17 h5py=2.7.0 gdal=2.1.3 descartes=1.1.0 pytest-cov=2.5.1 \
          pytest-mpl=0.7 configobj=5.0.6 pandas=0.20.2 sphinx=1.6.3 \
          scikit-learn=0.18.2 scikit-image=0.13.0 ipython=6.1.0 \
          branca=0.2.0 paramiko=2.1.2)

# turn off whatever other virtual environment user might be in
source activate root

#remove any previous virtual environments called gf
CWD=`pwd`
cd $HOME;
conda remove --name $VENV --all -y
cd $CWD

conda create --name $VENV --yes --channel conda-forge python=$PYVER ${DEPARRAY[*]} -y

if [ $? -eq 1 ];then
    echo "Environment creation failed - look at error message from conda create above."
    exit 1
fi

# activate the new environment
source activate $VENV

#install some items separately
conda install -y psutil=5.2.1

#do pip installs of those things that are not available via conda.
#we're using curl to fetch these zip files instead of pip because some
#of our systems fail on the pip command - using curl gives us more control
#over how long we wait for a download to complete, and how many tries before giving up
#on the download.

#download openquake, install it using pip locally, ignore specified dependencies,
#as these should be installed using conda above
curl --max-time 60 --retry 3 -L https://github.com/gem/oq-engine/archive/v2.5.0.zip -o openquake.zip
pip -v install --no-deps openquake.zip
rm openquake.zip

#download MapIO, install it using pip locally
curl --max-time 60 --retry 3 -L https://github.com/usgs/MapIO/archive/0.6.2.zip -o mapio.zip
pip install mapio.zip
rm mapio.zip

#download impactutils, install it using pip locally
curl --max-time 60 --retry 3 -L https://github.com/usgs/earthquake-impact-utils/archive/master.zip -o impact.zip
pip install impact.zip
rm impact.zip


# do pip installs of those things that are not available via conda.
pip install sphinx_rtd_theme
pip install folium
pip install pelican markdown


# tell the user they have to activate this environment
echo "Type 'source activate gf' to use this new virtual environment."
