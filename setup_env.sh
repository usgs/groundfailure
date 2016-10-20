#!/bin/sh
echo $PATH

VENV=gf
PYVER=3.5


DEPARRAY=(numpy scipy pip matplotlib rasterio fiona basemap basemap-data-hires shapely h5py gdal descartes pytest-cov pytest-mpl configobj pandas sphinx scikit-learn psutil ipython)

# turn off whatever other virtual environment user might be in
source deactivate

#remove any previous virtual environments called pager
CWD=`pwd`
cd $HOME;
conda remove --name $VENV --all -y
cd $CWD

conda create --name $VENV --yes --channel conda-forge python=$PYVER ${DEPARRAY[*]} -y

# activate the new environment
source activate $VENV

# do pip installs of those things that are not available via conda.
pip -v install https://github.com/usgs/MapIO/archive/master.zip
pip -v install https://github.com/usgs/earthquake-impact-utils/archive/master.zip
pip -v install https://github.com/gem/oq-hazardlib/archive/master.zip
pip install scikit-image
pip install sphinx_rtd_theme

# tell the user they have to activate this environment
echo "Type 'source activate gf' to use this new virtual environment."