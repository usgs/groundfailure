
#!/usr/bin/env python

#python 3 compatibility
from __future__ import print_function
import os.path
import sys
#stdlib imports
import abc
import textwrap
import glob
import os
import tempfile
from configobj import ConfigObj

#hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
logisticmodeldir = os.path.abspath(os.path.join(homedir, 'groundfailure'))
sys.path.insert(0, logisticmodeldir)  # put this at the front of the system path, ignoring any installed mapio stuff

# further imports
#stdlib imports
import numpy as np
import re
import groundfailure.logisticmodel as lm
from groundfailure.conf import correct_config_filepaths
import groundfailure.gdalfuncs as gdal
import groundfailure.savelayers as sl
import groundfailure.testmodels as tm
#third party imports
from mapio.shake import ShakeGrid, getHeaderData
from mapio.gmt import GMTGrid
from mapio.gdal import GDALGrid
from mapio.grid2d import Grid2D

# Make temporary files for each layer
slope = np.full((2, 2), 270., dtype=float)
rock = np.array([[-1.7010214, -0.7960331], [-0.6563501, -1.6414169]], dtype=float)
landcover = np.array([[1.0466108, 1.0706804], [0.9440571, 0.7256497]], dtype=float)
precip = np.full((2, 2), 3., dtype=float)
cti = np.full((2, 2), 4., dtype=float)
elev = np.full((2, 2), 5., dtype=float)
output = np.array([[1.34140885, 2.27046675], [2.28354445, 1.08005225]], dtype=float)

# save temp layers to files
test_temp_dir = os.mkdir('Users/kbiegel/Documents/GroundFailure/inputs/temp_test')
slopefile.save()
##########
######### Finish this section ##################
##########

# declare config/shakefiles file variable
configfile = 'Users/kbiegel/Documents/GroundFailure/config_test.ini'
shakefile = 'Users/kbiegel/Documents/GroundFailure/inputs/text.xml'

# validate config file
config = ConfigObj(configfile)
config = groundfailure.conf.correct_config_filepaths(config)

# Load in shakemap data
edict = ShakeGrid.load(shakefile, adjust='res').getEventDict()
temp = ShakeGrid.load(shakefile, adjust='res').getShakeDict()
edict['eventid'] = temp['shakemap_id']
edict['version'] = temp['shakemap_version']


# test getLogisticModelNames(config):
def test_getLogisticModelNames():
    print('Testing Name retrieval from config file:')
    data = ['nowicki_2015']
    names = lm.getLogisticModelNames(config)
    if data == names:
        print('Test passed.\n')
    else:
        print('Test not passed')
        print('Data: ', data, '\nNames: ', names, '\n')

    return names


def test_validateCoefficients(cmodel):
    print('Test Coefficient validation from config file:')
    data = {'b0': -8.3453199, 'b1': 1.737721, 'b2': 2.1773966, 'b3': 1.0, 'b4': 0.0484136, 'b5': 0.1634385, 'b6': 0.000949, 'b7': 1.0, 'b8': 0.0002273, 'b9': 0.477635}
    coeff = lm.validateCoefficients(cmodel)
    if data == coeff:
        print('Test passed.\n')
    else:
        print('Test not passed.')
        print('Data: ', data, '\nCoeff: ', coeff, '\n')

    return coeff


def test_validateLayers(cmodel):
    print('Test Layer validation from config file:')
    data = {'slope': 'gted_maxslope_30c.flt', 'rock': 'glim_copy_2.grd', 'landcover': 'modis_30c_copy_2.grd', 'precip': "", 'cti': "", 'elev': ""}
    layers = lm.validateLayers(cmodel)
    if data == layers:
        print('Test passed.\n')
    else:
        print('Test not passed.')
        print('Data: ', data, '\nLayers: ', layers, '\n')

    return layers


def test_validateTerms(cmodel, coeffs, layers):
    print('Test Term validation from config file:')
    data = {'b1': 'np.log(pgv)', 'b2': 'slope / 90.', 'b3': 'rock', 'b4': 'cti', 'b5': 'MW', 'b6': 'precipMONTH', 'b7': 'landcover', 'b8': 'elev', 'b9': 'np.log(PGV) * slope / 90.'}
    terms, time = lm.validateTerms(cmodel, coeffs, layers)
    if data == terms:
        print('Test passed. \n')
    else:
        print('Test not passed.')
        print('Data: ', data, '\nTerms: ', terms, '\n')

    return terms, time


def test_validateInterpolations(cmodel, layers):
    print('Test Interpolation validation from config file:')
    data = {'slope': 'linear', 'rock': 'nearest', 'landcover': 'nearest', 'precip': 'nearest', 'cti': 'linear', 'elev': 'linear'}
    interp = lm.validateInterpolations(cmodel, layers)
    if data == interp:
        print('Test passed. \n')
    else:
        print('Test not passed.')
        print('Data: ',data,'\nInterp: ',terms,'\n')

    return interp


def test_validateUnits(cmodel, layers):
    print('Test Unit validation from config file:')
    data = {'slope': 'unitless', 'rock': 'unitless', 'landcover': 'unitless', 'precip': 'mm/month', 'cti': 'unitless', 'elev': 'meters'}
    units = lm.validateUnites(cmodel, layers)
    if data == units:
        print('Test passed.\n')
    else:
        print('Test not passed.')
        print('Data: ',data,'\nUnits: ',units,'\n')

    return units


def test_LogisticModel_SelfEquation():


def test_LogisticModel_calculate(self):

















