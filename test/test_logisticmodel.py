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
from mapio.geodict import GeoDict


def makeData():
    # fake layers
    X = ['x1', 'x2', 'x3']
    fileLocations = {}
    config = {}
    config.setdefault('logistic_models', {}).setdefault('test_model_type', {})
    config['logistic_models']['test_model_type'].setdefault('coefficients', {})['b0'] = 3.5
    config['logistic_models']['test_model_type'].setdefault('terms', {})
    config['logistic_models']['test_model_type'].setdefault('layers', {})
    config['logistic_models']['test_model_type'].setdefault('interpolations', {})
    config['logistic_models']['test_model_type'].setdefault('units', {})

    for items in X:
        # Make fake data (cannot be integers for GDALGrid)
        data = np.arange(2,6.).reshape((2,2))
        #Make up a geodictionary with required fields that matches the data size
        geodict = GeoDict({'xmin':0.0,'xmax':1.0,
                        'ymin':0.0,'ymax':1.0,
                        'dx':1.0,'dy':1.0,
                        'ny':2,'nx':2})
        # Use these two pieces to make a GDALGrid object (which is based on the Grid2D class)
        testgrid = GDALGrid(data,geodict)
        print(testgrid)

        # Save the file
        filepath = 'layer_%s' % (items)
        testgrid.save(filepath, format='EHdr')

        # document file location
        filelocations = {}
        filelocations = {items: filepath}

        # fake config file
        config['logistic_models']['test_model_type']['layers'].update({items: {'file': filelocations[items]}})
        config['logistic_models']['test_model_type']['interpolations'].update({items: 'nearest'})
        config['logistic_models']['test_model_type']['units'].update({items: 'unitless'})
        config['logistic_models']['test_model_type']['baselayer'] = X[0]
        if items == X[0]:
            config['logistic_models']['test_model_type']['terms'].update({'b1': 'log(%s)' % items})
            config['logistic_models']['test_model_type']['coefficients'].update({'b1': 1.5})
        if items == X[1]:
            config['logistic_models']['test_model_type']['terms'].update({'b2': items})
            config['logistic_models']['test_model_type']['coefficients'].update({'b2': 2.5})
        if items == X[2]:
            config['logistic_models']['test_model_type']['terms'].update({'b3': 'log(%s) * %s / 90.' % (X[0], items)})
            config['logistic_models']['test_model_type']['coefficients'].update({'b3': 4.0})

    print(config)
    return config


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


def test_validateInterpolations(cmodel, layers):
    print('Test Interpolation validation from config file:')
    data = {'slope': 'linear', 'rock': 'nearest', 'landcover': 'nearest', 'precip': 'nearest', 'cti': 'linear', 'elev': 'linear'}
    interp = lm.validateInterpolations(cmodel, layers)
    if data == interp:
        print('Test passed. \n')
    else:
        print('Test not passed.')
        print('Data: ',data,'\nInterp: ',terms,'\n')


def test_validateUnits(cmodel, layers):
    print('Test Unit validation from config file:')
    data = {'slope': 'unitless', 'rock': 'unitless', 'landcover': 'unitless', 'precip': 'mm/month', 'cti': 'unitless', 'elev': 'meters'}
    units = lm.validateUnites(cmodel, layers)
    if data == units:
        print('Test passed.\n')
    else:
        print('Test not passed.')
        print('Data: ',data,'\nUnits: ',units,'\n')


def test_LogisticModel_SelfEquation():
    # Would test equation created in the logistic model class.  Not sure how to test for this however?
    pass


def test_LogisticModel_calculate(self):
    # would take the equation and then test the calculation of the actual grid values.  Not sure how to test this either.
    pass


if __name__ == '__main__':
    config = makeData()
    cmodel = config['logistic_models']['test_model_type']
    test_getLogisticModelNames()
    coeff = test_validateCoefficients(cmodel)
    layers = test_validateLayers(cmodel)
    test_validateTerms(cmodel, coeffs, layers)
    test_validateInterpolations(cmodel, layers)
    test_validateUnits(cmodel, layers)
