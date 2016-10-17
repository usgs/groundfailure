
#!/usr/bin/env python

# Needs some additional editing to work, also need to make a fake config file

#python 3 compatibility
from __future__ import print_function
import os.path
#import sys
#stdlib imports
import os
#import tempfile
from configobj import ConfigObj
import numpy as np
import groundfailure.logisticmodel as LM
from groundfailure.conf import correct_config_filepaths

#hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
datadir = os.path.abspath(os.path.join(homedir, 'data'))
#logisticmodeldir = os.path.abspath(os.path.join(homedir, 'groundfailure'))
#sys.path.insert(0, logisticmodeldir)  # put this at the front of the system path, ignoring any installed mapio stuff

configfile = os.path.join(datadir, 'data', 'test_config.ini')
config = ConfigObj(configfile)
# Test path correction (from conf.py)
config = correct_config_filepaths(config)  # Need to make it correct filepaths for homedir
cmodel = config['logistic_models']['TestModelLS']

shakefile = os.path.join(datadir, 'data', 'test_shakegrid.xml')
uncertfile = os.path.join(datadir, 'test_uncert.xml')
cofile = os.path.join(datadir, 'test_cohesion.bil')
slopefile = os.path.join(datadir, 'test_slope.bil')
vs30file = os.path.join(datadir, 'test_vs30.bil')
ctifile = os.path.join(datadir, 'test_cti.bil')
precipfolder = os.path.join(datadir, 'test_precip')


def test_logisticmodel():
    print('Making sure the logistic model runs with and without uncertainty and precipitation - these files are 4x4 cells')
    modelLQ = {'logistic_models': {'TestModelLQ': {'description': 'This is a test liquefaction model',
                                                'gfetype': 'liquefaction',
                                                'baselayer': 'cohesion',
                                                'layers': {'vs30': '%s' % vs30file,
                                                           'cti': '%s' % ctifile},
                                                'interpolations': {'vs30': 'nearest',
                                                                   'cti': 'linear'},
                                                'terms': {'b1': 'log((pga/100.0)*(power(MW,2.56)/power(10,2.24)))',
                                                          'b2': 'cti',
                                                          'b3': 'log(vs30)'},
                                                'coefficients': {'b0': 24.10,
                                                                 'b1': 2.067,
                                                                 'b2': 0.355,
                                                                 'b3': -4.784}}}}

    modelLS = {'logistic_models': {'TestModelLS': {'description': 'This is a test landslide model',
                                                    'gfetype': 'landslide',
                                                    'baselayer': 'cohesion',
                                                    'layers': {'cohesion': '%s' % cofile,
                                                               'slope': '%s' % slopefile,
                                                               'precip': '%s' % precipfolder},
                                                    'interpolations': {'cohesion': 'linear',
                                                                       'slope': 'linear',
                                                                       'precip': 'nearest'},
                                                    'terms': {'b1': 'pga',
                                                              'b2': 'slope',
                                                              'b3': 'precipMONTH',
                                                              'b4': 'pga*slope*MW'},
                                                    'coefficients': {'b0': -7.15,
                                                                     'b1': 0.0604,
                                                                     'b2': 0.000825,
                                                                     'b3': 0.0201,
                                                                     'b4': 1.45e-05}}}}

    ls = LM(modelLS, shakefile, 'nowicki_2014_global', uncertfile=None)
    print(ls.getEquation())
    ls.calculate()

    lsu = LM(modelLS, shakefile, 'nowicki_2014_global', uncertfile)
    print(lsu.getEquation())
    print(lsu.getEquations())
    print(lsu.getGeoDict())
    lsu.calculate()

    lq = LM(modelLQ, shakefile, 'zhu_2015', uncertfile=None)
    print(lq.getEquation())
    lq.calculate()

    # Check if results are as expected by manual calculation
    np.testing.assert_allclose(ls['model'].getData(), targetLS)
    np.testing.assert_allclose(lsu['model'].getData(), targetLSU)  # Need to check one of the uncertainties at least
    np.testing.assert_allclose(lq['model'].getData(), targetLQ)


# test getLogisticModelNames(config):
def test_getLogisticModelNames():
    print('Testing Name retrieval from config file:')
    # declare config/shakefiles file variable

    data = ['nowicki_2015']
    names = LM.getLogisticModelNames(config)
    if data == names:
        print('Test passed.\n')
    else:
        print('Test not passed')
        print('Data: ', data, '\nNames: ', names, '\n')

    return names


def test_validateCoefficients(cmodel):
    print('Test Coefficient validation from config file:')
    data = {'b0': -8.3453199, 'b1': 1.737721, 'b2': 2.1773966, 'b3': 1.0, 'b4': 0.0484136, 'b5': 0.1634385, 'b6': 0.000949, 'b7': 1.0, 'b8': 0.0002273, 'b9': 0.477635}
    coeff = LM.validateCoefficients(cmodel)
    if data == coeff:
        print('Test passed.\n')
    else:
        print('Test not passed.')
        print('Data: ', data, '\nCoeff: ', coeff, '\n')

    return coeff


def test_validateLayers(cmodel):
    print('Test Layer validation from config file:')
    data = {'slope': 'gted_maxslope_30c.flt', 'rock': 'glim_copy_2.grd', 'landcover': 'modis_30c_copy_2.grd', 'precip': "", 'cti': "", 'elev': ""}
    layers = LM.validateLayers(cmodel)
    if data == layers:
        print('Test passed.\n')
    else:
        print('Test not passed.')
        print('Data: ', data, '\nLayers: ', layers, '\n')

    return layers


def test_validateTerms(cmodel, coeffs, layers):
    print('Test Term validation from config file:')
    data = {'b1': 'np.log(pgv)', 'b2': 'slope / 90.', 'b3': 'rock', 'b4': 'cti', 'b5': 'MW', 'b6': 'precipMONTH', 'b7': 'landcover', 'b8': 'elev', 'b9': 'np.log(PGV) * slope / 90.'}
    terms, time = LM.validateTerms(cmodel, coeffs, layers)
    if data == terms:
        print('Test passed. \n')
    else:
        print('Test not passed.')
        print('Data: ', data, '\nTerms: ', terms, '\n')

    return terms, time


def test_validateInterpolations(cmodel, layers):
    print('Test Interpolation validation from config file:')
    data = {'slope': 'linear', 'rock': 'nearest', 'landcover': 'nearest', 'precip': 'nearest', 'cti': 'linear', 'elev': 'linear'}
    interp = LM.validateInterpolations(cmodel, layers)
    if data == interp:
        print('Test passed. \n')
    else:
        print('Test not passed.')
        print('Data: ', data, '\nInterp: ', interp, '\n')

    return interp


def test_validateUnits(cmodel, layers):
    print('Test Unit validation from config file:')
    data = {'slope': 'unitless', 'rock': 'unitless', 'landcover': 'unitless', 'precip': 'mm/month', 'cti': 'unitless', 'elev': 'meters'}
    units = LM.validateUnits(cmodel, layers)
    if data == units:
        print('Test passed.\n')
    else:
        print('Test not passed.')
        print('Data: ', data, '\nUnits: ', units, '\n')

    return units


def test_validateLogisticModels():
    pass


def test_validateRefs():
    pass


def test_checkTerm():
    pass
