
#!/usr/bin/env python

# Needs some additional editing to work, also need to make a fake config file

#python 3 compatibility
from __future__ import print_function
import os.path
import sys
#stdlib imports
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
import groundfailure.logisticmodel as LM
from groundfailure.conf import correct_config_filepaths

# Make temporary files for each layer
slope = np.full((2, 2), 270., dtype=float)
rock = np.array([[-1.7010214, -0.7960331], [-0.6563501, -1.6414169]], dtype=float)
landcover = np.array([[1.0466108, 1.0706804], [0.9440571, 0.7256497]], dtype=float)
precip = np.full((2, 2), 3., dtype=float)
cti = np.full((2, 2), 4., dtype=float)
elev = np.full((2, 2), 5., dtype=float)
output = np.array([[1.34140885, 2.27046675], [2.28354445, 1.08005225]], dtype=float)

# declare config/shakefiles file variable
configfile = os.path.join(logisticmodeldir, 'data', 'test_config.ini')
shakefile = os.path.join(logisticmodeldir, 'data', 'test_shakegrid.xml')

# validate config file
config = ConfigObj(configfile)
config = correct_config_filepaths(config)


def _test(shakefile, cofile, slopefile, precipfolder):
    model = {'logistic_models': {'nowicki_2014': {'description': 'This is the Nowicki Model of 2014, which uses cohesion and slope max as input.',
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

    lm = LM(model, shakefile, 'nowicki_2014_global')
    print(lm.getEquation())
    P = lm.calculate()


# test getLogisticModelNames(config):
def test_getLogisticModelNames():
    print('Testing Name retrieval from config file:')
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
    units = LM.validateUnites(cmodel, layers)
    if data == units:
        print('Test passed.\n')
    else:
        print('Test not passed.')
        print('Data: ', data, '\nUnits: ', units, '\n')

    return units


def test_LogisticModel_SelfEquation():
    pass


def test_LogisticModel_calculate(self):
    pass
