#!/usr/bin/env python

#python 3 compatibility
from __future__ import print_function
import os.path
import sys
#stdlib imports
import os
from collections import OrderedDict
import datetime
import numpy as np
#third party imports
from mapio.shake import ShakeGrid
from mapio.gdal import GDALGrid
from mapio.geodict import GeoDict
import groundfailure.logisticmodel as lm

################# Not sure this is necessary - Mike originally said to include, but my formatting is different than his.

#hack the path so functions can be debugged if needed
homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
logisticmodeldir = os.path.abspath(os.path.join(homedir, 'groundfailure'))
sys.path.insert(0, logisticmodeldir)  # put this at the front of the system path, ignoring any installed mapio stuff

#################


def makeData():
    """
    Function creates simplified fake shakemap data and config data (along with file locations for layer data) in order to run a test simulaiton of the logistic model class functions

    Takes no parameters (must be edited to change data outputs)
    """

    #Make up a geodictionary with required fields that matches the data size
    geodict = GeoDict({'xmin': 0.0, 'xmax': 1.0, 'ymin': 0.0, 'ymax': 1.0, 'dx': 1.0, 'dy': 1.0, 'ny': 2, 'nx': 2})

    ########
    # Fake shakemap file creation - hand edit to change values - requires layers, geodict, eventDict, shakeDict, uncertaintyDict
    ########
    new_date = datetime.datetime(2016, 1, 1, 1, 1, 1)
    # layer creation
    layers = OrderedDict()
    layers['lat'] = np.array([[0., 1.0], [0., 1.0]])
    layers['lon'] = np.array([[0., 0.], [1.0, 1.0]])
    layers['pga'] = np.array([[1.3, 1.6], [1.9, 2.3]], dtype=float)
    layers['pgv'] = np.array([[3.4, 3.7], [4.0, 4.3]], dtype=float)
    layers['mmi'] = np.array([[2.2, 2.3], [2.4, 2.5]], dtype=float)
    layers['psa03'] = np.arange(1, 5.).reshape((2, 2))
    layers['psa10'] = np.arange(1, 5.).reshape((2, 2))
    layers['psa30'] = np.arange(1, 5.).reshape((2, 2))
    layers['stdpga'] = np.full((2, 2), 0.5)
    layers['urat'] = np.full((2, 2), 1.0)
    layers['svel'] = np.full((2, 2), 425)
    # eventDict creation
    eventdict = {'event_id': 'TEST0001', 'magnitude': 6.1, 'depth': 11.1, 'lat': 0.5, 'lon': 0.5, 'event_timestamp': new_date, 'event_network': 'US', 'event_description': 'Test, US'}
    # shakedict creation
    shakedict = {'event_id': 'TEST0001', 'shakemap_id': '1', 'shakemap_version': 2, 'code_version': '2016', 'process_timestamp': new_date, 'shakemap_originator': 'us', 'map_status': 'Released', 'shakemap_event_type': 'SCENARIO'}
    # uncertaintyDict creation
    uncertaintydict = {'pga': (0.5, 4), 'pgv': (0.5, 4), 'mi': (0.5, 4), 'psa03': (0.5, 4), 'psa10': (0.5, 4), 'psa30': (0.5, 4)}
    # create and save shakegrid
    shakegrid = ShakeGrid(layers, geodict, eventdict, shakedict, uncertaintydict)
    shakegridpath = 'test_shakegrid.xml'
    shakegrid.save(shakegridpath)

    #######
    # fake config file and fake layer data creation
    #######
    # variable definitions
    X = ['x1', 'x2', 'x3', 'x4']
    filelocations = {}
    config = {}
    model_type = 'test_model_type'
    # preliminary config layout
    config['logistic_models'] = {}
    config['logistic_models'][model_type] = {}
    config['logistic_models'][model_type]['layers'] = {}
    config['logistic_models'][model_type]['interpolations'] = {}
    config['logistic_models'][model_type]['units'] = {}
    config['logistic_models'][model_type]['baselayer'] = X[0]

    # for loop - creates data for layers and adds information about them to the config
    for items in X:
        # Make fake data (cannot be integers for GDALGrid)
        if items == X[0]:
            data = np.arange(2, 6.).reshape((2, 2))
        if items == X[1]:
            data = np.array([[1.4, 1.6], [1.55, 1.34]], dtype=float)
        if items == X[2]:
            data = np.arange(-1, 3.).reshape((2, 2))
        # Use data and geodict to create layer grids
        testgrid = GDALGrid(data, geodict)
        # Save the file and add to config
        filepath = 'layer_%s.bil' % (items)
        testgrid.save(filepath, format='EHdr')
        # document file location
        filelocations.update({items: filepath})
        # fake config file
        config['logistic_models'][model_type]['layers'].update({items: {'file': filelocations[items]}})
        config['logistic_models'][model_type]['interpolations'].update({items: 'nearest'})
        config['logistic_models'][model_type]['units'].update({items: 'unitless'})

    # create coefficient and term OrderedDicts to add to config - make sure they are ordered or some tests will not pass
    # variable definitions
    coeffs = OrderedDict()
    terms = OrderedDict()
    b = ['b0', 'b1', 'b2', 'b3', 'b4']
    # loop over all coefficients
    for item in b:
        if item == 'b0':
            coeffs[item] = 3.5
        elif item == 'b1':
            terms[item] = 'log(x1)'
            coeffs[item] = 1.5
        elif item == 'b2':
            terms[item] = 'x2'
            coeffs[item] = -2.5
        elif item == 'b3':
            terms[item] = 'log(x1) * x3 / 90.'
            coeffs[item] = 4.0
        else:
            terms[item] = 'pgv'
            coeffs[item] = -1.0
    # add coefficients and terms to config dictionary
    config['logistic_models'][model_type]['coefficients'] = coeffs
    config['logistic_models'][model_type]['terms'] = terms

    # end makedata function
    return config, filelocations, model_type


###############
# Begin testing functions
###############

def test_getLogisticModelNames():
    """
    Runs getLogisticModelNames function and compares to expected values
    Pass test if values are the same
    """
    print('Testing Name retrieval from config file:')
    # retrieve expected values
    data = list(config['logistic_models'].keys())
    # run function being tested
    names = lm.getLogisticModelNames(config)
    # compare values
    if data == names:
        print('Test passed.\n')
    else:
        print('Test not passed')
        print('Data: ', data, '\nNames: ', names, '\n')


def test_validateCoefficients(cmodel):
    """
    Runs validateCoefficients function and compares to expected values
    Pass test if values are the same
    :param cmodel:
        dictionary - config['logistic_models'][model_type]
    """
    print('Test Coefficient validation from config file:')
    # retrieve expected values
    data = cmodel['coefficients']
    # run funciton being tested
    coeff = lm.validateCoefficients(cmodel)
    # compate values
    if data == coeff:
        print('Test passed.\n')
    else:
        print('Test not passed.')
        print('Data: ', data, '\nCoeff: ', coeff, '\n')

    # return coefficients for later test
    coeff_test = data
    return coeff, coeff_test


def test_validateLayers(cmodel):
    """
    Runs validatelayers function and compares to expected values
    Pass test if values are the same
    :param cmodel:
        dictionary - config['logistic_models'][model_type]
    """
    print('Test Layer validation from config file:')
    # retrieve expected values
    data = {}
    for keys in cmodel['layers'].keys():
        data.update({keys: cmodel['layers'][keys]['file']})
    # run function being tested
    layers = lm.validateLayers(cmodel)
    # compare values
    if data == layers:
        print('Test passed.\n')
    else:
        print('Test not passed.')
        print('Data: ', data, '\nLayers: ', layers, '\n')

    # return layers for later test
    return layers


def test_validateTerms(cmodel, coeffs, layers):
    """
    Runs validateTerms function and compares to expected values
    Pass test if values are the same
    :param cmodel:
        dictionary - config['logistic_models'][model_type]
    :param coeffs:
        dictionary - value returned from running validateCoefficients function from logistic model class
        imported from test_validateCoefficients function
    :param layers:
        dictionary - value returned from running validateLayers function from logistic model class
        imported from test_validateLayers
    """
    print('Test Term validation from config file:')
    # retrieve expected values
    data_precheck = cmodel['terms']
    data = {}
    for key, value in data_precheck.items():
        data_checked, rem, tTimefield = lm.checkTerm(value, layers)
        data.update({key: data_checked})
    # run function being tested
    terms, time = lm.validateTerms(cmodel, coeffs, layers)
    # compare values
    if data == terms:
        print('Test passed. \n')
    else:
        print('Test not passed.')
        print('Data: ', data, '\nTerms: ', terms, '\n')

    # return terms for later test
    term_test = data
    return term_test


def test_validateInterpolations(cmodel, layers):
    """
    Runs validateInterpolations function and compares to expected values
    Pass test if values are the same
    :param cmodel:
        dictionary - config['logistic_models'][model_type]
    :param layers:
        dictionary - value returned from running validateLayers function from logistic model class
        imported from test_validateLayers
    """
    print('Test Interpolation validation from config file:')
    # retrieve expected values
    data = cmodel['interpolations']
    # fun function being tested
    interp = lm.validateInterpolations(cmodel, layers)
    # compares values
    if data == interp:
        print('Test passed. \n')
    else:
        print('Test not passed.')
        print('Data: ', data, '\nInterp: ', interp, '\n')


def test_validateUnits(cmodel, layers):
    """
    Runs validateUnits function and compares to expected values
    Pass test if values are the same
    :param cmodel:
        dictionary - config['logistic_models'][model_type]
    :param layers:
        dictionary - value returned from running validateLayers function from logistic model class
        imported from test_validateLayers
    """
    print('Test Unit validation from config file:')
    # retrieve expected values
    data = cmodel['units']
    # run function being tested
    units = lm.validateUnits(cmodel, layers)
    # compare values
    if data == units:
        print('Test passed.\n')
    else:
        print('Test not passed.')
        print('Data: ', data, '\nUnits: ', units, '\n')


def test_LogisticModel_SelfEquation(config, shakefile, coeff, terms):
    """
    Creates logistic model class element equation and compares to expected values
    Pass test if values are the same
    :param config:
        dictionary - config dictionary created in makeData function
    :param shakefile:
        string - filepath location to ShakeGrid object
    :param coeff:
        dictionary - value returned from running validateCoefficients function from logistic model class
        imported from test_validateCoefficients
    :param terms:
        dictionary - values returned from running validateTerms function from logistic model class
        imported from test_validateTerms
    """
    print('Test lm class equation:')
    # retrieve expected values
    data_nuggets = []
    for items, values in coeff.items():
        if items == 'b0':
            data_nuggets.append(str(values))
        else:
            data_nuggets.append('(%g * %s)' % (values, terms[items]))
    data = ' + '.join(data_nuggets)
    # run function being tested
    L = lm.LogisticModel(config, shakefile, model_type)
    equation = L.getEquation()
    # compare values
    if data == equation:
        print('Test passed. \n')
    else:
        print('Test not passed.')
        print('Data: ', data, '\nEquation: ', equation, '\n')


def test_LogisticModel_calculate(config, shakefile, filelocations):
    """
    Creates logistic model class element equation and compares to expected values
    Pass test if values are the same
    :param config:
        dictionary - created in makeData function
    :param shakefile:
        string - filepath location to ShakeGrid object
    :param filelocations:
        dictionary - containing layer keys and location values
    """
    print('Test lm class calculate:')
    # retrieve expected values
    shakemap = ShakeGrid.load(shakefile)
    grid_data = {}
    pre_calculate_data = {}
    calculated_data = {}
    for items, values in filelocations.items():
        grid_data.update({items: GDALGrid.load(values).getData()})
    grid_data.update({'x4': shakemap.getLayer('pgv').getData()})
    for key in grid_data.keys():
        if key == 'x1':
            variable = np.log(grid_data[key])
            pre_calculate_data.update({'b1': variable})
        if key == 'x2':
            pre_calculate_data.update({'b2': grid_data[key]})
        if key == 'x3':
            variable = np.log(grid_data['x1']) * grid_data[key] / 90.
            pre_calculate_data.update({'b3': variable})
        if key == 'x4':
            pre_calculate_data.update({'b4': grid_data[key]})
    coeffs = config['logistic_models']['test_model_type']['coefficients']
    for key in coeffs.keys():
        if key == 'b0':
            calculated_data.update({key: coeff[key]})
        else:
            variable = coeff[key] * pre_calculate_data[key]
            calculated_data.update({key: variable})
    x = sum(calculated_data.values())
    data = 1 / (1 + np.exp(-x))

    # run function being tested
    L = lm.LogisticModel(config, shakefile, model_type)
    model_grid = L.calculate()

    # compare values
    try:
        np.testing.assert_allclose(data, model_grid['model']['grid'].getData(), rtol=1e-4)
        print('Test passed. \n')
    except:
        print('Test not passed.')
        print('Data: ', data, '\nModel grid: ', model_grid['model']['grid'].getData(), '\n')

    import pdb; pdb.set_trace()


# runs all test functions
if __name__ == '__main__':
    config, filelocations, model_type = makeData()
    shakefile = 'test_shakegrid.xml'
    cmodel = config['logistic_models'][model_type]
    test_getLogisticModelNames()
    coeff, coeff_test = test_validateCoefficients(cmodel)
    layers = test_validateLayers(cmodel)
    term_test = test_validateTerms(cmodel, coeff, layers)
    test_validateInterpolations(cmodel, layers)
    test_validateUnits(cmodel, layers)
    test_LogisticModel_SelfEquation(config, shakefile, coeff_test, term_test)
    test_LogisticModel_calculate(config, shakefile, filelocations)
