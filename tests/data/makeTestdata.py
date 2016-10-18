from mapio.geodict import GeoDict
from mapio.gdal import GDALGrid
from mapio.shake import ShakeGrid
import numpy as np
import os
from datetime import datetime
from collections import OrderedDict
from configobj import ConfigObj

# Make temporary files for each layer
friction = np.array([[5., 8.], [2., 9.]], dtype=float)  # kPa
slope = np.array([[1000., 3000.], [4700., 1500.]], dtype=float)  # degrees*100
vs30 = np.array([[244., 400.], [500., 1000.]], dtype=float)  # m/s
cti1 = np.array([[10., 5.], [3.5, 18.5]], dtype=float)
precip = np.array([[300., 876.], [180., 90.]], dtype=float)  # mm
units = ['degrees', 'degrees', 'm/s', 'unitless', 'mm']
terms = ['friction', 'slope/100.', 'log(vs30)', 'cti1', 'precipMONTH']
coefficients = [0.3, 2.1, -0.5, 0.01, 1.0]

pga = np.array([[24., 72.], [54., 11.]], dtype=float)  # %g
pgv = np.array([[50., 73.], [32., 84.]], dtype=float)  # cm/s
stdpga = np.array([[0.43, 0.49], [0.32, 0.16]], dtype=float)  # ln(%g)

geodict = GeoDict({'xmin': 0.0, 'xmax': 1.0,
                  'ymin': 0.0, 'ymax': 1.0,
                  'dx': 1.0, 'dy': 1.0,
                  'ny': 2, 'nx': 2})


def makeTestData():
    # make test layers
    X = ['friction', 'slope', 'vs30', 'cti1', 'precip']
    config = OrderedDict()
    config.setdefault('logistic_models', {}).setdefault('test_model', {})
    config['logistic_models']['test_model'].setdefault('shortref', 'Name et al. year')
    config['logistic_models']['test_model'].setdefault('longref', 'full reference')
    config['logistic_models']['test_model'].setdefault('layers', {})
    config['logistic_models']['test_model'].setdefault('interpolations', {})
    config['logistic_models']['test_model'].setdefault('terms', {})
    config['logistic_models']['test_model'].setdefault('coefficients', {})['b0'] = 3.5

    for k, items in enumerate(X):
        coef = 'b%1d' % (k+1)
        # make a GDALGrid object
        testgrid = GDALGrid(eval(items), geodict)
        # Save the file
        if items == 'precip':
            try:
                os.mkdir('test_precip')
            except:
                pass
            filename = 'test_precip/prec_Jan.bil'  # Only make January for testing
        else:
            filename = 'test_%s.bil' % (items)
        testgrid.save(filename, format='EHdr')

        # add to test config file
        config['logistic_models']['test_model']['layers'].update({items: {'file': filename.split('/')[0], 'units': units[k], 'longref': 'longref', 'shortref': 'shortref'}})
        config['logistic_models']['test_model']['interpolations'].update({items: 'nearest'})
        config['logistic_models']['test_model']['terms'].update({coef: terms[k]})
        config['logistic_models']['test_model']['coefficients'].update({coef: coefficients[k]})

    config['logistic_models']['test_model']['gfeype'] = 'landslide'
    config['logistic_models']['test_model']['baselayer'] = 'slope'
    config['logistic_models']['test_model']['slopemin'] = 5.
    config['logistic_models']['test_model']['slopemax'] = 90.

    # Make test_shakegrid and test_uncert
    eventDict = OrderedDict([('event_id', 'test'),
                            ('lon', 0.5),
                            ('lat', 0.5),
                            ('event_timestamp', datetime(2000, 1, 5, 0, 30, 55)),
                            ('event_network', 'na'),
                            ('magnitude', 6.0),
                            ('event_description', 'Test event'),
                            ('depth', 5.0)])
    shakeDict = OrderedDict([('process_timestamp',
                            datetime(2000, 1, 6, 20, 38, 19)),
                            ('event_id', 'test'),
                            ('shakemap_version', 2),
                            ('code_version', '1 billion'),
                            ('shakemap_event_type', 'TEST'),
                            ('map_status', 'TEST'),
                            ('shakemap_id', 'test'),
                            ('shakemap_originator', 'na')])
    uncertaintyDict = {}

    layers1 = {'pga': pga, 'pgv': pgv}
    shakegrid = ShakeGrid(layers1, geodict, eventDict, shakeDict, uncertaintyDict)
    shakegrid.save('test_shakegrid.xml')

    layers2 = {'stdpga': stdpga}
    uncertgrid = ShakeGrid(layers2, geodict, eventDict, shakeDict, uncertaintyDict)
    uncertgrid.save('test_uncert.xml')

    C = ConfigObj(config)
    C.filename = 'test.ini'
    C.write()

    return config

makeTestData()
