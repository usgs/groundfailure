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
slope = np.array([[10., 30.], [47., 15.]], dtype=float)  # degrees
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
    # make test layers and config files (model, and mapping)
    X = ['friction', 'slope', 'vs30', 'cti1', 'precip']

    configModel = OrderedDict()

    # Work on model configs
    configModel.setdefault('test_model', {})
    configModel['test_model'].setdefault('shortref', 'Name et al. year')
    configModel['test_model'].setdefault('longref', 'full reference')
    configModel['test_model'].setdefault('layers', {})
    configModel['test_model'].setdefault('interpolations', {})
    configModel['test_model'].setdefault('terms', {})
    configModel['test_model'].setdefault('coefficients', {})['b0'] = '3.5'
    configModel['test_model'].setdefault('display_options', {'lims': {'model': 'np.linspace(0., 0.4, 10)'},
                                                             'colors': {'default': 'cm.inferno', 'alpha': '0.7',
                                                                        'model': 'cm.jet'},
                                                             'logscale': {'model': 'False'},
                                                             'maskthresholds': {'model': '0.001'}})

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
        configModel['test_model']['layers'].update({items: {'file': filename.split('/')[0],
                                                            'units': units[k], 'longref': 'longref',
                                                            'shortref': 'shortref'}})
        configModel['test_model']['interpolations'].update({items: 'nearest'})
        configModel['test_model']['terms'].update({coef: terms[k]})
        configModel['test_model']['coefficients'].update({coef: coefficients[k]})
        configModel['test_model']['display_options']['lims'].update({items: 'None'})
        configModel['test_model']['display_options']['colors'].update({items: 'None'})
        configModel['test_model']['display_options']['logscale'].update({items: 'False'})
        configModel['test_model']['display_options']['maskthresholds'].update({items: 'None'})

    configModel['test_model']['gfeype'] = 'landslide'
    configModel['test_model']['baselayer'] = 'slope'
    configModel['test_model']['slopemin'] = 5.
    configModel['test_model']['slopemax'] = 90.
    configModel['test_model']['display_options'].setdefault('pga', {'lims': 'None', 'colors': 'cm.jet',
                                                                    'logscale': 'True'})
    configModel['test_model']['display_options'].setdefault('pgv', {'lims': 'None', 'colors': 'cm.jet',
                                                                    'logscale': 'False'})

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

    C = ConfigObj(configModel)
    C.filename = 'testconfig_logimodel.ini'
    C.write()

    configMap = OrderedDict({'dem': {'file': 'None'},
                             'roads': {'file': 'None'},
                             'cities': {'file': 'None'},
                             'ocean': {'file': 'None'},
                             'colors': {'roadcolor': '808080',
                                        'countrycolor': '474747',
                                        'watercolor': 'B8EEFF'}})
    C = ConfigObj(configMap)
    C.filename = 'testconfig_map.ini'
    C.write()


makeTestData()
