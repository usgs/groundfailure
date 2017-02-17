from mapio.geodict import GeoDict
from mapio.gdal import GDALGrid
from mapio.shake import ShakeGrid
import numpy as np
import os
from datetime import datetime
from collections import OrderedDict
from configobj import ConfigObj

# Make temporary files for each layer
friction = np.array([[25., 30.], [14., 20.]], dtype=float)  # degrees
cohesion = np.array([[5., 8.], [12., 9.]])  # kPa
slope = np.array([[10., 30.], [47., 15.]], dtype=float)  # degrees
slopeGodt = np.array([[1000., 3000.], [4700., 1500.]], dtype=float)  # degrees
slopequants = [0.15, 0.3, 0.5, 0.6, 0.7, 0.8, 1.0]  # multipliers for quants in godt model
slopequantnames = ['_min', '10', '30', '50', '70', '90', '_max']
vs30 = np.array([[244., 400.], [500., 1000.]], dtype=float)  # m/s
cti1 = np.array([[10., 5.], [3.5, 18.5]], dtype=float)
precip = np.array([[300., 876.], [180., 90.]], dtype=float)  # mm
susceptibility = np.array([[1., 1., 1., 1., 1.],
                          [1., 3., 4., 6., 1.],
                          [1., 3., 4., 6., 1.],
                          [8., 10., 10., 10., 1.],
                          [8., 10., 10., 10., 1.]])
watertable = np.array([[10., 0.1], [0.3, 2.]])
units = {'friction': 'degrees', 'cohesion': 'kPa', 'slope': 'degrees',
         'vs30': 'm/s', 'cti1': 'unitless', 'precip': 'mm', 'susceptibility': 'category', 'watertable': 'm',
         'pga': 'g', 'pgv': 'cm/s'}

pga = np.array([[24., 72.], [54., 11.]], dtype=float)  # %g
pgv = np.array([[50., 73.], [32., 84.]], dtype=float)  # cm/s
stdpga = np.array([[0.43, 0.49], [0.32, 0.16]], dtype=float)  # ln(%g)
stdpgv = np.array([[0.43, 0.49], [0.32, 0.16]], dtype=float)  # ln(%g)

geodict = GeoDict({'xmin': 0.5, 'xmax': 1.5,
                  'ymin': 0.5, 'ymax': 1.5,
                  'dx': 1., 'dy': 1.,
                  'ny': 2, 'nx': 2})

susgeodict = GeoDict({'xmin': 0., 'xmax': 2.0,
                      'ymin': 0., 'ymax': 2.0,
                      'dx': 0.5, 'dy': 0.5,
                      'ny': 5, 'nx': 5})


def makeTestData():
    # make test layers and config files (model, and mapping)
    Xmod = ['friction', 'slope', 'vs30', 'cti1', 'precip']
    terms = ['friction', 'slope/100.', 'log(vs30)', 'cti1', 'precipMONTH']
    coefficients = [0.3, 2.1, -0.5, 0.01, 1.0]
    Xhaz = ['susceptibility']
    Xgod = ['cohesion', 'friction', 'slope']
    Xcla = ['cohesion', 'friction', 'slope', 'watertable']

    configModel = OrderedDict()
    configHazus = OrderedDict()
    configGodt = OrderedDict()
    configClassic = OrderedDict()

    # Work on model configs

    #LOGISTIC
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
    configModel['test_model']['gfeype'] = 'landslide'
    configModel['test_model']['baselayer'] = 'slope'
    configModel['test_model']['slopemin'] = 5.
    configModel['test_model']['slopemax'] = 90.
    configModel['test_model']['display_options'].setdefault('pga', {'lims': 'None', 'colors': 'cm.jet',
                                                                    'logscale': 'True'})
    configModel['test_model']['display_options'].setdefault('pgv', {'lims': 'None', 'colors': 'cm.jet',
                                                                    'logscale': 'False'})

    # HAZUS
    configHazus.setdefault('hazus', {})
    configHazus['hazus'].setdefault('shortref', 'Name et al. year')
    configHazus['hazus'].setdefault('longref', 'full reference')
    configHazus['hazus']['gfeype'] = 'landslide'
    configHazus['hazus'].setdefault('layers', {})
    configHazus['hazus'].setdefault('parameters', {'dnthresh': 5.})
    configHazus['hazus'].setdefault('display_options', {'lims': {'model': 'np.linspace(0., 0.4, 10)'},
                                                        'colors': {'default': 'cm.inferno', 'alpha': '0.7',
                                                                   'model': 'cm.jet'},
                                                        'logscale': {'model': 'False'},
                                                        'maskthresholds': {'model': '0.001'}})
    configHazus['hazus']['display_options'].setdefault('susceptibility', {'lims': 'None', 'colors': 'cm.jet',
                                                                          'logscale': 'True'})

    # GODT
    configGodt.setdefault('godt_2008', {})
    configGodt['godt_2008'].setdefault('shortref', 'Name et al. year')
    configGodt['godt_2008'].setdefault('longref', 'full reference')
    configGodt['godt_2008']['gfeype'] = 'landslide'
    configGodt['godt_2008'].setdefault('layers', {})
    configGodt['godt_2008'].setdefault('parameters', {'thick': '2.4',
                                                      'uwt': '15.7',
                                                      'nodata_cohesion': '1.0',
                                                      'nodata_friction': '26.',
                                                      'dnthresh': '5.',
                                                      'fsthresh': '1.01',
                                                      'acthresh': '0.05'})
    configGodt['godt_2008'].setdefault('display_options', {'lims': {'model': 'np.linspace(0., 0.4, 10)'},
                                                           'colors': {'default': 'cm.inferno', 'alpha': '0.7',
                                                                      'model': 'cm.jet'},
                                                           'logscale': {'model': 'False'},
                                                           'maskthresholds': {'model': '0.001'}})
    configGodt['godt_2008']['display_options'].setdefault('pga', {'lims': 'None', 'colors': 'cm.jet',
                                                                  'logscale': 'True'})

    # NEWMARK
    configClassic.setdefault('classic_newmark', {})
    configClassic['classic_newmark'].setdefault('shortref', 'Name et al. year')
    configClassic['classic_newmark'].setdefault('longref', 'full reference')
    configClassic['classic_newmark']['gfeype'] = 'landslide'
    configClassic['classic_newmark'].setdefault('layers', {})
    configClassic['classic_newmark'].setdefault('parameters', {'thick': '2.4',
                                                               'uwt': '15.7',
                                                               'nodata_cohesion': '1.0',
                                                               'nodata_friction': '26.',
                                                               'dnthresh': '5.',
                                                               'fsthresh': '1.01',
                                                               'acthresh': '0.05',
                                                               'slopethresh': '5.',
                                                               'm': '0.5'})
    configClassic['classic_newmark'].setdefault('display_options', {'lims': {'model': 'np.linspace(0., 0.4, 10)'},
                                                                    'colors': {'default': 'cm.inferno', 'alpha': '0.7',
                                                                               'model': 'cm.jet'},
                                                                    'logscale': {'model': 'False'},
                                                                    'maskthresholds': {'model': '0.001'}})
    configClassic['classic_newmark']['display_options'].setdefault('pga', {'lims': 'None', 'colors': 'cm.jet',
                                                                           'logscale': 'True'})

    for k, items in enumerate(Xmod):
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
                                                            'units': units[items], 'longref': 'longref',
                                                            'shortref': 'shortref'}})
        configModel['test_model']['interpolations'].update({items: 'nearest'})
        configModel['test_model']['terms'].update({coef: terms[k]})
        configModel['test_model']['coefficients'].update({coef: coefficients[k]})
        configModel['test_model']['display_options']['lims'].update({items: 'None'})
        configModel['test_model']['display_options']['colors'].update({items: 'None'})
        configModel['test_model']['display_options']['logscale'].update({items: 'False'})
        configModel['test_model']['display_options']['maskthresholds'].update({items: 'None'})

    for k, items in enumerate(Xhaz):
        # Save the file if not already saved
        filename = 'test_%s.bil' % (items)
        if items == 'susceptibility':
            testgrid = GDALGrid(eval(items), susgeodict)
            testgrid.save(filename, format='EHdr')

        # add to test config file
        configHazus['hazus']['layers'].update({items: {'file': filename,
                                                       'units': units[items], 'longref': 'longref',
                                                       'shortref': 'shortref'}})
        configHazus['hazus']['display_options']['lims'].update({items: 'None'})
        configHazus['hazus']['display_options']['colors'].update({items: 'None'})
        configHazus['hazus']['display_options']['logscale'].update({items: 'False'})
        configHazus['hazus']['display_options']['maskthresholds'].update({items: 'None'})

    for k, items in enumerate(Xgod):
        # make a GDALGrid object
        if items == 'slope':
            testgrid = GDALGrid(eval('slopeGodt'), geodict)
        else:
            testgrid = GDALGrid(eval(items), geodict)
        # Save the file
        filename = 'test_%s.bil' % (items)
        word = 'file'
        if items == 'slope':
            word = 'filepath'
            try:
                os.mkdir('test_slope_quantiles_godt')
            except:
                pass
            for j, slp in enumerate(slopequants):
                filename = 'test_slope_quantiles_godt/slope%s.bil' % (slopequantnames[j])  # Only make January for testing
                temp = GDALGrid(eval('slopeGodt') * slopequants[j], geodict)
                temp.save(filename, format='EHdr')
            filename = 'test_slope_quantiles_godt'
        elif items == 'cohesion':
            testgrid.save(filename, format='EHdr')

        # add to test config file
        configGodt['godt_2008']['layers'].update({items: {word: filename,
                                                          'units': units[items], 'longref': 'longref',
                                                          'shortref': 'shortref'}})
        configGodt['godt_2008']['display_options']['lims'].update({items: 'None'})
        configGodt['godt_2008']['display_options']['colors'].update({items: 'None'})
        configGodt['godt_2008']['display_options']['logscale'].update({items: 'False'})
        configGodt['godt_2008']['display_options']['maskthresholds'].update({items: 'None'})

    for k, items in enumerate(Xcla):
        # make a GDALGrid object
        testgrid = GDALGrid(eval(items), geodict)
        # Save the file
        filename = 'test_%s.bil' % (items)
        if items == 'watertable':
            testgrid.save(filename, format='EHdr')

        # add to test config file
        configClassic['classic_newmark']['layers'].update({items: {'file': filename,
                                                           'units': units[items], 'longref': 'longref',
                                                           'shortref': 'shortref'}})
        configClassic['classic_newmark']['display_options']['lims'].update({items: 'None'})
        configClassic['classic_newmark']['display_options']['colors'].update({items: 'None'})
        configClassic['classic_newmark']['display_options']['logscale'].update({items: 'False'})
        configClassic['classic_newmark']['display_options']['maskthresholds'].update({items: 'None'})

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

    layers2 = {'stdpga': stdpga, 'stdpgv': stdpgv}
    uncertgrid = ShakeGrid(layers2, geodict, eventDict, shakeDict, uncertaintyDict)
    uncertgrid.save('test_uncert.xml')

    C = ConfigObj(configModel)
    C.filename = 'testconfig_logimodel.ini'
    C.write()

    C = ConfigObj(configHazus)
    C.filename = 'testconfig_hazus.ini'
    C.write()

    C = ConfigObj(configClassic)
    C.filename = 'testconfig_classic.ini'
    C.write()

    C = ConfigObj(configGodt)
    C.filename = 'testconfig_godt.ini'
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
