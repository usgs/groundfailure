from mapio.gdal import GDALGrid
from mapio.geodict import GeoDict
import numpy as np
import os

# Make temporary files for each layer
cohesion = np.array([[5., 8.], [2., 9.]], dtype=float)  # kPa
slope = np.array([[10., 30.], [47., 15.]], dtype=float)  # degrees
vs30 = np.array([[244., 400.], [500., 1000.]], dtype=float)  # m/s
cti = np.array([[10., 5.], [3.5, 18.5]], dtype=float)
precip = np.array([[300., 876.], [180., 900.]], dtype=float)  # mm
units = ['kPa', 'degrees', 'm/s', 'unitless', 'mm']
terms = 
coefficients = 

pga = np.array([[24., 72.], [54., 110.]], dtype=float)  # %g
pgv = np.array([[50., 73.], [32., 84.]], dtype=float)  # cm/s
stdpga = np.array([[0.43, 0.49], [0.32, 0.16]], dtype=float)  # ln(%g)
stdpgv = np.array([[0.46, 0.58], [0.53, 0.33]], dtype=float)  # ln(cm/s)

geodict = GeoDict({'xmin': 0.0, 'xmax': 1.0,
                  'ymin': 0.0, 'ymax': 1.0,
                  'dx': 1.0, 'dy': 1.0,
                  'ny': 2, 'nx': 2})


def makeTestData():
    # make test layers
    X = ['cohesion', 'slope', 'vs30', 'cti', 'precip']
    config = {}
    config.setdefault('logistic_models', {}).setdefault('test_model_type', {})
    config['logistic_models']['test_model_type'].setdefault('coefficients', {})['b0'] = 3.5
    config['logistic_models']['test_model_type'].setdefault('terms', {})
    config['logistic_models']['test_model_type'].setdefault('layers', {})
    config['logistic_models']['test_model_type'].setdefault('interpolations', {})
    config['logistic_models']['test_model_type'].setdefault('units', {})

    for k, items in enumerate(X):
        coef = 'b%1d' % (k)
        # make a GDALGrid object (which is based on the Grid2D class)
        testgrid = GDALGrid(eval(items), geodict)
        # Save the file
        if items == 'precip':
            try:
                os.mkdir('test_precip')
            except:
                pass
            filename = 'test_precip/prec_Apr.bil' % (items)  # Only make April for testing
        else:
            filename = 'test_%s.bil' % (items)
            testgrid.save(filename, format='EHdr')

        # add to test config file
        config['logistic_models']['test_model_type']['layers'].update({items: {'file': filename}})
        config['logistic_models']['test_model_type']['interpolations'].update({items: 'nearest'})
        config['logistic_models']['test_model_type']['units'].update({items: units[k]})
        config['logistic_models']['test_model_type']['terms'].update({coef: terms[k]})
        config['logistic_models']['test_model_type']['coefficients'].update({coef: coefficients[k]})

    # Add full file extension to config for testing correct_config_filepaths

    # test_shakegrid.xml
    # test_uncert.xml

    #config['logistic_models']['test_model_type']['baselayer'] = X[0]


    print(config)
    return config

makeTestData()

# fake shakemap file
#f = open('shakemap.xml','w')


## Load it back in
#testgrid2 = GDALGrid.load('test.bil')
#print testgrid

#header
#    ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>')
#shakemap_grid specs
#    ('<shakemap_grid xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns="http://earthquake.usgs.gov/eqcenter/shakemap" xsi:schemaLocation="http://earthquake.usgs.gov http://earthquake.usgs.gov/eqcenter/shakemap/xml/schemas/shakemap.xsd" event_id="1" shakemap_id="1" shakemap_version="2" code_version="3.5.1446" process_timestamp="2016-07-20T20:00:00Z" shakemap_originator="us" map_status="RELEASED" shakemap_event_type="ACTUAL">')
#event_data
#    ('<event event_id="1" magnitude="7" depth="20" lat="34.0" lon="-118.0" event_timestamp="1990-01-01T12:30:00UTC" event_network="us" event_description="Test, California" />')
#grid specs
#    ('''''<grid_specification lon_min="-119.0" lat_min="33.0" lon_max="-117.0" lat_max="35.0" nominal_lon_spacing="1.0" nominal_lat_spacing="1" nlon="2" nlat="2" /><event_specific_uncertainty name="pga" value="0.5" numsta="4" />
#        <event_specific_uncertainty name="pgv" value="0.5" numsta="4" />
 #       <event_specific_uncertainty name="mi" value="0.5" numsta="4" />
 #       <event_specific_uncertainty name="psa03" value="0.5" numsta="4" />
 #       <event_specific_uncertainty name="psa10" value="0.5" numsta="4" />
 #       <event_specific_uncertainty name="psa30" value="0.5" numsta="4" />
 #       <grid_field index="1" name="LON" units="dd" />
 #       <grid_field index="2" name="LAT" units="dd" />
 #       <grid_field index="3" name="PGA" units="pctg" />
 #       <grid_field index="4" name="PGV" units="cms" />
 #       <grid_field index="5" name="MMI" units="intensity" />
 #       <grid_field index="6" name="PSA03" units="pctg" />
 #       <grid_field index="7" name="PSA10" units="pctg" />
 #       <grid_field index="8" name="PSA30" units="pctg" />
 #       <grid_field index="9" name="STDPGA" units="ln(pctg)" />
 #       <grid_field index="10" name="URAT" units="" />
 #       <grid_field index="11" name="SVEL" units="ms" />''''')
#grid data
#    ('''<grid_data>
#        -117.5 33.5 2.0 6.0 4.0 5.0 4.0 1.0 0.5 1.0 425
#        -117.5 34.5 2.0 6.0 4.0 5.0 4.0 1.0 0.5 1.0 425
#        -118.5 33.5 2.0 6.0 4.0 5.0 4.0 1.0 0.5 1.0 555
#        -118.5 34.5 2.0 6.0 4.0 5.0 4.0 1.0 0.5 1.0 555
#        </grid_data>''')
#footer
#    ('</shakemap_grid>')
