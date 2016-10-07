from mapio.gdal import GDALGrid
from mapio.geodict import GeoDict
import numpy as np
import tempfile


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


makeData()

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
