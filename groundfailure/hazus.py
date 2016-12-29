#!/usr/bin/env python

#stdlib imports
import os.path
import warnings
import collections

#local imports
from mapio.shake import ShakeGrid
from mapio.shake import getHeaderData
from mapio.gdal import GDALGrid
from mapio.gmt import GMTGrid
from mapio.grid2d import Grid2D
from mapio.geodict import GeoDict

#third party imports
import numpy as np


def hazus_liq(shakefile, config, uncertfile=None, saveinputs=False,
              modeltype=None, displmodel=None,
              probtype=None, bounds=None):
    """
    Method for computing the probability of liquefaction using the Hazus method
    using the Wills et al. (2015) Vs30 map of California to define the
    susceptibility classes and the Fan et al. global water table model. 
    """
    layers = config['hazus_liq_cal']['layers']
    vs30_file = layers['vs30']['file']
    wtd_file = layers['watertable']['file']
    shkgdict = ShakeGrid.getFileGeoDict(shakefile)
    fgeodict = GMTGrid.getFileGeoDict(vs30_file)[0]

    
    #---------------------------------------------------------------------------
    # Loading
    #---------------------------------------------------------------------------
    shakemap = ShakeGrid.load(shakefile, fgeodict, resample=True,
                              method='linear', doPadding=True)
    PGA = shakemap.getLayer('pga').getData()/100 # convert to g
    griddict,eventdict,specdict,fields,uncertainties = getHeaderData(shakefile)
    mag = eventdict['magnitude']

    # Correction factor for moment magnitudes other than M=7.5
    k_m = 0.0027*mag**3 - 0.0267*mag**2 - 0.2055*mag + 2.9188

    #---------------------------------------------------------------------------
    # Susceptibility from Vs30
    #---------------------------------------------------------------------------
    vs30_grid = GMTGrid.load(vs30_file)

    vs30 = vs30_grid.getData()
    p_ml = np.zeros_like(vs30)
    a = np.zeros_like(vs30)
    b = np.zeros_like(vs30)
    for k,v in config['hazus_liq_cal']['parameters'].items():
        ind = np.where(vs30 == float(v[0]))
        if v[1] == "VH":
            p_ml[ind] = 0.25
            a[ind] = 9.09
            b[ind] = -0.82
        if v[1] == "H":
            p_ml[ind] = 0.2
            a[ind] = 7.67
            b[ind] = -0.92
        if v[1] == "M":
            p_ml[ind] = 0.1
            a[ind] = 6.67
            b[ind] = -1.0
        if v[1] == "L":
            p_ml[ind] = 0.05
            a[ind] = 5.57
            b[ind] = -1.18
        if v[1] == "VL":
            p_ml[ind] = 0.02
            a[ind] = 4.16
            b[ind] = -1.08

    # Conditional liquefaction probability for a given susceptibility category 
    # at a specified PGA
    p_liq_pga = a*PGA + b
    p_liq_pga = p_liq_pga.clip(min = 0, max = 1)

    #---------------------------------------------------------------------------
    # Water table
    #---------------------------------------------------------------------------
    wtd_grid = GMTGrid.load(wtd_file, fgeodict, resample=True, 
                            method=layers['watertable']['interpolation'], 
                            doPadding = True)
    tmp = wtd_grid._data
    tmp = np.nan_to_num(tmp)

    # Convert to ft
    wt_ft = tmp * 3.28084

    # Correction factor for groundwater depths other than five feet
    k_w = 0.022*wt_ft + 0.93
    
    #---------------------------------------------------------------------------
    # Combine to get conditional liquefaction probability
    #---------------------------------------------------------------------------
    p_liq_sc = p_liq_pga * p_ml / k_m / k_w

    #---------------------------------------------------------------------------
    # Turn output and inputs into into grids and put in maplayers dictionary
    #---------------------------------------------------------------------------
    maplayers = collections.OrderedDict()
    
    temp = shakemap.getShakeDict()
    shakedetail = '%s_ver%s' % (temp['shakemap_id'], temp['shakemap_version'])
    modelsref = config['hazus_liq_cal']['shortref']
    modellref = config['hazus_liq_cal']['longref']
    modeltype = 'Hazus/Wills'
    maplayers['model'] = {'grid': GDALGrid(p_liq_sc, fgeodict), 
                          'label': 'Probability', 
                          'type': 'output',
                          'description': {'name': modelsref, 
                                          'longref': modellref, 
                                          'units': 'coverage',
                                          'shakemap': shakedetail, 
                                          'parameters': {'modeltype': modeltype}
                                          }
                          }

    if saveinputs is True:
        maplayers['pga'] = {'grid': GDALGrid(PGA, fgeodict), 
                            'label': 'PGA (g)', 
                            'type': 'input',
                            'description': {'units': 'g', 'shakemap': shakedetail}}
        maplayers['vs30'] = {'grid': GDALGrid(vs30, fgeodict),
                             'label': 'Vs30 (m/s)', 
                             'type': 'input',
                             'description': {'units': 'm/s'}}
        maplayers['wtd'] = {'grid': GDALGrid(wtd_grid._data, fgeodict),
                            'label': 'wtd (m)', 
                            'type': 'input',
                            'description': {'units': 'm'}}
    return maplayers
