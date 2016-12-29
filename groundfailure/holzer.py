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


def holzer_liq(shakefile, config, uncertfile=None, saveinputs=False,
               modeltype=None, displmodel=None,
               probtype=None, bounds=None):
    """
    Method for computing the probability of liquefaction using the Holzer method
    using the Wills et al. (2015) Vs30 map of California to define the
    susceptibility classes and the Fan et al. global water table model. 
    """
    layers = config['holzer_liq_cal']['layers']
    vs30_file = layers['vs30']['file']
    wtd_file = layers['watertable']['file']
    shkgdict = ShakeGrid.getFileGeoDict(shakefile)
    fgeodict = GMTGrid.getFileGeoDict(vs30_file)[0]

    
    #---------------------------------------------------------------------------
    # Loading info
    #---------------------------------------------------------------------------
    shakemap = ShakeGrid.load(shakefile, fgeodict, resample=True,
                              method='linear', doPadding=True)
    PGA = shakemap.getLayer('pga').getData()/100 # convert to g
    griddict,eventdict,specdict,fields,uncertainties = getHeaderData(shakefile)
    mag = eventdict['magnitude']


    #---------------------------------------------------------------------------
    # Logistic funciton parameters from Vs30
    #---------------------------------------------------------------------------
    vs30_grid = GMTGrid.load(vs30_file)

    vs30 = vs30_grid.getData()
    a0 = np.zeros_like(vs30)
    b0 = np.zeros_like(vs30)
    c0 = np.zeros_like(vs30)
    a1 = np.zeros_like(vs30)
    b1 = np.zeros_like(vs30)
    c1 = np.zeros_like(vs30)
    for k,v in config['holzer_liq_cal']['parameters'].items():
        ind = np.where(vs30 == float(v[0]))
        a0[ind] = v[1]
        b0[ind] = v[2]
        c0[ind] = v[3]
        a1[ind] = v[4]
        b1[ind] = v[5]
        c1[ind] = v[6]


    #---------------------------------------------------------------------------
    # Water table
    #---------------------------------------------------------------------------
    wtd_grid = GMTGrid.load(wtd_file, fgeodict, resample=True, 
                            method=layers['watertable']['interpolation'], 
                            doPadding = True)
    tmp = wtd_grid._data
    tmp = np.nan_to_num(tmp)

    # Compute water weights
    w0, w1 = get_water_weights(tmp)
    
    #---------------------------------------------------------------------------
    # Compute probability of liquefaction
    #---------------------------------------------------------------------------
    prob0 = get_prob(PGA, a0, b0, c0, mag)
    prob1 = get_prob(PGA, a1, b1, c1, mag)
    prob = prob0*w0 + prob1*w1

    #---------------------------------------------------------------------------
    # Turn output and inputs into into grids and put in maplayers dictionary
    #---------------------------------------------------------------------------
    maplayers = collections.OrderedDict()
    
    temp = shakemap.getShakeDict()
    shakedetail = '%s_ver%s' % (temp['shakemap_id'], temp['shakemap_version'])
    modelsref = config['holzer_liq_cal']['shortref']
    modellref = config['holzer_liq_cal']['longref']
    modeltype = 'Holzer/Wills'
    maplayers['model'] = {'grid': GDALGrid(prob, fgeodict), 
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


def get_prob(pga, a, b, c, M):
    """Compute probability of liquefaction from logistic function and 
    magnitude scaling factor. 
    """
    msf = 10**2.24/(M**2.56)
    prob = np.nan_to_num(a/(1 + ((pga/msf)/b)**c))
    return prob


def get_water_weights(z):
    """Compute weights for the two different water table depths
    for an arbitrary depth z. 
    """
    w0 = np.ones_like(z)
    w0[(z > 1.5) & (z <= 5)] = 1 + 1.5/3.5 - z[(z > 1.5) & (z <= 5)]/3.5
    w0[z > 5] = 0
    w1 = np.zeros_like(z)
    w1[(z > 1.5) & (z <= 5)] = 1 - w0[(z > 1.5) & (z <= 5)]
    w1[(z > 5) & (z <= 20)] = 1 + 5/15 - z[(z > 5) & (z <= 20)]/15
    return w0, w1
