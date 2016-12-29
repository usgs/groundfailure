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


def slhrf_liq(shakefile, config, uncertfile=None, saveinputs=False,
              modeltype=None, displmodel=None,
              probtype=None, bounds=None):
    """
    Method for computing the probability of liquefaction using the SLHRF,
    primarily relying on the Wills et al. (2015) Vs30 map of California and
    Hydrosheds distance to rivers. 
    """
    layers = config['slhrf_liq_cal']['layers']
    vs30_file = layers['vs30']['file']
    elev_file = layers['elev']['file']
    dc_file = layers['dc']['file']
    dr_file = layers['dr']['file']
    fgeodict = GMTGrid.getFileGeoDict(vs30_file)[0]

    
    #---------------------------------------------------------------------------
    # Read in data layers
    #---------------------------------------------------------------------------
    shakemap = ShakeGrid.load(shakefile, fgeodict, resample=True,
                              method='linear', doPadding=True)
    PGA = shakemap.getLayer('pga').getData()/100 # convert to g
    griddict,eventdict,specdict,fields,uncertainties = getHeaderData(shakefile)
    mag = eventdict['magnitude']
    vs30_grid = GMTGrid.load(vs30_file)
    vs30 = vs30_grid.getData()
    elev = GDALGrid.load(elev_file, fgeodict, resample=True,
                        method=layers['elev']['interpolation'],
                        doPadding = True).getData()
    dc = GDALGrid.load(dc_file, fgeodict, resample=True,
                       method=layers['dc']['interpolation'],
                       doPadding = True).getData()
    dr = GDALGrid.load(dr_file, fgeodict, resample=True,
                       method=layers['dr']['interpolation'],
                       doPadding = True).getData()
    dw = np.minimum(dr, dc)


    #---------------------------------------------------------------------------
    # Evaluate the different factors
    #---------------------------------------------------------------------------
    Fgeo = np.zeros_like(vs30)
    for k,v in config['slhrf_liq_cal']['parameters'].items():
        ind = np.where(vs30 == float(v[0]))
        Fgeo[ind] = float(v[1])
    Fz = z_factor(elev)
    Fmag = mag_factor(mag)
    Fpga = pga_factor(PGA)
    Fdw = dw_factor(dw)

    
    #---------------------------------------------------------------------------
    # Combine factors
    #---------------------------------------------------------------------------
    SLHRF = Fz * Fmag * Fpga * Fdw * Fgeo

    # Transform into a 'probability'
    prob = 0.4 * (1 - np.exp(-0.2 * SLHRF**2) )

    #---------------------------------------------------------------------------
    # Turn output and inputs into into grids and put in maplayers dictionary
    #---------------------------------------------------------------------------
    maplayers = collections.OrderedDict()
    
    temp = shakemap.getShakeDict()
    shakedetail = '%s_ver%s' % (temp['shakemap_id'], temp['shakemap_version'])
    modelsref = config['slhrf_liq_cal']['shortref']
    modellref = config['slhrf_liq_cal']['longref']
    modeltype = 'SLHRF/Wills'
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
        maplayers['slhrf'] = {'grid': GDALGrid(SLHRF, fgeodict),
                              'label': 'SLHRF', 
                              'type': 'input',
                              'description': {'units': 'none'}}
        maplayers['pga'] = {'grid': GDALGrid(PGA, fgeodict), 
                            'label': 'PGA (g)', 
                            'type': 'input',
                            'description': {'units': 'g', 'shakemap': shakedetail}}
        maplayers['vs30'] = {'grid': GDALGrid(vs30, fgeodict),
                             'label': 'Vs30 (m/s)', 
                             'type': 'input',
                             'description': {'units': 'm/s'}}
        maplayers['dw'] = {'grid': GDALGrid(dw, fgeodict),
                           'label': 'dw (km)', 
                           'type': 'input',
                           'description': {'units': 'km'}}
        maplayers['elev'] = {'grid': GDALGrid(elev, fgeodict),
                             'label': 'elev (m)', 
                             'type': 'input',
                             'description': {'units': 'm'}}
        maplayers['FPGA'] = {'grid': GDALGrid(Fpga, fgeodict),
                             'label': 'Fpga', 
                             'type': 'input',
                             'description': {'units': 'none'}}
        maplayers['FDW'] = {'grid': GDALGrid(Fdw, fgeodict),
                            'label': 'Fdw', 
                            'type': 'input',
                            'description': {'units': 'none'}}
        maplayers['FGEO'] = {'grid': GDALGrid(Fgeo, fgeodict),
                             'label': 'Fgeo', 
                             'type': 'input',
                             'description': {'units': 'none'}}
        maplayers['FZ'] = {'grid': GDALGrid(Fz, fgeodict),
                           'label': 'Fz', 
                           'type': 'input',
                           'description': {'units': 'none'}}
    return maplayers


def z_factor(z):
    fac = np.ones_like(z) * 1.2
    fac[z > 100.0] = 1.0
    fac[z > 200.0] = 0.9
    return fac

def mag_factor(mag):
    mag = np.array([mag])
    fac = np.ones_like(mag) * 0.8
    fac[mag > 6] = 0.9
    fac[mag > 7] = 1.0
    fac[mag > 8] = 1.1
    fac[mag > 9] = 1.2
    return fac

def pga_factor(pga):
    # get rid of nans
    pga = np.nan_to_num(pga)
    fac = np.ones_like(pga) * 0.05
    fac[pga > 0.1] = 0.25
    fac[pga > 0.2] = 0.5
    fac[pga > 0.3] = 0.75
    fac[pga > 0.4] = 1.0
    return fac

def dw_factor(dw):
    fac = np.ones_like(dw) * 1.2
    fac[dw > 0.1] = 1.0
    fac[dw > 0.5] = 0.9
    return fac


