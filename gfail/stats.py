#!/usr/bin/env python
# -*- coding: utf-8 -*-

# stdlib imports
import numpy as np
import collections

# local imports
from mapio.shake import ShakeGrid

# Make fonts readable and recognizable by illustrator
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = ['Arial',
                                   'Bitstream Vera Serif',
                                   'sans-serif']


def computeStats(grid2D, probthresh=0.0, shakefile=None,
                shakethreshtype='pga', shakethresh=0.0):
    """
    Compute summary stats of a ground failure model output.

    Args:
        grid2D: grid2D object of model output.
        probthresh: Optional, Float or list of probability thresholds
            for use in Hagg2 computation.
        shakefile: Optional, path to shakemap file to use for ground motion
            threshold.
        shakethreshtype: Optional, Type of ground motion to use for
            shakethresh, 'pga', 'pgv', or 'mmi'.
        shakethresh: Optional, Float or list of shaking thresholds in %g for
            pga, cm/s for pgv, float for mmi. Used for Hagg computation

    Returns:
        Dictionary with the following keys:
            - Max
            - Median
            - Std
            - Hagg_# where # is the shaking threshold for each
            - Parea_# where # is the probability threshold
    """
    stats = collections.OrderedDict()
    grid = grid2D.getData()
    
    stats['Max'] = np.nanmax(grid)
    stats['Median'] = np.nanmedian(grid)
    stats['Std'] = np.nanstd(grid)
    Hagg = computeHagg(grid2D, probthresh=0.0, shakefile=shakefile,
                shakethreshtype=shakethreshtype, shakethresh=shakethresh)
    if type(Hagg) != list and type(Hagg) != list:
        shakethresh = [shakethresh]
        Hagg = [Hagg]
            
    for T, H in zip(shakethresh, Hagg):
        if T == 0.:
            stats['Hagg'] = H
        else:
            newkey = 'Hagg_%1.0f%%g' % T
            stats[newkey] = H

    Parea = computeParea(grid2D, probthresh=probthresh, shakefile=None,
                         shakethreshtype=shakethreshtype, shakethresh=0.0)
    if type(Parea) != list and type(Parea) != np.ndarray:
        probthresh = [probthresh]
        Parea = [Parea]
        
    for T, P in zip(probthresh, Parea):
        if T == 0.:
            stats['Parea'] = P
        else:
            newkey = 'Parea_%1.2f' % T
            stats[newkey] = P

    return stats


def computeHagg(grid2D, proj='moll', probthresh=0.0, shakefile=None,
                shakethreshtype='pga', shakethresh=0.0):
    """
    Computes the Aggregate Hazard (Hagg) which is equal to the
    probability * area of grid cell For models that compute areal coverage,
    this is equivalant to the total predicted area affected in km2.

    Args:
        grid2D: grid2D object of model output.
        proj: projection to use to obtain equal area, 'moll'  mollweide, or
            'laea' lambert equal area.
        probthresh: Probability threshold, any values less than this will not
            be included in aggregate hazard estimation.
        shakefile: Optional, path to shakemap file to use for ground motion
            threshold.
        shakethreshtype: Optional, Type of ground motion to use for
            shakethresh, 'pga', 'pgv', or 'mmi'.
        shakethresh: Optional, Float or list of shaking thresholds in %g for
            pga, cm/s for pgv, float for mmi.

    Returns:
        Float if no shakethresh defined or only one shakethresh defined,
        otherwise, a list of aggregate hazard for all shakethresh values.
    """
    Hagg = []
    bounds = grid2D.getBounds()
    lat0 = np.mean((bounds[2], bounds[3]))
    lon0 = np.mean((bounds[0], bounds[1]))
    projs = ('+proj=%s +lat_0=%f +lon_0=%f +x_0=0 +y_0=0 +ellps=WGS84 '
             '+units=km +no_defs' % (proj, lat0, lon0))
    geodict = grid2D.getGeoDict()

    if shakefile is not None:
        if type(shakethresh) != list and type(shakethresh) != np.ndarray:
            shakethresh = [shakethresh]
        for shaket in shakethresh:
            if shaket < 0.:
                raise Exception('shaking threshold must be equal or greater '
                                'than zero')
        # resample shakemap to grid2D
        temp = ShakeGrid.load(shakefile, samplegeodict=geodict, resample=True,
                              doPadding=True, adjust='res')
        shk = temp.getLayer(shakethreshtype)
        if shk.getGeoDict() != geodict:
            raise Exception('shakemap was not resampled to exactly the same '
                            'geodict as the model')

    if probthresh < 0.:
        raise Exception('probability threshold must be equal or greater '
                        'than zero')

    grid = grid2D.project(projection=projs)
    geodictRS = grid.getGeoDict()
    cell_area_km2 = geodictRS.dx * geodictRS.dy
    model = grid.getData()
    model[np.isnan(model)] = -1.
    if shakefile is not None:
        for shaket in shakethresh:
            modcop = model.copy()
            shkgrid = shk.project(projection=projs)
            shkdat = shkgrid.getData()
            # use -1 to avoid nan errors and warnings, will always be thrown
            # out because default is 0.
            shkdat[np.isnan(shkdat)] = -1.
            modcop[shkdat < shaket] = -1.
            Hagg.append(np.sum(modcop[modcop >= probthresh] * cell_area_km2))
    else:
        Hagg.append(np.sum(model[model >= probthresh] * cell_area_km2))
    if len(Hagg) == 1:
        Hagg = Hagg[0]
    return Hagg


def computeParea(grid2D, proj='moll', probthresh=0.0, shakefile=None,
                 shakethreshtype='pga', shakethresh=0.0):
    """
    Alternative to Aggregate Hazard (Hagg), which is equal to the
    the sum of the area of grid cells that exceeds a given probability.

    Args:
        grid2D: grid2D object of model output.
        proj: projection to use to obtain equal area, 'moll'  mollweide, or
            'laea' lambert equal area.
        probthresh: Optional, Float or list of probability thresholds.
        shakefile: Optional, path to shakemap file to use for ground motion
            threshold.
        shakethreshtype: Optional, Type of ground motion to use for
            shakethresh, 'pga', 'pgv', or 'mmi'.
        shakethresh: Optional, Float of shaking thresholds in %g for
            pga, cm/s for pgv, float for mmi.

    Returns:
        Float if no probthresh defined or only one probthresh defined,
        otherwise, a list of Parea for all probthresh values.
    """
    if type(probthresh) != list and type(probthresh) != np.ndarray:
        probthresh = [probthresh]

    Parea = []
    bounds = grid2D.getBounds()
    lat0 = np.mean((bounds[2], bounds[3]))
    lon0 = np.mean((bounds[0], bounds[1]))
    projs = ('+proj=%s +lat_0=%f +lon_0=%f +x_0=0 +y_0=0 +ellps=WGS84 '
             '+units=km +no_defs' % (proj, lat0, lon0))
    geodict = grid2D.getGeoDict()

    if shakefile is not None:
        if shakethresh < 0.:
            raise Exception('shaking threshold must be equal or greater '
                            'than zero')
        # resample shakemap to grid2D
        temp = ShakeGrid.load(shakefile, samplegeodict=geodict, resample=True,
                              doPadding=True, adjust='res')
        shk = temp.getLayer(shakethreshtype)
        if shk.getGeoDict() != geodict:
            raise Exception('shakemap was not resampled to exactly the same '
                            'geodict as the model')

    grid = grid2D.project(projection=projs)
    geodictRS = grid.getGeoDict()
    cell_area_km2 = geodictRS.dx * geodictRS.dy
    model = grid.getData()
    model[np.isnan(model)] = -1.
    for probt in probthresh:
        if probt < 0.:
            raise Exception('probability threshold must be equal or greater '
                            'than zero')
        modcop = model.copy()
        if shakefile is not None:
            shkgrid = shk.project(projection=projs)
            shkdat = shkgrid.getData()
            # use -1 to avoid nan errors and warnings, will always be thrown
            # out because default probthresh is 0 and must be positive.
            shkdat[np.isnan(shkdat)] = -1.
            modcop[shkdat < shakethresh] = -1.
        one_mat = np.ones_like(modcop)
        Parea.append(np.sum(one_mat[modcop >= probt] * cell_area_km2))

    if len(Parea) == 1:
        Parea = Parea[0]
    return Parea
