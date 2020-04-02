#!/usr/bin/env python
# -*- coding: utf-8 -*-

# stdlib imports
import numpy as np
import collections
import shutil
import tempfile
import os


# local imports
from mapio.shake import ShakeGrid
from mapio.gdal import GDALGrid
from mapio.geodict import GeoDict
from gfail.spatial import quickcut
from mapio.grid2d import Grid2D
from skimage.measure import block_reduce

from configobj import ConfigObj

# Make fonts readable and recognizable by illustrator
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = ['Arial',
                                   'Bitstream Vera Serif',
                                   'sans-serif']


def computeStats(grid2D, stdgrid2D=None, probthresh=None, shakefile=None,
                 shakethreshtype='pga', shakethresh=0.0,
                 statprobthresh=None, pop_file=None):
    """
    Compute summary stats of a ground failure model output.

    Args:
        grid2D: grid2D object of model output.
        stdgrid2D: grid2D object of model standard deviations (optional)
        probthresh: Optional, Float or list of probability thresholds
            for use in Parea computation.
        shakefile: Optional, path to shakemap file to use for ground motion
            threshold.
        shakethreshtype: Optional, Type of ground motion to use for
            shakethresh, 'pga', 'pgv', or 'mmi'.
        shakethresh: Optional, Float or list of shaking thresholds in %g for
            pga, cm/s for pgv, float for mmi. Used for Hagg and Exposure
            computation
        statprobthresh: Optional, None or float, exclude any cells with
            probabilities less than or equal to this value
        pop_file (str): File path to population file to use to compute exposure
            stats

    Returns:
        dict: Dictionary with the following keys:
            - Max
            - Median
            - Std
            - Hagg_# where # is the shaking threshold
            - Parea_# where # is the probability threshold
            - exp_pop_# where # is the shaking threshold (if pop_file specified)
    """
    stats = collections.OrderedDict()
    grid = grid2D.getData().copy()
    if statprobthresh is not None:
        grid = grid[grid > statprobthresh]
    else:
        statprobthresh = 0.0

    if len(grid) == 0:
        print('no probability values above statprobthresh')
        stats['Max'] = 0.  # float('nan')
        stats['Median'] = 0.  # float('nan')
        stats['Std'] = 0.  # float('nan')
    else:
        stats['Max'] = float(np.nanmax(grid))
        stats['Median'] = float(np.nanmedian(grid))
        stats['Std'] = float(np.nanstd(grid))
    Hagg = computeHagg(grid2D, probthresh=statprobthresh, shakefile=shakefile,
                       shakethreshtype=shakethreshtype,
                       shakethresh=shakethresh, stdgrid2D=stdgrid2D)
    if type(Hagg) != list:
        shakethresh = [shakethresh]
        Hagg = [Hagg]

    for T, H in zip(shakethresh, Hagg):
        if T == 0.:
            stats['Hagg'] = float(H)
        else:
            newkey = 'Hagg_%1.2fg' % (T/100.)
            stats[newkey] = float(H)

    if probthresh is not None:
        Parea = computeParea(grid2D, probthresh=probthresh, shakefile=None,
                             shakethreshtype=shakethreshtype, shakethresh=0.0)
        if type(Parea) != list and type(Parea) != np.ndarray:
            probthresh = [probthresh]
            Parea = [Parea]

        for T, P in zip(probthresh, Parea):
            if T == 0.:
                stats['Parea'] = float(P)
            else:
                newkey = 'Parea_%1.2f' % T
                stats[newkey] = float(P)

    if pop_file is None:
        try:
            # Try to find population file in .gfail_defaults
            default_file = os.path.join(
                os.path.expanduser('~'), '.gfail_defaults')
            defaults = ConfigObj(default_file)
            pop_file = defaults['popfile']
        except:
            print('No population file specified nor found in .gfail_defaults, '
                  'skipping exp_pop')

    if pop_file is not None:
        exp_dict = get_exposures(grid2D, pop_file, shakefile=shakefile,
                                 shakethreshtype=shakethreshtype,
                                 shakethresh=shakethresh,
                                 probthresh=statprobthresh)
        for k, v in exp_dict.items():
            stats[k] = v
    
    if stdgrid2D is not None:
        pass

    return stats


def computeHagg(grid2D, proj='moll', probthresh=0.0, shakefile=None,
                shakethreshtype='pga', shakethresh=0.0, stdgrid2D=None):
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
        stdgrid2D: grid2D object of model standard deviations (optional)

    Returns: Aggregate hazard (float) if no shakethresh or only one shakethresh
        was defined, otherwise, a list of floats of aggregate hazard for all
        shakethresh values.
    """
    Hagg = []
    StdH = None
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
        tmpdir = tempfile.mkdtemp()
        # resample shakemap to grid2D
        temp = ShakeGrid.load(shakefile)
        junkfile = os.path.join(tmpdir, 'temp.bil')
        GDALGrid.copyFromGrid(temp.getLayer(shakethreshtype)).save(junkfile)
        shk = quickcut(junkfile, geodict, precise=True, method='bilinear')
        shutil.rmtree(tmpdir)
        if shk.getGeoDict() != geodict:
            raise Exception('shakemap was not resampled to exactly the same '
                            'geodict as the model')

    if probthresh < 0.:
        raise Exception('probability threshold must be equal or greater '
                        'than zero')

    grid = grid2D.project(projection=projs, method='bilinear')
    geodictRS = grid.getGeoDict()
    cell_area_km2 = geodictRS.dx * geodictRS.dy
    model = grid.getData().copy()
    if stdgrid2D is not None:
        stdgrid = stdgrid2D.project(projection=projs, method='bilinear')
        StdH = []
        std = stdgrid.getData().copy()
        std[np.isnan(model)] = -1.
    model[np.isnan(model)] = -1.
    if shakefile is not None:
        shkgrid = shk.project(projection=projs)
        for shaket in shakethresh:
            shkdat = shkgrid.getData()
            # use -1 to avoid nan errors and warnings, will always be thrown
            # out because default probthresh is 0.
            model[np.isnan(shkdat)] = -1.
            Hagg.append(np.sum(model[model >= probthresh] * cell_area_km2))
            if stdgrid2D is not None:
                StdH.append(np.sum(model[model >= probthresh] * cell_area_km2))
    else:
        Hagg.append(np.sum(model[model >= probthresh] * cell_area_km2))
    if len(Hagg) == 1:
        Hagg = Hagg[0]
    return Hagg, StdH


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
        Parea (float) if no or only one probthresh defined,
        otherwise, a list of floats of Parea corresponding to all
        specified probthresh values.
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
        tmpdir = tempfile.mkdtemp()
        # resample shakemap to grid2D
        temp = ShakeGrid.load(shakefile)
        junkfile = os.path.join(tmpdir, 'temp.bil')
        GDALGrid.copyFromGrid(temp.getLayer(shakethreshtype)).save(junkfile)
        shk = quickcut(junkfile, geodict, precise=True,
                       method='bilinear')
        shutil.rmtree(tmpdir)
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


def get_exposures(grid, pop_file, shakefile=None, shakethreshtype=None,
                  shakethresh=None, probthresh=None):
    """
    Get exposure-based statistics.

    Args:
        grid: Model grid.
        pop_file (str):  Path to the landscan population grid.
        shakefile (str): Optional, path to shakemap file to use for ground
            motion threshold.
        shakethreshtype(str): Optional, Type of ground motion to use for
            shakethresh, 'pga', 'pgv', or 'mmi'.
        shakethresh: Optional, Float or list of shaking thresholds in %g for
            pga, cm/s for pgv, float for mmi.
        probthresh: Optional, None or float, exclude any cells with
            probabilities less than or equal to this value

    Returns:
        dict: Dictionary with keys named exp_pop_# where # is the shakethresh
    """

    # If probthresh defined, zero out any areas less than or equal to
    # probthresh before proceeding

    if probthresh is not None:
        origdata = grid.getData()
        moddat = origdata.copy()
        moddat[moddat <= probthresh] = 0.0
        moddat[np.isnan(origdata)] = float('nan')
    else:
        moddat = grid.getData()

    mdict = grid.getGeoDict()

    # Cut out area from population file
    popcut = quickcut(pop_file, mdict, precise=False,
                      extrasamp=2., method='nearest')
    popdat = popcut.getData()
    pdict = popcut.getGeoDict()

    # Pad grid with nans to beyond extent of pdict
    pad_dict = {}
    pad_dict['padleft'] = int(
        np.abs(np.ceil((mdict.xmin - pdict.xmin)/mdict.dx)))
    pad_dict['padright'] = int(
        np.abs(np.ceil((pdict.xmax - mdict.xmax)/mdict.dx)))
    pad_dict['padbottom'] = int(
        np.abs(np.ceil((mdict.ymin - pdict.ymin)/mdict.dy)))
    pad_dict['padtop'] = int(
        np.abs(np.ceil((pdict.ymax - mdict.ymax)/mdict.dy)))
    padgrid, mdict2 = Grid2D.padGrid(
        moddat, mdict, pad_dict)  # padds with inf
    padgrid[np.isinf(padgrid)] = float('nan')  # change to pad with nan
    padgrid = Grid2D(data=padgrid, geodict=mdict2)  # Turn into grid2d object

    # Resample model grid so as to be the nearest integer multiple of popdict
    factor = np.round(pdict.dx/mdict2.dx)

    # Create geodictionary that is a factor of X higher res but otherwise
    # identical
    ndict = GeoDict.createDictFromBox(
        pdict.xmin, pdict.xmax, pdict.ymin, pdict.ymax,
        pdict.dx/factor, pdict.dy/factor)

    # Resample
    grid2 = padgrid.interpolate2(ndict, method='linear')

    # Get proportion of each cell that has values (to account properly
    # for any nans)
    prop = block_reduce(~np.isnan(grid2.getData().copy()),
                        block_size=(int(factor), int(factor)),
                        cval=float('nan'), func=np.sum)/(factor**2.)

    # Now block reduce to same geodict as popfile
    modresamp = block_reduce(grid2.getData().copy(),
                             block_size=(int(factor), int(factor)),
                             cval=float('nan'), func=np.nanmean)

    exp_pop = {}
    if shakefile is not None:
        # Resample shakefile to population grid
        # , doPadding=True, padValue=0.)
        shakemap = ShakeGrid.load(shakefile, resample=False)
        shakemap = shakemap.getLayer(shakethreshtype)
        shakemap = shakemap.interpolate2(pdict)
        shkdat = shakemap.getData()
        for shaket in shakethresh:
            threshmult = shkdat > shaket
            threshmult = threshmult.astype(float)
            exp_pop['exp_pop_%1.2fg' % (shaket/100.,)] = np.nansum(
                popdat * prop * modresamp * threshmult)

    else:
        exp_pop['exp_pop_0.00g'] = np.nansum(popdat * prop * modresamp)

    return exp_pop
