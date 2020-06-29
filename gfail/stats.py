#!/usr/bin/env python
# -*- coding: utf-8 -*-

# stdlib imports
import numpy as np
import collections
import shutil
import tempfile
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import convolve
from numpy import matlib

# local imports
from mapio.shake import ShakeGrid
from mapio.gdal import GDALGrid
from mapio.geodict import GeoDict
from gfail.spatial import quickcut
from mapio.grid2d import Grid2D
from skimage.measure import block_reduce

from configobj import ConfigObj

# Make fonts readable and recognizable by illustrator
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = ['Arial',
                                   'Bitstream Vera Serif',
                                   'sans-serif']


def computeStats(grid2D, stdgrid2D=None, shakefile=None,
                 shakethreshtype='pga', shakethresh=0.0,
                 probthresh=None, pop_file=None, stdtype='full',
                 maxP=1., proj='moll'):
    """
    Compute summary stats of a ground failure model output.

    Args:
        grid2D: grid2D object of model output.
        stdgrid2D: grid2D object of model standard deviations (optional)
        shakefile: Optional, path to shakemap file to use for ground motion
            threshold.
        shakethreshtype: Optional, Type of ground motion to use for
            shakethresh, 'pga', 'pgv', or 'mmi'.
        shakethresh: Optional, Float in %g for
            pga, cm/s for pgv, float for mmi. Used for Hagg and Exposure
            computation
        probthresh: Optional, None or float, exclude any cells with
            probabilities less than or equal to this value
        pop_file (str): File path to population file to use to compute exposure
            stats
        stdtype (str): assumption of spatial correlation used to compute
            the stdev of the statistics, 'max', 'min' or 'mean' of max and min,
            or full (default) estimates std considering covariance
        maxP = maximum possible value of P (1 default, but coverge models
            have smaller values, 0.487 and 0.256 for LQ and LS)

    Returns:
        dict: Dictionary with all or some of the following keys
            (depending on input options):
            - Max
            - Median
            - Std
            - hagg_# where # is the shaking threshold input
            - hagg_std_# standard deviation of hagg_#
            - hlim_# maximum possible value of hagg (for given Pmax value)
            - p_hagg_# beta function shapefile p for hagg
            - q_hagg_# beta function shapefile q for hagg
            - exp_pop_# where # is the shaking threshold (if pop_file specified)
            - exp_std_# (optional) standard deviation of exp_pop_#
            - elim_# maximum possible value of E (for given Pmax value)
            - p_exp_# beta function shapefile p for exp_pop
            - q_exp_# beta function shapefile q for exp_pop

    """
    stats = collections.OrderedDict()
    grid = grid2D.getData().copy()
    if probthresh is not None:
        grid = grid[grid > probthresh]
    else:
        probthresh = 0.0

    if len(grid) == 0:
        print('no probability values above probthresh')
        stats['Max'] = 0.  # float('nan')
        stats['Median'] = 0.  # float('nan')
        stats['Std'] = 0.  # float('nan')
    else:
        stats['Max'] = float(np.nanmax(grid))
        stats['Median'] = float(np.nanmedian(grid))
        stats['Std'] = float(np.nanstd(grid))

    output = computeHagg(grid2D, probthresh=probthresh, shakefile=shakefile,
                         shakethreshtype=shakethreshtype, stdtype=stdtype,
                         shakethresh=shakethresh, stdgrid2D=stdgrid2D,
                         maxP=maxP, sill1=None, range1=None, proj=proj)

    hagg_dict, sill1, range1km = output

    for k, v in hagg_dict.items():
        stats[k] = v

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
        # convert range1km to range1 in pixels of size of population grid
        if stdtype == 'full' and range1km is not None:
            range1 = range1km/(GDALGrid.getFileGeoDict(pop_file)[0].dy * 111.)
        else:
            range1 = None
        exp_dict = computePexp(grid2D, pop_file, shakefile=shakefile,
                               shakethreshtype=shakethreshtype,
                               shakethresh=shakethresh,
                               probthresh=probthresh,
                               stdgrid2D=stdgrid2D,
                               stdtype=stdtype, maxP=maxP,
                               sill1=sill1, range1=range1, proj=proj)
        for k, v in exp_dict.items():
            stats[k] = v

    return stats


def computeHagg(grid2D, proj='moll', probthresh=0.0, shakefile=None,
                shakethreshtype='pga', shakethresh=0.0, stdgrid2D=None,
                stdtype='full', maxP=1., sill1=None, range1=None):
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
        stdtype (str): assumption of spatial correlation used to compute
            the stdev of the statistics, 'max', 'min', 'mean' of max and min,
            or 'full' (default) which estimates the range of correlation and
            accounts for covariance
        maxP (float): the maximum possible probability of the model
        sill1 (float):
        range1 (float):

    Returns:
        dict: Dictionary with keys:
            hagg_#g where # is the shakethresh
            std_# if stdgrid2D is supplied (stdev of exp_pop)
            hlim_#, the maximum exposure value possible with the
            applied thresholds and given maxP value
            N_# the number of cells exceeding that value (in projected coords)
            cell_area_km2 grid cell area
            p_hagg_# beta distribution shape factor p (sometimes called alpha)
            q_hagg_# beta distribution shape factor q (sometimes called beta)
    """
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
        temp = ShakeGrid.load(shakefile)
        shk = temp.getLayer(shakethreshtype)
        shk = shk.interpolate2(geodict)
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
        std = stdgrid.getData().copy()

    Hagg = {}

    if shakefile is not None:
        shkgrid = shk.project(projection=projs)
        shkdat = shkgrid.getData()
        # use -1 to avoid nan errors and warnings, will always be thrown
        # out because default probthresh is 0.
        model[shkdat < shakethresh] = float('nan')
    else:
        shakethresh == 0.
        shkdat = None

    mu = np.nansum(model[model >= probthresh] * cell_area_km2)
    Hagg['hagg_%1.2fg' % (shakethresh/100.,)] = mu
    Hagg['cell_area_km2'] = cell_area_km2
    N = np.nansum([model >= probthresh])
    Hagg['N_%1.2fg' % (shakethresh/100.,)] = N
    hlim = cell_area_km2*N*maxP
    Hagg['hlim_%1.2fg' % (shakethresh/100.,)] = hlim

    sill1 = None
    range1km = None
    if stdgrid2D is not None:
        if np.nanmax(std) > 0. and np.nanmax(model) >= probthresh:
            totalmin = cell_area_km2 * np.sqrt(np.nansum((std[model >= probthresh])**2.))
            totalmax = np.nansum(std[model >= probthresh] * cell_area_km2)
            if stdtype == 'max':
                Hagg['hagg_std_%1.2fg' % (shakethresh/100.,)] = totalmax
            elif stdtype == 'min':
                Hagg['hagg_std_%1.2fg' % (shakethresh/100.,)] = totalmin
            elif stdtype == 'mean':
                Hagg['hagg_std_%1.2fg' % (shakethresh/100.,)] = (totalmax+totalmin)/2.
            else:
                if sill1 is None or range1 is None:
                    modelfresh = grid.getData().copy()
                    range1, sill1 = semivario(modelfresh, probthresh,
                                              shakethresh=shakethresh,
                                              shakegrid=shkdat)
                    
                range1km = range1 * geodictRS.dx
                stdz = std.copy()
                stdz[model < probthresh] = 0.
                Hagg['hagg_std_%1.2fg' % (shakethresh/100.,)] = np.sqrt(svar(
                        stdz, range1, sill1, scale=cell_area_km2))
            var = Hagg['hagg_std_%1.2fg' % (shakethresh/100.,)]**2.
            # Beta distribution shape factors
            Hagg['p_hagg_%1.2fg' % (shakethresh/100.,)] = (mu/hlim)*((hlim*mu-mu**2)/var-1)
            Hagg['q_hagg_%1.2fg' % (shakethresh/100.,)] = (1-mu/hlim)*((hlim*mu-mu**2)/var-1)
        else:
            print('No std values above threshold, skipping uncertainty '
                  'and filling with zeros')
            Hagg['hagg_std_%1.2fg' % (shakethresh/100.,)] = 0.
            Hagg['p_hagg_%1.2fg' % (shakethresh/100.,)] = 0.
            Hagg['q_hagg_%1.2fg' % (shakethresh/100.,)] = 0.

    return Hagg, sill1, range1km


def computePexp(grid, pop_file, shakefile=None, shakethreshtype=None,
                shakethresh=0.0, probthresh=None, stdgrid2D=None,
                stdtype='full', maxP=1., sill1=None, range1=None,
                proj='moll'):
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
        stdgrid2D: grid2D object of model standard deviations (optional)
        stdtype (str): assumption of spatial correlation used to compute
            the stdev of the statistics, 'max', 'min' or 'mean' of max and min,
            'full' estimates the std with covariance

    Returns:
        dict: Dictionary with keys named exp_pop_# where # is the shakethresh
            and exp_std_# if stdgrid2D is supplied (stdev of exp_pop)
            and elim_#, the maximum exposure value possible with the
            applied thresholds and given maxP value
            p_exp_# beta distribution shape factor p (sometimes called alpha)
            q_exp_# beta distribution shape factor q (sometimes called beta)
    """

    # If probthresh defined, zero out any areas less than or equal to
    # probthresh before proceeding
    moddat = grid.getData().copy()
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

    if stdgrid2D is not None:
        stddat = stdgrid2D.getData().copy()
        if stdtype=='full':
            # Get range and sill from original data if not provided
            if sill1 is None or range1 is None:
                # Compute variogram of data in distance preserving space
                bounds = grid.getBounds()
                lat0 = np.mean((bounds[2], bounds[3]))
                lon0 = np.mean((bounds[0], bounds[1]))
                projs = ('+proj=%s +lat_0=%f +lon_0=%f +x_0=0 +y_0=0 +ellps=WGS84 '
                         '+units=km +no_defs' % (proj, lat0, lon0))
                gridP = grid.project(projection=projs)
                geodictRS = gridP.getGeoDict()
                if shakefile is not None:
                    # resample shakemap to grid2D
                    temp = ShakeGrid.load(shakefile)
                    shk = temp.getLayer(shakethreshtype)
                    shk = shk.interpolate2(grid.getGeoDict())
                    shkgrid = shk.project(projection=projs)
                    shkdat = shkgrid.getData()
                else:
                    shkdat = None
                    shakethresh = 0.
                model = gridP.getData().copy()
                if np.nanmax(stddat) > 0. and np.nanmax(model) >= probthresh:
                    range1, sill1 = semivario(model, probthresh,
                                              shakethresh=shakethresh,
                                              shakegrid=shkdat)
                range1km = range1 * geodictRS.dx
                # Convert to approx # of gridpts relative to pop_file
                range1 = range1km/(pdict.dy * 111.)

        padstdgrid, mdict3 = Grid2D.padGrid(stddat, mdict, pad_dict)  # padds with inf
        padstdgrid[np.isinf(padstdgrid)] = float('nan')  # change to pad with nan
        padstdgrid = Grid2D(data=padstdgrid, geodict=mdict3)  # Turn into grid2d object

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

    if stdgrid2D is not None:
        grid2std = padstdgrid.interpolate2(ndict, method='linear')
        modresampstd = block_reduce(grid2std.getData().copy(),
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
        threshmult = (shkdat >= shakethresh) & (modresamp >= probthresh)
    else:
        shakethresh = 0.0
        shkdat = None
        threshmult = modresamp >= probthresh
    threshmult = threshmult.astype(float)

    dat2 = popdat * prop * modresamp * threshmult

    mu = np.nansum(dat2)
    exp_pop['exp_pop_%1.2fg' % (shakethresh/100.,)] = mu
    elim = maxP*np.nansum(popdat * prop * threshmult)
    exp_pop['elim_%1.2fg' % (shakethresh/100.,)] = elim
    if stdgrid2D is not None:
        if np.nanmax(modresampstd) > 0. and np.sum(threshmult) > 0:
            if stdtype=='full':
                exp_pop['exp_std_%1.2fg' % (shakethresh/100.,)] = \
                    np.sqrt(svar(np.copy(modresampstd), range1, sill1, scale=popdat))
            else:                
                datstd2 = modresampstd * threshmult
                totalmax = np.nansum(popdat * prop * datstd2)
                totalmin = np.sqrt(np.nansum(popdat * prop * datstd2**2.))
                if stdtype == 'max':
                    exp_pop['exp_std_%1.2fg' % (shakethresh/100.,)] = totalmax
                elif stdtype == 'min':
                    exp_pop['exp_std_%1.2fg' % (shakethresh/100.,)] = totalmin
                elif stdtype == 'mean':
                    exp_pop['exp_std_%1.2fg' % (shakethresh/100.,)] = (totalmax+totalmin)/2.

            # Beta distribution shape factors
            var = exp_pop['exp_std_%1.2fg' % (shakethresh/100.,)]**2.
            exp_pop['p_exp_%1.2fg' % (shakethresh/100.,)] = (mu/elim)*((elim*mu-mu**2)/var-1)
            exp_pop['q_exp_%1.2fg' % (shakethresh/100.,)] = (1-mu/elim)*((elim*mu-mu**2)/var-1)
        else:
            print('no std values above zero, filling with zeros')
            exp_pop['exp_std_%1.2fg' % (shakethresh/100.,)] = 0.
            exp_pop['p_exp_%1.2fg' % (shakethresh/100.,)] = 0.
            exp_pop['q_exp_%1.2fg' % (shakethresh/100.,)] = 0.

    return exp_pop


def semivario(model, threshold=0., maxlag=100, npts=1000, ndists=200,
              nvbins=20, makeplots=False, shakegrid=None, shakethresh=0.):
    """
    Quickly estimate semivariogram with by selecting seed points and then
    computing semivariogram between each of those points and ndists random
    locations around it that are within maxlag of the seed point. Will result 
    in npts x ndists total distance pairs. Uses spherical model.

    Args:
        model: array of raster to estimate semivariogram for
        threshold: 
        maxlag: in pixels
        npts: number of seed points to sample from
        ndists: number of points to sample at random distances from each seed point
        nvbins: number of semivariogram bins

    Returns:
        range, sill

    """
    if threshold is None:
        threshold = 0.

    if shakegrid is None or shakethresh == 0.:
        shakegrid == np.zeros(np.shape(model))
    
    if np.shape(shakegrid) != np.shape(model):
        raise Exception('Shakegrid is not the same shape as the model')

    # prepare data
    nrows, ncols = np.shape(model)
    rows = matlib.repmat(np.arange(nrows), ncols, 1).T
    cols = matlib.repmat(np.arange(ncols), nrows, 1)
    values = model.flatten()
    shkvals = shakegrid.flatten()
    rowsf = rows.flatten()
    colsf = cols.flatten()

    # Select npts seed points
    indx = np.where((values >= threshold) & (shkvals >= shakethresh))[0]
    if len(indx) == 0:
        print('No values above thresholds in model. Returning empty results')
        return None, None
    
    np.random.seed(47)  # Always use same seed so results are repeatable
    # just select all points if there aren't npts above the threshold
    picks = np.random.choice(len(indx), size=np.min((npts, len(indx))),
                             replace=False)
    seedpts = indx[picks]
    
    # Get lags and differences for seed point vs. ndists other points around it
    #TODO vectorize this to remove loop
    lags = np.array([])
    diffs = np.array([])
    #vals = np.array([])
    seednums = range(len(seedpts))
    seednums2 = reversed(range(len(seedpts)))
    for seed, seednum1, seednum2 in zip(seedpts, seednums, seednums2):
        row1 = rowsf[seed] 
        col1 = colsf[seed]
        np.random.seed(seednum1)
        addr = np.random.randint(np.max((-maxlag, -row1)), np.min((maxlag, nrows-row1)), size=ndists)
        np.random.seed(seednum2)
        addc = np.random.randint(np.max((-maxlag, -col1)), np.min((maxlag, ncols-col1)), size=ndists)
        indr = row1 + addr
        indc = col1 + addc
        newvalues = model[indr, indc]
        ptval = model[row1, col1]
        dists = np.sqrt(addr**2 + addc**2)  # distance in pixels
        # Remove any invalid samples
        dists = dists[np.isfinite(newvalues)]
        newvalues = newvalues[np.isfinite(newvalues)]
        difvals = np.abs(newvalues - ptval)
        diffs = np.hstack((diffs, difvals))
        lags = np.hstack((lags, dists))
        #vals = np.hstack((vals, newvalues, ptval))

    diffs2 = diffs**2

    # % Make variogram out of these samples
    binedges = np.linspace(0, maxlag, num=nvbins+1, endpoint=True)
    binmid = (binedges[:-1] + binedges[1:])/2
    
    subs = np.zeros(nvbins)
    N = np.zeros(nvbins)
    for b in range(nvbins):
        inrange = diffs2[(lags > binedges[b]) & (lags <= binedges[b+1])]
        subs[b] = np.nansum(inrange)
        N[b] = len(inrange)
    semiv = 1./(2*N)*subs
    # Fit model using weighting by 1/N to weigh bins with more samples more highly
    popt, pcov = curve_fit(spherical, binmid[np.isfinite(semiv)],
                           semiv[np.isfinite(semiv)],
                           sigma=1./N[np.isfinite(semiv)],
                           absolute_sigma=False, bounds=(0, [maxlag, 1.]))
    if makeplots:
        plt.figure()    
        plt.plot(binmid, semiv, 'ob')
        plt.xlabel('Lags (pixels)')
        plt.ylabel('Semivariance')
        plt.plot(binmid, spherical(binmid, *popt), '-b')

    range2, sill2 = popt
    
    return range2, sill2


def spherical(lag, range1, sill):#, nugget=0):
    """
    https://github.com/mmaelicke/scikit-gstat/blob/master/skgstat/models.py#L23
    
    nugget = value of independent variable at distance of zero
    """
    nugget = 0.
    range1 = range1 / 1.

    out = nugget + sill * ((1.5 * (lag / range1)) - (0.5 * ((lag / range1) ** 3.0)))
    if isinstance(out, float):
        if lag > range1:
            out = nugget + sill
    else:
        out[lag > range1] = nugget + sill
    return out


def svar(stds, range1, sill1, scale=1.):
    """
    Estimate variance of aggregate statistic using correlation from 
    semivariogram and std values for each pair of cells that are within range
    of each other, add up quickly by creating kernal of the correlations and
    convolving with the image, then multiply by 2*std to equal sum of 2*std1*std2*corr
    over each valid cell
    
    Args:
        stds:
        range1:
        sill1:
        scale: float or array same size as std, factor to multiply by
            (area or population) and include in conv
    
    Returns:
        variance of aggregate statistic
    
    """
    range5 = int(np.ceil(range1))
    # Prepare kernal that is size of range of spherical equation
    nrows = 2*range5 + 1
    ncols = 2*range5 + 1
    # get distance in row and col from center point
    rows = matlib.repmat(np.arange(nrows), ncols, 1).T - (range5)
    cols = matlib.repmat(np.arange(ncols), nrows, 1) - (range5)
    dists = np.sqrt(rows**2 + cols**2)
    # Convert from semivariance to correlation and build kernal
    kernal = (sill1-spherical(dists, range1, sill1))/sill1
    
    # convolve with stds, equivalent to sum of corr * std at each pixel for within range1
    #out = convolve2d(stds, kernal, mode='same')
    # Replace all nans with zeros so can use fft convolve
    stdzeros = stds.copy()
    stdzeros[np.isnan(stds)] = 0.
    if not isinstance(scale, float):
        scale[np.isnan(scale)] = 0.
    stdzeros *= scale
    out = convolve(stdzeros, kernal, mode='same')
    # multiply by stds and scale again
    full1 = scale * out * stds
    # add up
    var2 = np.nansum(full1)
    return var2
