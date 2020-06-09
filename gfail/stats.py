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
                 statprobthresh=None, pop_file=None, stdtype='mean',
                 maxP=1.):
    """
    Compute summary stats of a ground failure model output.

    Args:
        grid2D: grid2D object of model output.
        stdgrid2D: grid2D object of model standard deviations (optional)
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
        stdtype (str): assumption of spatial correlation used to compute
            the stdev of the statistics, 'max', 'min' or 'mean' of max and min
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

    hagg_dict = computeHagg(grid2D, probthresh=statprobthresh, shakefile=shakefile,
                            shakethreshtype=shakethreshtype, stdtype=stdtype,
                            shakethresh=shakethresh, stdgrid2D=stdgrid2D,
                            maxP=maxP)

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
        exp_dict = get_exposures(grid2D, pop_file, shakefile=shakefile,
                                 shakethreshtype=shakethreshtype,
                                 shakethresh=shakethresh,
                                 probthresh=statprobthresh,
                                 stdgrid2D=stdgrid2D,
                                 stdtype=stdtype, maxP=maxP
                                 )
        for k, v in exp_dict.items():
            stats[k] = v

    return stats


def computeHagg(grid2D, proj='moll', probthresh=0.0, shakefile=None,
                shakethreshtype='pga', shakethresh=0.0, stdgrid2D=None,
                stdtype='mean', maxP=1.):
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
            or 'full' (default) which estimates the range and accounts for
            covariance
        maxP (float): the maximum possible probability of the model

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
        std = stdgrid.getData().copy()
        std[np.isnan(model)] = -1.

    Hagg = {}
    model[np.isnan(model)] = -1.
    if shakefile is not None:
        shkgrid = shk.project(projection=projs)
        shkdat = shkgrid.getData()
        for shaket in shakethresh:
            # use -1 to avoid nan errors and warnings, will always be thrown
            # out because default probthresh is 0.
            model[np.isnan(shkdat)] = -1.
            model[shkdat < shaket] = -1.
            mu = np.sum(model[model >= probthresh] * cell_area_km2)
            Hagg['hagg_%1.2fg' % (shaket/100.,)] = mu
            Hagg['cell_area_km2'] = cell_area_km2
            N = np.sum([model >= probthresh])
            Hagg['N_%1.2fg' % (shaket/100.,)] = N
            hlim = cell_area_km2*N*maxP
            Hagg['hlim_%1.2fg' % (shaket/100.,)] = hlim
            if stdgrid2D is not None:
                totalmin = cell_area_km2 * np.sqrt(np.nansum((std[model >= probthresh])**2.))
                totalmax = np.nansum(std[model >= probthresh] * cell_area_km2)
                if stdtype == 'max':
                    Hagg['hagg_std_%1.2fg' % (shaket/100.,)] = totalmax
                elif stdtype == 'min':
                    Hagg['hagg_std_%1.2fg' % (shaket/100.,)] = totalmin
                elif stdtype == 'mean':
                    Hagg['hagg_std_%1.2fg' % (shaket/100.,)] = (totalmax+totalmin)/2.
                else:
                    range1, sill1 = semivario(model, probthresh)
                    Hagg['hagg_std_%1.2fg' % (shaket/100.,)] = svar(
                            model, std, range1, sill1)
                var = Hagg['hagg_std_%1.2fg' % (shaket/100.,)]**2.
                # Beta distribution shape factors
                Hagg['p_hagg_%1.2fg' % (shaket/100.,)] = (mu/hlim)*((hlim*mu-mu**2)/var-1)
                Hagg['q_hagg_%1.2fg' % (shaket/100.,)] = (1-mu/hlim)*((hlim*mu-mu**2)/var-1)
    else:
        mu = np.sum(model[model >= probthresh] * cell_area_km2)
        Hagg['hagg_0.00g'] = mu
        Hagg['cell_area_km2'] = cell_area_km2
        N = np.sum([model >= probthresh])
        Hagg['N_0.00g'] = N
        hlim = cell_area_km2*N*maxP
        Hagg['hlim_0.00g'] = hlim
        if stdgrid2D is not None:
            totalmax = np.nansum(std[model >= probthresh] * cell_area_km2)
            totalmin = cell_area_km2 * np.sqrt(np.nansum((std[model >= probthresh])**2.))
            if stdtype == 'max':
                Hagg['hagg_std_0.00g'] = totalmax
            elif stdtype == 'min':
                Hagg['hagg_std_0.00g'] = totalmin
            elif stdtype == 'mean':
                Hagg['std_0.00g'] = (totalmax+totalmin)/2.
            else:
                range1, sill1 = semivario(model, probthresh)
                Hagg['std_0.00g'] = svar(model, std, range1, sill1)

            var = Hagg['hagg_std_0.00g']
            # Beta distribution shape factors
            Hagg['p_hagg_0.00g'] = (mu/hlim)*((hlim*mu-mu**2)/var-1)
            Hagg['q_hagg_0.00g'] = (1-mu/hlim)*((hlim*mu-mu**2)/var-1)

    return Hagg


def get_exposures(grid, pop_file, shakefile=None, shakethreshtype=None,
                  shakethresh=0.0, probthresh=None, stdgrid2D=None,
                  stdtype='mean', maxP=1.):
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
            the stdev of the statistics, 'max', 'min' or 'mean' of max and min

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
    if type(shakethresh) != list and type(shakethresh) != np.ndarray:
        shakethresh = [shakethresh]
    if probthresh is not None:
        origdata = grid.getData()
        moddat = origdata.copy()
        moddat[moddat <= probthresh] = 0.0
        moddat[np.isnan(origdata)] = float('nan')
        if stdgrid2D is not None:
            stddat = stdgrid2D.getData().copy()
            stddat[moddat <= probthresh] = 0.0
            stddat[np.isnan(origdata)] = 0.0
    else:
        moddat = grid.getData().copy()
        if stdgrid2D is not None:
            stddat = stdgrid2D.getData().copy()

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

    if stdgrid2D is not None:
        padstdgrid, mdict3 = Grid2D.padGrid(
            stddat, mdict, pad_dict)  # padds with inf
        padstdgrid[np.isinf(padstdgrid)] = float('nan')  # change to pad with nan
        padstdgrid = Grid2D(data=padstdgrid, geodict=mdict3)  # Turn into grid2d object

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

    if stdgrid2D is not None:
        grid2std = padstdgrid.interpolate2(ndict, method='linear')
        propstd = block_reduce(~np.isnan(grid2std.getData().copy()),
                               block_size=(int(factor), int(factor)),
                               cval=float('nan'), func=np.sum)/(factor**2.)
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
        for shaket in shakethresh:
            threshmult = shkdat > shaket
            threshmult = threshmult.astype(float)
            dat2 = popdat * prop * modresamp * threshmult
            mu = np.nansum(dat2)
            exp_pop['exp_pop_%1.2fg' % (shaket/100.,)] = mu
            elim = maxP*np.nansum(popdat * prop * threshmult)
            exp_pop['elim_%1.2fg' % (shaket/100.,)] = elim
            if stdgrid2D is not None:
                datstd2 = popdat * propstd * modresampstd * threshmult
                totalmax = np.nansum(datstd2)
                totalmin = np.sqrt(np.nansum(datstd2**2.))
                if stdtype == 'max':
                    exp_pop['exp_std_%1.2fg' % (shaket/100.,)] = totalmax
                elif stdtype == 'min':
                    exp_pop['exp_std_%1.2fg' % (shaket/100.,)] = totalmin
                elif stdtype == 'mean':
                    exp_pop['exp_std_%1.2fg' % (shaket/100.,)]=(totalmax+totalmin)/2.
                else:
                    range1, sill1 = semivario(dat2)
                    exp_pop['exp_std_%1.2fg' % (shaket/100.,)] = svar(dat2,
                            datstd2, range1, sill1)

                # Beta distribution shape factors
                var = exp_pop['exp_std_%1.2fg' % (shaket/100.,)]**2.
                exp_pop['p_exp_%1.2fg' % (shaket/100.,)] = (mu/elim)*((elim*mu-mu**2)/var-1)
                exp_pop['q_exp_%1.2fg' % (shaket/100.,)] = (1-mu/elim)*((elim*mu-mu**2)/var-1)

    else:
        dat2 = popdat * prop * modresamp
        mu = np.nansum(dat2)
        exp_pop['exp_pop_0.00g'] = mu
        elim = maxP*np.nansum(popdat * prop)
        exp_pop['elim_0.00g'] = elim
        if stdgrid2D is not None:
            datstd2 = popdat * propstd * modresampstd
            totalmax = np.nansum(datstd2)
            totalmin = np.sqrt(np.nansum(datstd2**2.))
            if stdtype == 'max':
                exp_pop['exp_std_0.00g'] = totalmax
            elif stdtype == 'min':
                exp_pop['exp_std_0.00g'] = totalmin
            elif stdtype == 'mean':
                exp_pop['exp_std_0.00g'] = (totalmax+totalmin)/2.
            else:
                range1, sill1 = semivario(dat2)
                exp_pop['exp_std_%1.2fg' % (shaket/100.,)] = svar(dat2,
                        datstd2, range1, sill1)
                exp_pop['exp_std_%1.2fg' % (shaket/100.,)] = semivario(propstd,
                        popdat, threshmult)
                    
            # Beta distribution shape factors
            var = exp_pop['exp_std_0.00g']**2.
            exp_pop['exp_std_0.00g'] = (mu/elim)*((elim*mu-mu**2)/var-1)
            exp_pop['exp_std_0.00g'] = (1-mu/elim)*((elim*mu-mu**2)/var-1)
            #exp_pop['exp_std_0.00g'] = np.nansum(popdat * propstd * modresampstd)

    return exp_pop


def semivario(model, threshold=0., maxlag=100, nsep=5, nptspb=200, ndists=100,
              nvbins=15, spacing='log', variomodel='spherical', maxrange=100.,
              makeplots=False):
    """
    Quickly estimate semivariogram with emphasis on sampling from full
    range of model values above defined threshold by sampling equal number
    of seed points from each range of values and computing a range of random
    distances pairs from those seed points that are within a maxlag range of
    pixels. If nsep is set to 1, then nptspb seed points will be selected
    randomly

    Args:
        model: 2D array of raster to estimate semivariogram for
        threshold: 
        maxlag: in pixels
        nsep: number of bins to sample from equally, spaced from min to max
            value of model
        nptspb: number of seed points to sample from in each separation bin
        ndists: number of points to sample at random distances from each seed point
        nvbins: number of semivariogram bins
        spacing: spacing of nsep, 'log' for logarithmic spacing, or 'linear'
            for linear spacing
        model: model to fit to semivariogram

    Returns:
        range, sill

    """
    #TODO deal with cases where there are no values above threshold
    if np.nanmax(model) < threshold:
        raise Exception('No values above threshold in model')

    # prepare data
    nrows, ncols = np.shape(model)
    rows = np.matlib.repmat(np.arange(nrows), ncols, 1).T
    cols = np.matlib.repmat(np.arange(ncols), nrows, 1)
    values = model.flatten()
    rowsf = rows.flatten()
    colsf = cols.flatten()
    # Divide range of values into X equally spaced bins and ensure sufficient samples are collected from each
    if spacing=='log':
        bins = np.logspace(np.log10(threshold), np.log10(np.nanmax(values)),
                           num=nsep+1, endpoint=True)
    else:
        bins = np.linspace(threshold, np.nanmax(values), num=nsep+1,
                           endpoint=True)

    # Select nptspb points in each bin
    seedpts = np.array([], dtype=int)
    for i in range(nsep):
        indx = np.where((values>bins[i]) & (values<bins[i+1]))[0]
        if len(indx) == 0:
            continue
        elif len(indx) < nptspb:
            # keep all points and pick some more than once to equal nptspb
            picks = np.random.choice(len(indx), size=nptspb, replace=True)
        elif len(indx) >= nptspb:
            # randomly select without replacement
            picks = np.random.choice(len(indx), size=nptspb, replace=False)
        seedpts = np.hstack((seedpts, indx[picks]))
    
    # Get lags and differences for seed point vs. ndists other points around it
    #TODO vectorize this to remove loop
    lags = np.array([])
    diffs = np.array([])
    vals = np.array([])
    for seed in seedpts:
        row1 = rowsf[seed] 
        col1 = colsf[seed]
        addr = np.random.randint(np.max((-maxlag, -row1)), np.min((maxlag, nrows-row1)), size=ndists)
        addc = np.random.randint(np.max((-maxlag, -col1)), np.min((maxlag, ncols-col1)), size=ndists)
        indr = row1 + addr
        indc = col1 + addc
        newvalues = model[indr, indc]
        ptval = model[row1, col1]
        dists = np.sqrt(addr**2 + addc**2)  # distance in pixels
        difvals = np.abs(newvalues - ptval)
        diffs = np.hstack((diffs, difvals))
        lags = np.hstack((lags, dists))
        vals = np.hstack((vals, newvalues, ptval))

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

    model1 = eval(variomodel)
    
    # Fit model using weighting by 1/N to weigh bins with more samples more highly
    popt, pcov = curve_fit(model1, binmid, semiv, sigma=1./N,
                           absolute_sigma=False, bounds=(0, [maxrange, 1.]))
    if makeplots:
        plt.figure()    
        plt.plot(binmid, semiv, 'ob')
        plt.xlabel('Lags (pixels)')
        plt.ylabel('Semivariance')
        plt.plot(binmid, spherical(binmid, *popt), '-b')

    range2, sill2 = popt
    
    return range2, sill2


def spherical(lag, range1, sill, nugget=0):
    """
    https://github.com/mmaelicke/scikit-gstat/blob/master/skgstat/models.py#L23
    
    nugget = value of independent variable at distance of zero
    """
    range1 = range1 / 1.

    out = nugget + sill * ((1.5 * (lag / range1)) - (0.5 * ((lag / range1) ** 3.0)))
    out[lag > range1] = nugget + sill
    return out


def exponential(lag, range1, sill, nugget=0):
    """
    https://github.com/mmaelicke/scikit-gstat/blob/master/skgstat/models.py#L87
    
    """
    a = range1 / 3.
    return nugget + sill * (1. - np.exp(-(lag / a)))


def gaussian(lag, range1, sill, nugget=0):
    a = range1 / 2.
    return nugget + sill * (1. - np.exp(- (lag ** 2 / a ** 2)))


def svar(model, stds, range1, sill1, variomodel='spherical'):
    """
    Estimate variance of aggregate statistic using correlation from 
    semivariogram and std values for each pair of cells that are within range
    of each other
    
    
    """
    # prepare data
    nrows, ncols = np.shape(model)
    rows = np.matlib.repmat(np.arange(nrows), ncols, 1).T
    cols = np.matlib.repmat(np.arange(ncols), nrows, 1)
    values = model.flatten()
    rowsf = rows.flatten()
    colsf = cols.flatten()
    
    model1 = eval(variomodel)

    var2 = 0 # Simulated using correlation estimated from semivariogram and
    #std for each cell, only include if two cells are within range of each other

    #TODO vectorize/make more efficient
    for i in range(len(values)):
        for j in range(len(values)):
            if i == j:
                var2 += stds[i]**2
            elif i>j:
                dist = np.sqrt((rowsf[i]-rowsf[j])**2 + (colsf[i]-colsf[j]))
                sigij = sill1-model1(dist, range1, sill1)
                corr = sigij/sill1
                if dist < range1:
                    var2 += 2 * corr*stds[i]*stds[j]
    return var2
