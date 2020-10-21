#!/usr/bin/env python
# -*- coding: utf-8 -*-

# stdlib imports
import numpy as np
import collections
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import convolve
from numpy import matlib

# local imports
from mapio.shake import ShakeGrid
from mapio.gdal import GDALGrid
from gfail.spatial import quickcut


from configobj import ConfigObj

# Turn off warnings that will pop up regarding nan's in greater than operations
np.warnings.filterwarnings('ignore')


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
        maxP (float): maximum possible value of P (1 default, but coverage
            models  have smaller values, 0.487 and 0.256 for LQ and LS)
        proj (str): projection string to use when computing stats

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

    hagg_dict = output

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
        exp_dict = computePexp(grid2D, pop_file, shakefile=shakefile,
                               shakethreshtype=shakethreshtype,
                               shakethresh=shakethresh,
                               probthresh=probthresh,
                               stdgrid2D=stdgrid2D,
                               stdtype=stdtype, maxP=maxP,
                               sill1=None, range1=None)
        for k, v in exp_dict.items():
            stats[k] = v

    return stats


def computeHagg(grid2D, proj='moll', probthresh=0., shakefile=None,
                shakethreshtype='pga', shakethresh=0., stdgrid2D=None,
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
            accounts for covariance. Will return 'mean' if
            ridge and sill cannot be estimated.
        maxP (float): the maximum possible probability of the model
        sill1 (float, None): If known, the sill of the variogram of grid2D,
            will be estimated if None and stdtype='full'
        range1 (float, None): If known, the range of the variogram of grid2D,
            will be estimated if None and stdtype='full'

    Returns:
        dict: Dictionary with keys:
            hagg_#g where # is the shakethresh
            std_# if stdgrid2D is supplied (stdev of exp_pop)
            hlim_#, the maximum exposure value possible with the
            applied thresholds and given maxP value
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

    Hagg = {}

    if shakefile is not None:
        shkgrid = shk.project(projection=projs)
        shkdat = shkgrid.getData()
        model[shkdat < shakethresh] = float('nan')
    else:
        shakethresh = 0.
        shkdat = None

    mu = np.nansum(model[model >= probthresh] * cell_area_km2)
    Hagg['hagg_%1.2fg' % (shakethresh/100.,)] = mu
    Hagg['cell_area_km2'] = cell_area_km2
    N = np.nansum([model >= probthresh])
    #Hagg['N_%1.2fg' % (shakethresh/100.,)] = N
    hlim = cell_area_km2*N*maxP
    Hagg['hlim_%1.2fg' % (shakethresh/100.,)] = hlim

    if stdgrid2D is not None:
        stdgrid = GDALGrid.copyFromGrid(stdgrid2D)  # Make a copy
        stdgrid = stdgrid.project(projection=projs, method='bilinear')
        std = stdgrid.getData().copy()
        if np.nanmax(std) > 0. and np.nanmax(model) >= probthresh:
            totalmin = cell_area_km2 * np.sqrt(np.nansum((std[model >= probthresh])**2.))
            totalmax = np.nansum(std[model >= probthresh] * cell_area_km2)
            if stdtype == 'full':
                if sill1 is None or range1 is None:
                    range1, sill1 = semivario(grid.getData().copy(),
                                              probthresh,
                                              shakethresh=shakethresh,
                                              shakegrid=shkdat)
                if range1 is None:
                    # Use mean
                    Hagg['hagg_std_%1.2fg' % (shakethresh/100.,)] = (totalmax+totalmin)/2.
                else:
                    # Zero out std at cells where the model probability was
                    # below the threshold because we aren't including those
                    # cells in Hagg
                    stdz = std.copy()
                    stdz[model < probthresh] = 0.
                    svar1 = svar(stdz, range1, sill1, scale=cell_area_km2)
                    Hagg['hagg_std_%1.2fg' % (shakethresh/100.,)] = np.sqrt(svar1)
                    #Hagg['hagg_range_%1.2fg' % (shakethresh/100.,)] = range1
                    #Hagg['hagg_sill_%1.2fg' % (shakethresh/100.,)] = sill1 
            elif stdtype == 'max':
                Hagg['hagg_std_%1.2fg' % (shakethresh/100.,)] = totalmax
            elif stdtype == 'min':
                Hagg['hagg_std_%1.2fg' % (shakethresh/100.,)] = totalmin
            else:
                Hagg['hagg_std_%1.2fg' % (shakethresh/100.,)] = (totalmax+totalmin)/2.

            var = Hagg['hagg_std_%1.2fg' % (shakethresh/100.,)]**2.
            # Beta distribution shape factors
            Hagg['p_hagg_%1.2fg' % (shakethresh/100.,)] = (mu/hlim)*((hlim*mu-mu**2)/var-1)
            Hagg['q_hagg_%1.2fg' % (shakethresh/100.,)] = (1-mu/hlim)*((hlim*mu-mu**2)/var-1)
        else:
            print('No model values above threshold, skipping uncertainty '
                  'and filling with zeros')
            Hagg['hagg_std_%1.2fg' % (shakethresh/100.,)] = 0.
            Hagg['p_hagg_%1.2fg' % (shakethresh/100.,)] = 0.
            Hagg['q_hagg_%1.2fg' % (shakethresh/100.,)] = 0.
    else:
        print('No uncertainty provided, filling with zeros')
        Hagg['hagg_std_%1.2fg' % (shakethresh/100.,)] = 0.
        Hagg['p_hagg_%1.2fg' % (shakethresh/100.,)] = 0.
        Hagg['q_hagg_%1.2fg' % (shakethresh/100.,)] = 0.

    return Hagg


def computePexp(grid, pop_file, shakefile=None, shakethreshtype='pga',
                shakethresh=0., probthresh=0., stdgrid2D=None,
                stdtype='full', maxP=1., sill1=None, range1=None):
    """
    Get exposure-based statistics.

    Args:
        grid: Model grid.
        pop_file (str):  Path to the landscan population grid.
        shakefile (str): Optional, path to shakemap file to use for ground
            motion threshold.
        shakethreshtype(str): Optional, Type of ground motion to use for
            shakethresh, 'pga', 'pgv', or 'mmi'.
        shakethresh: Float or list of shaking thresholds in %g for
            pga, cm/s for pgv, float for mmi.
        probthresh: Float, exclude any cells with
            probabilities less than or equal to this value
        stdgrid2D: grid2D object of model standard deviations (optional)
        stdtype (str): assumption of spatial correlation used to compute
            the stdev of the statistics, 'max', 'min', 'mean' of max and min,
            or 'full' (default) which estimates the range of correlation and
            accounts for covariance. Will return 'mean' if
            ridge and sill cannot be estimated.
        maxP (float): the maximum possible probability of the model
        sill1 (float, None): If known, the sill of the variogram of grid2D,
            will be estimated if None and stdtype='full'
        range1 (float, None): If known, the range of the variogram of grid2D,
            will be estimated if None and stdtype='full'

    Returns:
        dict: Dictionary with keys named exp_pop_# where # is the shakethresh
            and exp_std_# if stdgrid2D is supplied (stdev of exp_pop)
            and elim_#, the maximum exposure value possible with the
            applied thresholds and given maxP value
            p_exp_# beta distribution shape factor p (sometimes called alpha)
            q_exp_# beta distribution shape factor q (sometimes called beta)
    """

    model = grid.getData().copy()
    mdict = grid.getGeoDict()

    # Figure out difference in resolution of popfile to shakefile
    ptemp, J = GDALGrid.getFileGeoDict(pop_file)
    factor = ptemp.dx/mdict.dx

    # Cut out area from population file
    popcut1 = quickcut(pop_file, mdict, precise=False,
                       extrasamp=2, method='nearest')
    #tot1 = np.sum(popcut1.getData())
    # Adjust for factor to prepare for upsampling to avoid creating new people
    popcut1.setData(popcut1.getData()/factor**2)

    # Upsample to mdict
    popcut = popcut1.interpolate2(mdict, method='nearest')
    popdat = popcut.getData()
    exp_pop = {}

    if shakefile is not None:
        if shakethresh < 0.:
            raise Exception('shaking threshold must be equal or greater '
                            'than zero')
        # resample shakemap to grid2D
        temp = ShakeGrid.load(shakefile)
        shk = temp.getLayer(shakethreshtype)
        shk = shk.interpolate2(mdict)
        if shk.getGeoDict() != mdict:
            raise Exception('shakemap was not resampled to exactly the same '
                            'geodict as the model')
        shkdat = shk.getData()
        model[shkdat < shakethresh] = float('nan')
    else:
        shakethresh = 0.
        shkdat = None

    mu = np.nansum(model[model >= probthresh] * popdat[model >= probthresh])
    exp_pop['exp_pop_%1.2fg' % (shakethresh/100.,)] = mu
    #N = np.nansum([model >= probthresh])
    #exp_pop['N_%1.2fg' % (shakethresh/100.,)] = N
    elim = np.nansum(popdat[model >= probthresh])*maxP
    exp_pop['elim_%1.2fg' % (shakethresh/100.,)] = elim

    if stdgrid2D is not None:
        std = stdgrid2D.getData().copy()
        if np.nanmax(std) > 0. and np.nanmax(model) >= probthresh:
            totalmin = np.sqrt(np.nansum((popdat[model >= probthresh]*std[model >= probthresh])**2.))
            totalmax = np.nansum(std[model >= probthresh] * popdat[model >= probthresh])
            if stdtype == 'full':
                if sill1 is None or range1 is None:
                    modelfresh = grid.getData().copy()
                    range1, sill1 = semivario(modelfresh, probthresh,
                                              shakethresh=shakethresh,
                                              shakegrid=shkdat)
                if range1 is None:
                    # Use mean
                    exp_pop['exp_std_%1.2fg' % (shakethresh/100.,)] = \
                        (totalmax+totalmin)/2.
                else:
                    # Zero out std at cells where the model probability was
                    # below the threshold because we aren't including those
                    # cells in Hagg
                    stdz = std.copy()
                    stdz[model < probthresh] = 0.
                    svar1 = svar(stdz, range1, sill1, scale=popdat)
                    exp_pop['exp_std_%1.2fg' % (shakethresh/100.,)] = \
                        np.sqrt(svar1)
                    #exp_pop['exp_range_%1.2fg' % (shakethresh/100.,)] = range1
                    #exp_pop['exp_sill_%1.2fg' % (shakethresh/100.,)] = sill1

            elif stdtype == 'max':
                exp_pop['exp_std_%1.2fg' % (shakethresh/100.,)] = totalmax
            elif stdtype == 'min':
                exp_pop['exp_std_%1.2fg' % (shakethresh/100.,)] = totalmin
            else:
                exp_pop['exp_std_%1.2fg' % (shakethresh/100.,)] = \
                    (totalmax+totalmin)/2.
            # Beta distribution shape factors
            var = exp_pop['exp_std_%1.2fg' % (shakethresh/100.,)]**2.
            exp_pop['p_exp_%1.2fg' % (shakethresh/100.,)] = \
                (mu/elim)*((elim*mu-mu**2)/var-1)
            exp_pop['q_exp_%1.2fg' % (shakethresh/100.,)] = \
                (1-mu/elim)*((elim*mu-mu**2)/var-1)
        else:
            print('no std values above zero, filling with zeros')
            exp_pop['exp_std_%1.2fg' % (shakethresh/100.,)] = 0.
            exp_pop['p_exp_%1.2fg' % (shakethresh/100.,)] = 0.
            exp_pop['q_exp_%1.2fg' % (shakethresh/100.,)] = 0.
    else:
        exp_pop['exp_std_%1.2fg' % (shakethresh/100.,)] = 0.
        exp_pop['p_exp_%1.2fg' % (shakethresh/100.,)] = 0.
        exp_pop['q_exp_%1.2fg' % (shakethresh/100.,)] = 0.

    return exp_pop


def semivario(model, threshold=0., maxlag=100, npts=1000, ndists=200,
              nvbins=20, makeplots=False, shakegrid=None, shakethresh=0.,
              minpts=50):
    """
    Quickly estimate semivariogram with by selecting seed points and then
    computing semivariogram between each of those points and ndists random
    locations around it that are within maxlag of the seed point. Will result 
    in npts x ndists total distance pairs. Uses spherical model.

    Args:
        model (array): array of raster to estimate semivariogram for
        threshold: 
        maxlag (int): in pixels
        npts (int): number of seed points to sample from
        ndists (int): number of points to sample at random distances from each
            seed point
        nvbins (int): number of semivariogram bins
        minpts (float): minimum number of samples above threshold required
        makeplots (bool): create semivariogram plots
        shakegrid (array): array of shaking that is the same size as model
        shakethresh (float): Shaking threshold for seed point selection

    Returns:
        range, sill

    """
    model = model.copy()
    if threshold is None:
        threshold = 0.

    if shakegrid is None or shakethresh == 0.:
        shakegrid = np.zeros(np.shape(model))
    
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
    if len(indx) < minpts:
        print('Not enough values above thresholds in model. '
              'Returning empty results')
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
        addr = np.random.randint(np.max((-maxlag, -row1)),
                                 np.min((maxlag, nrows-row1)), size=ndists)
        np.random.seed(seednum2)
        addc = np.random.randint(np.max((-maxlag, -col1)),
                                 np.min((maxlag, ncols-col1)), size=ndists)
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
    # Fit model using weighting by 1/N to weigh bins with more samples higher
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


def spherical(lag, range1, sill):  # , nugget=0):
    """
    Spherical variogram model assuming nugget = 0
    
    Args:
        lag: float or array of lags as # of pixels/cells
        range1 (float): range of spherical model
        sill (float): sill of spherical model
        
    Returns:
        semivariance as float or array, depending on type(lag)
    """
    nugget = 0.
    range1 = range1 / 1.

    out = nugget + sill * ((1.5 * (lag / range1)) -
                           (0.5 * ((lag / range1) ** 3.0)))
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
    convolving with the image, then multiply by std to equal sum of 
    std1*std2*corr*scale1*scale2 over each valid cell
    
    Args:
        stds (array): grid of standard deviation of model
        range1 (float): range of empirical variogram used to estimate
            correlation model
        sill1 (float): sill of empirical variogram used to estimate
            correlation model
        scale: float or array same size as std, factor to multiply by
            (area or population) and include in convolution
    
    Returns:
        variance of aggregate statistic
    
    """
    range5 = int(range1)
    # Prepare kernal that is size of range of spherical equation
    nrows = 2*range5 + 1
    ncols = 2*range5 + 1
    # get distance in row and col from center point
    rows = matlib.repmat(np.arange(nrows), ncols, 1).T - range5
    cols = matlib.repmat(np.arange(ncols), nrows, 1) - range5
    dists = np.sqrt(rows**2 + cols**2)
    # Convert from semivariance to correlation and build kernal
    kernal = (sill1-spherical(dists, range1, sill1))/sill1
    
    # convolve with stds, equivalent to sum of corr * std at each pixel for
    # within range1
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
