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
                    Hagg['hagg_std_%1.2fg' % (shaket/100.,)] = fullvario(std, model, probthresh)
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
                Hagg['std_0.00g'] = fullvario(std, model, probthresh)

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
            mu = np.nansum(popdat * prop * modresamp * threshmult)
            exp_pop['exp_pop_%1.2fg' % (shaket/100.,)] = mu
            elim = maxP*np.nansum(popdat * prop * threshmult)
            exp_pop['elim_%1.2fg' % (shaket/100.,)] = elim
            if stdgrid2D is not None:
                totalmax = np.nansum(popdat * propstd * modresampstd * threshmult)
                totalmin = np.sqrt(np.nansum((popdat * propstd * modresampstd * threshmult)**2.))
                if stdtype == 'max':
                    exp_pop['exp_std_%1.2fg' % (shaket/100.,)] = totalmax
                elif stdtype == 'min':
                    exp_pop['exp_std_%1.2fg' % (shaket/100.,)] = totalmin
                elif stdtype == 'mean':
                    exp_pop['exp_std_%1.2fg' % (shaket/100.,)]=(totalmax+totalmin)/2.
                else:
                    #TODO
                    exp_pop['exp_std_%1.2fg' % (shaket/100.,)] = fullvario(propstd, popdat, threshmult)
                    
                    
                    
                    
                # Beta distribution shape factors
                var = exp_pop['exp_std_%1.2fg' % (shaket/100.,)]**2.
                exp_pop['p_exp_%1.2fg' % (shaket/100.,)] = (mu/elim)*((elim*mu-mu**2)/var-1)
                exp_pop['q_exp_%1.2fg' % (shaket/100.,)] = (1-mu/elim)*((elim*mu-mu**2)/var-1)

    else:
        mu = np.nansum(popdat * prop * modresamp)
        exp_pop['exp_pop_0.00g'] = mu
        elim = maxP*np.nansum(popdat * prop)
        exp_pop['elim_0.00g'] = elim
        if stdgrid2D is not None:
            totalmax = np.nansum(popdat * propstd * modresampstd)
            totalmin = np.sqrt(np.nansum((popdat * propstd * modresampstd)**2.))
            if stdtype == 'max':
                exp_pop['exp_std_0.00g'] = totalmax
            elif stdtype == 'min':
                exp_pop['exp_std_0.00g'] = totalmin
            elif stdtype == 'mean':
                exp_pop['exp_std_0.00g'] = (totalmax+totalmin)/2.
            else:
                #TODO
                exp_pop['exp_std_%1.2fg' % (shaket/100.,)] = fullvario(propstd, popdat, threshmult)
                    
            # Beta distribution shape factors
            var = exp_pop['exp_std_0.00g']**2.
            exp_pop['exp_std_0.00g'] = (mu/elim)*((elim*mu-mu**2)/var-1)
            exp_pop['exp_std_0.00g'] = (1-mu/elim)*((elim*mu-mu**2)/var-1)
            #exp_pop['exp_std_0.00g'] = np.nansum(popdat * propstd * modresampstd)

    return exp_pop


def fullvario():
    """
    
    """
    pass
