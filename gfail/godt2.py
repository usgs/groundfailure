#!/usr/bin/env python
"""
This module contains functions that can be used to run Newmark-based
mechanistic landslide models.
"""

# stdlib imports
import os.path
# import warnings
import collections
import tempfile
import shutil

# local imports
from mapio.shake import ShakeGrid
from mapio.gdal import GDALGrid
from mapio.reader import read, get_file_geodict
from mapio.geodict import GeoDict
from gfail.spatial import trim_ocean2

# third party imports
import numpy as np


def godt2008_2(shakefile, config, uncertfile=None, saveinputs=False,
               displmodel=None, bounds=None, slopediv=100.,
               codiv=10., numstd=None, trimfile=None):
    """
    This function runs the Godt and others (2008) global method for a given
    ShakeMap. The Factor of Safety is calculated using infinite slope analysis
    assumuing dry conditions. The method uses threshold newmark displacement
    and estimates areal coverage by doing the calculations for each slope
    quantile.

    Args:
        shakefile (str): Path to shakemap xml file.
        config (ConfigObj): ConfigObj of config file containing inputs required
            for running the model
        uncertfile (str): Path to shakemap uncertainty xml file (optional).
        saveinputs (bool): Whether or not to return the model input layers,
            False (default) returns only the model output (one layer).
        displmodel (str): Newmark displacement regression model to use

            * ``'J_PGA'`` (default) -- PGA-based model, equation 6 from
              Jibson (2007).
            * ``'J_PGA_M'`` -- PGA and M-based model, equation 7 from
              Jibson (2007).
            * ``'RS_PGA_M'`` -- PGA and M-based model from from Rathje and
              Saygili (2009).
            * ``'RS_PGA_PGV'`` -- PGA and PGV-based model, equation 6
              from Saygili and Rathje (2008).

        bounds (dict): Optional dictionary with keys 'xmin', 'xmax', 'ymin',
            'ymax' that defines a subset of the shakemap area to compute.
        slopediv (float): Divide slope by this number to get slope in degrees
            (Verdin datasets need to be divided by 100).
        codiv (float): Divide cohesion input layer by this number
            (For Godt method, need to divide by 10 because that is how it was
            calibrated).
        numstd (float): Number of (+/-) standard deviations to use if
            uncertainty is computed (uncertfile must be supplied).
        trimfile (str): shapefile of earth's land masses to trim offshore areas
            of model

    Returns:
        dict: Dictionary containing output and input layers (if
        saveinputs=True):

        .. code-block:: python

            {
                'grid': mapio grid2D object,
                'label': 'label for colorbar and top line of subtitle',
                'type': 'output or input to model',
                'description': {'name': 'short reference of model',
                                'longref': 'full model reference',
                                'units': 'units of output',
                                'shakemap': 'information about shakemap used',
                                'event_id': 'shakemap event id',
                                'parameters': 'dictionary of model parameters
                                               used'

                }
            }

    Raises:
         NameError: when unable to parse the config correctly (probably a
             formatting issue in the configfile) or when unable to find the
             shakefile (Shakemap filepath) -- these cause program to end.

    """
    # TODO:
    #    - Add 'all' -- averages Dn from all four equations, add term to
    #      convert PGA and PGV to Ia and use other equations, add Ambraseys and
    #      Menu (1988) option.

    # Empty refs
    slopesref = 'unknown'
    slopelref = 'unknown'
    cohesionlref = 'unknown'
    cohesionsref = 'unknown'
    frictionsref = 'unknown'
    frictionlref = 'unknown'
    modellref = 'unknown'
    modelsref = 'unknown'

    # See if trimfile exists
    if trimfile is not None:
        if not os.path.exists(trimfile):
            print('trimfile defined does not exist: %s\n'
                  'Ocean will not be trimmed' % trimfile)
            trimfile = None
        if os.path.splitext(trimfile)[1] != '.shp':
            print('trimfile must be a shapefile, ocean will not be trimmed')
            trimfile = None

    # Parse config
    try:    # May want to add error handling so if refs aren't given, just
        # includes unknown
        slopefilepath = config['godt_2008']['layers']['slope']['filepath']
        cohesionfile = config['godt_2008']['layers']['cohesion']['file']
        frictionfile = config['godt_2008']['layers']['friction']['file']

        thick = float(config['godt_2008']['parameters']['thick'])
        uwt = float(config['godt_2008']['parameters']['uwt'])
        nodata_cohesion = \
            float(config['godt_2008']['parameters']['nodata_cohesion'])
        nodata_friction = \
            float(config['godt_2008']['parameters']['nodata_friction'])
        dnthresh = float(config['godt_2008']['parameters']['dnthresh'])
        fsthresh = float(config['godt_2008']['parameters']['fsthresh'])
        acthresh = float(config['godt_2008']['parameters']['acthresh'])
        try:
            slopemin = float(config['godt_2008']['parameters']['slopemin'])
        except Exception:
            slopemin = 0.01
            print('No slopemin found in config file, using 0.01 deg '
                  'for slope minimum')
    except Exception as e:
        raise NameError('Could not parse configfile, %s' % e)

    if displmodel is None:
        try:
            displmodel = config['godt_2008']['parameters']['displmodel']
        except Exception:
            print('No regression model specified, using default of J_PGA_M')
            displmodel = 'J_PGA_M'

    # TO DO: ADD ERROR CATCHING ON UNITS, MAKE SURE THEY ARE WHAT THEY SHOULD
    #        BE FOR THIS MODEL

    try:  # Try to fetch source information from config
        modelsref = config['godt_2008']['shortref']
        modellref = config['godt_2008']['longref']
        slopesref = config['godt_2008']['layers']['slope']['shortref']
        slopelref = config['godt_2008']['layers']['slope']['longref']
        cohesionsref = config['godt_2008']['layers']['cohesion']['shortref']
        cohesionlref = config['godt_2008']['layers']['cohesion']['longref']
        frictionsref = config['godt_2008']['layers']['friction']['shortref']
        frictionlref = config['godt_2008']['layers']['friction']['longref']
    except Exception:
        print('Was not able to retrieve all references from config file. '
              'Continuing')

    # Figure out how/if need to cut anything
    geodict = ShakeGrid.getFileGeoDict(shakefile)  # , adjust='res')
    if bounds is not None:  # Make sure bounds are within ShakeMap Grid
        if geodict.xmin < geodict.xmax:  # only if signs are not opposite
            if (geodict.xmin > bounds['xmin'] or
                    geodict.xmax < bounds['xmax'] or
                    geodict.ymin > bounds['ymin'] or
                    geodict.ymax < bounds['ymax']):
                print('Specified bounds are outside shakemap area, using '
                      'ShakeMap bounds instead.')
                bounds = None

    if bounds is not None:
        tempgdict = GeoDict.createDictFromBox(
            bounds['xmin'], bounds['xmax'],
            bounds['ymin'], bounds['ymax'],
            geodict.dx, geodict.dy, inside=False)
        # If Shakemap geodict crosses 180/-180 line, fix geodict so things
        # don't break
        if geodict.xmin > geodict.xmax:
            if tempgdict.xmin < 0:
                geodict._xmin -= 360.
            else:
                geodict._xmax += 360.
        geodict = geodict.getBoundsWithin(tempgdict)

    slpfile = os.path.join(slopefilepath, 'slope_min.bil')
    basegeodict = get_file_geodict(slpfile)
    if basegeodict == geodict:
        sampledict = geodict
    else:
        sampledict = basegeodict.getBoundsWithin(geodict)

    # Do we need to subdivide baselayer?
    if 'divfactor' in config['godt_2008'].keys():
        divfactor = float(config['godt_2008']['divfactor'])
        if divfactor != 1.:
            # adjust sampledict so everything will be resampled (cut one cell
            # of each edge so will be inside bounds)
            newxmin = sampledict.xmin - sampledict.dx / 2. + \
                sampledict.dx / (2. * divfactor) + sampledict.dx
            newymin = sampledict.ymin - sampledict.dy / 2. + \
                sampledict.dy / (2. * divfactor) + sampledict.dy
            newxmax = sampledict.xmax + sampledict.dx / 2. - \
                sampledict.dx / (2. * divfactor) - sampledict.dx
            newymax = sampledict.ymax + sampledict.dy / 2. - \
                sampledict.dy / (2. * divfactor) - sampledict.dy
            newdx = sampledict.dx / divfactor
            newdy = sampledict.dy / divfactor

            sampledict = GeoDict.createDictFromBox(
                newxmin, newxmax, newymin,
                newymax, newdx, newdy, inside=True)

    tmpdir = tempfile.mkdtemp()

    # Load in ShakeMap and get new geodictionary
    temp = ShakeGrid.load(shakefile)  # , adjust='res')
    pga = temp.getLayer('pga')
    pga = pga.interpolate2(sampledict, method='linear')
    pgv = temp.getLayer('pgv')
    pgv = pgv.interpolate2(sampledict, method='linear')
    # sampledict = pga.getGeoDict()

    t2 = temp.getEventDict()
    M = t2['magnitude']
    event_id = t2['event_id']
    shakedict = temp.getShakeDict()
    del(temp)

    # read in uncertainty if present
    if uncertfile is not None:
        try:
            temp = ShakeGrid.load(uncertfile)  # , adjust='res')
            uncertpga = temp.getLayer('stdpga')
            uncertpga = uncertpga.interpolate2(sampledict, method='linear')
            uncertpgv = uncertpgv = temp.getLayer('stdpgv')
            uncertpgv.interpolate2(sampledict, method='linear')
        except Exception:
            print('Could not read uncertainty file, ignoring uncertainties')
            uncertfile = None

        if numstd is None:
            numstd = 1.

    # Read in all the slope files, divide all by 100 to get to slope in
    # degrees (because input files are multiplied by 100.)
    slopes = []
    quantiles = ['slope_min.bil', 'slope10.bil', 'slope30.bil', 'slope50.bil',
                 'slope70.bil', 'slope90.bil', 'slope_max.bil']
    for quant in quantiles:
        tmpslp = read(os.path.join(slopefilepath, quant),
                      samplegeodict=sampledict)
        tgd = tmpslp.getGeoDict()
        if tgd != sampledict:
            raise Exception('Input layers are not aligned to same geodict')
        else:
            slopes.append(tmpslp.getData() / slopediv)

    slopestack = np.dstack(slopes)

    # Change any zero slopes to a very small number to avoid dividing by
    # zero later
    slopestack[slopestack == 0] = 1e-8

    # Read in the cohesion and friction files and duplicate layers so they
    # are same shape as slope structure

    tempco = read(cohesionfile, samplegeodict=sampledict,
                  resample=True, method='nearest',
                  doPadding=True, padValue=np.nan)
    tempco = tempco.getData()[:, :, np.newaxis] / codiv
    cohesion = np.repeat(tempco, 7, axis=2)
    cohesion[cohesion == -999.9] = nodata_cohesion
    cohesion = np.nan_to_num(cohesion)
    cohesion[cohesion == 0] = nodata_cohesion

    tempfric = read(frictionfile, samplegeodict=sampledict,
                    resample=True, method='nearest',
                    doPadding=True, padValue=np.nan)
    tempfric = tempfric.getData().astype(float)[:, :, np.newaxis]
    friction = np.repeat(tempfric, 7, axis=2)
    friction[friction == -9999] = nodata_friction
    friction = np.nan_to_num(friction)
    friction[friction == 0] = nodata_friction

    # Do the calculations using Jibson (2007) PGA only model for Dn
    FS = (cohesion / (uwt * thick * np.sin(slopestack * (np.pi / 180.))) +
          np.tan(friction * (np.pi / 180.)) / np.tan(slopestack * (np.pi / 180.)))
    FS[FS < fsthresh] = fsthresh

    # Compute critical acceleration, in g
    # This gives ac in g, equations that multiply by g give ac in m/s2
    Ac = (FS - 1) * np.sin(slopestack * (np.pi / 180.)).astype(float)
    Ac[Ac < acthresh] = acthresh

    # Get PGA in g (PGA is %g in ShakeMap, convert to g)
    PGA = np.repeat(pga.getData()[:, :, np.newaxis] / 100., 7,
                    axis=2).astype(float)
    if 'PGV' in displmodel:  # Load in PGV also, in cm/sec
        PGV = np.repeat(pgv.getData()[:, :, np.newaxis], 7,
                        axis=2).astype(float)
    else:
        PGV = None

    if uncertfile is not None:
        stdpga = np.repeat(uncertpga.getData()[:, :, np.newaxis], 7,
                           axis=2).astype(float)
        stdpgv = np.repeat(uncertpgv.getData()[:, :, np.newaxis], 7,
                           axis=2).astype(float)
        # estimate PGA +- 1std
        PGAmin = np.exp(np.log(PGA * 100) - numstd * stdpga) / 100
        PGAmax = np.exp(np.log(PGA * 100) + numstd * stdpga) / 100
        if 'PGV' in displmodel:
            PGVmin = np.exp(np.log(PGV) - numstd * stdpgv)
            PGVmax = np.exp(np.log(PGV) + numstd * stdpgv)
        else:
            PGVmin = None
            PGVmax = None

    # Ignore errors so still runs when Ac > PGA, just leaves nan instead
    # of crashing.
    np.seterr(invalid='ignore')

    Dn, logDnstd, logtype = NMdisp(Ac, PGA, model=displmodel, M=M, PGV=PGV)
    if uncertfile is not None:
        Dnmin, logDnstdmin, logtype = NMdisp(
            Ac, PGAmin, model=displmodel, M=M, PGV=PGVmin)
        Dnmax, logDnstdmax, logtype = NMdisp(
            Ac, PGAmax, model=displmodel, M=M, PGV=PGVmax)

    PROB = Dn.copy()
    PROB[PROB < dnthresh] = 0.
    PROB[PROB >= dnthresh] = 1.
    PROB = np.sum(PROB, axis=2)
    if uncertfile is not None:
        PROBmin = Dnmin.copy()
        PROBmin[PROBmin <= dnthresh] = 0.
        PROBmin[PROBmin > dnthresh] = 1.
        PROBmin = np.sum(PROBmin, axis=2)
        PROBmax = Dnmax.copy()
        PROBmax[PROBmax <= dnthresh] = 0.
        PROBmax[PROBmax > dnthresh] = 1.
        PROBmax = np.sum(PROBmax, axis=2)

    PROB[PROB == 1.] = 0.01
    PROB[PROB == 2.] = 0.10
    PROB[PROB == 3.] = 0.30
    PROB[PROB == 4.] = 0.50
    PROB[PROB == 5.] = 0.70
    PROB[PROB == 6.] = 0.90
    PROB[PROB == 7.] = 0.99

    if uncertfile is not None:
        PROBmin[PROBmin == 1.] = 0.01
        PROBmin[PROBmin == 2.] = 0.10
        PROBmin[PROBmin == 3.] = 0.30
        PROBmin[PROBmin == 4.] = 0.50
        PROBmin[PROBmin == 5.] = 0.70
        PROBmin[PROBmin == 6.] = 0.90
        PROBmin[PROBmin == 7.] = 0.99
        PROBmax[PROBmax == 1.] = 0.01
        PROBmax[PROBmax == 2.] = 0.10
        PROBmax[PROBmax == 3.] = 0.30
        PROBmax[PROBmax == 4.] = 0.50
        PROBmax[PROBmax == 5.] = 0.70
        PROBmax[PROBmax == 6.] = 0.90
        PROBmax[PROBmax == 7.] = 0.99

    if slopemin is not None:
        PROB[slopestack[:, :, 6] <= slopemin] = 0.
        # uncert too
        if uncertfile is not None:
            PROBmin[slopestack[:, :, 6] <= slopemin] = 0.
            PROBmax[slopestack[:, :, 6] <= slopemin] = 0.

    # Turn output and inputs into into grids and put in mapLayers dictionary
    maplayers = collections.OrderedDict()

    shakedetail = '%s_ver%s' % (
        shakedict['shakemap_id'], shakedict['shakemap_version'])

    description = {
        'name': modelsref,
        'longref': modellref,
        'units': 'Proportion of Area Affected',
        'shakemap': shakedetail,
        'event_id': event_id,
        'parameters': {
            'displmodel': displmodel,
            'thickness_m': thick,
            'unitwt_kNm3': uwt,
            'dnthresh_cm': dnthresh,
            'acthresh_g': acthresh,
            'fsthresh': fsthresh,
            'modeltype': 'Landslide'
        }
    }
    PROBgrid = GDALGrid(PROB, sampledict)
    if trimfile is not None:
        PROBgrid = trim_ocean2(PROBgrid, trimfile)

    maplayers['model'] = {
        'grid': PROBgrid,
        'label': 'Landslide - Proportion of Area Affected',
        'type': 'output',
        'description': description
    }

    if uncertfile is not None:
        PROBmingrid = GDALGrid(PROBmin, sampledict)
        PROBmaxgrid = GDALGrid(PROBmax, sampledict)
        if trimfile is not None:
            PROBmingrid = trim_ocean2(PROBmingrid, trimfile)
            PROBmaxgrid = trim_ocean2(PROBmaxgrid, trimfile)
        maplayers['modelmin'] = {
            'grid': PROBmingrid,
            'label': 'Landslide Probability-%1.2fstd' % numstd,
            'type': 'output',
            'description': description
        }
        maplayers['modelmax'] = {
            'grid': PROBmaxgrid,
            'label': 'Landslide Probability+%1.2fstd' % numstd,
            'type': 'output',
            'description': description
        }

    if saveinputs is True:
        maplayers['pga'] = {
            'grid': GDALGrid(PGA[:, :, 0], sampledict),
            'label': 'PGA (g)',
            'type': 'input',
            'description': {
                'units': 'g',
                'shakemap': shakedetail}
        }
        if 'PGV' in displmodel:
            maplayers['pgv'] = {
                'grid': GDALGrid(PGV[:, :, 0], sampledict),
                'label': 'PGV (cm/s)',
                'type': 'input',
                'description': {
                    'units': 'cm/s',
                    'shakemap': shakedetail}
            }
        maplayers['minFS'] = {
            'grid': GDALGrid(np.min(FS, axis=2), sampledict),
            'label': 'Min Factor of Safety',
            'type': 'input',
            'description': {
                'units': 'unitless'}
        }
        maplayers['max slope'] = {
            'grid': GDALGrid(slopestack[:, :, -1], sampledict),
            'label': r'Maximum slope ($^\circ$)',
            'type': 'input',
            'description': {
                'units': 'degrees',
                'name': slopesref,
                'longref': slopelref}
        }
        maplayers['cohesion'] = {
            'grid': GDALGrid(cohesion[:, :, 0], sampledict),
            'label': 'Cohesion (kPa)',
            'type': 'input',
            'description': {
                'units': 'kPa (adjusted)',
                'name': cohesionsref,
                'longref': cohesionlref}
        }
        maplayers['friction angle'] = {
            'grid': GDALGrid(friction[:, :, 0], sampledict),
            'label': r'Friction angle ($^\circ$)',
            'type': 'input',
            'description': {
                'units': 'degrees',
                'name': frictionsref,
                'longref': frictionlref}
        }
        if uncertfile is not None:
            maplayers['pgamin'] = {
                'grid': GDALGrid(PGAmin[:, :, 0], sampledict),
                'label': 'PGA - %1.2fstd (g)' % numstd,
                'type': 'input',
                'description': {
                    'units': 'g',
                    'shakemap': shakedetail}
            }
            maplayers['pgamax'] = {
                'grid': GDALGrid(PGAmax[:, :, 0], sampledict),
                'label': 'PGA + %1.2fstd (g)' % numstd,
                'type': 'input',
                'description': {
                    'units': 'g',
                    'shakemap': shakedetail}
            }
        if 'PGV' in displmodel:
            if uncertfile is not None:
                maplayers['pgvmin'] = {
                    'grid': GDALGrid(PGVmin[:, :, 0], sampledict),
                    'label': 'PGV - %1.2fstd (cm/s)' % numstd,
                    'type': 'input',
                    'description': {
                        'units': 'cm/s',
                        'shakemap': shakedetail}
                }
                maplayers['pgvmax'] = {
                    'grid': GDALGrid(PGVmax[:, :, 0], sampledict),
                    'label': 'PGV + %1.2fstd (cm/s)' % numstd,
                    'type': 'input',
                    'description': {
                        'units': 'cm/s',
                        'shakemap': shakedetail}
                }

    shutil.rmtree(tmpdir)

    return maplayers


def NMdisp(Ac, PGA, model='J_PGA', M=None, PGV=None):
    """
    PGA-based Newmark Displacement model

    Args:
        Ac (array): NxM array of critical accelerations in units of g.
        PGA (array): NxM Array of PGA values in units of g.
        model (str):

            * ``'J_PGA'`` -- PGA only model from Jibson (2007), equation 6.
              Applicable for Magnitude range of dataset (5.3-7.6).
            * ``'J_PGA_M'`` -- PGA-and M- based Newmark Displacement model
              from Jibson(2007), equation 7. Applicable for Magnitude range
              of dataset (5.3-7.6).
            * ``'RS_PGA_M'`` -- PGA and M-based Newmark displacement model
              from Rathje and Saygili (2009).
            * ``'RS_PGA_PGV'`` -- PGA and PGV-based model from Saygili and
              Rathje (2008) -- eq 6.
            * ``'BT_PGA_M'`` -- PGA and M-based model from Bray and
              Travasarou (2007) assuming natural fundamental period of
              sliding mass Ts = 0 (eq 6).

        M (float): Magnitude -- only needed for models with M in the name.
        PGV (float): NxM Array of PGV values in units of cm/sec -- only needed
            for models with PGV in the name.

    Returns:
        tuple: (Dn, logDnstd, logtype) where:
            * Dn: Newmark displacement in cm
            * logDnstd: Log of standard deviation of Dn
            * logtype: Type of log used in logDnstd (log10 or ln)
    """
    # Deal with non-array inputs
    if isinstance(Ac, float) or isinstance(Ac, int):
        flag = 1
    else:
        flag = 0
    # Ignore errors so still runs when Ac > PGA, just leaves nan instead of
    # crashing
    np.seterr(invalid='ignore')

    if model == 'J_PGA':
        C1 = 0.215  # additive constant in newmark displacement calculation
        C2 = 2.341  # first exponential constant
        C3 = -1.438  # second exponential constant
        Dn = np.array(
            10.**(C1 + np.log10(((1 - Ac / PGA)**C2) * (Ac / PGA)**C3)))
        Dn[np.isnan(Dn)] = 0.
        logDnstd = np.ones(np.shape(Dn)) * 0.51
        logtype = 'log10'

    elif model == 'J_PGA_M':
        if M is None:
            raise Exception('M (magnitude) not found, cannot use RS_PGA_M '
                            'model')
        else:
            C1 = -2.71  # additive constant in newmark displacement calculation
            C2 = 2.335  # first exponential constant
            C3 = -1.478  # second exponential constant
            C4 = 0.424
            Dn = np.array(10.**(C1 + np.log10(((1 - Ac / PGA)**C2) * (Ac / PGA)**C3) +
                                C4 * M))
            Dn[np.isnan(Dn)] = 0.
            logDnstd = np.ones(np.shape(Dn)) * 0.454
            logtype = 'log10'

    elif model == 'RS_PGA_M':
        if M is None:
            raise Exception('You must enter a value for M to use the '
                            'RS_PGA_M model')
        C1 = 4.89
        C2 = -4.85
        C3 = -19.64
        C4 = 42.49
        C5 = -29.06
        C6 = 0.72
        C7 = 0.89
        # Equation from Saygili and Rathje (2008)/Rathje and Saygili (2009)
        Dn = np.array(np.exp(C1 + C2 * (Ac / PGA) + C3 * (Ac / PGA)**2 +
                             C4 * (Ac / PGA)**3 + C5 * (Ac / PGA)**4 +
                             C6 * np.log(PGA) + C7 * (M - 6)))
        Dn[np.isnan(Dn)] = 0.
        logDnstd = 0.732 + 0.789 * (Ac / PGA) - 0.539 * (Ac / PGA)**2
        logtype = 'ln'

    elif model == 'RS_PGA_PGV':
        if PGV is None:
            raise Exception('You must enter a value for M to use the '
                            'RS_PGA_PGV model')
        C1 = -1.56
        C2 = -4.58
        C3 = -20.84
        C4 = 44.75
        C5 = -30.50
        C6 = -0.64
        C7 = 1.55
        # Equation from Saygili and Rathje (2008)/Rathje and Saygili (2009)
        Dn = np.array(np.exp(C1 + C2 * (Ac / PGA) + C3 * (Ac / PGA)**2 +
                             C4 * (Ac / PGA)**3 + C5 * (Ac / PGA)**4 +
                             C6 * np.log(PGA) + C7 * np.log(PGV)))
        Dn[np.isnan(Dn)] = 0.
        logDnstd = 0.405 + 0.524 * (Ac / PGA)
        logtype = 'ln'

    elif model == 'BT_PGA_M':
        if M is None:
            raise Exception('You must enter a value for M to use the '
                            'BT_PGA_M model')
        Dn = np.array(
            np.exp(-0.22 - 2.83 * np.log(Ac) - 0.333 * (np.log(Ac))**2 +
                   0.566 * np.log(Ac) * np.log(PGA) +
                   3.04 * np.log(PGA) - 0.244 * (np.log(PGA))**2 + 0.278 * (M - 7.)))
        Dn[np.isnan(Dn)] = 0.
        logDnstd = np.ones(np.shape(Dn)) * 0.66
        logtype = 'log10'

    if flag == 1:
        Dn = float(Dn)
        logDnstd = float(logDnstd)

    return Dn, logDnstd, logtype
