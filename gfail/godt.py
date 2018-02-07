#!/usr/bin/env python
"""
This module contains functions that can be used to run Newmark-based
mechanistic landslide models.
"""

# stdlib imports
import os.path
import warnings
import collections

# local imports
from mapio.shake import ShakeGrid
from mapio.gdal import GDALGrid
from mapio.geodict import GeoDict

# third party imports
import numpy as np

# turn off all warnings...
warnings.filterwarnings('ignore')


def godt2008(shakefile, config, uncertfile=None, saveinputs=False,
             displmodel=None, bounds=None, slopediv=100.,
             codiv=10., numstd=None):
    """
    This function runs the Godt et al. (2008) global method for a given
    ShakeMap. The Factor of Safety is calculated using infinite slope analysis
    assumuing dry conditions. The method uses threshold newmark displacement
    and estimates areal coverage by doing the calculations for each slope
    quantile.

    TODO:
        - Add 'all' -- averages Dn from all four equations, add term to
          convert PGA and PGV to Ia and use other equations, add Ambraseys and
          Menu (1988) option.

    Args:
        shakefile (str): Path to shakemap xml file.
        config (ConfigObj): ConfigObj of config file containing inputs required
            for running the model
        uncertfile (str): Path to shakemap uncertainty xml file.
        saveinputs (bool): Whether or not to return the model input layers,
            False (defeault) returns only the model output (one layer).
        displmodel (str): Newmark displacement regression model to use

            * ``'J_PGA'`` (default) -- PGA-based model, equation 6 from
              Jibson (2007).
            * ``'J_PGA_M'`` -- PGA and M-based model, equation 7 from
              Jibson (2007).
            * ``'RS_PGA_M'`` -- PGA and M-based model from from Rathje and
              Saygili (2009).
            * ``'RS_PGA_PGV'`` -- PGA and PGV-based model, equation 6
              from Saygili and Rathje (2008).

        probtype (str): Method used to estimate probability.

            * ``'jibson2000'`` uses equation 5 from Jibson et al. (2000) to
              estimate probability from Newmark displacement.
            * ``'threshold'`` uses a specified threshold of Newmark displacement
              (defined in config file) and assumes anything greater than this
              threshold fails.

        slopediv (float): Divide slope by this number to get slope in degrees
            (Verdin datasets need to be divided by 100).
        codiv (float): Divide cohesion by this number to get reasonable numbers
            (For Godt method, need to divide by 10 because that is how it was
            calibrated, but values are reasonable without multiplying for
            regular analysis).
        numstd (float): Number of +/- standard deviations to use if uncertainty
            is computed (uncertfile is not None).

    Returns:
        dict: Dictionary containing output and input layers (if
        saveinputs=True) along with metadata formatted like:

        .. code-block:: python

            {
                'grid': mapio grid2D object,
                'label': 'label for colorbar and top line of subtitle',
                'type': 'output or input to model',
                'description': 'detailed description'
            }

    Raises:
         NameError: when unable to parse the config correctly (probably a
             formatting issue in the configfile) or when unable to find the
             shakefile (Shakemap filepath) -- these cause program to end.

    """

    # Empty refs
    slopesref = 'unknown'
    slopelref = 'unknown'
    cohesionlref = 'unknown'
    cohesionsref = 'unknown'
    frictionsref = 'unknown'
    frictionlref = 'unknown'
    modellref = 'unknown'
    modelsref = 'unknown'

    # Parse config
    try:    # May want to add error handling so if refs aren't given, just
            # includes unknown
        slopefilepath = config['godt_2008']['layers']['slope']['filepath']
        slopeunits = config['godt_2008']['layers']['slope']['units']
        cohesionfile = config['godt_2008']['layers']['cohesion']['file']
        cohesionunits = config['godt_2008']['layers']['cohesion']['units']
        frictionfile = config['godt_2008']['layers']['friction']['file']
        frictionunits = config['godt_2008']['layers']['friction']['units']

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
        except:
            slopemin = 0.01
            print('No slopemin found in config file, using 0.01 deg '
                  'for slope minimum')
    except Exception as e:
        raise NameError('Could not parse configfile, %s' % e)

    if displmodel is None:
        try:
            displmodel = config['godt_2008']['parameters']['displmodel']
        except:
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
    except:
        print('Was not able to retrieve all references from config file. '
              'Continuing')

    sampledict = ShakeGrid.getFileGeoDict(shakefile, adjust='res')
    
    # Do we need to subdivide baselayer?
    resample = False
    if 'divfactor' in config['godt_2008'].keys():
        divfactor = float(config['godt_2008']['divfactor'])
        if divfactor != 1.:
            # adjust sampledict so everything will be resampled (cut one cell of each edge so will be inside bounds)
            newxmin = sampledict.xmin - sampledict.dx/2. + sampledict.dx/(2.*divfactor) + sampledict.dx
            newymin = sampledict.ymin - sampledict.dy/2. + sampledict.dy/(2.*divfactor) + sampledict.dy
            newxmax = sampledict.xmax + sampledict.dx/2. - sampledict.dx/(2.*divfactor) - sampledict.dx
            newymax = sampledict.ymax + sampledict.dy/2. - sampledict.dy/(2.*divfactor) - sampledict.dy
            newdx = sampledict.dx/divfactor
            newdy = sampledict.dy/divfactor
            
            sampledict = GeoDict.createDictFromBox(newxmin, newxmax, newymin,
                                                   newymax, newdx, newdy, inside=True)
            resample = True
    
    if bounds is not None:  # Make sure bounds are within ShakeMap Grid
        if (sampledict.xmin > bounds['xmin'] or
                sampledict.xmax < bounds['xmax'] or
                sampledict.ymin > bounds['ymin'] or
                sampledict.ymax < bounds['ymax']):
            print('Specified bounds are outside shakemap area, using '
                  'ShakeMap bounds instead')
            bounds = None
    if bounds is not None:
        tempgdict = GeoDict.createDictFromBox(
            bounds['xmin'], bounds['xmax'], bounds['ymin'], bounds['ymax'],
            sampledict.dx, sampledict.dy, inside=False)
        if tempgdict == sampledict:
            gdict = tempgdict
        else:
            gdict = sampledict.getBoundsWithin(tempgdict)
        shakemap = ShakeGrid.load(shakefile, samplegeodict=gdict,
                                  resample=True, method='linear',
                                  adjust='bounds')
    elif resample:
        shakemap = ShakeGrid.load(shakefile, samplegeodict=sampledict, 
                                  resample=True, method='linear',
                                  adjust='bounds')
    else:
        shakemap = ShakeGrid.load(shakefile, adjust='res', resample=False)

    shkgdict = shakemap.getGeoDict()  # Get updated geodict
    t2 = shakemap.getEventDict()
    M = t2['magnitude']
    event_id = t2['event_id']

    # read in uncertainty if present
    if uncertfile is not None:
        try:
            uncert = ShakeGrid.load(uncertfile, samplegeodict=shkgdict,
                                    resample=True, method='linear',
                                    adjust='res')
        except:
            print('Could not read uncertainty file, ignoring uncertainties')
            uncertfile = None
        if numstd is None:
            numstd = 1.

    # Read in all the slope files, divide all by 100 to get to slope in
    # degrees (because input files are multiplied by 100.)
    slopes = []
    slopes.append(
        GDALGrid.load(
            os.path.join(slopefilepath, 'slope_min.bil'),
            samplegeodict=shkgdict,
            resample=True, method='linear').getData()/slopediv)
    slopes.append(
        GDALGrid.load(
            os.path.join(slopefilepath, 'slope10.bil'),
            samplegeodict=shkgdict,
            resample=True, method='linear').getData()/slopediv)
    slopes.append(
        GDALGrid.load(
            os.path.join(slopefilepath, 'slope30.bil'),
            samplegeodict=shkgdict,
            resample=True, method='linear').getData()/slopediv)
    slopes.append(
        GDALGrid.load(
            os.path.join(slopefilepath, 'slope50.bil'),
            samplegeodict=shkgdict,
            resample=True, method='linear').getData()/slopediv)
    slopes.append(
        GDALGrid.load(
            os.path.join(slopefilepath, 'slope70.bil'),
            samplegeodict=shkgdict,
            resample=True, method='linear').getData()/slopediv)
    slopes.append(
        GDALGrid.load(
            os.path.join(slopefilepath, 'slope90.bil'),
            samplegeodict=shkgdict,
            resample=True, method='linear').getData()/slopediv)
    slopes.append(
        GDALGrid.load(
            os.path.join(slopefilepath, 'slope_max.bil'),
            samplegeodict=shkgdict,
            resample=True, method='linear').getData()/slopediv)
    slopestack = np.dstack(slopes)

    # Change any zero slopes to a very small number to avoid dividing by
    # zero later
    slopestack[slopestack == 0] = 1e-8

    # Read in the cohesion and friction files and duplicate layers so they
    # are same shape as slope structure
    cohesion = np.repeat(
        GDALGrid.load(cohesionfile,
                      samplegeodict=shkgdict,
                      resample=True,
                      method='nearest').getData()[:, :, np.newaxis]/codiv,
        7,
        axis=2)
    cohesion[cohesion == -999.9] = nodata_cohesion
    cohesion = np.nan_to_num(cohesion)
    cohesion[cohesion == 0] = nodata_cohesion
    friction = np.repeat(
        GDALGrid.load(frictionfile,
                      samplegeodict=shkgdict,
                      resample=True,
                      method='nearest').getData().astype(float)[:, :,
                                                                np.newaxis],
        7,
        axis=2)
    friction[friction == -9999] = nodata_friction
    friction = np.nan_to_num(friction)
    friction[friction == 0] = nodata_friction

    # Do the calculations using Jibson (2007) PGA only model for Dn
    FS = (cohesion/(uwt*thick*np.sin(slopestack*(np.pi/180.))) +
          np.tan(friction*(np.pi/180.))/np.tan(slopestack*(np.pi/180.)))
    FS[FS < fsthresh] = fsthresh

    # Compute critical acceleration, in g
    # This gives ac in g, equations that multiply by g give ac in m/s2
    Ac = (FS-1)*np.sin(slopestack*(np.pi/180.)).astype(float)
    Ac[Ac < acthresh] = acthresh

    # Get PGA in g (PGA is %g in ShakeMap, convert to g)
    PGA = np.repeat(
        shakemap.getLayer('pga').getData()[:, :, np.newaxis]/100.,
        7,
        axis=2).astype(float)
    if 'PGV' in displmodel:  # Load in PGV also, in cm/sec
        PGV = np.repeat(
            shakemap.getLayer('pgv').getData()[:, :, np.newaxis],
            7,
            axis=2).astype(float)
    else:
        PGV = None

    if uncertfile is not None:
        stdpga = np.repeat(
            uncert.getLayer('stdpga').getData()[:, :, np.newaxis],
            7,
            axis=2).astype(float)
        stdpgv = np.repeat(
            uncert.getLayer('stdpgv').getData()[:, :, np.newaxis],
            7,
            axis=2).astype(float)
        # estimate PGA +- 1std
        PGAmin = np.exp(np.log(PGA*100) - numstd*stdpga)/100
        PGAmax = np.exp(np.log(PGA*100) + numstd*stdpga)/100
        if 'PGV' in displmodel:
            PGVmin = np.exp(np.log(PGV) - numstd*stdpgv)
            PGVmax = np.exp(np.log(PGV) + numstd*stdpgv)
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

    temp = shakemap.getShakeDict()
    shakedetail = '%s_ver%s' % (temp['shakemap_id'], temp['shakemap_version'])

    description = {
        'name': modelsref,
        'longref': modellref,
        'units': 'coverage',
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

    maplayers['model'] = {
        'grid': GDALGrid(PROB, shakemap.getGeoDict()),
        'label': 'Landslide Areal coverage',
        'type': 'output',
        'description': description
    }

    if uncertfile is not None:
        maplayers['modelmin'] = {
            'grid': GDALGrid(PROBmin, shkgdict),
            'label': 'Probability-%1.2fstd' % numstd,
            'type': 'output',
            'description': description
        }
        maplayers['modelmax'] = {
            'grid': GDALGrid(PROBmax, shkgdict),
            'label': 'Probability+%1.2fstd' % numstd,
            'type': 'output',
            'description': description
        }

    if saveinputs is True:
        maplayers['pga'] = {
            'grid': GDALGrid(PGA[:, :, 0], shakemap.getGeoDict()),
            'label': 'PGA (g)',
            'type': 'input',
            'description': {
                'units': 'g',
                'shakemap': shakedetail}
        }
        if 'PGV' in displmodel:
            maplayers['pgv'] = {
                'grid': GDALGrid(PGV[:, :, 0], shakemap.getGeoDict()),
                'label': 'PGV (cm/s)',
                'type': 'input',
                'description': {
                    'units': 'cm/s',
                    'shakemap': shakedetail}
            }
        maplayers['minFS'] = {
            'grid': GDALGrid(np.min(FS, axis=2), shakemap.getGeoDict()),
            'label': 'Min Factor of Safety',
            'type': 'input',
            'description': {
                'units': 'unitless'}
        }
        maplayers['max slope'] = {
            'grid': GDALGrid(slopestack[:, :, -1], shakemap.getGeoDict()),
            'label': 'Maximum slope ($^\circ$)',
            'type': 'input',
            'description': {
                'units': 'degrees',
                'name': slopesref,
                'longref': slopelref}
        }
        maplayers['cohesion'] = {
            'grid': GDALGrid(cohesion[:, :, 0], shakemap.getGeoDict()),
            'label': 'Cohesion (kPa)',
            'type': 'input',
            'description': {
                'units': 'kPa (adjusted)',
                'name': cohesionsref,
                'longref': cohesionlref}
        }
        maplayers['friction angle'] = {
            'grid': GDALGrid(friction[:, :, 0], shakemap.getGeoDict()),
            'label': 'Friction angle ($^\circ$)',
            'type': 'input',
            'description': {
                'units': 'degrees',
                'name': frictionsref,
                'longref': frictionlref}
        }
        if uncertfile is not None:
            maplayers['pgamin'] = {
                'grid': GDALGrid(PGAmin[:, :, 0], shakemap.getGeoDict()),
                'label': 'PGA - %1.2fstd (g)' % numstd,
                'type': 'input',
                'description': {
                    'units': 'g',
                    'shakemap': shakedetail}
            }
            maplayers['pgamax'] = {
                'grid': GDALGrid(PGAmax[:, :, 0], shakemap.getGeoDict()),
                'label': 'PGA + %1.2fstd (g)' % numstd,
                'type': 'input',
                'description': {
                    'units': 'g',
                    'shakemap': shakedetail}
            }
        if 'PGV' in displmodel:
            if uncertfile is not None:
                maplayers['pgvmin'] = {
                    'grid': GDALGrid(PGVmin[:, :, 0], shakemap.getGeoDict()),
                    'label': 'PGV - %1.2fstd (cm/s)' % numstd,
                    'type': 'input',
                    'description': {
                        'units': 'cm/s',
                        'shakemap': shakedetail}
                }
                maplayers['pgvmax'] = {
                    'grid': GDALGrid(PGVmax[:, :, 0], shakemap.getGeoDict()),
                    'label': 'PGV + %1.2fstd (cm/s)' % numstd,
                    'type': 'input',
                    'description': {
                        'units': 'cm/s',
                        'shakemap': shakedetail}
                }

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
        Tuple of Dn, logDnstd, and logtype.
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
        Dn = np.array(10.**(C1 + np.log10(((1-Ac/PGA)**C2)*(Ac/PGA)**C3)))
        Dn[np.isnan(Dn)] = 0.
        logDnstd = np.ones(np.shape(Dn))*0.51
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
            Dn = np.array(10.**(C1 + np.log10(((1-Ac/PGA)**C2)*(Ac/PGA)**C3) +
                                C4*M))
            Dn[np.isnan(Dn)] = 0.
            logDnstd = np.ones(np.shape(Dn))*0.454
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
        Dn = np.array(np.exp(C1 + C2*(Ac/PGA) + C3*(Ac/PGA)**2 +
                             C4*(Ac/PGA)**3 + C5*(Ac/PGA)**4 +
                             C6*np.log(PGA)+C7*(M-6)))
        Dn[np.isnan(Dn)] = 0.
        logDnstd = 0.732 + 0.789*(Ac/PGA) - 0.539*(Ac/PGA)**2
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
        Dn = np.array(np.exp(C1 + C2*(Ac/PGA) + C3*(Ac/PGA)**2 +
                             C4*(Ac/PGA)**3 + C5*(Ac/PGA)**4 +
                             C6*np.log(PGA)+C7*np.log(PGV)))
        Dn[np.isnan(Dn)] = 0.
        logDnstd = 0.405 + 0.524*(Ac/PGA)
        logtype = 'ln'

    elif model == 'BT_PGA_M':
        if M is None:
            raise Exception('You must enter a value for M to use the '
                            'BT_PGA_M model')
        Dn = np.array(
            np.exp(-0.22 - 2.83*np.log(Ac) - 0.333*(np.log(Ac))**2 +
                   0.566*np.log(Ac)*np.log(PGA) +
                   3.04*np.log(PGA) - 0.244*(np.log(PGA))**2 + 0.278*(M-7.)))
        Dn[np.isnan(Dn)] = 0.
        logDnstd = np.ones(np.shape(Dn))*0.66
        logtype = 'log10'

    if flag == 1:
        Dn = float(Dn)
        logDnstd = float(logDnstd)

    return Dn, logDnstd, logtype
