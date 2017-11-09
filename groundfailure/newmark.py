#!/usr/bin/env python
"""
This module contains functions that can be used to run Newmark-based mechanistic landslide models
"""

#stdlib imports
import os.path
import warnings
import collections

#turn off all warnings...
warnings.filterwarnings('ignore')

#local imports
from mapio.shake import ShakeGrid
from mapio.gdal import GDALGrid
from mapio.geodict import GeoDict

#third party imports
import numpy as np


def hazus(shakefile, config, uncertfile=None, saveinputs=False, modeltype=None, displmodel=None,
          probtype=None, bounds=None):

    """This function runs the HAZUS landslide procedure (FEMA, 2003, Chapter 4) using susceptiblity categories (I-X)
    defined by the HAZUS manual.

    TODO add numstd
    :param shakefile: complete file path to the location of the Shakemap to use as input
    :type shakefile: string:
    :param config: Model configuration file object containing locations of input files and other input values config =
      ConfigObj(configfilepath)
    :type config: ConfigObj
    :param saveinputs: Whether or not to return the model input layers, False (defeault) returns only the model output
      (one layer)
    :type saveinputs: boolean
    :param modeltype: OVERWRITES VALUE IN CONFIG FILE IF SPECIFIED
        * 'coverage' (default) if critical acceleration is exceeded by pga, this gives the  estimated areal coverage of landsliding for that cell.
        * 'dn_hazus' - Outputs Newmark displacement using HAZUS methods without relating to probability of failure.
        * 'dn_prob' - Estimates Newmark displacement using HAZUS methods and relates to probability of failure using param probtype.
        * 'ac_classic_dn' - Uses the critical acceleration defined by HAZUS methodology and uses regression model defined by displmodel param to get Newmark displacement without relating to probability of failure.
        * 'ac_classic_prob' - Uses the critical acceleration defined by HAZUS methodology and uses regression model defined by displmodel param to get Newmark displacement and probability defined by probtype method.
    :type modeltype: string
    :param displmodel: Newmark displacement regression model to use OVERWRITES VALUE IN CONFIG FILE IF SPECIFIED\n
        * 'J_PGA' (default) - PGA-based model from Jibson (2007) - equation 6.
        * 'J_PGA_M' - PGA and M-based model from Jibson (2007) - equation 7.
        * 'RS_PGA_M' - PGA and M-based model from from Rathje and Saygili (2009).
        * 'RS_PGA_PGV' - PGA and PGV-based model from Saygili and Rathje (2008) - equation 6.
    :type displmodel: string
    :param probtype: Method used to estimate probability. OVERWRITES VALUE IN CONFIG FILE IF SPECIFIED\n
        * 'jibson2000' (default) uses equation 5 from Jibson et al. (2000) to estimate probability from Newmark displacement.
        * 'threshold' uses a specified threshold of Newmark displacement (defined in config file) and assumes anything greater than this threshold fails
    :type probtype: string
    :param bounds: Boundaries to compute over if different from ShakeMap boundaries as dictionary with keys 'xmin',
      'xmax', 'ymin', 'ymax'
    :type bounds: dictionary
    :returns:
        maplayers(OrderedDict): Dictionary containing output and input layers (if saveinputs=True) along with metadata formatted like:

        maplayers['layer name']={'grid': mapio grid2D object, 'label': 'label for colorbar and top line of subtitle',
        'type': 'output or input to model', 'description': 'detailed description of layer for subtitle,
        potentially including source information'}

    """

    # Empty refs
    suslref = 'unknown'
    sussref = 'unknown'
    modellref = 'unknown'
    modelsref = 'unknown'

    # Parse config and read in files
    sus = None
    susdat = None

    # Read in susceptiblity file
    #try:
    susfile = config['hazus']['layers']['susceptibility']['file']
    shkgdict = ShakeGrid.getFileGeoDict(shakefile, adjust='res')
    susdict, first_column_duplicated = GDALGrid.getFileGeoDict(susfile)
    if bounds is not None:  # Make sure bounds are within ShakeMap Grid
        if shkgdict.xmin > bounds['xmin'] or shkgdict.xmax < bounds['xmax'] or shkgdict.ymin > bounds['ymin'] or\
           shkgdict.ymax < bounds['ymax']:
            print('Specified bounds are outside shakemap area, using ShakeMap bounds instead')
            bounds = None
    if bounds is not None:
        tempgdict1 = GeoDict({'xmin': bounds['xmin'], 'ymin': bounds['ymin'], 'xmax': bounds['xmax'],
                             'ymax': bounds['ymax'], 'dx': 100., 'dy': 100., 'nx': 100., 'ny': 100.}, adjust='res')
        tempgdict = susdict.getBoundsWithin(tempgdict1)
    else:
        tempgdict = susdict.getBoundsWithin(shkgdict)
    #sus = GDALGrid.load(susfile, samplegeodict=tempgdict, resample=False)
    sus = GDALGrid.load(susfile, samplegeodict=tempgdict, resample=True, method='linear')
    gdict = sus.getGeoDict()
    susdat = sus.getData()
    #except Exception as e:
    #    raise IOError('Unable to read in susceptibility category file specified in config, %s,' % e)
    #    return

    try:  # Try to fetch source information from config
        modelsref = config['hazus']['shortref']
        modellref = config['hazus']['longref']
        sussref = config['hazus']['layers']['susceptibility']['shortref']
        suslref = config['hazus']['layers']['susceptibility']['longref']
    except:
        print('Was not able to retrieve all references from config file. Continuing')

    try:
        dnthresh = float(config['hazus']['parameters']['dnthresh'])
    except:
        if probtype == 'threshold':
            dnthresh = 5.
            print('Unable to find dnthresh in config, using 5cm')

    if modeltype is None:
        try:
            modeltype = config['hazus']['parameters']['modeltype']
        except:
            print('No modeltype specified, using default of coverage')
            modeltype = 'coverage'
    if displmodel is None:
        try:
            displmodel = config['hazus']['parameters']['displmodel']
        except:
            print('No regression model specified, using default of J_PGA_M')
            displmodel = 'J_PGA_M'
    if probtype is None:
        try:
            probtype = config['hazus']['parameters']['probtype']
        except:
            print('No probability type (probtype) specified, using default of jibson2000')
            probtype = 'jibson2000'

    # Load in shakemap, resample to susceptibility file
    shakemap = ShakeGrid.load(shakefile, adjust='res')
    if uncertfile is not None:
        try:
            uncert = ShakeGrid.load(uncertfile, adjust='res')
        except:
            print('Could not read uncertainty file, ignoring uncertainties')
            uncertfile = None
    PGA = shakemap.getLayer('pga').subdivide(gdict).getData().astype(float)/100.  # in units of g
    PGV = shakemap.getLayer('pgv').subdivide(gdict).getData().astype(float)  # cm/sec
    if uncertfile is not None:
        stdpga = uncert.getLayer('stdpga').subdivide(gdict).getData().astype(float)
        stdpgv = uncert.getLayer('stdpgv').subdivide(gdict).getData().astype(float)
        # estimate PGA +- 1std
        PGAmin = np.exp(np.log(PGA*100.) - stdpga)/100.
        PGAmax = np.exp(np.log(PGA*100.) + stdpga)/100.
        PGVmin = np.exp(np.log(PGV) - stdpgv)
        PGVmax = np.exp(np.log(PGV) + stdpgv)
    M = shakemap.getEventDict()['magnitude']

    # Get critical accelerations in g
    Ac = np.empty(np.shape(susdat))
    Ac[(susdat < 1) & (susdat > 10)] = 9999.
    Ac[susdat == 1] = 0.6
    Ac[susdat == 2] = 0.5
    Ac[susdat == 3] = 0.4
    Ac[susdat == 4] = 0.35
    Ac[susdat == 5] = 0.3
    Ac[susdat == 6] = 0.25
    Ac[susdat == 7] = 0.2
    Ac[susdat == 8] = 0.15
    Ac[susdat == 9] = 0.1
    Ac[susdat == 10] = 0.05

    # can delete sus and susdat now, if don't need to output it, to free up memory
    if saveinputs is False:
        del susdat, sus

    if modeltype == 'coverage':
        areal = np.zeros(np.shape(PGA))
        # This seems to be slow for large matrices
        areal[(PGA >= Ac) & (Ac == 0.6)] = 0.01
        areal[(PGA >= Ac) & (Ac == 0.5)] = 0.02
        areal[(PGA >= Ac) & (Ac == 0.4)] = 0.03
        areal[(PGA >= Ac) & (Ac == 0.35)] = 0.05
        areal[(PGA >= Ac) & (Ac == 0.3)] = 0.08
        areal[(PGA >= Ac) & (Ac == 0.25)] = 0.1
        areal[(PGA >= Ac) & (Ac == 0.2)] = 0.15
        areal[(PGA >= Ac) & (Ac == 0.15)] = 0.2
        areal[(PGA >= Ac) & (Ac == 0.1)] = 0.25
        areal[(PGA >= Ac) & (Ac == 0.05)] = 0.3
        if uncertfile is not None:
            # areal coverage minimum
            arealmin = np.zeros(np.shape(PGA))
            # This seems to be slow for large matrices
            arealmin[(PGAmin >= Ac) & (Ac == 0.6)] = 0.01
            arealmin[(PGAmin >= Ac) & (Ac == 0.5)] = 0.02
            arealmin[(PGAmin >= Ac) & (Ac == 0.4)] = 0.03
            arealmin[(PGAmin >= Ac) & (Ac == 0.35)] = 0.05
            arealmin[(PGAmin >= Ac) & (Ac == 0.3)] = 0.08
            arealmin[(PGAmin >= Ac) & (Ac == 0.25)] = 0.1
            arealmin[(PGAmin >= Ac) & (Ac == 0.2)] = 0.15
            arealmin[(PGAmin >= Ac) & (Ac == 0.15)] = 0.2
            arealmin[(PGAmin >= Ac) & (Ac == 0.1)] = 0.25
            arealmin[(PGAmin >= Ac) & (Ac == 0.05)] = 0.3
            # areal coverage maximum
            arealmax = np.zeros(np.shape(PGA))
            # This seems to be slow for large matrices
            arealmax[(PGAmax >= Ac) & (Ac == 0.6)] = 0.01
            arealmax[(PGAmax >= Ac) & (Ac == 0.5)] = 0.02
            arealmax[(PGAmax >= Ac) & (Ac == 0.4)] = 0.03
            arealmax[(PGAmax >= Ac) & (Ac == 0.35)] = 0.05
            arealmax[(PGAmax >= Ac) & (Ac == 0.3)] = 0.08
            arealmax[(PGAmax >= Ac) & (Ac == 0.25)] = 0.1
            arealmax[(PGAmax >= Ac) & (Ac == 0.2)] = 0.15
            arealmax[(PGAmax >= Ac) & (Ac == 0.15)] = 0.2
            arealmax[(PGAmax >= Ac) & (Ac == 0.1)] = 0.25
            arealmax[(PGAmax >= Ac) & (Ac == 0.05)] = 0.3
    elif modeltype == 'dn_hazus' or modeltype == 'dn_prob':
        ed_low, ed_high = est_disp(Ac, PGA)
        ed_mean = np.mean((np.dstack((ed_low, ed_high))), axis=2)  # Get mean estimated displacements
        Dn = ed_mean * numcycles(M) * PGA
        if uncertfile is not None:
            Dnmin = ed_mean * numcycles(M) * PGAmin
            Dnmax = ed_mean * numcycles(M) * PGAmax
    else:  # Calculate newmark displacement using a regression model
        Dn, logDnstd, logtype = NMdisp(Ac, PGA, model=displmodel, M=M, PGV=PGV)
        if uncertfile is not None:
            Dnmin, logDnstdmin, logtype = NMdisp(Ac, PGAmin, model=displmodel, M=M, PGV=PGVmin)
            Dnmax, logDnstdmax, logtype = NMdisp(Ac, PGAmax, model=displmodel, M=M, PGV=PGVmax)

    # Calculate probability from Dn, if necessary for selected model
    if modeltype == 'ac_classic_prob' or modeltype == 'dn_prob':
        if probtype.lower() in 'jibson2000':
            PROB = 0.335*(1-np.exp(-0.048*Dn**1.565))
            dnthresh = None
            if uncertfile is not None:
                PROBmin = 0.335*(1-np.exp(-0.048*Dnmin**1.565))
                PROBmax = 0.335*(1-np.exp(-0.048*Dnmax**1.565))
        elif probtype.lower() in 'threshold':
            PROB = Dn.copy()
            PROB[PROB <= dnthresh] = 0
            PROB[PROB > dnthresh] = 1
            units = 'prediction'
            label = 'Predicted Landslides'
            if uncertfile is not None:
                PROBmin = Dnmin.copy()
                PROBmin[PROBmin <= dnthresh] = 0
                PROBmin[PROBmin > dnthresh] = 1
                PROBmax = Dnmax.copy()
                PROBmax[PROBmax <= dnthresh] = 0
                PROBmax[PROBmax > dnthresh] = 1
        else:
            print('invalid probtype, assuming jibson2000')
            PROB = 0.335*(1-np.exp(-0.048*Dn**1.565))
            dnthresh = None
            if uncertfile is not None:
                PROBmin = 0.335*(1-np.exp(-0.048*Dnmin**1.565))
                PROBmax = 0.335*(1-np.exp(-0.048*Dnmax**1.565))

    # Turn output and inputs into into grids and put in maplayers dictionary
    maplayers = collections.OrderedDict()

    temp = shakemap.getShakeDict()
    shakedetail = '%s_ver%s' % (temp['shakemap_id'], temp['shakemap_version'])

    if modeltype == 'coverage':
        maplayers['model'] = {'grid': GDALGrid(areal, gdict), 'label': 'Areal coverage', 'type': 'output',
                              'description': {'name': modelsref, 'longref': modellref, 'units': 'coverage',
                                              'shakemap': shakedetail, 'parameters': {'modeltype': 'Landslide - ' + modeltype}}}
        if uncertfile is not None:
            maplayers['modelmin'] = {'grid': GDALGrid(arealmin, gdict), 'label': 'Areal coverage', 'type': 'output',
                                     'description': {'name': modelsref, 'longref': modellref, 'units': 'coverage',
                                                     'shakemap': shakedetail, 'parameters': {'modeltype': 'Landslide - ' + modeltype}}}
            maplayers['modelmax'] = {'grid': GDALGrid(arealmax, gdict), 'label': 'Areal coverage', 'type': 'output',
                                     'description': {'name': modelsref, 'longref': modellref, 'units': 'coverage',
                                                     'shakemap': shakedetail, 'parameters': {'modeltype': 'Landslide - ' + modeltype}}}
    elif modeltype == 'dn_hazus':
        maplayers['model'] = {'grid': GDALGrid(Dn, gdict), 'label': 'Dn (cm)', 'type': 'output',
                              'description': {'name': modelsref, 'longref': modellref, 'units': 'displacement',
                                              'shakemap': shakedetail,
                                              'parameters': {'displmodel': displmodel, 'modeltype': 'Landslide - ' + modeltype}}}
        if uncertfile is not None:
            maplayers['modelmin'] = {'grid': GDALGrid(Dnmin, gdict), 'label': 'Dn (cm)', 'type': 'output',
                                     'description': {'name': modelsref, 'longref': modellref, 'units': 'displacement',
                                                     'shakemap': shakedetail, 'parameters': {'displmodel': displmodel,
                                                                                             'modeltype': 'Landslide - ' + modeltype}}}
            maplayers['modelmax'] = {'grid': GDALGrid(Dnmax, gdict), 'label': 'Dn (cm)', 'type': 'output',
                                     'description': {'name': modelsref, 'longref': modellref, 'units': 'displacement',
                                                     'shakemap': shakedetail,
                                                     'parameters': {'displmodel': displmodel, 'modeltype': 'Landslide - ' + modeltype}}}
    elif modeltype == 'dn_prob':
        maplayers['model'] = {'grid': GDALGrid(PROB, gdict), 'label': 'Landslide Probability', 'type': 'output',
                              'description': {'name': modelsref, 'longref': modellref, 'units': 'probability',
                                              'shakemap': shakedetail, 'parameters': {'displmodel': displmodel,
                                                                                      'dnthresh_cm': dnthresh,
                                                                                      'modeltype': 'Landslide - ' + modeltype,
                                                                                      'probtype': probtype}}}
    elif modeltype == 'ac_classic_prob':
        maplayers['model'] = {'grid': GDALGrid(PROB, gdict), 'label': 'Landslide Probability', 'type': 'output',
                              'description': {'name': modelsref, 'longref': modellref, 'units': 'probability',
                                              'shakemap': shakedetail, 'parameters': {'displmodel': displmodel,
                                                                                      'dnthresh_cm': dnthresh,
                                                                                      'modeltype': 'Landslide - ' + modeltype,
                                                                                      'probtype': probtype}}}
    elif modeltype == 'ac_classic_dn':
        maplayers['model'] = {'grid': GDALGrid(Dn, gdict), 'label': 'Landslide Probability', 'type': 'output',
                              'description': {'name': modelsref, 'longref': modellref, 'units': 'probability',
                                              'shakemap': shakedetail, 'parameters': {'displmodel': displmodel,
                                                                                      'dnthresh_cm': dnthresh,
                                                                                      'modeltype': 'Landslide - ' + modeltype,
                                                                                      'probtype': probtype}}}
        if uncertfile is not None:
            maplayers['modelmin'] = {'grid': GDALGrid(Dnmin, gdict), 'label': 'Dn (cm)', 'type': 'output',
                                     'description': {'name': modelsref, 'longref': modellref, 'units': 'displacement',
                                                     'shakemap': shakedetail, 'parameters': {'displmodel': displmodel,
                                                                                             'modeltype': 'Landslide - ' + modeltype}}}
            maplayers['modelmax'] = {'grid': GDALGrid(Dnmax, gdict), 'label': 'Dn (cm)', 'type': 'output',
                                     'description': {'name': modelsref, 'longref': modellref, 'units': 'displacement',
                                                     'shakemap': shakedetail,
                                                     'parameters': {'displmodel': displmodel, 'modeltype': 'Landslide - ' + modeltype}}}

    label = 'Probability'
    if modeltype != 'coverage' and modeltype != 'dn_hazus' and modeltype != 'ac_classic_dn':
        if uncertfile is not None:
            maplayers['modelmin'] = {'grid': GDALGrid(PROBmin, gdict), 'label': label+' -1std', 'type': 'output',
                                     'description': {}}
            maplayers['modelmax'] = {'grid': GDALGrid(PROBmax, gdict), 'label': label+' +1std', 'type': 'output',
                                     'description': {}}

    if saveinputs is True:
        maplayers['suscat'] = {'grid': sus, 'label': 'Susceptibility Category', 'type': 'input',
                               'description': {'name': sussref, 'longref': suslref, 'units': 'Category'}}
        maplayers['Ac'] = {'grid': GDALGrid(Ac, gdict), 'label': 'Ac (g)', 'type': 'output',
                           'description': {'units': 'g', 'shakemap': shakedetail}}
        maplayers['pga'] = {'grid': GDALGrid(PGA, gdict), 'label': 'PGA (g)', 'type': 'input',
                            'description': {'units': 'g', 'shakemap': shakedetail}}
        if 'pgv' in displmodel.lower():
            maplayers['pgv'] = {'grid': GDALGrid(PGV, gdict), 'label': 'PGV (cm/s)', 'type': 'input',
                                'description': {'units': 'cm/s', 'shakemap': shakedetail}}
        if 'dn' not in modeltype.lower() and modeltype != 'coverage':
            maplayers['dn'] = {'grid': GDALGrid(Dn, gdict), 'label': 'Dn (cm)', 'type': 'output',
                               'description': {'units': 'displacement', 'shakemap': shakedetail,
                                               'parameters': {'displmodel': displmodel,
                                                              'modeltype': 'Landslide - ' + modeltype}}}

    return maplayers


def est_disp(Ac, PGA):
    """est_disp retrieved the estimated displacement factor from HAZUS. This was digitized from figure in HAZUS manual
    pg. 4-37, which is based on Makdisi and Seed (1978) according to HAZUS, but the equations don't appear in that
    document. It also assumes that PGA is equal to the ais (induced acceleration) and any ac/ais ratios lower than 0.1,
    which is the minimum shown on the graph, have estimated displacements of 20 and 40 cm for low and high, respectively.
    This is based solely on a horizontal projection out from the end of the lines shown on Fig 4.14.

    :param Ac: Critical acceleration in the same units as PGA (this is ratio based)
    :type Ac: numpy array
    :param PGA: Peak ground acceleration in the same units as Ac
    :type PGA: numpy array
    :returns:
        * ed_low: low estimate of expected displacement factor
        * ed_high: high estimate of expected displacement factor

    """
    from scipy.interpolate import interp1d

    acais = Ac/PGA  # Get "Expected displacement factor"
    xlow = np.array((0.0994941, 0.198426, 0.278246, 0.370433, 0.450253, 0.53457, 0.594154, 0.640247, 0.684092, 0.72344,
                    0.770658, 0.802136, 0.851602, 0.897695))
    ylow = np.array((22.0186, 11.4578, 7.09682, 3.85734, 2.28738, 1.24326, 0.822041, 0.543531, 0.367292, 0.237622,
                    0.137873, 0.101646, 0.044438, 0.0211954))
    xhigh = np.array((0.100618, 0.179314, 0.251265, 0.319843, 0.406408, 0.48398, 0.567173, 0.634626, 0.700956,
                     0.751546, 0.80326, 0.857223, 0.901068))
    yhigh = np.array((39.6377, 21.5443, 13.638, 9.21593, 4.90126, 2.84381, 1.6145, 1.02201, 0.592993, 0.351641,
                      0.208521, 0.0931681, 0.0495491))
    f_low = interp1d(xlow, ylow, bounds_error=False, fill_value=0.)
    f_high = interp1d(xhigh, yhigh, bounds_error=False, fill_value=0.)
    ed_low = f_low(acais)
    ed_high = f_high(acais)

    # Fix sections with ac/ais ratios < 0.1
    ed_low[acais < 0.1] = 20.
    ed_high[acais < 0.1] = 40.
    return ed_low, ed_high


def numcycles(M):
    """
    Estimate the number of earthquake cycles using Seed and Idriss (1982) relationship as presented in HAZUS manual, pg. 4-35 and Fig. 4.13

    :param M: earthquake magnitude
    :type M: float
    :returns:
        n(float): number of cycles

    """
    n = 0.3419*M**3 - 5.5214*M**2 + 33.6154*M - 70.7692
    return n


def classic(shakefile, config, uncertfile=None, saveinputs=False, displmodel=None, probtype=None,
            slopediv=1., codiv=1., bounds=None):
    """This function uses the Newmark method to estimate probability of failure at each grid cell.
    Factor of Safety and critcal accelerations are calculated following Jibson et al. (2000) and the
    Newmark displacement is estimated using PGA, PGV, and/or Magnitude (depending on equation used)
    from Shakemap with regression equations from Jibson (2007), Rathje and Saygili (2008) and
    Saygili and Rathje (2009)

    TODO add numstd
    :param shakefile: complete file path to the location of the Shakemap to use as input
    :type shakefile: string:
    :param config: Model configuration file object containing locations of input files and other input values
        config = ConfigObj(configfilepath)
    :type config: ConfigObj
    :param uncertfile: complete file path to the location of the uncertainty.xml for the shakefile,
        if this is not None, it will compute the model for +-std in addition to the best estimate
    :param saveinputs: Whether or not to return the model input layers, False (defeault) returns only the model output (one layer)
    :type saveinputs: boolean
    :param displmodel: Newmark displacement regression model to use OVERWRITES CONFIG FILE CHOICES\n
        * 'J_PGA' (default) - PGA-based model from Jibson (2007) - equation 6.
        * 'J_PGA_M' - PGA and M-based model from Jibson (2007) - equation 7.
        * 'RS_PGA_M' - PGA and M-based model from from Rathje and Saygili (2009).
        * 'RS_PGA_PGV' - PGA and PGV-based model from Saygili and Rathje (2008) - equation 6.
    :type displmodel: string
    :param probtype: Method used to estimate probability. OVERWRITES CONFIG FILE CHOICES\n
        * 'jibson2000' (default) uses equation 5 from Jibson et al. (2000) to estimate probability from Newmark displacement.
        * 'threshold' uses a specified threshold of Newmark displacement (defined in config file) and assumes
            anything greather than this threshold fails
    :type probtype: string
    :param slopediv: Divide slope by this number to get slope in degrees (Verdin datasets need to be divided by 100)
    :type slopediv: float
    :param codiv: Divide cohesion by this number to get reasonable numbers (For Godt method, need to divide by 10
        because it was calibrated downwards to adjust for lower resolution slope estimates)
    :type codiv: float

    :returns:
        maplayers(OrderedDict): Dictionary containing output and input layers (if saveinputs=True) along with metadata
        formatted like:\n

        maplayers['layer name']={'grid': mapio grid2D object, 'label': 'label for colorbar and top line
        of subtitle', 'type': 'output or input to model', 'description': 'detailed description of layer for subtitle,
        potentially including source information'}

    :raises:
        NameError: when unable to parse the config correctly (probably a formatting issue in the configfile) or when
        unable to find the shakefile (Shakemap filepath) - these cause program to end

        NameError: when probtype does not match a predifined probability type, will cause to default to 'jibson2000'

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
    try:
        slopefile = config['classic_newmark']['layers']['slope']['file']
        slopeunits = config['classic_newmark']['layers']['slope']['units']
        cohesionfile = config['classic_newmark']['layers']['cohesion']['file']
        cohesionunits = config['classic_newmark']['layers']['cohesion']['units']
        frictionfile = config['classic_newmark']['layers']['friction']['file']
        frictionunits = config['classic_newmark']['layers']['friction']['units']
        if displmodel is None:
            try:
                displmodel = config['classic_newmark']['parameters']['displmodel']
            except:
                print('No regression model specified, using default of J_PGA_M')
                displmodel = 'J_PGA_M'
        if probtype is None:
            try:
                probtype = config['classic_newmark']['parameters']['probtype']
            except:
                print('No probability type (probtype) specified, using default of jibson2000')
                probtype = 'jibson2000'

        thick = float(config['classic_newmark']['parameters']['thick'])
        uwt = float(config['classic_newmark']['parameters']['uwt'])
        try:
            uwtw = float(config['classic_newmark']['parameters']['uwtw'])
        except:
            print('Could not read soil wet unit weight, using 18.8 kN/m3')
            uwtw = 18.8
        nodata_cohesion = float(config['classic_newmark']['parameters']['nodata_cohesion'])
        nodata_friction = float(config['classic_newmark']['parameters']['nodata_friction'])
        try:
            dnthresh = float(config['classic_newmark']['parameters']['dnthresh'])
        except:
            if probtype == 'threshold':
                dnthresh = 5.
                print('Unable to find dnthresh in config, using 5cm')
            else:
                dnthresh = None
        fsthresh = float(config['classic_newmark']['parameters']['fsthresh'])
        acthresh = float(config['classic_newmark']['parameters']['acthresh'])
        slopethresh = float(config['classic_newmark']['parameters']['slopethresh'])
        if config['classic_newmark']['parameters']['m'] == 'file':
            wtfile = 1
        else:
            wtfile = 0
            try:
                m = float(config['classic_newmark']['parameters']['m'])
            except:
                print('no constant saturated thickness specified, setting m=0')
                m = 0.
    except Exception as e:
        raise NameError('Could not parse configfile, %s' % e)
        return

    try:  # Try to fetch source information from config
        modelsref = config['classic_newmark']['shortref']
        modellref = config['classic_newmark']['longref']
        slopesref = config['classic_newmark']['layers']['slope']['shortref']
        slopelref = config['classic_newmark']['layers']['slope']['longref']
        cohesionsref = config['classic_newmark']['layers']['cohesion']['shortref']
        cohesionlref = config['classic_newmark']['layers']['cohesion']['longref']
        frictionsref = config['classic_newmark']['layers']['friction']['shortref']
        frictionlref = config['classic_newmark']['layers']['friction']['longref']
    except:
        print('Was not able to retrieve all references from config file. Continuing')

    # Cut and resample all files
    shkgdict = ShakeGrid.getFileGeoDict(shakefile, adjust='res')
    slpdict, first_column_duplicated = GDALGrid.getFileGeoDict(slopefile)
    if bounds is not None:  # Make sure bounds are within ShakeMap Grid
        if shkgdict.xmin > bounds['xmin'] or shkgdict.xmax < bounds['xmax'] or shkgdict.ymin > bounds['ymin'] or\
           shkgdict.ymax < bounds['ymax']:
            print('Specified bounds are outside shakemap area, using ShakeMap bounds instead')
            bounds = None
    if bounds is not None:
        tempgdict = GeoDict({'xmin': bounds['xmin'], 'ymin': bounds['ymin'], 'xmax': bounds['xmax'],
                            'ymax': bounds['ymax'], 'dx': 100., 'dy': 100., 'nx': 100., 'ny': 100.}, adjust='res')
        gdict = slpdict.getBoundsWithin(tempgdict)
    else:  # Get boundaries from shakemap if not specified
        shkgdict = ShakeGrid.getFileGeoDict(shakefile, adjust='res')
        slpdict, first_column_duplicated = GDALGrid.getFileGeoDict(slopefile)
        gdict = slpdict.getBoundsWithin(shkgdict)

    # Load in slope file
    slopegrid = GDALGrid.load(slopefile, samplegeodict=gdict, resample=True, method='linear')
    #slopegrid = GDALGrid.load(slopefile, samplegeodict=gdict, resample=False)
    gdict = slopegrid.getGeoDict()  # Get this again just in case it changed
    slope = slopegrid.getData().astype(float)/slopediv  # Adjust slope to degrees, if needed
    # Change any zero slopes to a very small number to avoid dividing by zero later
    slope[slope == 0] = 1e-8

    # Load in shakemap, resample to slope file (this will be important when go to higher res)
    shakemap = ShakeGrid.load(shakefile, samplegeodict=gdict, resample=True, method='linear', adjust='res')
    M = shakemap.getEventDict()['magnitude']
    # Read in uncertainty layer, if present
    if uncertfile is not None:
        try:
            uncert = ShakeGrid.load(uncertfile, samplegeodict=gdict, resample=True, method='linear', adjust='res')
        except:
            print('Could not read uncertainty file, ignoring uncertainties')
            uncertfile = None

    # Read in the cohesion and friction files, resampled to slope grid
    cohesion = GDALGrid.load(cohesionfile, samplegeodict=gdict, resample=True, method='nearest').getData().astype(float)/codiv
    cohesion[np.isnan(cohesion)] = nodata_cohesion
    friction = GDALGrid.load(frictionfile, samplegeodict=gdict, resample=True, method='nearest').getData().astype(float)
    friction[np.isnan(friction)] = nodata_friction

    # See if there is a water table depth file and read it in if there is
    if wtfile:
        try:
            waterfile = config['classic_newmark']['layers']['watertable']['file']
            watertable = GDALGrid.load(waterfile, samplegeodict=gdict, resample=True, method='linear').getData()  # Needs to be in meters!
            try:
                watersref = config['classic_newmark']['layers']['watertable']['shortref']
                waterlref = config['classic_newmark']['layers']['watertable']['longref']
            except:
                print('Was not able to retrieve water table references from config file. Continuing')

        except:
            print('Water table file not specified or readable, assuming constant saturated thickness proportion of 0.')
            wtfile = 0

    # Factor of safety
    if wtfile:
        watertable[watertable > thick] = thick
        m = (thick - watertable)/thick
        m[np.isnan(m)] = 0.
    FS = cohesion/(uwt*thick*np.sin(slope*(np.pi/180.))) + np.tan(friction*(np.pi/180.))/np.tan(slope*(np.pi/180.)) - (m*uwtw*np.tan(friction*(np.pi/180.)))/(uwt*np.tan(slope*(np.pi/180.)))
    FS[FS < fsthresh] = fsthresh

    # Compute critical acceleration, in g
    Ac = (FS-1.)*np.sin(slope*(np.pi/180.))  # This gives ac in g, equations that multiply by g give ac in m/s2
    Ac[Ac < acthresh] = acthresh
    Ac[slope < slopethresh] = float('nan')

    # Get PGA in g (PGA is %g in ShakeMap, convert to g)
    PGA = shakemap.getLayer('pga').getData().astype(float)/100.
    PGV = shakemap.getLayer('pgv').getData().astype(float)
    if uncertfile is not None:
        stdpga = uncert.getLayer('stdpga')
        stdpgv = uncert.getLayer('stdpgv')
        # Estimate PGA +- 1std
        PGAmin = np.exp(np.log(PGA*100.) - stdpga.getData())/100.
        PGAmax = np.exp(np.log(PGA*100.) + stdpga.getData())/100.
        PGVmin = np.exp(np.log(PGV) - stdpgv.getData())
        PGVmax = np.exp(np.log(PGV) + stdpgv.getData())

    np.seterr(invalid='ignore')  # Ignore errors so still runs when Ac > PGA, just leaves nan instead of crashing

    Dn, logDnstd, logtype = NMdisp(Ac, PGA, model=displmodel, M=M, PGV=PGV)
    if uncertfile is not None:
        Dnmin, logDnstdmin, logtype = NMdisp(Ac, PGAmin, model=displmodel, M=M, PGV=PGVmin)
        Dnmax, logDnstdmax, logtype = NMdisp(Ac, PGAmax, model=displmodel, M=M, PGV=PGVmax)

    units = 'probability'
    label = 'Landslide Probability'
    if probtype.lower() in 'jibson2000':
        PROB = 0.335*(1-np.exp(-0.048*Dn**1.565))
        dnthresh = None
        if uncertfile is not None:
            PROBmin = 0.335*(1-np.exp(-0.048*Dnmin**1.565))
            PROBmax = 0.335*(1-np.exp(-0.048*Dnmax**1.565))
    elif probtype.lower() in 'threshold':
        PROB = Dn.copy()
        PROB[PROB <= dnthresh] = 0
        PROB[PROB > dnthresh] = 1
        units = 'prediction'
        label = 'Predicted Landslides'
        if uncertfile is not None:
            PROBmin = Dnmin.copy()
            PROBmin[PROBmin <= dnthresh] = 0
            PROBmin[PROBmin > dnthresh] = 1
            PROBmax = Dnmax.copy()
            PROBmax[PROBmax <= dnthresh] = 0
            PROBmax[PROBmax > dnthresh] = 1
    else:
        print('invalid probtype, assuming jibson2000')
        PROB = 0.335*(1-np.exp(-0.048*Dn**1.565))
        dnthresh = None
        if uncertfile is not None:
            PROBmin = 0.335*(1-np.exp(-0.048*Dnmin**1.565))
            PROBmax = 0.335*(1-np.exp(-0.048*Dnmax**1.565))

    # Turn output and inputs into into grids and put in mapLayers dictionary
    maplayers = collections.OrderedDict()

    temp = shakemap.getShakeDict()
    shakedetail = '%s_ver%s' % (temp['shakemap_id'], temp['shakemap_version'])

    if wtfile:
        des = 'variable'
    else:
        des = m
    description = {'name': modelsref, 'longref': modellref, 'units': units, 'shakemap': shakedetail,
                   'parameters': {'displmodel': displmodel, 'thickness_m': thick, 'unitwt_kNm3': uwt,
                                  'dnthresh_cm': dnthresh, 'acthresh_g': acthresh, 'fsthresh': fsthresh,
                                  'slopethresh': slopethresh, 'sat_proportion': des,
                                  'modeltype': 'Landslide'}}

    maplayers['model'] = {'grid': GDALGrid(PROB, gdict), 'label': label, 'type': 'output', 'description': description}
    if uncertfile is not None:
        maplayers['modelmin'] = {'grid': GDALGrid(PROBmin, gdict), 'label': label+' -1std', 'type': 'output',
                                 'description': description}
        maplayers['modelmax'] = {'grid': GDALGrid(PROBmax, gdict), 'label': label+' +1std', 'type': 'output',
                                 'description': description}

    if saveinputs is True:
        maplayers['pga'] = {'grid': GDALGrid(PGA, gdict), 'label': 'PGA (g)', 'type': 'input',
                            'description': {'units': 'g', 'shakemap': shakedetail}}
        maplayers['FS'] = {'grid': GDALGrid(FS, gdict), 'label': 'Factor of Safety', 'type': 'input',
                           'description': {'units': 'unitless'}}
        maplayers['Ac'] = {'grid': GDALGrid(Ac, gdict), 'label': 'Critical acceleration (g)', 'type': 'input'}
        maplayers['Dn'] = {'grid': GDALGrid(Dn, gdict), 'label': 'Newmark Displacement (cm)', 'type': 'input'}
        maplayers['slope'] = {'grid': GDALGrid(slope, gdict), 'label': 'Max slope ($^\circ$)', 'type': 'input',
                              'description': {'units': 'degrees', 'name': slopesref, 'longref': slopelref}}
        maplayers['cohesion'] = {'grid': GDALGrid(cohesion, gdict), 'label': 'Cohesion (kPa)', 'type': 'input',
                                 'description': {'units': 'kPa (adjusted)', 'name': cohesionsref, 'longref': cohesionlref}}
        maplayers['friction angle'] = {'grid': GDALGrid(friction, gdict), 'label': 'Friction angle ($^\circ$)',
                                       'type': 'input', 'description': {'units': 'degrees', 'name': frictionsref,
                                                                        'longref': frictionlref}}
        if uncertfile is not None:
            maplayers['pgamin'] = {'grid': GDALGrid(PGAmin, gdict), 'label': 'PGA - 1std (g)', 'type': 'input',
                                   'description': {'units': 'g', 'shakemap': shakedetail}}
            maplayers['pgamax'] = {'grid': GDALGrid(PGAmax, gdict), 'label': 'PGA + 1std (g)', 'type': 'input',
                                   'description': {'units': 'g', 'shakemap': shakedetail}}
        if 'PGV' in displmodel:
            maplayers['pgv'] = {'grid': GDALGrid(PGV, gdict), 'label': 'PGV (cm/s)', 'type': 'input',
                                'description': {'units': 'cm/s', 'shakemap': shakedetail}}
            if uncertfile is not None:
                maplayers['pgvmin'] = {'grid': GDALGrid(PGVmin, gdict), 'label': 'PGV - 1std (cm/s)', 'type': 'input',
                                       'description': {'units': 'cm/s', 'shakemap': shakedetail}}
                maplayers['pgvmax'] = {'grid': GDALGrid(PGVmax, gdict), 'label': 'PGV + 1std (cm/s)', 'type': 'input',
                                       'description': {'units': 'cm/s', 'shakemap': shakedetail}}
        if wtfile:
            maplayers['sat thick prop'] = {'grid': GDALGrid(m, gdict), 'label': 'Saturated thickness proprtion [0,1]',
                                           'type': 'input', 'description': {'units': 'meters', 'name': watersref,
                                                                            'longref': waterlref}}

    return maplayers


def godt2008(shakefile, config, uncertfile=None, saveinputs=False, displmodel=None, bounds=None, slopediv=100.,
             codiv=10., numstd=None):
    """ This function runs the Godt et al. (2008) global method for a given ShakeMap. The Factor of Safety
    is calculated using infinite slope analysis assumuing dry conditions. The method uses threshold newmark
    displacement and estimates areal coverage by doing the calculations for each slope quantile
    TO DO - add 'all' - averages Dn from all four equations, add term to convert PGA and PGV to Ia and use other
    equations, add Ambraseys and Menu (1988) option

    :param shakefile: filepath to shakemap xml file
    :type shakefile: string
    :param config: ConfigObj of config file containing inputs required for running the model
    :type config: ConfigObj
    :param saveinputs: Whether or not to return the model input layers, False (defeault) returns only the model output (one layer)
    :type saveinputs: boolean
    :param displmodel: Newmark displacement regression model to use\n
        * 'J_PGA' (default) - PGA-based model from Jibson (2007) - equation 6.
        * 'J_PGA_M' - PGA and M-based model from Jibson (2007) - equation 7.
        * 'RS_PGA_M' - PGA and M-based model from from Rathje and Saygili (2009).
        * 'RS_PGA_PGV' - PGA and PGV-based model from Saygili and Rathje (2008) - equation 6.
    :type displmodel: string
    :param probtype: Method used to estimate probability.\n
        * 'jibson2000' uses equation 5 from Jibson et al. (2000) to estimate probability from Newmark displacement.
        * 'threshold' uses a specified threshold of Newmark displacement (defined in config file) and assumes anything greater than this threshold fails
    :type probtype: string
    :param slopediv: Divide slope by this number to get slope in degrees (Verdin datasets need to be divided by 100)
    :type slopediv: float
    :param codiv: Divide cohesion by this number to get reasonable numbers (For Godt method, need to divide by 10
                  because that is how it was calibrated, but values are reasonable without multiplying for regular
                  analysis)
    :type codiv: float
    :param numstd: number of +/- standard deviations to use if uncertainty is computed (uncertfile is not None)
    :returns: maplayers(OrderedDict): Dictionary containing output and input layers (if saveinputs=True) along with metadata formatted like:\n
                maplayers['layer name']={'grid': mapio grid2D object, 'label': 'label for
                colorbar and top line of subtitle', 'type': 'output or input to model',
                'description': 'detailed description of layer for subtitle, potentially including source information'}

    :raises:
         NameError: when unable to parse the config correctly (probably a formatting issue in the configfile) or when
         unable to find the shakefile (Shakemap filepath) - these cause program to end

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
    try:  # May want to add error handling so if refs aren't given, just includes unknown
        slopefilepath = config['godt_2008']['layers']['slope']['filepath']
        slopeunits = config['godt_2008']['layers']['slope']['units']
        cohesionfile = config['godt_2008']['layers']['cohesion']['file']
        cohesionunits = config['godt_2008']['layers']['cohesion']['units']
        frictionfile = config['godt_2008']['layers']['friction']['file']
        frictionunits = config['godt_2008']['layers']['friction']['units']

        thick = float(config['godt_2008']['parameters']['thick'])
        uwt = float(config['godt_2008']['parameters']['uwt'])
        nodata_cohesion = float(config['godt_2008']['parameters']['nodata_cohesion'])
        nodata_friction = float(config['godt_2008']['parameters']['nodata_friction'])
        dnthresh = float(config['godt_2008']['parameters']['dnthresh'])
        fsthresh = float(config['godt_2008']['parameters']['fsthresh'])
        acthresh = float(config['godt_2008']['parameters']['acthresh'])
        try:
            slopemin = float(config['godt_2008']['parameters']['slopemin'])
        except:
            slopemin = 0.01
            print('No slopemin found in config file, using 0.01 deg for slope minimum')
    except Exception as e:
        raise NameError('Could not parse configfile, %s' % e)
        #return

    if displmodel is None:
        try:
            displmodel = config['godt_2008']['parameters']['displmodel']
        except:
            print('No regression model specified, using default of J_PGA_M')
            displmodel = 'J_PGA_M'

    # TO DO, ADD ERROR CATCHING ON UNITS, MAKE SURE THEY ARE WHAT THEY SHOULD BE FOR THIS MODEL

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
        print('Was not able to retrieve all references from config file. Continuing')

    shkgdict = ShakeGrid.getFileGeoDict(shakefile, adjust='res')
    if bounds is not None:  # Make sure bounds are within ShakeMap Grid
        if shkgdict.xmin > bounds['xmin'] or shkgdict.xmax < bounds['xmax'] or shkgdict.ymin > bounds['ymin'] or\
           shkgdict.ymax < bounds['ymax']:
            print('Specified bounds are outside shakemap area, using ShakeMap bounds instead')
            bounds = None
    if bounds is not None:
        tempgdict = GeoDict.createDictFromBox(bounds['xmin'], bounds['xmax'], bounds['ymin'], bounds['ymax'],
                                              shkgdict.dx, shkgdict.dy, inside=False)
        gdict = shkgdict.getBoundsWithin(tempgdict)
        #shakemap = ShakeGrid.load(shakefile, samplegeodict=gdict, adjust='bounds')
        shakemap = ShakeGrid.load(shakefile, samplegeodict=gdict, resample=True, method='linear', adjust='bounds')
    else:
        shakemap = ShakeGrid.load(shakefile, adjust='res')
    shkgdict = shakemap.getGeoDict()  # Get updated geodict
    M = shakemap.getEventDict()['magnitude']

    # read in uncertainty if present
    if uncertfile is not None:
        try:
            uncert = ShakeGrid.load(uncertfile, samplegeodict=shkgdict, resample=True, method='linear', adjust='res')
        except:
            print('Could not read uncertainty file, ignoring uncertainties')
            uncertfile = None
        if numstd is None:
            numstd = 1.

    # Read in all the slope files, divide all by 100 to get to slope in degrees (because input files are multiplied by 100.)
    slopes = []
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope_min.bil'), samplegeodict=shkgdict,
                  resample=True, method='linear').getData()/slopediv)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope10.bil'), samplegeodict=shkgdict,
                  resample=True, method='linear').getData()/slopediv)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope30.bil'), samplegeodict=shkgdict,
                  resample=True, method='linear').getData()/slopediv)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope50.bil'), samplegeodict=shkgdict,
                  resample=True, method='linear').getData()/slopediv)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope70.bil'), samplegeodict=shkgdict,
                  resample=True, method='linear').getData()/slopediv)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope90.bil'), samplegeodict=shkgdict,
                  resample=True, method='linear').getData()/slopediv)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope_max.bil'), samplegeodict=shkgdict,
                  resample=True, method='linear').getData()/slopediv)
    slopestack = np.dstack(slopes)

    # Change any zero slopes to a very small number to avoid dividing by zero later
    slopestack[slopestack == 0] = 1e-8

    # Read in the cohesion and friction files and duplicate layers so they are same shape as slope structure
    cohesion = np.repeat(GDALGrid.load(cohesionfile, samplegeodict=shkgdict, resample=True,
                         method='nearest').getData()[:, :, np.newaxis]/codiv, 7, axis=2)
    cohesion[cohesion == -999.9] = nodata_cohesion
    cohesion = np.nan_to_num(cohesion)
    cohesion[cohesion == 0] = nodata_cohesion
    friction = np.repeat(GDALGrid.load(frictionfile, samplegeodict=shkgdict, resample=True,
                         method='nearest').getData().astype(float)[:, :, np.newaxis], 7, axis=2)
    friction[friction == -9999] = nodata_friction
    friction = np.nan_to_num(friction)
    friction[friction == 0] = nodata_friction

    # Do the calculations using Jibson (2007) PGA only model for Dn
    FS = cohesion/(uwt*thick*np.sin(slopestack*(np.pi/180.))) + np.tan(friction*(np.pi/180.))/np.tan(slopestack*(np.pi/180.))
    FS[FS < fsthresh] = fsthresh

    # Compute critical acceleration, in g
    Ac = (FS-1)*np.sin(slopestack*(np.pi/180.)).astype(float)  # This gives ac in g, equations that multiply by g give ac in m/s2
    Ac[Ac < acthresh] = acthresh

    # Get PGA in g (PGA is %g in ShakeMap, convert to g)
    PGA = np.repeat(shakemap.getLayer('pga').getData()[:, :, np.newaxis]/100., 7, axis=2).astype(float)
    if 'PGV' in displmodel:  # Load in PGV also, in cm/sec
        PGV = np.repeat(shakemap.getLayer('pgv').getData()[:, :, np.newaxis], 7, axis=2).astype(float)
    else:
        PGV = None

    if uncertfile is not None:
        stdpga = np.repeat(uncert.getLayer('stdpga').getData()[:, :, np.newaxis], 7, axis=2).astype(float)
        stdpgv = np.repeat(uncert.getLayer('stdpgv').getData()[:, :, np.newaxis], 7, axis=2).astype(float)
        # estimate PGA +- 1std
        PGAmin = np.exp(np.log(PGA*100) - numstd*stdpga)/100
        PGAmax = np.exp(np.log(PGA*100) + numstd*stdpga)/100
        if 'PGV' in displmodel:
            PGVmin = np.exp(np.log(PGV) - numstd*stdpgv)
            PGVmax = np.exp(np.log(PGV) + numstd*stdpgv)
        else:
            PGVmin = None
            PGVmax = None

    np.seterr(invalid='ignore')  # Ignore errors so still runs when Ac > PGA, just leaves nan instead of crashing

    Dn, logDnstd, logtype = NMdisp(Ac, PGA, model=displmodel, M=M, PGV=PGV)
    if uncertfile is not None:
        Dnmin, logDnstdmin, logtype = NMdisp(Ac, PGAmin, model=displmodel, M=M, PGV=PGVmin)
        Dnmax, logDnstdmax, logtype = NMdisp(Ac, PGAmax, model=displmodel, M=M, PGV=PGVmax)

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
        #uncert too
        if uncertfile is not None:
            PROBmin[slopestack[:, :, 6] <= slopemin] = 0.
            PROBmax[slopestack[:, :, 6] <= slopemin] = 0.

    # Turn output and inputs into into grids and put in mapLayers dictionary
    maplayers = collections.OrderedDict()

    temp = shakemap.getShakeDict()
    shakedetail = '%s_ver%s' % (temp['shakemap_id'], temp['shakemap_version'])

    description = {'name': modelsref, 'longref': modellref, 'units': 'coverage', 'shakemap': shakedetail,
                   'parameters': {'displmodel': displmodel, 'thickness_m': thick, 'unitwt_kNm3': uwt,
                                  'dnthresh_cm': dnthresh, 'acthresh_g': acthresh, 'fsthresh': fsthresh,
                                  'modeltype': 'Landslide'}}

    maplayers['model'] = {'grid': GDALGrid(PROB, shakemap.getGeoDict()), 'label': 'Areal coverage', 'type': 'output',
                          'description': description}
    if uncertfile is not None:
        maplayers['modelmin'] = {'grid': GDALGrid(PROBmin, shkgdict), 'label': 'Probability-%1.2fstd' % numstd, 'type': 'output',
                                 'description': {}}
        maplayers['modelmax'] = {'grid': GDALGrid(PROBmax, shkgdict), 'label': 'Probability+%1.2fstd' % numstd, 'type': 'output',
                                 'description': {}}

    if saveinputs is True:
        maplayers['pga'] = {'grid': GDALGrid(PGA[:, :, 0], shakemap.getGeoDict()), 'label': 'PGA (g)', 'type': 'input',
                            'description': {'units': 'g', 'shakemap': shakedetail}}
        if 'PGV' in displmodel:
            maplayers['pgv'] = {'grid': GDALGrid(PGV[:, :, 0], shakemap.getGeoDict()), 'label': 'PGV (cm/s)',
                                'type': 'input', 'description': {'units': 'cm/s', 'shakemap': shakedetail}}
        maplayers['minFS'] = {'grid': GDALGrid(np.min(FS, axis=2), shakemap.getGeoDict()),
                              'label': 'Min Factor of Safety', 'type': 'input', 'description': {'units': 'unitless'}}
        maplayers['max slope'] = {'grid': GDALGrid(slopestack[:, :, -1], shakemap.getGeoDict()),
                                  'label': 'Maximum slope ($^\circ$)', 'type': 'input',
                                  'description': {'units': 'degrees', 'name': slopesref, 'longref': slopelref}}
        maplayers['cohesion'] = {'grid': GDALGrid(cohesion[:, :, 0], shakemap.getGeoDict()), 'label': 'Cohesion (kPa)',
                                 'type': 'input', 'description': {'units': 'kPa (adjusted)', 'name': cohesionsref,
                                                                  'longref': cohesionlref}}
        maplayers['friction angle'] = {'grid': GDALGrid(friction[:, :, 0], shakemap.getGeoDict()),
                                       'label': 'Friction angle ($^\circ$)', 'type': 'input',
                                       'description': {'units': 'degrees', 'name': frictionsref, 'longref': frictionlref}}
        if uncertfile is not None:
            maplayers['pgamin'] = {'grid': GDALGrid(PGAmin[:, :, 0], shakemap.getGeoDict()), 'label': 'PGA - %1.2fstd (g)' % numstd,
                                   'type': 'input', 'description': {'units': 'g', 'shakemap': shakedetail}}
            maplayers['pgamax'] = {'grid': GDALGrid(PGAmax[:, :, 0], shakemap.getGeoDict()), 'label': 'PGA + %1.2fstd (g)' % numstd,
                                   'type': 'input', 'description': {'units': 'g', 'shakemap': shakedetail}}
        if 'PGV' in displmodel:
            if uncertfile is not None:
                maplayers['pgvmin'] = {'grid': GDALGrid(PGVmin[:, :, 0], shakemap.getGeoDict()), 'label': 'PGV - %1.2fstd (cm/s)' % numstd,
                                       'type': 'input', 'description': {'units': 'cm/s', 'shakemap': shakedetail}}
                maplayers['pgvmax'] = {'grid': GDALGrid(PGVmax[:, :, 0], shakemap.getGeoDict()), 'label': 'PGV + %1.2fstd (cm/s)' % numstd,
                                       'type': 'input', 'description': {'units': 'cm/s', 'shakemap': shakedetail}}

    return maplayers


def NMdisp(Ac, PGA, model='J_PGA', M=None, PGV=None):
    """
    PGA-based Newmark Displacement model

    :param Ac: NxM Array of critical accelerations in units of g
    :type Ac: numpy Array
    :param PGA: NxM Array of PGA values in units of g
    :type PGA: numpy Array

    :param model:\n
      *'J_PGA' - PGA only model from Jibson (2007), equation 6 Applicable for Magnitude range of dataset (5.3-7.6)
      *'J_PGA_M' - PGA-and M- based Newmark Displacement model from Jibson (2007), equation 7 Applicable for Magnitude range of dataset (5.3-7.6)
      *'RS_PGA_M' - PGA and M-based Newmark displacement model from Rathje and Saygili (2009)
      *'RS_PGA_PGV' - PGA and PGV-based model from Saygili and Rathje (2008) - equation 6
      *'BT_PGA_M' - PGA and M-based model from Bray and Travasarou, 2007 assuming natural fundamental
        period of sliding mass Ts = 0 (equation 6)
    :type model: string

    :param M: Magnitude - only needed for models with M in the name
    :type M: float
    :param PGV: NxM Array of PGV values in units of cm/sec - only needed for models with PGV in the name
    :type PGV: numpy array or float

    :returns:
        Dn(array): NxM array of Newmark displacements in cm\n
        logDnstd(array): NxM array of sigma Dn in log units
    """
    # Deal with non-array inputs
    if isinstance(Ac, float) or isinstance(Ac, int):
        flag = 1
    else:
        flag = 0
    np.seterr(invalid='ignore')  # Ignore errors so still runs when Ac > PGA, just leaves nan instead of crashing

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
            raise Exception('M (magnitude) not found, cannot use RS_PGA_M model')
        else:
            C1 = -2.71  # additive constant in newmark displacement calculation
            C2 = 2.335  # first exponential constant
            C3 = -1.478  # second exponential constant
            C4 = 0.424
            #Dn = np.exp(C1 + np.log(((1-Ac/PGA)**C2)*(Ac/PGA)**C3) + C4*M)
            Dn = np.array(10.**(C1 + np.log10(((1-Ac/PGA)**C2)*(Ac/PGA)**C3) + C4*M))
            Dn[np.isnan(Dn)] = 0.
            logDnstd = np.ones(np.shape(Dn))*0.454
            logtype = 'log10'

    elif model == 'RS_PGA_M':
        if M is None:
            raise Exception('You must enter a value for M to use the RS_PGA_M model')
        C1 = 4.89
        C2 = -4.85
        C3 = -19.64
        C4 = 42.49
        C5 = -29.06
        C6 = 0.72
        C7 = 0.89
        Dn = np.array(np.exp(C1 + C2*(Ac/PGA) + C3*(Ac/PGA)**2 + C4*(Ac/PGA)**3 + C5*(Ac/PGA)**4 +
                      C6*np.log(PGA)+C7*(M-6)))  # Equation from Saygili and Rathje (2008)/Rathje and Saygili (2009)
        Dn[np.isnan(Dn)] = 0.
        logDnstd = 0.732 + 0.789*(Ac/PGA) - 0.539*(Ac/PGA)**2
        logtype = 'ln'

    elif model == 'RS_PGA_PGV':
        if PGV is None:
            raise Exception('You must enter a value for M to use the RS_PGA_PGV model')
        C1 = -1.56
        C2 = -4.58
        C3 = -20.84
        C4 = 44.75
        C5 = -30.50
        C6 = -0.64
        C7 = 1.55
        Dn = np.array(np.exp(C1 + C2*(Ac/PGA) + C3*(Ac/PGA)**2 + C4*(Ac/PGA)**3 + C5*(Ac/PGA)**4
                      + C6*np.log(PGA)+C7*np.log(PGV)))  # Equation from Saygili and Rathje (2008)/Rathje and Saygili (2009)
        Dn[np.isnan(Dn)] = 0.
        logDnstd = 0.405 + 0.524*(Ac/PGA)
        logtype = 'ln'

    elif model == 'BT_PGA_M':
        if M is None:
            raise Exception('You must enter a value for M to use the BT_PGA_M model')
        Dn = np.array(np.exp(-0.22 - 2.83*np.log(Ac) - 0.333*(np.log(Ac))**2 + 0.566*np.log(Ac)*np.log(PGA)
                      + 3.04*np.log(PGA) - 0.244*(np.log(PGA))**2 + 0.278*(M-7.)))
        Dn[np.isnan(Dn)] = 0.
        logDnstd = np.ones(np.shape(Dn))*0.66
        logtype = 'log10'

    if flag == 1:
        Dn = float(Dn)
        logDnstd = float(logDnstd)

    return Dn, logDnstd, logtype
