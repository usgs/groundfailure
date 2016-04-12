#!/usr/bin/env python

"""
Newmark based landslide mechanistic_models
"""

#stdlib imports
import os.path
import warnings
import urllib2
import tempfile
import collections

#turn off all warnings...
warnings.filterwarnings('ignore')

#local imports
from mapio.shake import ShakeGrid
from mapio.gdal import GDALGrid
from mapio.geodict import GeoDict

#third party imports
import numpy as np


def getGridURL(gridurl):
    """This function downloads the url and returns the corresponding file object

    :param gridurl: string defining a url
    :type gridurl: string

    :returns gridfile: file object corresponding to the url
    """
    gridfile = None
    try:
        fh = urllib2.urlopen(gridurl)
        data = fh.read()
        fd, gridfile = tempfile.mkstemp()
        os.close(fd)
        f = open(gridfile, 'wt')
        f.write(data)
        f.close()
        fh.close()
    except:
        raise IOError('Could not retrieve data from %s' % gridurl)
    return gridfile


def isURL(gridurl):
    """This function determines if a string is a valid url

    :param gridurl: string defining the potential url
    :type gridurl: string

    :returns isURL: True if griurl is a valid url, False otherwise
    """
    isURL = False
    try:
        urllib2.urlopen(gridurl)
        isURL = True
    except:
        pass
    return isURL


def HAZUS(shakefile, config, saveinputs=False, modeltype='coverage', regressionmodel='J_PGA', probtype='jibson2000', bounds=None):
    """
    Runs HAZUS landslide procedure (FEMA, 2003, Chapter 4) using susceptiblity categories from defined by HAZUS manual (I-X)

    :param shakefile: URL or complete file path to the location of the Shakemap to use as input
    :type shakefile: string:
    :param config: Model configuration file object containing locations of input files and other input values config = ConfigObj(configfilepath)
    :type config: ConfigObj
    :param saveinputs: Whether or not to return the model input layers, False (defeault) returns only the model output (one layer)
    :type saveinputs: boolean
    :param modeltype: 'coverage' if critical acceleration is exceeded by pga, this gives the estimated areal coverage of landsliding for that cell
        'dn_hazus' - Outputs Newmark displacement using HAZUS methods without relating to probability of failure
        'dn_prob' - Estimates Newmark displacement using HAZUS methods and relates to probability of failure using param probtype
        'ac_classic_dn' - Uses the critical acceleration defined by HAZUS methodology and uses regression model defined by regressionmodel param to get Newmark displacement without relating to probability of failure
        'ac_classic_prob' - Uses the critical acceleration defined by HAZUS methodology and uses regression model defined by regressionmodel param to get Newmark displacement and probability defined by probtype method
    :type modeltype: string
    :param regressionmodel:
        Newmark displacement regression model to use
        'J_PGA' (default) - PGA-based model from Jibson (2007) - equation 6
        'J_PGA_M' - PGA and M-based model from Jibson (2007) - equation 7
        'RS_PGA_M' - PGA and M-based model from from Rathje and Saygili (2009)
        'RS_PGA_PGV' - PGA and PGV-based model from Saygili and Rathje (2008) - equation 6
    :type regressionmodel: string
    :param probtype: Method used to estimate probability. Entering 'jibson2000' uses equation 5 from Jibson et al. (2000) to estimate probability from Newmark displacement. 'threshold' uses a specified threshold of Newmark displacement (defined in config file) and assumes anything greather than this threshold fails
    :type probtype: string
    :param bounds: Boundaries to compute over if different from ShakeMap boundaries as tuple (xmin, ymin, xmax, ymax)

    :returns maplayers:  Dictionary containing output and input layers (if saveinputs=True) along with metadata formatted like maplayers['layer name']={'grid': mapio grid2D object, 'label': 'label for colorbar and top line of subtitle', 'type': 'output or input to model', 'description': 'detailed description of layer for subtitle, potentially including source information'}
    :type maplayers: OrderedDict
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
    try:
        susfile = config['mechanistic_models']['hazus']['layers']['susceptibility']['file']
        shkgdict = ShakeGrid.getFileGeoDict(shakefile, adjust='res')
        susdict = GDALGrid.getFileGeoDict(susfile)
        if bounds is not None:  # Make sure bounds are within ShakeMap Grid
            if shkgdict.xmin > bounds[0] or shkgdict.xmax < bounds[2] or shkgdict.ymin > bounds[1] or shkgdict.ymax < bounds[3]:
                print('Specified bounds are outside shakemap area, using ShakeMap bounds instead')
                bounds = None
        if bounds is not None:
            tempgdict1 = GeoDict({'xmin': bounds[0], 'ymin': bounds[1], 'xmax': bounds[2], 'ymax': bounds[3], 'dx': 100., 'dy': 100., 'nx': 100., 'ny': 100.}, adjust='res')
            tempgdict = susdict.getBoundsWithin(tempgdict1)
        else:
            tempgdict = susdict.getBoundsWithin(shkgdict)
        sus = GDALGrid.load(susfile, samplegeodict=tempgdict, resample=False)
        gdict = sus.getGeoDict()
        susdat = sus.getData()
    except Exception as e:
        raise IOError('Unable to read in susceptibility category file specified in config, %s,' % e)
        return

    try:  # Try to fetch source information from config
        modelsref = config['mechanistic_models']['hazus']['shortref']
        modellref = config['mechanistic_models']['hazus']['longref']
        sussref = config['mechanistic_models']['hazus']['layers']['susceptibility']['shortref']
        suslref = config['mechanistic_models']['hazus']['layers']['susceptibility']['longref']
    except:
        print('Was not able to retrieve all references from config file. Continuing')

    try:
        dnthresh = float(config['mechanistic_models']['hazus']['values']['dnthresh'])
    except:
        if probtype == 'threshold':
            dnthresh = 5.
            print('Unable to find dnthresh in config, using 5cm')

    # Load in shakemap, resample to susceptibility file
    shakemap = ShakeGrid.load(shakefile, adjust='res')

    PGA = shakemap.getLayer('pga').subdivide(gdict).getData().astype(float)/100.  # in units of g
    PGV = shakemap.getLayer('pgv').subdivide(gdict).getData().astype(float)  # cm/sec
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
        # # But this way is even slower, takes 2x as long
        # numrows, numcols = np.shape(areal)
        # for j in np.arange(numrows):
        #     for k in np.arange(numcols):
        #         acval = Ac[j, k]
        #         if PGA[j, k] >= acval:
        #             if acval == 0.6:
        #                 areal[j, k] = 0.01
        #             elif acval == 0.5:
        #                 areal[j, k] = 0.02
        #             elif acval == 0.4:
        #                 areal[j, k] = 0.03
        #             elif acval == 0.35:
        #                 areal[j, k] = 0.05
        #             elif acval == 0.3:
        #                 areal[j, k] = 0.08
        #             elif acval == 0.25:
        #                 areal[j, k] = 0.1
        #             elif acval == 0.2:
        #                 areal[j, k] = 0.15
        #             elif acval == 0.15:
        #                 areal[j, k] = 0.2
        #             elif acval == 0.1:
        #                 areal[j, k] = 0.25
        #             elif acval == 0.05:
        #                 areal[j, k] = 0.3

    elif modeltype == 'dn_hazus' or modeltype == 'dn_prob':
        ed_low, ed_high = est_disp(Ac, PGA)
        ed_mean = np.mean((np.dstack((ed_low, ed_high))), axis=2)  # Get mean estimated displacements
        dn = ed_mean * numcycles(M) * PGA
    else:  # Calculate newmark displacement using a regression model
        if regressionmodel is 'J_PGA':
            dn = J_PGA(Ac, PGA)
        elif regressionmodel is 'J_PGA_M':
            dn = J_PGA_M(Ac, PGA, M)
        elif regressionmodel is 'RS_PGA_M':
            dn = RS_PGA_M(Ac, PGA, M)
        elif regressionmodel is 'RS_PGA_PGV':
            dn = RS_PGA_PGV(Ac, PGA, PGV)
        else:
            print('Unrecognized model, using J_PGA\n')
            dn = J_PGA(Ac, PGA)

    # Calculate probability from dn, if necessary for selected model
    if modeltype == 'ac_classic_prob' or modeltype == 'dn_prob':
        if probtype.lower() in 'jibson2000':
            PROB = 0.335*(1-np.exp(-0.048*dn**1.565))
            dnthresh = None
        elif probtype.lower() in 'threshold':
            PROB = dn.copy()
            PROB[PROB <= dnthresh] = 0
            PROB[PROB > dnthresh] = 1
        else:
            raise NameError('invalid probtype, assuming jibson2000')
            PROB = 0.335*(1-np.exp(-0.048*dn**1.565))
            dnthresh = None

    # Turn output and inputs into into grids and put in maplayers dictionary
    maplayers = collections.OrderedDict()

    temp = shakemap.getShakeDict()
    shakedetail = '%s_ver%s' % (temp['shakemap_id'], temp['shakemap_version'])

    if modeltype == 'coverage':
        maplayers['model'] = {'grid': GDALGrid(areal, gdict), 'label': 'Areal coverage', 'type': 'output', 'description': {'name': modelsref, 'longref': modellref, 'units': 'coverage', 'shakemap': shakedetail, 'parameters': {'modeltype': modeltype}}}
    elif modeltype == 'dn_hazus':
        maplayers['model'] = {'grid': GDALGrid(dn, gdict), 'label': 'Dn (cm)', 'type': 'output', 'description': {'name': modelsref, 'longref': modellref, 'units': 'displacement', 'shakemap': shakedetail, 'parameters': {'regressionmodel': regressionmodel, 'modeltype': modeltype}}}
    elif modeltype == 'ac_classic_dn':
        maplayers['model'] = {'grid': GDALGrid(dn, gdict), 'label': 'Dn (cm)', 'type': 'output', 'description': {'name': modelsref, 'longref': modellref, 'units': 'displacement', 'shakemap': shakedetail, 'parameters': {'regressionmodel': regressionmodel, 'modeltype': modeltype}}}
    elif modeltype == 'dn_prob':
        maplayers['model'] = {'grid': GDALGrid(PROB, gdict), 'label': 'Landslide Probability', 'type': 'output', 'description': {'name': modelsref, 'longref': modellref, 'units': 'probability', 'shakemap': shakedetail, 'parameters': {'regressionmodel': regressionmodel, 'dnthresh_cm': dnthresh, 'modeltype': modeltype, 'probtype': probtype}}}
    elif modeltype == 'ac_classic_prob':
        maplayers['model'] = {'grid': GDALGrid(PROB, gdict), 'label': 'Landslide Probability', 'type': 'output', 'description': {'name': modelsref, 'longref': modellref, 'units': 'probability', 'shakemap': shakedetail, 'parameters': {'regressionmodel': regressionmodel, 'dnthresh_cm': dnthresh, 'modeltype': modeltype, 'probtype': probtype}}}

    if saveinputs is True:
        maplayers['suscat'] = {'grid': sus, 'label': 'Susceptibility Category', 'type': 'input', 'description': {'name': sussref, 'longref': suslref, 'units': 'Category'}}
        maplayers['Ac'] = {'grid': GDALGrid(Ac, gdict), 'label': 'Ac (g)', 'type': 'output', 'description': {'units': 'g', 'shakemap': shakedetail}}
        maplayers['pga'] = {'grid': GDALGrid(PGA, gdict), 'label': 'PGA (g)', 'type': 'input', 'description': {'units': 'g', 'shakemap': shakedetail}}
        if 'pgv' in regressionmodel.lower():
            maplayers['pgv'] = {'grid': GDALGrid(PGV, gdict), 'label': 'PGV (cm/s)', 'type': 'input', 'description': {'units': 'cm/s', 'shakemap': shakedetail}}
        if 'dn' not in modeltype.lower() and modeltype != 'coverage':
            maplayers['dn'] = {'grid': GDALGrid(dn, gdict), 'label': 'Dn (cm)', 'type': 'output', 'description': {'units': 'displacement', 'shakemap': shakedetail, 'parameters': {'regressionmodel': regressionmodel, 'modeltype': modeltype}}}

    return maplayers


def est_disp(Ac, PGA):
    """
    Get estimated displacement factor, digitized from figure in HAZUS manual pg. 4-37, which is based on Makdisi and Seed (1978) according to HAZUS, but the equations don't appear in that document
    Assumes PGA is equal to ais (induced acceleration)
    """
    from scipy.interpolate import interp1d

    acais = Ac/PGA  # Get "Expected displacement factor"
    xlow = np.array((0.0994941, 0.198426, 0.278246, 0.370433, 0.450253, 0.53457, 0.594154, 0.640247, 0.684092, 0.72344, 0.770658, 0.802136, 0.851602, 0.897695))
    ylow = np.array((22.0186, 11.4578, 7.09682, 3.85734, 2.28738, 1.24326, 0.822041, 0.543531, 0.367292, 0.237622, 0.137873, 0.101646, 0.044438, 0.0211954))
    xhigh = np.array((0.100618, 0.179314, 0.251265, 0.319843, 0.406408, 0.48398, 0.567173, 0.634626, 0.700956, 0.751546, 0.80326, 0.857223, 0.901068))
    yhigh = np.array((39.6377, 21.5443, 13.638, 9.21593, 4.90126, 2.84381, 1.6145, 1.02201, 0.592993, 0.351641, 0.208521, 0.0931681, 0.0495491))
    f_low = interp1d(xlow, ylow, bounds_error=False, fill_value=0.)
    f_high = interp1d(xhigh, yhigh, bounds_error=False, fill_value=0.)
    ed_low = f_low(acais)
    ed_high = f_high(acais)
    return ed_low, ed_high


def numcycles(M):
    """
    Estimate the number of earthquake cycles using Seed and Idriss (1982) relationship as presented in HAZUS manual, pg. 4-35 and Fig. 4.13

    :param M: earthquake magnitude
    :returns n: number of cycles
    """
    n = 0.3419*M**3 - 5.5214*M**2 + 33.6154*M - 70.7692
    return n


def classic(shakefile, config, saveinputs=False, regressionmodel='J_PGA', probtype='jibson2000', slopediv=1., codiv=1., bounds=None):
    """This function uses the Newmark method to estimate probability of failure at each grid cell.
    Factor of Safety and critcal accelerations are calculated following Jibson et al. (2000) and the
    Newmark displacement is estimated using PGA, PGV, and/or Magnitude (depending on equation used)
    from Shakemap with regression equations from Jibson (2007), Rathje and Saygili (2008) and
    Saygili and Rathje (2009)

    :param shakefile: URL or complete file path to the location of the Shakemap to use as input
    :type shakefile: string:
    :param config: Model configuration file object containing locations of input files and other input values config = ConfigObj(configfilepath)
    :type config: ConfigObj
    :param saveinputs: Whether or not to return the model input layers, False (defeault) returns only the model output (one layer)
    :type saveinputs: boolean
    :param regressionmodel:
        Newmark displacement regression model to use
        'J_PGA' (default) - PGA-based model from Jibson (2007) - equation 6
        'J_PGA_M' - PGA and M-based model from Jibson (2007) - equation 7
        'RS_PGA_M' - PGA and M-based model from from Rathje and Saygili (2009)
        'RS_PGA_PGV' - PGA and PGV-based model from Saygili and Rathje (2008) - equation 6
    :type regressionmodel: string
    :param probtype: Method used to estimate probability. Entering 'jibson2000' uses equation 5 from Jibson et al. (2000) to estimate probability from Newmark displacement. 'threshold' uses a specified threshold of Newmark displacement (defined in config file) and assumes anything greather than this threshold fails
    :type probtype: string
    :param slopediv: Divide slope by this number to get slope in degrees (Verdin datasets need to be divided by 100)
    :type slopediv: float
    :param codiv: Divide cohesion by this number to get reasonable numbers (For Godt method, need to divide by 10 because that is how it was calibrated, but values are reasonable without multiplying for regular analysis)
    :type codiv: float

    :returns maplayers:  Dictionary containing output and input layers (if saveinputs=True) along with metadata formatted like maplayers['layer name']={'grid': mapio grid2D object, 'label': 'label for colorbar and top line of subtitle', 'type': 'output or input to model', 'description': 'detailed description of layer for subtitle, potentially including source information'}
    :type maplayers: OrderedDict

    :raises NameError: when unable to parse the config correctly (probably a formatting issue in the configfile) or when unable to find the shakefile (Shakemap URL or filepath) - these cause program to end
    :raises NameError: when probtype does not match a predifined probability type, will cause to default to 'jibson2000'

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

    # Parse config - should make it so it uses defaults if any are missing...
    try:
        slopefile = config['mechanistic_models']['classic_newmark']['layers']['slope']['file']
        slopeunits = config['mechanistic_models']['classic_newmark']['layers']['slope']['units']
        cohesionfile = config['mechanistic_models']['classic_newmark']['layers']['cohesion']['file']
        cohesionunits = config['mechanistic_models']['classic_newmark']['layers']['cohesion']['units']
        frictionfile = config['mechanistic_models']['classic_newmark']['layers']['friction']['file']
        frictionunits = config['mechanistic_models']['classic_newmark']['layers']['friction']['units']

        thick = float(config['mechanistic_models']['classic_newmark']['parameters']['thick'])
        uwt = float(config['mechanistic_models']['classic_newmark']['parameters']['uwt'])
        nodata_cohesion = float(config['mechanistic_models']['classic_newmark']['parameters']['nodata_cohesion'])
        nodata_friction = float(config['mechanistic_models']['classic_newmark']['parameters']['nodata_friction'])
        try:
            dnthresh = float(config['mechanistic_models']['classic_newmark']['parameters']['dnthresh'])
        except:
            if probtype == 'threshold':
                dnthresh = 5.
                print('Unable to find dnthresh in config, using 5cm')
            else:
                dnthresh = None
        fsthresh = float(config['mechanistic_models']['classic_newmark']['parameters']['fsthresh'])
        acthresh = float(config['mechanistic_models']['classic_newmark']['parameters']['acthresh'])
        slopethresh = float(config['mechanistic_models']['classic_newmark']['parameters']['slopethresh'])
        try:
            m = float(config['mechanistic_models']['classic_newmark']['parameters']['m'])
        except:
            print('no constant saturated thickness specified, m=0 if no watertable file is found')
            m = 0.
    except Exception as e:
        raise NameError('Could not parse configfile, %s' % e)
        return

    try:  # Try to fetch source information from config
        modelsref = config['mechanistic_models']['classic_newmark']['shortref']
        modellref = config['mechanistic_models']['classic_newmark']['longref']
        slopesref = config['mechanistic_models']['classic_newmark']['layers']['slope']['shortref']
        slopelref = config['mechanistic_models']['classic_newmark']['layers']['slope']['longref']
        cohesionsref = config['mechanistic_models']['classic_newmark']['layers']['cohesion']['shortref']
        cohesionlref = config['mechanistic_models']['classic_newmark']['layers']['cohesion']['longref']
        frictionsref = config['mechanistic_models']['classic_newmark']['layers']['friction']['shortref']
        frictionlref = config['mechanistic_models']['classic_newmark']['layers']['friction']['longref']
    except:
        print('Was not able to retrieve all references from config file. Continuing')

    # Cut and resample all files
    shkgdict = ShakeGrid.getFileGeoDict(shakefile, adjust='res')
    slpdict = GDALGrid.getFileGeoDict(slopefile)
    if bounds is not None:  # Make sure bounds are within ShakeMap Grid
        if shkgdict.xmin > bounds[0] or shkgdict.xmax < bounds[2] or shkgdict.ymin > bounds[1] or shkgdict.ymax < bounds[3]:
            print('Specified bounds are outside shakemap area, using ShakeMap bounds instead')
            bounds = None
    if bounds is not None:
        tempgdict = GeoDict({'xmin': bounds[0], 'ymin': bounds[1], 'xmax': bounds[2], 'ymax': bounds[3], 'dx': 100., 'dy': 100., 'nx': 100., 'ny': 100.}, adjust='res')
        gdict = slpdict.getBoundsWithin(tempgdict)
    else:  # Get boundaries from shakemap if not specified
        shkgdict = ShakeGrid.getFileGeoDict(shakefile, adjust='res')
        slpdict = GDALGrid.getFileGeoDict(slopefile)
        gdict = slpdict.getBoundsWithin(shkgdict)

    # Load in slope file
    slopegrid = GDALGrid.load(slopefile, samplegeodict=gdict, resample=False)
    gdict = slopegrid.getGeoDict()  # Get this again just in case it changed
    slope = slopegrid.getData()/slopediv  # Adjust slope to degrees, if needed
    # Change any zero slopes to a very small number to avoid dividing by zero later
    slope[slope == 0] = 1e-8

    # Load in shakefile
    if not os.path.isfile(shakefile):
        if isURL(shakefile):
            shakefile = getGridURL(shakefile)  # returns a file object
        else:
            raise NameError('Could not find "%s" as a file or a valid url' % (shakefile))
            return

    # Load in shakemap, resample to slope file (this will be important when go to higher res)
    shakemap = ShakeGrid.load(shakefile, samplegeodict=gdict, resample=True, method='linear', adjust='res')
    M = shakemap.getEventDict()['magnitude']

    # Read in the cohesion and friction files, resampled to slope grid
    cohesion = GDALGrid.load(cohesionfile, samplegeodict=gdict, resample=True, method='nearest').getData()/codiv
    cohesion[np.isnan(cohesion)] = nodata_cohesion
    friction = GDALGrid.load(frictionfile, samplegeodict=gdict, resample=True, method='nearest').getData()
    friction[np.isnan(friction)] = nodata_friction

    # See if there is a water table depth file and read it in if there is
    try:
        waterfile = config['mechanistic_models']['classic_newmark']['layers']['watertable']['file']
        watertable = GDALGrid.load(waterfile, samplegeodict=gdict, resample=True, method='linear').getData()  # Needs to be in meters!
        uwtw = float(config['mechanistic_models']['classic_newmark']['parameters']['uwtw'])
        try:
            watersref = config['mechanistic_models']['classic_newmark']['layers']['watertable']['shortref']
            waterlref = config['mechanistic_models']['classic_newmark']['layers']['watertable']['longref']
        except:
            print('Was not able to retrieve water table references from config file. Continuing')

    except:
        print('Water table file not specified or readable, assuming constant saturated thickness proportion of %0.1f' % m)
        watertable = None
        try:
            uwtw = float(config['mechanistic_models']['classic_newmark']['parameters']['uwtw'])
        except:
            print('Could not read soil wet unit weight, using 18.8 kN/m3')
            uwtw = 18.8

    # Factor of safety
    if watertable is not None:
        watertable[watertable > thick] = thick
        m = (thick - watertable)/thick
    FS = cohesion/(uwt*thick*np.sin(slope*(np.pi/180.))) + np.tan(friction*(np.pi/180.))/np.tan(slope*(np.pi/180.)) - (m*uwtw*np.tan(friction*(np.pi/180.)))/(uwt*np.tan(slope*(np.pi/180.)))
    FS[FS < fsthresh] = fsthresh

    # Compute critical acceleration, in g
    Ac = (FS-1)*np.sin(slope*(np.pi/180.)).astype(float)  # This gives ac in g, equations that multiply by g give ac in m/s2
    Ac[Ac < acthresh] = acthresh
    Ac[slope < slopethresh] = float('nan')

    # Get PGA in g (PGA is %g in ShakeMap, convert to g)
    PGA = shakemap.getLayer('pga').getData().astype(float)/100.
    PGV = shakemap.getLayer('pgv').getData().astype(float)

    np.seterr(invalid='ignore')  # Ignore errors so still runs when Ac > PGA, just leaves nan instead of crashing

    if regressionmodel is 'J_PGA':
        Dn = J_PGA(Ac, PGA)

    if regressionmodel is 'J_PGA_M':
        Dn = J_PGA_M(Ac, PGA, M)

    if regressionmodel is 'RS_PGA_M':
        Dn = RS_PGA_M(Ac, PGA, M)

    if regressionmodel is 'RS_PGA_PGV':
        Dn = RS_PGA_PGV(Ac, PGA, PGV)

    units = 'probability'
    label = 'Landslide Probability'
    if probtype.lower() in 'jibson2000':
        PROB = 0.335*(1-np.exp(-0.048*Dn**1.565))
        dnthresh = None
    elif probtype.lower() in 'threshold':
        PROB = Dn.copy()
        PROB[PROB <= dnthresh] = 0
        PROB[PROB > dnthresh] = 1
        units = 'prediction'
        label = 'Predicted Landslides'
    else:
        raise NameError('invalid probtype, assuming jibson2000')
        PROB = 0.335*(1-np.exp(-0.048*Dn**1.565))
        dnthresh = None

    # Turn output and inputs into into grids and put in mapLayers dictionary
    maplayers = collections.OrderedDict()

    temp = shakemap.getShakeDict()
    shakedetail = '%s_ver%s' % (temp['shakemap_id'], temp['shakemap_version'])

    if watertable is not None:
        des = 'variable'
    else:
        des = m
    description = {'name': modelsref, 'longref': modellref, 'units': units, 'shakemap': shakedetail, 'parameters': {'regressionmodel': regressionmodel, 'thickness_m': thick, 'unitwt_kNm3': uwt, 'dnthresh_cm': dnthresh, 'acthresh_g': acthresh, 'fsthresh': fsthresh, 'slopethresh': slopethresh, 'sat_proportion': des}}

    maplayers['model'] = {'grid': GDALGrid(PROB, gdict), 'label': label, 'type': 'output', 'description': description}

    if saveinputs is True:
        maplayers['pga'] = {'grid': GDALGrid(PGA, gdict), 'label': 'PGA (g)', 'type': 'input', 'description': {'units': 'g', 'shakemap': shakedetail}}
        maplayers['FS'] = {'grid': GDALGrid(FS, gdict), 'label': 'Factor of Safety', 'type': 'input', 'description': {'units': 'unitless'}}
        maplayers['Ac'] = {'grid': GDALGrid(Ac, gdict), 'label': 'Critical acceleration (g)', 'type': 'input'}
        maplayers['Dn'] = {'grid': GDALGrid(Dn, gdict), 'label': 'Newmark Displacement (cm)', 'type': 'input'}
        maplayers['slope'] = {'grid': GDALGrid(slope, gdict), 'label': 'Max slope ($^\circ$)', 'type': 'input', 'description': {'units': 'degrees', 'name': slopesref, 'longref': slopelref}}
        maplayers['cohesion'] = {'grid': GDALGrid(cohesion, gdict), 'label': 'Cohesion (kPa)', 'type': 'input', 'description': {'units': 'kPa (adjusted)', 'name': cohesionsref, 'longref': cohesionlref}}
        maplayers['friction angle'] = {'grid': GDALGrid(friction, gdict), 'label': 'Friction angle ($^\circ$)', 'type': 'input', 'description': {'units': 'degrees', 'name': frictionsref, 'longref': frictionlref}}
        if 'PGV' in regressionmodel:
            maplayers['pgv'] = {'grid': GDALGrid(PGV, gdict), 'label': 'PGV (cm/s)', 'type': 'input', 'description': {'units': 'cm/s', 'shakemap': shakedetail}}
        if watertable is not None:
            maplayers['sat thick prop'] = {'grid': GDALGrid(m, gdict), 'label': 'Saturated thickness proprtion [0,1]', 'type': 'input', 'description': {'units': 'meters', 'name': watersref, 'longref': waterlref}}

    return maplayers


def godt2008(shakefile, config, saveinputs=False, regressionmodel='J_PGA', bounds=None, slopediv=100., codiv=10.):
    """ This function runs the Godt et al. (2008) global method for a given ShakeMap. The Factor of Safety
    is calculated using infinite slope analysis assumuing dry conditions. The method uses threshold newmark
    displacement and estimates areal coverage by doing the calculations for each slope quantile
    TO DO - add 'all' - averages Dn from all four equations, add term to convert PGA and PGV to Ia and use other equations, add Ambraseys and Menu (1988) option

    :param shakefile: url or filepath to shakemap xml file
    :type shakefile: string
    :param config: ConfigObj of config file containing inputs required for running the model
    :type config: ConfigObj
    :param saveinputs: Whether or not to return the model input layers, False (defeault) returns only the model output (one layer)
    :type saveinputs: boolean
    :param regressionmodel:
        Newmark displacement regression model to use
        'J_PGA' (default) - PGA-based model from Jibson (2007) - equation 6
        'J_PGA_M' - PGA and M-based model from Jibson (2007) - equation 7
        'RS_PGA_M' - PGA and M-based model from from Rathje and Saygili (2009)
        'RS_PGA_PGV' - PGA and PGV-based model from Saygili and Rathje (2008) - equation 6
    :type regressionmodel: string
    :param probtype: Method used to estimate probability. Entering 'jibson2000' uses equation 5 from Jibson et al. (2000) to estimate probability from Newmark displacement. 'threshold' uses a specified threshold of Newmark displacement (defined in config file) and assumes anything greather than this threshold fails
    :type probtype: string
    :param slopediv: Divide slope by this number to get slope in degrees (Verdin datasets need to be divided by 100)
    :type slopediv: float
    :param codiv: Divide cohesion by this number to get reasonable numbers (For Godt method, need to divide by 10 because that is how it was calibrated, but values are reasonable without multiplying for regular analysis)
    :type codiv: float

    :returns maplayers:  Dictionary containing output and input layers (if saveinputs=True) along with metadata formatted like maplayers['layer name']={'grid': mapio grid2D object, 'label': 'label for colorbar and top line of subtitle', 'type': 'output or input to model', 'description': 'detailed description of layer for subtitle, potentially including source information'}
    :type maplayers: OrderedDict

    :raises NameError: when unable to parse the config correctly (probably a formatting issue in the configfile) or when unable to find the shakefile (Shakemap URL or filepath) - these cause program to end
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
        slopefilepath = config['mechanistic_models']['godt_2008']['layers']['slope']['filepath']
        slopeunits = config['mechanistic_models']['godt_2008']['layers']['slope']['units']
        cohesionfile = config['mechanistic_models']['godt_2008']['layers']['cohesion']['file']
        cohesionunits = config['mechanistic_models']['godt_2008']['layers']['cohesion']['units']
        frictionfile = config['mechanistic_models']['godt_2008']['layers']['friction']['file']
        frictionunits = config['mechanistic_models']['godt_2008']['layers']['friction']['units']

        thick = float(config['mechanistic_models']['godt_2008']['parameters']['thick'])
        uwt = float(config['mechanistic_models']['godt_2008']['parameters']['uwt'])
        nodata_cohesion = float(config['mechanistic_models']['godt_2008']['parameters']['nodata_cohesion'])
        nodata_friction = float(config['mechanistic_models']['godt_2008']['parameters']['nodata_friction'])
        dnthresh = float(config['mechanistic_models']['godt_2008']['parameters']['dnthresh'])
        fsthresh = float(config['mechanistic_models']['godt_2008']['parameters']['fsthresh'])
        acthresh = float(config['mechanistic_models']['godt_2008']['parameters']['acthresh'])
    except Exception as e:
        raise NameError('Could not parse configfile, %s' % e)
        return

    # TO DO, ADD ERROR CATCHING ON UNITS, MAKE SURE THEY ARE WHAT THEY SHOULD BE FOR THIS MODEL

    try:  # Try to fetch source information from config
        modelsref = config['mechanistic_models']['godt_2008']['shortref']
        modellref = config['mechanistic_models']['godt_2008']['longref']
        slopesref = config['mechanistic_models']['godt_2008']['layers']['slope']['shortref']
        slopelref = config['mechanistic_models']['godt_2008']['layers']['slope']['longref']
        cohesionsref = config['mechanistic_models']['godt_2008']['layers']['cohesion']['shortref']
        cohesionlref = config['mechanistic_models']['godt_2008']['layers']['cohesion']['longref']
        frictionsref = config['mechanistic_models']['godt_2008']['layers']['friction']['shortref']
        frictionlref = config['mechanistic_models']['godt_2008']['layers']['friction']['longref']
    except:
        print('Was not able to retrieve all references from config file. Continuing')

    # Load in shakefile
    if not os.path.isfile(shakefile):
        if isURL(shakefile):
            shakefile = getGridURL(shakefile)  # returns a file object
        else:
            raise NameError('Could not find "%s" as a file or a valid url' % (shakefile))
            return

    shkgdict = ShakeGrid.getFileGeoDict(shakefile, adjust='res')
    if bounds is not None:  # Make sure bounds are within ShakeMap Grid
        if shkgdict.xmin > bounds[0] or shkgdict.xmax < bounds[2] or shkgdict.ymin > bounds[1] or shkgdict.ymax < bounds[3]:
            print('Specified bounds are outside shakemap area, using ShakeMap bounds instead')
            bounds = None
    if bounds is not None:
        tempgdict = GeoDict({'xmin': bounds[0], 'ymin': bounds[1], 'xmax': bounds[2], 'ymax': bounds[3], 'dx': shkgdict.dx, 'dy': shkgdict.dy, 'nx': shkgdict.nx, 'ny': shkgdict.ny}, adjust='res')
        gdict = shkgdict.getBoundsWithin(tempgdict)
        shakemap = ShakeGrid.load(shakefile, samplegeodict=gdict, adjust='bounds')
    else:
        shakemap = ShakeGrid.load(shakefile, adjust='res')
    shkgdict = shakemap.getGeoDict()  # Get updated geodict
    M = shakemap.getEventDict()['magnitude']

    # Read in all the slope files, divide all by 100 to get to slope in degrees (because input files are multiplied by 100.)
    slopes = []
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope_min.bil'), samplegeodict=shkgdict, resample=True, method='linear').getData()/slopediv)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope10.bil'), samplegeodict=shkgdict, resample=True, method='linear').getData()/slopediv)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope30.bil'), samplegeodict=shkgdict, resample=True, method='linear').getData()/slopediv)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope50.bil'), samplegeodict=shkgdict, resample=True, method='linear').getData()/slopediv)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope70.bil'), samplegeodict=shkgdict, resample=True, method='linear').getData()/slopediv)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope90.bil'), samplegeodict=shkgdict, resample=True, method='linear').getData()/slopediv)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope_max.bil'), samplegeodict=shkgdict, resample=True, method='linear').getData()/slopediv)
    slopestack = np.dstack(slopes)

    # Change any zero slopes to a very small number to avoid dividing by zero later
    slopestack[slopestack == 0] = 1e-8

    # Read in the cohesion and friction files and duplicate layers so they are same shape as slope structure
    #import pdb; pdb.set_trace()
    cohesion = np.repeat(GDALGrid.load(cohesionfile, samplegeodict=shakemap.getGeoDict(), resample=True, method='nearest').getData()[:, :, np.newaxis]/codiv, 7, axis=2)
    cohesion[cohesion == -999.9] = nodata_cohesion
    cohesion[cohesion == 0] = nodata_cohesion
    friction = np.repeat(GDALGrid.load(frictionfile, samplegeodict=shakemap.getGeoDict(), resample=True, method='nearest').getData()[:, :, np.newaxis], 7, axis=2)
    friction[friction == -9999] = nodata_friction
    friction[friction == 0] = nodata_friction

    # Do the calculations using Jibson (2007) PGA only model for Dn
    FS = cohesion/(uwt*thick*np.sin(slopestack*(np.pi/180.))) + np.tan(friction*(np.pi/180.))/np.tan(slopestack*(np.pi/180.))
    FS[FS < fsthresh] = fsthresh

    # Compute critical acceleration, in g
    Ac = (FS-1)*np.sin(slopestack*(np.pi/180.)).astype(float)  # This gives ac in g, equations that multiply by g give ac in m/s2
    Ac[Ac < acthresh] = acthresh

    # Get PGA in g (PGA is %g in ShakeMap, convert to g)
    PGA = np.repeat(shakemap.getLayer('pga').getData()[:, :, np.newaxis]/100., 7, axis=2).astype(float)

    if 'PGV' in regressionmodel:  # Load in PGV also, in cm/sec
        PGV = np.repeat(shakemap.getLayer('pgv').getData()[:, :, np.newaxis], 7, axis=2).astype(float)

    np.seterr(invalid='ignore')  # Ignore errors so still runs when Ac > PGA, just leaves nan instead of crashing

    if regressionmodel is 'J_PGA':
        Dn = J_PGA(Ac, PGA)

    if regressionmodel is 'J_PGA_M':
        Dn = J_PGA_M(Ac, PGA, M)

    if regressionmodel is 'RS_PGA_M':
        Dn = RS_PGA_M(Ac, PGA, M)

    if regressionmodel is 'RS_PGA_PGV':
        Dn = RS_PGA_PGV(Ac, PGA, PGV)

    PROB = Dn.copy()
    PROB[PROB < dnthresh] = 0.
    PROB[PROB >= dnthresh] = 1.
    PROB = np.sum(PROB, axis=2)
    PROB[PROB == 1.] = 0.01
    PROB[PROB == 2.] = 0.10
    PROB[PROB == 3.] = 0.30
    PROB[PROB == 4.] = 0.50
    PROB[PROB == 5.] = 0.70
    PROB[PROB == 6.] = 0.90
    PROB[PROB == 7.] = 0.99

    # Turn output and inputs into into grids and put in mapLayers dictionary
    maplayers = collections.OrderedDict()

    temp = shakemap.getShakeDict()
    shakedetail = '%s_ver%s' % (temp['shakemap_id'], temp['shakemap_version'])

    description = {'name': modelsref, 'longref': modellref, 'units': 'coverage', 'shakemap': shakedetail, 'parameters': {'regressionmodel': regressionmodel, 'thickness_m': thick, 'unitwt_kNm3': uwt, 'dnthresh_cm': dnthresh, 'acthresh_g': acthresh, 'fsthresh': fsthresh}}

    maplayers['model'] = {'grid': GDALGrid(PROB, shakemap.getGeoDict()), 'label': 'Areal coverage', 'type': 'output', 'description': description}

    if saveinputs is True:
        maplayers['pga'] = {'grid': GDALGrid(PGA[:, :, 0], shakemap.getGeoDict()), 'label': 'PGA (g)', 'type': 'input', 'description': {'units': 'g', 'shakemap': shakedetail}}
        if 'PGV' in regressionmodel:
            maplayers['pgv'] = {'grid': GDALGrid(PGV[:, :, 0], shakemap.getGeoDict()), 'label': 'PGV (cm/s)', 'type': 'input', 'description': {'units': 'cm/s', 'shakemap': shakedetail}}
        maplayers['minFS'] = {'grid': GDALGrid(np.min(FS, axis=2), shakemap.getGeoDict()), 'label': 'Min Factor of Safety', 'type': 'input', 'description': {'units': 'unitless'}}
        maplayers['max slope'] = {'grid': GDALGrid(slopestack[:, :, -1], shakemap.getGeoDict()), 'label': 'Maximum slope ($^\circ$)', 'type': 'input', 'description': {'units': 'degrees', 'name': slopesref, 'longref': slopelref}}
        maplayers['cohesion'] = {'grid': GDALGrid(cohesion[:, :, 0], shakemap.getGeoDict()), 'label': 'Cohesion (kPa)', 'type': 'input', 'description': {'units': 'kPa (adjusted)', 'name': cohesionsref, 'longref': cohesionlref}}
        maplayers['friction angle'] = {'grid': GDALGrid(friction[:, :, 0], shakemap.getGeoDict()), 'label': 'Friction angle ($^\circ$)', 'type': 'input', 'description': {'units': 'degrees', 'name': frictionsref, 'longref': frictionlref}}

    return maplayers


def Saade2016():
    """
    Limit equilibrium approach combining mohr-coulomb for shallower slopes and GSI for steeper. No assumption of failure depth required (this could be moved to a different module since it doesn't exactly use Newmark)
    """
    print('Saade2016 not implemented yet')


def multiNewmark():
    """
    Run Classic or Godt model for set of different thicknesses, cell sizes, and unit weights to simulate different landslide sizes
    (borrow from )
    """
    print('multiNewmark not implemented yet')


def J_PGA(Ac, PGA):
    """
    PGA-based Newmark Displacement model from Jibson (2007), equation 6

    :param Ac: NxM Array of critical accelerations in units of g
    :type Ac: numpy Array
    :param PGA: NxM Array of PGA values in units of g
    :type PGA: numpy Array

    :returns Dn: NxM array of Newmark displacements in cm
    """
    np.seterr(invalid='ignore')  # Ignore errors so still runs when Ac > PGA, just leaves nan instead of crashing
    C1 = 0.215  # additive constant in newmark displacement calculation
    C2 = 2.341  # first exponential constant
    C3 = -1.438  # second exponential constant
    Dn = np.exp(C1 + np.log(((1-Ac/PGA)**C2)*(Ac/PGA)**C3))
    Dn[np.isnan(Dn)] = 0.
    return Dn


def J_PGA_M(Ac, PGA, M):
    """
    PGA-and M- based Newmark Displacement model from Jibson (2007), equation 7

    :param Ac: NxM Array of critical accelerations in units of g
    :type Ac: numpy Array
    :param PGA: NxM Array of PGA values in units of g
    :type PGA: numpy Array
    :param M: Magnitude
    :type M: float

    :returns Dn: NxM array of Newmark displacements in cm
    """
    np.seterr(invalid='ignore')  # Ignore errors so still runs when Ac > PGA, just leaves nan instead of crashing
    C1 = -2.71  # additive constant in newmark displacement calculation
    C2 = 2.335  # first exponential constant
    C3 = -1.478  # second exponential constant
    C4 = 0.424
    Dn = np.exp(C1 + np.log(((1-Ac/PGA)**C2)*(Ac/PGA)**C3) + C4*M)
    Dn[np.isnan(Dn)] = 0.
    return Dn


def RS_PGA_M(Ac, PGA, M):
    """
    PGA and M-based Newmark displacement model from Rathje and Saygili (2009)

    :param Ac: NxM Array of critical accelerations in units of g
    :type Ac: numpy Array
    :param PGA: NxM Array of PGA values in units of g
    :type PGA: numpy Array
    :param M: Magnitude
    :type M: float

    :returns Dn: NxM array of Newmark displacements in cm
    """
    np.seterr(invalid='ignore')
    C1 = 4.89
    C2 = -4.85
    C3 = -19.64
    C4 = 42.49
    C5 = -29.06
    C6 = 0.72
    C7 = 0.89
    Dn = np.exp(C1 + C2*(Ac/PGA) + C3*(Ac/PGA)**2 + C4*(Ac/PGA)**3 + C5*(Ac/PGA)**4 + C6*np.log(PGA)+C7*(M-6))  # Equation from Saygili and Rathje (2008)/Rathje and Saygili (2009)
    Dn[np.isnan(Dn)] = 0.
    return Dn


def RS_PGA_PGV(Ac, PGA, PGV):
    """
    PGA and PGV-based model from Saygili and Rathje (2008) - equation 6
    :param Ac: NxM Array of critical accelerations in units of g
    :type Ac: numpy Array
    :param PGA: NxM Array of PGA values in units of g
    :type PGA: numpy Array
    :param PGV: NxM Array of PGV values in units of cm/sec
    :type PGV: numpy array

    :returns Dn: NxM array of Newmark displacements in cm
    """
    np.seterr(invalid='ignore')
    C1 = -1.56
    C2 = -4.58
    C3 = -20.84
    C4 = 44.75
    C5 = -30.50
    C6 = -0.64
    C7 = 1.55
    Dn = np.exp(C1 + C2*(Ac/PGA) + C3*(Ac/PGA)**2 + C4*(Ac/PGA)**3 + C5*(Ac/PGA)**4 + C6*np.log(PGA)+C7*np.log(PGV))  # Equation from Saygili and Rathje (2008)/Rathje and Saygili (2009)
    Dn[np.isnan(Dn)] = 0.
    return Dn
