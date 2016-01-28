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


def classic(shakefile, config, saveinputs=False, regressionmodel='J_PGA', probtype='jibson2000'):
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
    :type saveinputs: Binary
    :param regressionmodel:
        Newmark displacement regression model to use
        'J_PGA' (default) - PGA-based model from Jibson (2007) - equation 6
        'J_PGA_M' - PGA and M-based model from Jibson (2007) - equation 7
        'RS_PGA_M' - PGA and M-based model from from Rathje and Saygili (2009)
        'RS_PGA_PGV' - PGA and PGV-based model from Saygili and Rathje (2008) - equation 6
    :type regressionmodel: string
    :param probtype: Method used to estimate probability. Entering 'jibson2000' uses equation 5 from Jibson et al. (2000) to estimate probability from Newmark displacement. 'threshold' uses a specified threshold of Newmark displacement (defined in config file) and assumes anything greather than this threshold fails
    :type probtype: string

    :returns maplayers:  Dictionary containing output and input layers (if saveinputs=True) along with metadata formatted like maplayers['layer name']={'grid': mapio grid2D object, 'label': 'label for colorbar and top line of subtitle', 'type': 'output or input to model', 'description': 'detailed description of layer for subtitle, potentially including source information'}
    :type maplayers: OrderedDict

    :raises NameError: when unable to parse the config correctly (probably a formatting issue in the configfile) or when unable to find the shakefile (Shakemap URL or filepath) - these cause program to end
    :raises NameError: when probtype does not match a predifined probability type, will cause to default to 'jibson2000'

    """
    # Parse config - should make it so it uses defaults if any are missing...
    try:
        slopefile = config['mechanistic_models']['classic_newmark']['layers']['slopefile']
        cohesionfile = config['mechanistic_models']['classic_newmark']['layers']['cohesionfile']
        frictionfile = config['mechanistic_models']['classic_newmark']['layers']['frictionfile']
        thick = float(config['mechanistic_models']['classic_newmark']['values']['thick'])
        uwt = float(config['mechanistic_models']['classic_newmark']['values']['uwt'])
        nodata_cohesion = float(config['mechanistic_models']['classic_newmark']['values']['nodata_cohesion'])
        nodata_friction = float(config['mechanistic_models']['classic_newmark']['values']['nodata_friction'])
        try:
            dnthresh = float(config['mechanistic_models']['classic_newmark']['values']['dnthresh'])
        except:
            if probtype == 'threshold':
                dnthresh = 5.
                print('Unable to find dnthresh in config, using 5cm')
        fsthresh = float(config['mechanistic_models']['classic_newmark']['values']['fsthresh'])
        acthresh = float(config['mechanistic_models']['classic_newmark']['values']['acthresh'])
    except Exception as e:
        raise NameError('Could not parse configfile, %s' % e)
        return

    # Get boundaries from shakemap
    shkgdict = ShakeGrid.load(shakefile).getGeoDict()
    gdict = GDALGrid.getBoundsWithin(slopefile, shkgdict)

    # Load in slope file
    slopegrid = GDALGrid.load(slopefile, samplegeodict=gdict, resample=False, preserve='dims')  # Need to divide values by 100 if using slope_max.bil
    slope = slopegrid.getData()/100.
    # Change any zero slopes to a very small number to avoid dividing by zero later
    slope[slope == 0] = 0.0000001

    # Load in shakefile
    if not os.path.isfile(shakefile):
        if isURL(shakefile):
            shakefile = getGridURL(shakefile)  # returns a file object
        else:
            raise NameError('Could not find "%s" as a file or a valid url' % (shakefile))
            return

    # Load in shakemap, resample to slope file (this will be important when go to higher res)
    shakemap = ShakeGrid.load(shakefile, samplegeodict=slopegrid.getGeoDict(), resample=True, preserve='shape', method='linear')
    M = shakemap.getEventDict()['magnitude']

    # Read in the cohesion and friction files, resampled to slope grid
    cohesion = GDALGrid.load(cohesionfile, samplegeodict=slopegrid.getGeoDict(), resample=True, preserve='shape', method='nearest').getData()/10.
    cohesion[np.isnan(cohesion)] = nodata_cohesion
    friction = GDALGrid.load(frictionfile, samplegeodict=shakemap.getGeoDict(), resample=True, preserve='shape', method='nearest').getData()
    friction[np.isnan(friction)] = nodata_friction

    # See if there is a water table depth file and read it in if there is
    try:
        waterfile = config['mechanistic_models']['classic_newmark']['layers']['watertable']
        watertable = GDALGrid.load(waterfile, samplegeodict=slopegrid.getGeoDict(), resample=True, preserve='shape', method='linear').getData()  # Needs to be in meters!
        uwtw = float(config['mechanistic_models']['classic_newmark']['values']['uwtw'])
    except:
        watertable = None

    # Factor of safety
    if watertable is not None:
        FS = cohesion/(uwt*thick*np.sin(slope*(np.pi/180.))) + np.tan(friction*(np.pi/180.))/np.tan(slope*(np.pi/180.)) - (watertable*uwtw*np.tan(friction*(np.pi/180.)))/(uwt*np.tan(slope*(np.pi/180.)))
    else:
        FS = cohesion/(uwt*thick*np.sin(slope*(np.pi/180.))) + np.tan(friction*(np.pi/180.))/np.tan(slope*(np.pi/180.))
    FS[FS < fsthresh] = fsthresh

    # Compute critical acceleration, in g
    Ac = (FS-1)*np.sin(slope*(np.pi/180.)).astype(float)  # This gives ac in g, equations that multiply by g give ac in m/s2
    Ac[Ac < acthresh] = acthresh

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

    if probtype.lower() in 'jibson2000':
        PROB = 0.335*(1-np.exp(-0.048*Dn)**1.565)
    elif probtype.lower() in 'threshold':
        PROB = Dn.copy()
        PROB[PROB <= dnthresh] = 0
        PROB[PROB > dnthresh] = 1
    else:
        raise NameError('invalid probtype, assuming jibson2000')
        PROB = 0.335*(1-np.exp(-0.048*Dn)**1.565)

    # Turn output and inputs into into grids and put in mapLayers dictionary
    maplayers = collections.OrderedDict()

    maplayers['probability'] = {'grid': GDALGrid(PROB, slopegrid.getGeoDict()), 'label': 'Probability of > one landslide', 'type': 'output'}

    if saveinputs is True:
        maplayers['pga'] = {'grid': GDALGrid(PGA, slopegrid.getGeoDict()), 'label': 'PGA (g)', 'type': 'input'}
        maplayers['FS'] = {'grid': GDALGrid(FS, slopegrid.getGeoDict()), 'label': 'Factor of Safety', 'type': 'input'}
        #maplayers['Ac'] = {'grid': GDALGrid(FS, slopegrid.getGeoDict()), 'label': 'Critical acceleration (g)', 'type': 'input'}
        maplayers['slope'] = {'grid': GDALGrid(slope, slopegrid.getGeoDict()), 'label': 'Slope ($^\circ$)', 'type': 'input'}
        maplayers['cohesion'] = {'grid': GDALGrid(cohesion, slopegrid.getGeoDict()), 'label': 'Cohesion (kPa)', 'type': 'input'}
        maplayers['friction angle'] = {'grid': GDALGrid(friction, slopegrid.getGeoDict()), 'label': 'Friction angle ($^\circ$)', 'type': 'input'}
        if watertable is not None:
            maplayers['water depth'] = {'grid': GDALGrid(watertable, slopegrid.getGeoDict()), 'label': 'Water table depth (m)', 'type': 'input'}

    return maplayers


def godt2008(shakefile, config, saveinputs=False, regressionmodel='J_PGA'):
    """ This function runs the Godt et al. (2008) global method for a given ShakeMap. The Factor of Safety
    is calculated using infinite slope analysis assumuing dry conditions. The method uses threshold newmark
    displacement and estimates areal/spatial 'probabilities' by doing the calculations for each slope quantile
    TO DO - add 'all' - averages Dn from all four equations, add term to convert PGA and PGV to Ia and use other equations, add Ambraseys and Menu (1988) option

    :param shakefile: url or filepath to shakemap xml file
    :type shakefile: string
    :param config: ConfigObj of config file containing inputs required for running the model
    :type config: ConfigObj
    :param saveinputs: Whether or not to return the model input layers, False (defeault) returns only the model output (one layer)
    :type saveinputs: Binary
    :param regressionmodel:
        Newmark displacement regression model to use
        'J_PGA' (default) - PGA-based model from Jibson (2007) - equation 6
        'J_PGA_M' - PGA and M-based model from Jibson (2007) - equation 7
        'RS_PGA_M' - PGA and M-based model from from Rathje and Saygili (2009)
        'RS_PGA_PGV' - PGA and PGV-based model from Saygili and Rathje (2008) - equation 6
    :type regressionmodel: string
    :param probtype: Method used to estimate probability. Entering 'jibson2000' uses equation 5 from Jibson et al. (2000) to estimate probability from Newmark displacement. 'threshold' uses a specified threshold of Newmark displacement (defined in config file) and assumes anything greather than this threshold fails
    :type probtype: string

    :returns maplayers:  Dictionary containing output and input layers (if saveinputs=True) along with metadata formatted like maplayers['layer name']={'grid': mapio grid2D object, 'label': 'label for colorbar and top line of subtitle', 'type': 'output or input to model', 'description': 'detailed description of layer for subtitle, potentially including source information'}
    :type maplayers: OrderedDict

    :raises NameError: when unable to parse the config correctly (probably a formatting issue in the configfile) or when unable to find the shakefile (Shakemap URL or filepath) - these cause program to end
    """

    # Parse config
    try:
        slopefilepath = config['mechanistic_models']['godt_2008']['layers']['slopefilepath']
        cohesionfile = config['mechanistic_models']['godt_2008']['layers']['cohesionfile']
        frictionfile = config['mechanistic_models']['godt_2008']['layers']['frictionfile']
        thick = float(config['mechanistic_models']['godt_2008']['values']['thick'])
        uwt = float(config['mechanistic_models']['godt_2008']['values']['uwt'])
        nodata_cohesion = float(config['mechanistic_models']['godt_2008']['values']['nodata_cohesion'])
        nodata_friction = float(config['mechanistic_models']['godt_2008']['values']['nodata_friction'])
        dnthresh = float(config['mechanistic_models']['godt_2008']['values']['dnthresh'])
        fsthresh = float(config['mechanistic_models']['godt_2008']['values']['fsthresh'])
        acthresh = float(config['mechanistic_models']['godt_2008']['values']['acthresh'])
    except Exception as e:
        raise NameError('Could not parse configfile, %s' % e)
        return

    # Load in shakefile
    if not os.path.isfile(shakefile):
        if isURL(shakefile):
            shakefile = getGridURL(shakefile)  # returns a file object
        else:
            raise NameError('Could not find "%s" as a file or a valid url' % (shakefile))
            return

    shakemap = ShakeGrid.load(shakefile)
    M = shakemap.getEventDict()['magnitude']

    # Read in all the slope files, divide all by 100 to get to slope in degrees (because input files are multiplied by 100.)
    slopes = []
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope_min.bil'), samplegeodict=shakemap.getGeoDict(), resample=True, preserve='shape', method='linear').getData()/100.)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope10.bil'), samplegeodict=shakemap.getGeoDict(), resample=True, preserve='shape', method='linear').getData()/100.)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope30.bil'), samplegeodict=shakemap.getGeoDict(), resample=True, preserve='shape', method='linear').getData()/100.)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope50.bil'), samplegeodict=shakemap.getGeoDict(), resample=True, preserve='shape', method='linear').getData()/100.)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope70.bil'), samplegeodict=shakemap.getGeoDict(), resample=True, preserve='shape', method='linear').getData()/100.)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope90.bil'), samplegeodict=shakemap.getGeoDict(), resample=True, preserve='shape', method='linear').getData()/100.)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope_max.bil'), samplegeodict=shakemap.getGeoDict(), resample=True, preserve='shape', method='linear').getData()/100.)
    slopestack = np.dstack(slopes)

    # Change any zero slopes to a very small number to avoid dividing by zero later
    slopestack[slopestack == 0] = 0.0000001

    # Read in the cohesion and friction files and duplicate layers so they are same shape as slope structure
    #import pdb; pdb.set_trace()
    cohesion = np.repeat(GDALGrid.load(cohesionfile, samplegeodict=shakemap.getGeoDict(), resample=True, preserve='shape', method='nearest').getData()[:, :, np.newaxis]/10., 7, axis=2)
    cohesion[cohesion == -999.9] = nodata_cohesion
    cohesion[cohesion == 0] = nodata_cohesion
    friction = np.repeat(GDALGrid.load(frictionfile, samplegeodict=shakemap.getGeoDict(), resample=True, preserve='shape', method='nearest').getData()[:, :, np.newaxis], 7, axis=2)
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

    maplayers['probability'] = {'grid': GDALGrid(PROB, shakemap.getGeoDict()), 'label': 'Areal probability', 'type': 'output'}

    if saveinputs is True:
        maplayers['pga'] = {'grid': GDALGrid(PGA[:, :, 0], shakemap.getGeoDict()), 'label': 'PGA (g)', 'type': 'input'}
        maplayers['minFS'] = {'grid': GDALGrid(np.min(FS, axis=2), shakemap.getGeoDict()), 'label': 'Min Factor of Safety', 'type': 'input'}
        maplayers['max slope'] = {'grid': GDALGrid(slopestack[:, :, -1], shakemap.getGeoDict()), 'label': 'Maximum slope ($^\circ$)', 'type': 'input'}
        maplayers['cohesion'] = {'grid': GDALGrid(cohesion[:, :, 0], shakemap.getGeoDict()), 'label': 'Cohesion (kPa)', 'type': 'input'}
        maplayers['friction angle'] = {'grid': GDALGrid(friction[:, :, 0], shakemap.getGeoDict()), 'label': 'Friction angle ($^\circ$)', 'type': 'input'}

    return maplayers


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
