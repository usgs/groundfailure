#!/usr/bin/env python

"""
All Newmark based landslide mechanistic_models
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
    isURL = False
    try:
        fh = urllib2.urlopen(gridurl)
        isURL = True
    except:
        pass
    return isURL


def classic(shakefile, config, saveinputs=False, regressionmodel='J_PGA', threshold=False):
    pass


def godt2008(shakefile, config, saveinputs=False, regressionmodel='J_PGA'):
    """
    Applies the Godt et al. 2008 global method to a given ShakeMap. Factor of Safety calculated using infinite slope analysis assumuing dry conditions. Uses threshold newmark displacement and estimates 'probabilities' by doing the calculations for each slope quartile (which basically means probabilities are areal and 100 percent failure after threshold)
    INPUTS
    shakefile = url or filepath to shakemap xml file
    config = ConfigObj of config file
    OUTPUTS
    maps, .png and .pdf
    bounds = 'shakemap' will show entire area computed for shakemap, 'zoom' will only show areas predicted to be affected by landslides
    regressionmodel - what regression model to use
    'J_PGA' - PGA model from Jibson (2007)
    'J_PGA_M' - PGA and M model from Jibson (2007)
    'RS_PGA_M' - PGA and M model from Saygili and Rathje 2008 and Rathje and Saygili (2009)
    'RS_PGA_PGV' - PGA PGV model from Saygili and Rathje 2008 and Rathje and Saygili (2009)
    STILL NEED TO ADD 'all' - averages Dn from all three...have to figure out how to do this
    TO DO - add term to convert PGA and PGV to Ia and use other equations...
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
    cohesion = np.repeat(GDALGrid.load(cohesionfile, samplegeodict=shakemap.getGeoDict(), resample=True, preserve='shape', method='nearest').getData()[:, :, np.newaxis], 7, axis=2)
    cohesion[np.isnan(cohesion)] = nodata_cohesion
    friction = np.repeat(GDALGrid.load(frictionfile, samplegeodict=shakemap.getGeoDict(), resample=True, preserve='shape', method='nearest').getData()[:, :, np.newaxis], 7, axis=2)
    friction[np.isnan(friction)] = nodata_friction

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
        maplayers['max slope'] = {'grid': GDALGrid(slopestack[:, :, -1], shakemap.getGeoDict()), 'label': 'Maximum slope ($^\circ$)', 'type': 'input'}
        maplayers['cohesion'] = {'grid': GDALGrid(cohesion[:, :, 0], shakemap.getGeoDict()), 'label': 'Cohesion (kPa)', 'type': 'input'}
        maplayers['friction angle'] = {'grid': GDALGrid(friction[:, :, 0], shakemap.getGeoDict()), 'label': 'Friction angle ($^\circ$)', 'type': 'input'}

    return maplayers


def J_PGA(Ac, PGA):
    np.seterr(invalid='ignore')  # Ignore errors so still runs when Ac > PGA, just leaves nan instead of crashing
    C1 = 0.215  # additive constant in newmark displacement calculation
    C2 = 2.341  # first exponential constant
    C3 = -1.438  # second exponential constant
    Dn = np.exp(C1 + np.log(((1-Ac/PGA)**C2)*(Ac/PGA)**C3))
    Dn[np.isnan(Dn)] = 0.
    return Dn


def J_PGA_M(Ac, PGA, M):
    np.seterr(invalid='ignore')  # Ignore errors so still runs when Ac > PGA, just leaves nan instead of crashing
    C1 = -2.71  # additive constant in newmark displacement calculation
    C2 = 2.335  # first exponential constant
    C3 = -1.478  # second exponential constant
    C4 = 0.424
    Dn = np.exp(C1 + np.log(((1-Ac/PGA)**C2)*(Ac/PGA)**C3) + C4*M)
    Dn[np.isnan(Dn)] = 0.
    return Dn


def RS_PGA_M(Ac, PGA, M):
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
