#!/usr/bin/env python

"""
godt2008.py applies the Godt et al. 2008 method to a given ShakeMap as written at 1km resolution, applicable globally. Factor of Safety calculated using infinite slope analysis for 2.4m thickness, no water. This applies Jibson (2007) regression equation relating Newmark displacement to critical acceleration and PGA only. Uses threshold displacement of 5cm and estimates 'probabilities' by doing the calculations for each slope quartile (which basically means probabilities are areal and 100 percent failure after threshold)
"""
#stdlib imports
from configobj import ConfigObj
import os.path
import warnings
import urllib2
import tempfile
import matplotlib.cm as cm

#turn off all warnings...
warnings.filterwarnings('ignore')

#local imports
from mapio.shake import ShakeGrid
from mapio.gdal import GDALGrid
from secondary.mapnew import makeMap#, saveMap

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
        raise IOError, 'Could not retrieve data from %s' % gridurl
    return gridfile


def isURL(gridurl):
    isURL = False
    try:
        fh = urllib2.urlopen(gridurl)
        isURL = True
    except:
        pass
    return isURL


def runmodel(shakefile, configfile, maproads=True, mapcities=True, mapProb=True, mapshaking=False, mapcohesion=False, mapfriction=False, mapmaxslope=False, mapmedslope=False, saveProb=True, saveinputs=False, bounds='shakemap'):
    """
    main()
    INPUTS
    shakefile = url or filepath to shakemap xml file
    configfile = full path to config file containing all necessary info
    OUTPUTS
    maps, .png and .pdf
    bounds = 'shakemap' will show entire area computed for shakemap, 'zoom' will only show areas predicted to be affected by landslides
    """
    # Parse config file
    try:
        config = ConfigObj(configfile)['mechanistic_models']['godt_2008']
        modelname = 'godt2008'
        slopefilepath = config['layers']['slopefilepath']
        cohesionfile = config['layers']['cohesionfile']
        frictionfile = config['layers']['frictionfile']
        thick = float(config['values']['thick'])
        uwt = float(config['values']['uwt'])
        nodata_cohesion = float(config['values']['nodata_cohesion'])
        nodata_friction = float(config['values']['nodata_friction'])
        dnthresh = float(config['values']['dnthresh'])
        fsthresh = float(config['values']['fsthresh'])
        acthresh = float(config['values']['acthresh'])
    except Exception as e:
        print e
        print 'Could not parse configfile, missing necessary values or incorrect structure'
        return

    # Load in shakefile
    if not os.path.isfile(shakefile):
        if isURL(shakefile):
            shakefile = getGridURL(shakefile)  # returns a file object
        else:
            print 'Could not find "%s" as a file or a url.  Returning.' % (shakefile)
            return

    shakemap = ShakeGrid.load(shakefile)
    #
    shakeevent = shakemap.getEventDict()
    shakeshake = shakemap.getShakeDict()
    edict = {'mag': shakeevent['magnitude'],
             'time': shakeevent['event_timestamp'],
             'loc': shakeevent['event_description'],
             'epicenter': (shakeevent['lat'], shakeevent['lon']),
             'version': int(shakeshake['shakemap_version']),
             'eventid': shakeshake['event_id']}
    network = shakeshake['shakemap_originator']
    eventcode = shakeshake['shakemap_id']
    if eventcode.startswith(network):
        eventid = eventcode
    else:
        eventid = network + eventcode
    edict['eventid'] = eventid

    # Read in all the slope files, divide all by 100 to get to slope in degrees (input files are multiplied by 100.)
    slopes = []
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope_min.bil'), samplegeodict=shakemap.getGeoDict(), resample=True, method='linear').getData()/100.)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope10.bil'), samplegeodict=shakemap.getGeoDict(), resample=True, method='linear').getData()/100.)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope30.bil'), samplegeodict=shakemap.getGeoDict(), resample=True, method='linear').getData()/100.)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope50.bil'), samplegeodict=shakemap.getGeoDict(), resample=True, method='linear').getData()/100.)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope70.bil'), samplegeodict=shakemap.getGeoDict(), resample=True, method='linear').getData()/100.)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope90.bil'), samplegeodict=shakemap.getGeoDict(), resample=True, method='linear').getData()/100.)
    slopes.append(GDALGrid.load(os.path.join(slopefilepath, 'slope_max.bil'), samplegeodict=shakemap.getGeoDict(), resample=True, method='linear').getData()/100.)
    slopestack = np.dstack(slopes)

    # Change any zero slopes to a very small number to avoid dividing by zero later
    slopestack[slopestack == 0] = 0.00001

    # Read in the cohesion and friction files and duplicate layers so they are same shape as slope structure
    cohesion = np.repeat(GDALGrid.load(cohesionfile, samplegeodict=shakemap.getGeoDict(), resample=True, method='nearest').getData()[:, :, np.newaxis], 7, axis=2)
    cohesion[np.isnan(cohesion)] = nodata_cohesion
    friction = np.repeat(GDALGrid.load(frictionfile, samplegeodict=shakemap.getGeoDict(), resample=True, method='nearest').getData()[:, :, np.newaxis], 7, axis=2)
    friction[np.isnan(friction)] = nodata_friction

    # Do the calculations using Jibson (2007) PGA only model for Dn
    FS = cohesion/(uwt*thick*np.sin(slopestack*(np.pi/180.))) + np.tan(friction*(np.pi/180.))/np.tan(slopestack*(np.pi/180.))
    FS[FS < fsthresh] = fsthresh

    # Compute critical acceleration, in g
    Ac = (FS-1)*np.sin(slopestack*(np.pi/180.)).astype(float)  # This gives ac in g, equations that multiply by g give ac in m/s2
    Ac[Ac < acthresh] = acthresh

    # Get PGA in g (PGA is %g in ShakeMap)
    PGA = np.repeat(shakemap.getLayer('pga').getData()[:, :, np.newaxis]/100., 7, axis=2).astype(float)

    # Estimate Newmark displacement
    C1 = 0.215  # additive constant in newmark displacement calculation
    C2 = 2.341  # first exponential constant
    C3 = -1.438  # second exponential constant
    np.seterr(invalid='ignore')  # Ignore errors so still runs when Ac > PGA, just leaves nan instead of crashing
    Dn = np.exp(C1 + np.log(((1-Ac/PGA)**C2)*(Ac/PGA)**C3))
    Dn[np.isnan(Dn)] = 0.
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

    lsgrid = GDALGrid(PROB, shakemap.getGeoDict())

    if bounds == 'shakemap':
        boundaries1 = shakemap.getGeoDict()
    else:
        # get lat lons of areas affected and add, if no areas affected, switch to shakemap boundaries
        xmin, xmax, ymin, ymax = shakemap.getBounds()
        lons = np.arange(xmin, xmax, shakemap.getGeoDict()['xdim'])
        lons = lons[:shakemap.getGeoDict()['ncols']]  # make sure right length
        lats = np.arange(ymax, ymin, -shakemap.getGeoDict()['ydim'])  # backwards so it plots right
        lats = lats[:shakemap.getGeoDict()['nrows']]
        llons, llats = np.meshgrid(lons, lats)  # make meshgrid
        llons1 = llons[PROB > 0]
        llats1 = llats[PROB > 0]
        boundaries1 = {}
        boundaries1['xmin'] = llons1.min()-0.2*(llons1.max()-llons1.min())
        boundaries1['xmax'] = llons1.max()+0.2*(llons1.max()-llons1.min())
        boundaries1['ymin'] = llats1.min()-0.2*(llats1.max()-llats1.min())
        boundaries1['ymax'] = llats1.max()+0.2*(llats1.max()-llats1.min())

    # Turn other into into grids and put in a dictionary if mapping of them is desired
    flag = 0.
    maplayers = {}
    plotorder = []
    colormaps = []
    if mapProb is True:
        maplayers['areal probability'] = lsgrid
        plotorder.append('areal probability')
        colormaps.append(cm.jet)  # cm.autumn_r)
        # Always make a single panel version of this one
        makeMap(maplayers, edict, configfile, boundaries=boundaries1, modelname=modelname+'_prob', maproads=maproads, mapcities=mapcities, colormaps=colormaps)
        #saveMap(maplayers, edict, configfile, modelname+'_prob')
    if mapshaking is True:
        maplayers['pga (g)'] = GDALGrid(PGA[:, :, 0], shakemap.getGeoDict())
        flag = 1.
        plotorder.append('pga (g)')
        colormaps.append(cm.jet)
    if mapmaxslope is True:
        maplayers['max slope'] = GDALGrid(slopestack[:, :, -1], shakemap.getGeoDict())
        flag = 1.
        plotorder.append('max slope')
        colormaps.append(cm.gnuplot2)
    if mapmedslope is True:
        maplayers['median slope'] = GDALGrid(slopestack[:, :, 3], shakemap.getGeoDict())
        flag = 1.
        plotorder.append('median slope')
        colormaps.append(cm.gnuplot2)
    if mapcohesion is True:
        maplayers['cohesion (kPa)'] = GDALGrid(cohesion[:, :, 0], shakemap.getGeoDict())
        flag = 1.
        plotorder.append('cohesion (kPa)')
        colormaps.append(cm.jet)
    if mapfriction is True:
        maplayers['friction angle'] = GDALGrid(friction[:, :, 0], shakemap.getGeoDict())
        flag = 1.
        plotorder.append('friction angle')
        colormaps.append(cm.jet)

    if flag == 1.:  # Only call multipanel map if more than probabilities should be mapped
        makeMap(maplayers, edict, configfile, modelname=modelname, maproads=maproads, mapcities=mapcities, colormaps=colormaps, plotorder=plotorder, boundaries=boundaries1)
        #saveMap(maplayers, edict, modelname)

    return maplayers, edict
