#!/usr/bin/env python

#stdlib imports
import os.path
import gc
import math
import glob
from mpl_toolkits.basemap import maskoceans
import copy
import datetime
import matplotlib as mpl
from matplotlib.colors import LightSource
#from matplotlib.colorbar import ColorbarBase

#third party imports
import matplotlib.cm as cm
import numpy as np
import matplotlib.pyplot as plt
import fiona
from shapely.geometry import mapping, shape
from shapely.geometry import Polygon as PolygonSH
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import cm as cm2
from matplotlib.patches import Polygon, Rectangle
#from matplotlib.collections import PatchCollection
from skimage.measure import block_reduce
import collections
from descartes import PolygonPatch

#local imports
from mapio.gmt import GMTGrid
from mapio.gdal import GDALGrid
from mapio.geodict import GeoDict
from mapio.grid2d import Grid2D
from neicmap.city import PagerCity
from neicutil.text import ceilToNearest, floorToNearest, roundToNearest
#from mapio.basemapcity import BasemapCities

# Make fonts readable and recognizable by illustrator
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = ['Arial', 'Bitstream Vera Serif', 'sans-serif']


def parseMapConfig(config):
    # Parse config object
    # ADD PLOTORDER TO CONFIG? OTHER THINGS LIKE COLORMAPS?
    topofile = None
    roadfolder = None
    cityfile = None
    roadcolor = '6E6E6E'
    countrycolor = '177F10'
    watercolor = 'B8EEFF'
    ALPHA = 0.7
    outputdir = None
    oceanfile = None

    try:
        config1 = config['mapdata']
        if 'dem' in config1:
            topofile = config1['dem']['file']
            if os.path.exists(topofile) is False:
                print('DEM not valid - hillshade will not be possible\n')
        if 'ocean' in config1:
            oceanfile = config1['ocean']['file']
            try:
                oceanref = config1['ocean']['shortref']
            except:
                oceanref = 'unknown'
        if 'roads' in config1:
            roadfolder = config1['roads']['folder']
            if os.path.exists(roadfolder) is False:
                print('roadfolder not valid - roads will not be displayed\n')
                roadfolder = None
            try:
                roadref = config1['roads']['shortref']
            except:
                roadref = 'unknown'
        if 'cities' in config1:
            cityfile = config1['cities']['file']
            try:
                cityref = config1['cities']['shortref']
            except:
                cityref = 'unknown'
            if os.path.exists(cityfile):
                try:
                    PagerCity(cityfile)
                    #BasemapCities.loadFromGeoNames(cityfile=cityfile)
                except Exception as e:
                    print e
                    print('cities file not valid - cities will not be displayed\n')
                    cityfile = None
            else:
                print('cities file not valid - cities will not be displayed\n')
                cityfile = False
        if 'roadcolor' in config1['colors']:
            roadcolor = config1['colors']['roadcolor']
        if 'countrycolor' in config1['colors']:
            countrycolor = config1['colors']['countrycolor']
        if 'watercolor' in config1['colors']:
            watercolor = config1['colors']['watercolor']
        if 'alpha' in config1['colors']:
            ALPHA = float(config1['colors']['alpha'])
        try:
            outputdir = config['output']['folder']
        except:
            outputdir = None
    except Exception as e:
        print('%s - mapdata missing from or misformatted in config' % e)

    countrycolor = '#'+countrycolor
    watercolor = '#'+watercolor
    roadcolor = '#'+roadcolor

    mapin = {'topofile': topofile, 'roadfolder': roadfolder, 'cityfile': cityfile, 'roadcolor': roadcolor, 'countrycolor': countrycolor, 'watercolor': watercolor, 'ALPHA': ALPHA, 'outputdir': outputdir, 'roadref': roadref, 'cityref': cityref, 'oceanfile': oceanfile, 'oceanref': oceanref}

    return mapin


def parseConfigLayers(maplayers, config):
    """
    Parse things that need to coodinate with each layer (like lims, logscale, colormaps etc.) from config file, in right order - takes orders from maplayers
    # TO DO, ADD ABILITY TO INTERPRET CUSTOM COLOR MAPS
    """
    # get all key names, create a plotorder list in case maplayers is not an ordered dict, making sure that anything called 'model' is first
    keys = maplayers.keys()
    plotorder = []

    try:
        limits = config['mapdata']['lims']
        lims = []
    except:
        lims = None
        limits = None

    try:
        colors = config['mapdata']['colors']
        colormaps = []
    except:
        colormaps = None
        colors = None

    try:
        logs = config['mapdata']['logscale']
        logscale = []
    except:
        logscale = False
        logs = None

    try:
        masks = config['mapdata']['maskthresholds']
        maskthreshes = []
    except:
        maskthreshes = None
        masks = None

    try:
        default = config['mapdata']['colors']['default']
        default = eval(default)
    except:
        default = None

    for i, key in enumerate(keys):
        plotorder += [key]
        if limits is not None:
            found = False
            for l in limits:
                getlim = None
                if l in key:
                    if type(limits[l]) is list:
                        getlim = np.array(limits[l]).astype(np.float)
                    else:
                        try:
                            getlim = eval(limits[l])
                        except:
                            getlim = None
                    lims.append(getlim)
                    found = True
            if not found:
                lims.append(None)

        if colors is not None:
            found = False
            for c in colors:
                #colorobject = default
                if c in key:
                    getcol = colors[c]
                    colorobject = eval(getcol)
                    if colorobject is None:
                        colorobject = default
                    colormaps.append(colorobject)
                    found = True
            if not found:
                colormaps.append(default)

        if logs is not None:
            found = False
            for g in logs:
                getlog = False
                if g in key:
                    if logs[g].lower() == 'true':
                        getlog = True
                    logscale.append(getlog)
                    found = True
            if not found:
                logscale.append(False)

        if masks is not None:
            found = False
            for m in masks:
                #getmask = None
                if m in key:
                    getmask = eval(masks[m])
                    maskthreshes.append(getmask)
                    found = True
            if not found:
                maskthreshes.append(None)

    # Reorder everything so model is first, if it's not already
    if plotorder[0] != 'model':
        indx = [idx for idx, key in enumerate(plotorder) if key == 'model']
        if len(indx) == 1:
            indx = indx[0]
            firstpo = plotorder.pop(indx)
            plotorder = [firstpo] + plotorder
            firstlog = logscale.pop(indx)
            logscale = [firstlog] + logscale
            firstlim = lims.pop(indx)
            lims = [firstlim] + lims
            firstcol = colormaps.pop(indx)
            colormaps = [firstcol] + colormaps

    return plotorder, logscale, lims, colormaps, maskthreshes


def modelMap(grids, edict=None, suptitle=None, inventory_shapefile=None, plotorder=None, maskthreshes=None, colormaps=None, boundaries=None, zthresh=0, scaletype='continuous', lims=None, logscale=False, ALPHA=0.7, maproads=True, mapcities=True, isScenario=False, roadfolder=None, topofile=None, cityfile=None, oceanfile=None, roadcolor='#6E6E6E', watercolor='#B8EEFF', countrycolor='#177F10', outputdir=None, savepdf=True, savepng=True, showplots=False, roadref='unknown', cityref='unknown', oceanref='unknown', printparam=False, ds=True, dstype='mean', upsample=False):

    """
    This function creates maps of mapio grid layers (e.g. liquefaction or landslide models with their input layers)
    All grids must use the same bounds
    TO DO change so that all input layers do not have to have the same bounds, test plotting multiple probability layers, and add option so that if PDF and PNG aren't output, opens plot on screen using plt.show()

    :param grids: Dictionary of N layers and metadata formatted like maplayers['layer name']={'grid': mapio grid2D object, 'label': 'label for colorbar and top line of subtitle', 'type': 'output or input to model', 'description': 'detailed description of layer for subtitle'}. Layer names must be unique.
    :type name: Dictionary or Ordered dictionary - import collections; grids = collections.OrderedDict()
    :param edict: Optional event dictionary for ShakeMap, used for title and naming output folder
    :type edit: Shakemap Event Dictionary
    :param suptitle: This will be displayed at the top of the plots and in the figure names
    :type suptitle: string
    :param plotorder: List of keys describing the order to plot the grids, if None and grids is an ordered dictionary, it will use the order of the dictionary, otherwise it will choose order which may be somewhat random but it will always put a probability grid first
    :type plotorder: list
    :param maskthreshes: N x 1 array or list of lower thresholds for masking corresponding to order in plotorder or order of OrderedDict if plotorder is None. If grids is not an ordered dict and plotorder is not specified, this will not work right. If None (default), nothing will be masked
    :param colormaps: List of strings of matplotlib colormaps (e.g. cm.autumn_r) corresponding to plotorder or order of dictionary if plotorder is None. The list can contain both strings and None e.g. colormaps = ['cm.autumn', None, None, 'cm.jet'] and None's will default to default colormap
    :param boundaries: None to show entire study area, 'zoom' to zoom in on the area of action (only works if there is a probability layer) using zthresh as a threshold, or a dictionary defining lats and lons in the form of boundaries.xmin = minlon, boundaries.xmax = maxlon, boundaries.ymin = min lat, boundaries.ymax = max lat
    :param zthresh: threshold for computing zooming bounds, only used if boundaries = 'zoom'
    :type zthresh: float
    :param scaletype: Type of scale for plotting, 'continuous' or 'binned' - will be reflected in colorbar
    :type scaletype: string
    :param lims: None or Nx1 list of tuples or numpy arrays corresponding to plotorder defining the limits for saturating the colorbar (vmin, vmax) if scaletype is continuous or the bins to use (clev) if scaletype if binned. The list can contain tuples, arrays, and Nones, e.g. lims = [(0., 10.), None, (0.1, 1.5), np.linspace(0., 1.5, 15)]. When None is specified, the program will estimate the limits, when an array is specified but the scaletype is continuous, vmin will be set to min(array) and vmax will be set to max(array)
    :param lims: None or Nx1 list of Trues and Falses corresponding to plotorder defining whether to use a linear or log scale (log10) for plotting the layer. This will be reflected in the labels
    :param ALPHA: Transparency for mapping, if there is a hillshade that will plot below each layer, it is recommended to set this to at least 0.7
    :type ALPHA: float
    :param maproads: Whether to show roads or not, default True, but requires that roadfile is specified and valid to work
    :type maproads: boolean
    :param mapcities: Whether to show cities or not, default True, but requires that cityfile is specified and valid to work
    :type mapcities: boolean
    :param isScenario: Whether this is a scenario (True) or a real event (False) (default False)
    :type isScenario: boolean
    :param roadfolder: Full file path to folder containing road shapefiles
    :type roadfolder: string
    :param topofile: Full file path to topography grid (GDAL compatible) - this is only needed to make a hillshade if a premade hillshade is not specified
    :type topofile: string
    :param cityfile: Full file path to Pager file containing city & population information
    :type cityfile: string
    :param roadcolor: Color to use for roads, if plotted, default #6E6E6E
    :type roadcolor: Hex color or other matplotlib compatible way of defining color
    :param watercolor: Color to use for oceans, lakes, and rivers, default #B8EEFF
    :type watercolor: Hex color or other matplotlib compatible way of defining color
    :param countrycolor: Color for country borders, default #177F10
    :type countrycolor: Hex color or other matplotlib compatible way of defining color
    :param outputdir: File path for outputting figures, if edict is defined, a subfolder based on the event id will be created in this folder. If None, will use current directory
    :param savepdf: True to save pdf figure, False to not
    :param savepng: True to save png figure, False to not
    :param ds: True to allow downsampling for display (necessary when arrays are quite large, False to not allow)
    :param dstype: What function to use in downsampling, options are 'min', 'max', 'median', or 'mean'
    :param upsample: True to upsample the layer to the DEM resolution for better looking hillshades

    :returns:  PDF and/or PNG of map
    :returns newgrids: Downsampled and trimmed version of input grids. If no modification was needed for plotting, this will be identical to grids but without the metadata
    """
    if suptitle is None:
        suptitle = ' '

    plt.ioff()

    defaultcolormap = cm.jet

    # Get output file location
    if outputdir is None:
        print('No output location given, using current directory for outputs\n')
        outputdir = os.getcwd()
    if edict is not None:
        outfolder = os.path.join(outputdir, edict['eventid'])
    else:
        outfolder = outputdir
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)

    # Get plotting order, if not specified
    if plotorder is None:
        plotorder = grids.keys()

    # Get boundaries to use for all plots
    cut = True
    if boundaries is None:
        cut = False
        keytemp = grids.keys()
        boundaries = grids[keytemp[0]]['grid'].getGeoDict()
    elif boundaries == 'zoom':
        # Find probability layer (will just take the maximum bounds if there is more than one)
        keytemp = grids.keys()
        key1 = [key for key in keytemp if 'model' in key.lower()]
        if len(key1) == 0:
            print('Could not find model layer to use for zoom, using default boundaries')
            keytemp = grids.keys()
            boundaries = grids[keytemp[0]]['grid'].getGeoDict()
        else:
            lonmax = -1.e10
            lonmin = 1.e10
            latmax = -1.e10
            latmin = 1.e10
            for key in key1:
                # get lat lons of areas affected and add, if no areas affected, switch to shakemap boundaries
                temp = grids[key]['grid']
                xmin, xmax, ymin, ymax = temp.getBounds()
                lons = np.linspace(xmin, xmax, temp.getGeoDict().nx)
                lats = np.linspace(ymax, ymin, temp.getGeoDict().ny)  # backwards so it plots right
                row, col = np.where(temp.getData() > float(zthresh))
                lonmin = lons[col].min()
                lonmax = lons[col].max()
                latmin = lats[row].min()
                latmax = lats[row].max()
                # llons, llats = np.meshgrid(lons, lats)  # make meshgrid
                # llons1 = llons[temp.getData() > float(zthresh)]
                # llats1 = llats[temp.getData() > float(zthresh)]
                # if llons1.min() < lonmin:
                #     lonmin = llons1.min()
                # if llons1.max() > lonmax:
                #     lonmax = llons1.max()
                # if llats1.min() < latmin:
                #     latmin = llats1.min()
                # if llats1.max() > latmax:
                #     latmax = llats1.max()
            boundaries1 = {'dx': 100, 'dy': 100., 'nx': 100., 'ny': 100}  # dummy fillers, only really care about bounds
            if xmin < lonmin-0.15*(lonmax-lonmin):
                boundaries1['xmin'] = lonmin-0.1*(lonmax-lonmin)
            else:
                boundaries1['xmin'] = xmin
            if xmax > lonmax+0.15*(lonmax-lonmin):
                boundaries1['xmax'] = lonmax+0.1*(lonmax-lonmin)
            else:
                boundaries1['xmax'] = xmax
            if ymin < latmin-0.15*(latmax-latmin):
                boundaries1['ymin'] = latmin-0.1*(latmax-latmin)
            else:
                boundaries1['ymin'] = ymin
            if ymax > latmax+0.15*(latmax-latmin):
                boundaries1['ymax'] = latmax+0.1*(latmax-latmin)
            else:
                boundaries1['ymax'] = ymax
            boundaries = GeoDict(boundaries1, adjust='res')
    else:
        # SEE IF BOUNDARIES ARE SAME AS BOUNDARIES OF LAYERS
        keytemp = grids.keys()
        tempgdict = grids[keytemp[0]]['grid'].getGeoDict()
        if np.abs(tempgdict.xmin-boundaries['xmin']) < 0.05 and np.abs(tempgdict.ymin-boundaries['ymin']) < 0.05 and np.abs(tempgdict.xmax-boundaries['xmax'] < 0.05 and np.abs(tempgdict.ymax - boundaries['ymax']) < 0.05):
            print('Input boundaries are almost the same as specified boundaries, no cutting needed')
            boundaries = tempgdict
            cut = False
        else:
            try:
                if boundaries['xmin'] > boundaries['xmax'] or boundaries['ymin'] > boundaries['ymax']:
                    print('Input boundaries are not usable, using default boundaries')
                    keytemp = grids.keys()
                    boundaries = grids[keytemp[0]]['grid'].getGeoDict()
                    cut = False
                else:
                    # Build dummy GeoDict
                    boundaries = GeoDict({'xmin': boundaries['xmin'], 'xmax': boundaries['xmax'], 'ymin': boundaries['ymin'], 'ymax': boundaries['ymax'], 'dx': 100., 'dy': 100., 'ny': 100., 'nx': 100.}, adjust='res')
            except:
                print('Input boundaries are not usable, using default boundaries')
                keytemp = grids.keys()
                boundaries = grids[keytemp[0]]['grid'].getGeoDict()
                cut = False

    # Pull out bounds for various uses
    bxmin, bxmax, bymin, bymax = boundaries.xmin, boundaries.xmax, boundaries.ymin, boundaries.ymax

    # Determine if need a single panel or multi-panel plot and if multi-panel, how many and how it will be arranged
    fig = plt.figure()
    numpanels = len(grids)
    if numpanels == 1:
        rowpan = 1
        colpan = 1
        # create the figure and axes instances.
        fig.set_figwidth(5)
    elif numpanels == 2 or numpanels == 4:
        rowpan = np.ceil(numpanels/2.)
        colpan = 2
        fig.set_figwidth(10)
    else:
        rowpan = np.ceil(numpanels/3.)
        colpan = 3
        fig.set_figwidth(15)
    if rowpan == 1:
        fig.set_figheight(rowpan*6.0)
    else:
        fig.set_figheight(rowpan*5.3)

    # Need to update naming to reflect the shakemap version once can get getHeaderData to work, add edict['version'] back into title, maybe shakemap id also?
    fontsizemain = 14.
    fontsizesub = 12.
    fontsizesmallest = 10.
    if rowpan == 1.:
        fontsizemain = 12.
        fontsizesub = 10.
        fontsizesmallest = 8.
    if edict is not None:
        if isScenario:
            title = edict['event_description']
        else:
            timestr = edict['event_timestamp'].strftime('%b %d %Y')
            title = 'M%.1f %s v%i - %s' % (edict['magnitude'], timestr, edict['version'], edict['event_description'])
        plt.suptitle(title+'\n'+suptitle, fontsize=fontsizemain)
    else:
        plt.suptitle(suptitle, fontsize=fontsizemain)

    clear_color = [0, 0, 0, 0.0]

    # Cut all of them and release extra memory

    xbuff = (bxmax-bxmin)/10.
    ybuff = (bymax-bymin)/10.
    cutxmin = bxmin-xbuff
    cutymin = bymin-ybuff
    cutxmax = bxmax+xbuff
    cutymax = bymax+ybuff
    if cut is True:
        newgrids = collections.OrderedDict()
        for k, layer in enumerate(plotorder):
            templayer = grids[layer]['grid']
            try:
                newgrids[layer] = {'grid': templayer.cut(cutxmin, cutxmax, cutymin, cutymax, align=True)}
            except Exception as e:
                print('Cutting failed, %s, continuing with full layers' % e)
                newgrids = grids
                continue
        del templayer
        gc.collect()
    else:
        newgrids = grids
    tempgdict = newgrids[grids.keys()[0]]['grid'].getGeoDict()

    # Upsample layers to same as topofile if desired for better looking hillshades
    if upsample is True and topofile is not None:
        try:
            topodict = GDALGrid.getFileGeoDict(topofile)
            if topodict.dx >= tempgdict.dx or topodict.dy >= tempgdict.dy:
                print('Upsampling not possible, resolution of results already smaller than DEM')
                pass
            else:
                tempgdict1 = GeoDict({'xmin': tempgdict.xmin-xbuff, 'ymin': tempgdict.ymin-ybuff, 'xmax': tempgdict.xmax+xbuff, 'ymax': tempgdict.ymax+ybuff, 'dx': topodict.dx, 'dy': topodict.dy, 'nx': topodict.nx, 'ny': topodict.ny}, adjust='res')
                tempgdict2 = tempgdict1.getBoundsWithin(tempgdict)
                for k, layer in enumerate(plotorder):
                    newgrids[layer]['grid'] = newgrids[layer]['grid'].subdivide(tempgdict2)
        except:
            print('Upsampling failed, continuing')

    # Downsample all of them for plotting, if needed, and replace them in grids (to save memory)
    tempgrid = newgrids[grids.keys()[0]]['grid']
    xsize = tempgrid.getGeoDict().nx
    ysize = tempgrid.getGeoDict().ny
    inchesx, inchesy = fig.get_size_inches()
    divx = int(np.round(xsize/(500.*inchesx)))
    divy = int(np.round(ysize/(500.*inchesy)))
    xmin, xmax, ymin, ymax = tempgrid.getBounds()
    gdict = tempgrid.getGeoDict()  # Will be replaced if downsampled
    del tempgrid
    gc.collect()

    if divx <= 1:
        divx = 1
    if divy <= 1:
        divy = 1
    if (divx > 1. or divy > 1.) and ds:
        if dstype == 'max':
            func = np.nanmax
        elif dstype == 'min':
            func = np.nanmin
        elif dstype == 'med':
            func = np.nanmedian
        else:
            func = np.nanmean
        for k, layer in enumerate(plotorder):
            layergrid = newgrids[layer]['grid']
            dat = block_reduce(layergrid.getData().copy(), block_size=(divy, divx), cval=float('nan'), func=func)
            if k == 0:
                lons = block_reduce(np.linspace(xmin, xmax, layergrid.getGeoDict().nx), block_size=(divx,), func=np.mean, cval=float('nan'))
                if math.isnan(lons[-1]):
                    lons[-1] = lons[-2] + (lons[1]-lons[0])
                lats = block_reduce(np.linspace(ymax, ymin, layergrid.getGeoDict().ny), block_size=(divy,), func=np.mean, cval=float('nan'))
                if math.isnan(lats[-1]):
                    lats[-1] = lats[-2] + (lats[1]-lats[0])
                gdict = GeoDict({'xmin': lons.min(), 'xmax': lons.max(), 'ymin': lats.min(), 'ymax': lats.max(), 'dx': np.abs(lons[1]-lons[0]), 'dy': np.abs(lats[1]-lats[0]), 'nx': len(lons), 'ny': len(lats)}, adjust='res')
            newgrids[layer]['grid'] = Grid2D(dat, gdict)
        del layergrid, dat
    else:
        lons = np.linspace(xmin, xmax, xsize)
        lats = np.linspace(ymax, ymin, ysize)  # backwards so it plots right side up

    #make meshgrid
    llons1, llats1 = np.meshgrid(lons, lats)

    # See if there is an oceanfile for masking
    bbox = PolygonSH(((cutxmin, cutymin), (cutxmin, cutymax), (cutxmax, cutymax), (cutxmax, cutymin)))
    if oceanfile is not None:
        try:
            f = fiona.open(oceanfile)
            oc = f.next()
            f.close
            shapes = shape(oc['geometry'])
            # make boundaries into a shape
            ocean = shapes.intersection(bbox)
        except:
            print('Not able to read specified ocean file, will use default ocean masking')
            oceanfile = None
    if inventory_shapefile is not None:
        try:
            f = fiona.open(inventory_shapefile)
            invshp = list(f.items(bbox=(bxmin, bymin, bxmax, bymax)))
            f.close()
            inventory = [shape(inv[1]['geometry']) for inv in invshp]
        except:
            print('unable to read inventory shapefile specified, will not plot inventory')
            inventory_shapefile = None

    # # Find cities that will be plotted
    # if mapcities is True and cityfile is not None:
    #     try:
    #         mycity = BasemapCities.loadFromGeoNames(cityfile=cityfile)
    #         bcities = mycity.limitByBounds((bxmin, bxmax, bymin, bymax))
    #         #bcities = bcities.limitByPopulation(40000)
    #         bcities = bcities.limitByGrid(nx=4, ny=4, cities_per_grid=2)
    #     except:
    #         print('Could not read in cityfile, not plotting cities')
    #         mapcities = False
    #         cityfile = None

    # Load in topofile
    if topofile is not None:
        try:
            topomap = GDALGrid.load(topofile, resample=True, method='linear', samplegeodict=gdict)
        except:
            topomap = GMTGrid.load(topofile, resample=True, method='linear', samplegeodict=gdict)
        topodata = topomap.getData().copy()
        # mask oceans if don't have ocean shapefile
        if oceanfile is None:
            topodata = maskoceans(llons1, llats1, topodata, resolution='h', grid=1.25, inlands=True)
    else:
        print('no hillshade is possible\n')
        topomap = None
        topodata = None

    # Load in roads, if needed
    if maproads is True and roadfolder is not None:
        try:
            roadslist = []
            for folder in os.listdir(roadfolder):
                road1 = os.path.join(roadfolder, folder)
                shpfiles = glob.glob(os.path.join(road1, '*.shp'))
                if len(shpfiles):
                    shpfile = shpfiles[0]
                    f = fiona.open(shpfile)
                    shapes = list(f.items(bbox=(bxmin, bymin, bxmax, bymax)))
                    for shapeid, shapedict in shapes:
                        roadslist.append(shapedict)
                    f.close()
        except:
            print('Not able to plot roads')
            roadslist = None

    val = 1
    for k, layer in enumerate(plotorder):
        layergrid = newgrids[layer]['grid']
        if 'label' in grids[layer].keys():
            label1 = grids[layer]['label']
        else:
            label1 = layer
        try:
            sref = grids[layer]['description']['name']
        except:
            sref = None
        ax = fig.add_subplot(rowpan, colpan, val)
        val += 1
        clat = bymin + (bymax-bymin)/2.0
        clon = bxmin + (bxmax-bxmin)/2.0
        # setup of basemap ('lcc' = lambert conformal conic).
        # use major and minor sphere radii from WGS84 ellipsoid.
        m = Basemap(llcrnrlon=bxmin, llcrnrlat=bymin, urcrnrlon=bxmax, urcrnrlat=bymax,
                    rsphere=(6378137.00, 6356752.3142),
                    resolution='i', area_thresh=1000., projection='lcc',
                    lat_1=clat, lon_0=clon, ax=ax)

        x1, y1 = m(llons1, llats1)  # get projection coordinates
        axsize = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        if k == 0:
            wid, ht = axsize.width, axsize.height
        if colormaps is not None and len(colormaps) == len(newgrids) and colormaps[k] is not None:
            palette = eval(colormaps[k])
        else:  # Find preferred default color map for each type of layer
            if 'prob' in layer.lower() or 'pga' in layer.lower() or 'pgv' in layer.lower() or 'cohesion' in layer.lower() or 'friction' in layer.lower() or 'fs' in layer.lower():
                palette = cm.jet
            elif 'slope' in layer.lower():
                palette = cm.gnuplot2
            elif 'precip' in layer.lower():
                palette = cm2.s3pcpn
            else:
                palette = defaultcolormap

        if topodata is not None:
            if k == 0:
                ptopo = m.transform_scalar(np.flipud(topodata), lons+0.5*gdict.dx, lats[::-1]-0.5*gdict.dy, np.round(300.*wid), np.round(300.*ht), returnxy=False, checkbounds=False, order=1, masked=False)
                #use lightsource class to make our shaded topography
                ls = LightSource(azdeg=135, altdeg=45)
                ls1 = LightSource(azdeg=120, altdeg=45)
                ls2 = LightSource(azdeg=225, altdeg=45)
                intensity1 = ls1.hillshade(ptopo, fraction=0.25, vert_exag=1.)
                intensity2 = ls2.hillshade(ptopo, fraction=0.25, vert_exag=1.)
                intensity = intensity1*0.5 + intensity2*0.5
                #hillshm_im = m.transform_scalar(np.flipud(hillshm), lons, lats[::-1], np.round(300.*wid), np.round(300.*ht), returnxy=False, checkbounds=False, order=0, masked=False)
            #m.imshow(hillshm_im, cmap='Greys', vmin=0., vmax=3., zorder=1, interpolation='none')  # vmax = 3 to soften colors to light gray
            #m.pcolormesh(x1, y1, hillshm, cmap='Greys', linewidth=0., rasterized=True, vmin=0., vmax=3., edgecolors='none', zorder=1);
            # plt.draw()

        # Get the data
        dat = layergrid.getData().copy()

        # mask out anything below any specified thresholds
        # Might need to move this up to before downsampling...might give illusion of no hazard in places where there is some that just got averaged out
        if maskthreshes is not None and len(maskthreshes) == len(newgrids):
            if maskthreshes[k] is not None:
                dat[dat <= maskthreshes[k]] = float('NaN')
                dat = np.ma.array(dat, mask=np.isnan(dat))

        if logscale is not False and len(logscale) == len(newgrids):
            if logscale[k] is True:
                dat = np.log10(dat)
                label1 = r'$log_{10}$(' + label1 + ')'

        if scaletype.lower() == 'binned':
            # Find order of range to know how to scale
            order = np.round(np.log(np.nanmax(dat) - np.nanmin(dat)))
            if order < 1.:
                scal = 10**-order
            else:
                scal = 1.
            if lims is None or len(lims) != len(newgrids):
                clev = (np.linspace(np.floor(scal*np.nanmin(dat)), np.ceil(scal*np.nanmax(dat)), 10))/scal
            else:
                if lims[k] is None:
                    clev = (np.linspace(np.floor(scal*np.nanmin(dat)), np.ceil(scal*np.nanmax(dat)), 10))/scal
                else:
                    clev = lims[k]
            # Adjust to colorbar levels
            dat[dat < clev[0]] = clev[0]
            for j, level in enumerate(clev[:-1]):
                dat[(dat >= clev[j]) & (dat < clev[j+1])] = clev[j]
            # So colorbar saturates at top
            dat[dat > clev[-1]] = clev[-1]
            #panelhandle = m.contourf(x1, y1, datm, clev, cmap=palette, linewidth=0., alpha=ALPHA, rasterized=True)
            vmin = clev[0]
            vmax = clev[-1]
        else:
            if lims is not None and len(lims) == len(newgrids):
                if lims[k] is None:
                    vmin = None
                    vmax = None
                else:
                    vmin = lims[k][0]
                    vmax = lims[k][-1]
            else:
                vmin = np.nanmin(dat)
                vmax = np.nanmax(dat)

        # Mask out cells overlying oceans or block with a shapefile if available
        if oceanfile is None:
            dat = maskoceans(llons1, llats1, dat, resolution='h', grid=1.25, inlands=True)
        else:
            #patches = []
            if type(ocean) is PolygonSH:
                ocean = [ocean]
            for oc in ocean:
                patch = getProjectedPatch(oc, m, edgecolor="#006280", facecolor=watercolor, lw=0.5, zorder=4.)
                #x, y = m(oc.exterior.xy[0], oc.exterior.xy[1])
                #xy = zip(x, y)
                #patch = Polygon(xy, facecolor=watercolor, edgecolor="#006280", lw=0.5, zorder=4.)
                ##patches.append(Polygon(xy, facecolor=watercolor, edgecolor=watercolor, zorder=500.))
                ax.add_patch(patch)
            ##ax.add_collection(PatchCollection(patches))

        if inventory_shapefile is not None:
            for in1 in inventory:
                x, y = m(in1.exterior.xy[0], in1.exterior.xy[1])
                xy = zip(x, y)
                patch = Polygon(xy, facecolor='none', edgecolor='k', lw=0.5, zorder=10.)
                #patches.append(Polygon(xy, facecolor=watercolor, edgecolor=watercolor, zorder=500.))
                ax.add_patch(patch)
        palette.set_bad(clear_color, alpha=0.0)
        # Plot it up
        dat_im = m.transform_scalar(np.flipud(dat), lons+0.5*gdict.dx, lats[::-1]-0.5*gdict.dy, np.round(300.*wid), np.round(300.*ht), returnxy=False, checkbounds=False, order=0, masked=True)
        if topodata is not None:  # Drape over hillshade
            #turn data into an RGBA image
            cmap = palette
            #adjust data so scaled between vmin and vmax and between 0 and 1
            dat1 = dat_im.copy()
            dat1[dat1 < vmin] = vmin
            dat1[dat1 > vmax] = vmax
            dat1 = (dat1 - vmin)/(vmax-vmin)
            rgba_img = cmap(dat1)
            maskvals = np.dstack((dat1.mask, dat1.mask, dat1.mask))
            rgb = np.squeeze(rgba_img[:, :, 0:3])
            rgb[maskvals] = 1.
            draped_hsv = ls.blend_hsv(rgb, np.expand_dims(intensity, 2))
            m.imshow(draped_hsv, zorder=3., interpolation='none')
            panelhandle = m.imshow(dat_im, cmap=palette, zorder=0., vmin=vmin, vmax=vmax)  # This is just a dummy layer that will be deleted to make the colorbar look right
        else:
            panelhandle = m.imshow(dat_im, cmap=palette, zorder=3., vmin=vmin, vmax=vmax, interpolation='none')
        #panelhandle = m.pcolormesh(x1, y1, dat, linewidth=0., cmap=palette, vmin=vmin, vmax=vmax, alpha=ALPHA, rasterized=True, zorder=2.);
        #panelhandle.set_edgecolors('face')
        # add colorbar
        cbfmt = '%1.1f'
        if vmax is not None and vmin is not None:
            if (vmax - vmin) < 1.:
                cbfmt = '%1.2f'
            elif vmax > 5.:  # (vmax - vmin) > len(clev):
                cbfmt = '%1.0f'

        #norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        if scaletype.lower() == 'binned':
            cbar = fig.colorbar(panelhandle, spacing='proportional', ticks=clev, boundaries=clev, fraction=0.036, pad=0.04, format=cbfmt, extend='both')
            #cbar1 = ColorbarBase(cbar.ax, cmap=palette, norm=norm, spacing='proportional', ticks=clev, boundaries=clev, fraction=0.036, pad=0.04, format=cbfmt, extend='both', extendfrac='auto')

        else:
            cbar = fig.colorbar(panelhandle, fraction=0.036, pad=0.04, extend='both', format=cbfmt)
            #cbar1 = ColorbarBase(cbar.ax, cmap=palette, norm=norm, fraction=0.036, pad=0.04, extend='both', extendfrac='auto', format=cbfmt)

        if topodata is not None:
            panelhandle.remove()

        cbar.set_label(label1, fontsize=10)
        cbar.ax.tick_params(labelsize=8)

        parallels = m.drawparallels(getMapLines(bymin, bymax, 3), labels=[1, 0, 0, 0], linewidth=0.5, labelstyle='+/-', fontsize=9, xoffset=-0.8, color='gray', zorder=100.)
        m.drawmeridians(getMapLines(bxmin, bxmax, 3), labels=[0, 0, 0, 1], linewidth=0.5, labelstyle='+/-', fontsize=9, color='gray', zorder=100.)
        for par in parallels:
            try:
                parallels[par][1][0].set_rotation(90)
            except:
                pass

        #draw roads on the map, if they were provided to us
        if maproads is True and roadslist is not None:
            try:
                for road in roadslist:
                    xy = list(road['geometry']['coordinates'])
                    roadx, roady = zip(*xy)
                    mapx, mapy = m(roadx, roady)
                    m.plot(mapx, mapy, roadcolor, lw=0.5, zorder=9)
            except Exception as e:
                print('Failed to plot roads, %s' % e)

        #add city names to map
        # if mapcities is True and cityfile is not None:
        #     try:
        #         fontname = 'Arial'
        #         fontsize = 8
        #         if k == 0:  # Only need to choose cities first time and then apply to rest
        #             fcities = bcities.limitByMapCollision(m, fontname=fontname, fontsize=fontsize)
        #             ctlats, ctlons, names = fcities.getCities()
        #             cxis, cyis = m(ctlons, ctlats)
        #         for ctlat, ctlon, cxi, cyi, name in zip(ctlats, ctlons, cxis, cyis, names):
        #             m.scatter(ctlon, ctlat, c='k', latlon=True, marker='.', zorder=100000)
        #             ax.text(cxi, cyi, name, fontname=fontname, fontsize=fontsize, zorder=100000)
        #     except Exception as e:
        #         print('Failed to plot cities, %s' % e)

        if mapcities is True and cityfile is not None:
            try:
                dmin = 0.1*(m.ymax-m.ymin)
                xyplotted = []
                cities = PagerCity(cityfile)
                #Find cities within bounding box
                boundcity = cities.findCitiesByRectangle(bounds=(boundaries.xmin, boundaries.xmax, boundaries.ymin, boundaries.ymax))
                #Just keep 5 biggest cities
                if len(boundcity) < 5:
                    value = len(boundcity)
                else:
                    value = 5
                thresh = sorted([cit['pop'] for cit in boundcity])[-value]
                plotcity = [cit for cit in boundcity if cit['pop'] >= thresh]
                #For cities that are more than one xth of the xwidth apart, keep only the larger one
                pass  # do later
                #Plot cities
                for cit in plotcity:  # should sort so it plots them in order of population so larger cities are preferentially plotted - do later
                    xi, yi = m(cit['lon'], cit['lat'])
                    dist = [np.sqrt((xi-x0)**2+(yi-y0)**2) for x0, y0 in xyplotted]
                    xdist = [np.abs(xi-x0) for x0, y0 in xyplotted]
                    ydist = [np.abs(yi-y0) for x0, y0 in xyplotted]
                    if not dist or np.min(dist) > dmin:
                        if len(dist) > 0:
                            if np.min(xdist) < 0.2*(m.xmax-m.xmin) and np.min(ydist) < 0.1*(m.ymax-m.ymin):
                                pass
                            else:
                                m.scatter(cit['lon'], cit['lat'], c='k', latlon=True, marker='.', zorder=100000)
                                ax.text(xi, yi, cit['name'], ha='right', va='top', fontsize=10, zorder=100000)
                                xyplotted.append((xi, yi))
                        elif len(dist) == 0:
                            m.scatter(cit['lon'], cit['lat'], c='k', latlon=True, marker='.', zorder=100000)
                            ax.text(xi, yi, cit['name'], ha='right', va='top', fontsize=10, zorder=100000)
                            xyplotted.append((xi, yi))
            except Exception as e:
                print('Failed to plot cities, %s' % e)

        #draw star at epicenter
        plt.sca(ax)
        if edict is not None:
            elat, elon = edict['lat'], edict['lon']
            ex, ey = m(elon, elat)
            plt.plot(ex, ey, '*', markeredgecolor='k', mfc='None', mew=1.0, ms=15, zorder=10000.)

        m.drawmapboundary(fill_color=watercolor)

        m.fillcontinents(color=clear_color, lake_color=watercolor)
        m.drawrivers(color=watercolor)
        #m.drawcoastlines()

        #draw country boundaries
        m.drawcountries(color=countrycolor, linewidth=1.0)

        #add map scale
        m.drawmapscale((bxmax+bxmin)/2., (bymin+(bymax-bymin)/9.), clon, clat, np.round((((bxmax-bxmin)*111)/5)/10.)*10, barstyle='fancy', zorder=10)

        # Add border
        autoAxis = ax.axis()
        rec = Rectangle((autoAxis[0]-0.7, autoAxis[2]-0.2), (autoAxis[1]-autoAxis[0])+1, (autoAxis[3]-autoAxis[2])+0.4, fill=False, lw=1, zorder=1e8)
        rec = ax.add_patch(rec)
        rec.set_clip_on(False)

        plt.draw()

        if sref is not None:
            label2 = '%s\nsource: %s' % (label1, sref)  # '%s\n' % label1 + r'{\fontsize{10pt}{3em}\selectfont{}%s}' % sref  #
        else:
            label2 = label1
        plt.title(label2, axes=ax, fontsize=fontsizesub)

        #draw scenario watermark, if scenario
        if isScenario:
            plt.sca(ax)
            cx, cy = m(clon, clat)
            plt.text(cx, cy, 'SCENARIO', rotation=45, alpha=0.10, size=72, ha='center', va='center', color='red')

        #if ds: # Could add this to print "downsampled" on map
        #    plt.text()

        if k == 1 and rowpan == 1:
            # adjust single level plot
            axsize = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            ht2 = axsize.height
            fig.set_figheight(ht2*1.6)
        else:
            plt.tight_layout()

        # Make room for suptitle - tight layout doesn't account for it
        plt.subplots_adjust(top=0.92)

    if printparam is True:
        try:
            fig = plt.gcf()
            dictionary = grids['model']['description']['parameters']
            paramstring = 'Model parameters: '
            halfway = np.ceil(len(dictionary)/2.)
            for i, key in enumerate(dictionary):
                if i == halfway and colpan == 1:
                    paramstring += '\n'
                paramstring += ('%s = %s; ' % (key, dictionary[key]))
            print paramstring
            fig.text(0.01, 0.015, paramstring, fontsize=fontsizesmallest)
            plt.draw()
        except:
            print('Could not display model parameters')

    if edict is not None:
        eventid = edict['eventid']
    else:
        eventid = ''

    time1 = datetime.datetime.utcnow().strftime('%d%b%Y_%H%M')
    outfile = os.path.join(outfolder, '%s_%s_%s.pdf' % (eventid, suptitle, time1))
    pngfile = os.path.join(outfolder, '%s_%s_%s.png' % (eventid, suptitle, time1))

    if savepdf is True:
        print 'Saving map output to %s' % outfile
        plt.savefig(outfile, dpi=300)
    if savepng is True:
        print 'Saving map output to %s' % pngfile
        plt.savefig(pngfile)
    if showplots is True:
        plt.show()
    else:
        plt.close(fig)

    return newgrids


def make_hillshade(topogrid, azimuth=315., angle_altitude=50.):
    """
    Computes a hillshade from a digital elevation model. Most of this script borrowed from http://geoexamples.blogspot.com/2014/03/shaded-relief-images-using-gdal-python.html last accessed 9/2/2015

    :param topogrid: Digital elevation model
    :type topogrid: Grid2D object
    :param azimuth: azimuth of illumination in degrees
    :type azimuth: float
    :param angle_altitude: altitude angle of illumination
    :type angle_altitude: float

    :returns hillshade: Hillshade map layer
    :rtype hillshade: Grid2D object

    """
    topotmp = topogrid.getData().copy()
    #make a masked array
    topotmp = np.ma.array(topotmp)
    topodat = np.ma.masked_where(np.isnan(topotmp), topotmp)
    #topodat = np.ma.masked_where(topodat == SEA_LEVEL, topodat)
    x, y = np.gradient(topodat)
    slope = np.pi/2. - np.arctan(np.sqrt(x*x + y*y))
    aspect = np.arctan2(-x, y)
    azimuthrad = azimuth*np.pi / 180.
    altituderad = angle_altitude*np.pi / 180.
    shaded = np.sin(altituderad) * np.sin(slope) + np.cos(altituderad) * np.cos(slope) * np.cos(azimuthrad - aspect)
    hillshade = copy.deepcopy(topogrid)
    hillshade.setData(255*(shaded + 1)/2)

    return hillshade


def getMapLines(dmin, dmax, nlines):
    drange = dmax-dmin
    if drange > 4:
        near = 1
    else:
        if drange >= 0.5:
            near = 0.25
        else:
            near = 0.125
    inc = roundToNearest(drange/nlines, near)
    if inc == 0:
        near = np.power(10, round(math.log10(drange)))  # make the increment the closest power of 10
        inc = ceilToNearest(drange/nlines, near)
        newdmin = floorToNearest(dmin, near)
        newdmax = ceilToNearest(dmax, near)
    else:
        newdmin = ceilToNearest(dmin, near)
        newdmax = floorToNearest(dmax, near)
    darray = np.arange(newdmin, newdmax+inc, inc)
    if darray[-1] > dmax:
        darray = darray[0:-1]
    return darray


def getProjectedPatch(polygon, m, edgecolor, facecolor, lw=1., zorder=10):
    polyjson = mapping(polygon)
    tlist = []
    for sequence in polyjson['coordinates']:
        lon, lat = zip(*sequence)
        x, y = m(lon, lat)
        tlist.append(tuple(zip(x, y)))
    polyjson['coordinates'] = tuple(tlist)
    ppolygon = shape(polyjson)
    patch = PolygonPatch(ppolygon, facecolor=facecolor, edgecolor=edgecolor,
                         zorder=zorder, linewidth=lw, fill=True, visible=True)
    return patch


if __name__ == '__main__':
    pass
