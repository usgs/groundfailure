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
mpl.use('Agg')  # so figures will still be created even without display
from matplotlib.colors import LightSource, LogNorm
import re
#from matplotlib.colorbar import ColorbarBase

#third party imports
import matplotlib.cm as cm
import branca.colormap as cmb
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
import shapefile
#from json import dumps
import folium
from folium import plugins
from folium.features import GeoJson, RectangleMarker


#local imports
from mapio.gmt import GMTGrid
from mapio.gdal import GDALGrid
from mapio.geodict import GeoDict
from mapio.grid2d import Grid2D
#from neicutil.text import ceilToNearest, floorToNearest, roundToNearest
from mapio.basemapcity import BasemapCities
from mapio.shake import ShakeGrid

# Make fonts readable and recognizable by illustrator
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'Bitstream Vera Serif', 'sans-serif']


def parseMapConfig(config, fileext=None):
    # Parse config object
    # fileext will be prepended to any file paths in config
    # ADD PLOTORDER TO CONFIG? OTHER THINGS LIKE COLORMAPS?
    topofile = None
    roadfolder = None
    cityfile = None
    roadcolor = '6E6E6E'
    countrycolor = '177F10'
    watercolor = 'B8EEFF'
    ALPHA = 0.7
    oceanfile = None
    oceanref = None
    roadref = None
    cityref = None

    #try:
    if fileext is None:
        fileext = '.'
    if 'dem' in config:
        topofile = os.path.join(fileext, config['dem']['file'])
        if os.path.exists(topofile) is False:
            print('DEM not valid - hillshade will not be possible\n')
    if 'ocean' in config:
        oceanfile = os.path.join(fileext, config['ocean']['file'])
        try:
            oceanref = config['ocean']['shortref']
        except:
            oceanref = 'unknown'
    if 'roads' in config:
        roadfolder = os.path.join(fileext, config['roads']['file'])
        if os.path.exists(roadfolder) is False:
            print('roadfolder not valid - roads will not be displayed\n')
            roadfolder = None
        try:
            roadref = config['roads']['shortref']
        except:
            roadref = 'unknown'
    if 'cities' in config:
        cityfile = os.path.join(fileext, config['cities']['file'])
        try:
            cityref = config['cities']['shortref']
        except:
            cityref = 'unknown'
        if os.path.exists(cityfile):
            try:
                #PagerCity(cityfile)
                BasemapCities.loadFromGeoNames(cityfile=cityfile)
            except Exception as e:
                print(e)
                print('cities file not valid - cities will not be displayed\n')
                cityfile = None
        else:
            print('cities file not valid - cities will not be displayed\n')
            cityfile = False
    if 'roadcolor' in config['colors']:
        roadcolor = config['colors']['roadcolor']
    if 'countrycolor' in config['colors']:
        countrycolor = config['colors']['countrycolor']
    if 'watercolor' in config['colors']:
        watercolor = config['colors']['watercolor']
    if 'alpha' in config['colors']:
        ALPHA = float(config['colors']['alpha'])
    #except Exception as e:
    #    print(('%s - mapping options missing from or misformatted in config' % e))

    countrycolor = '#'+countrycolor
    watercolor = '#'+watercolor
    roadcolor = '#'+roadcolor

    mapin = {'topofile': topofile, 'roadfolder': roadfolder,
             'cityfile': cityfile, 'roadcolor': roadcolor,
             'countrycolor': countrycolor, 'watercolor': watercolor,
             'ALPHA': ALPHA, 'roadref': roadref,
             'cityref': cityref, 'oceanfile': oceanfile, 'oceanref': oceanref}

    return mapin


def parseConfigLayers(maplayers, config, keys=None):
    """
    Parse things that need to coodinate with each layer (like lims, logscale,
    colormaps etc.) from config file, in right order - takes orders from
    maplayers
    # TO DO, ADD ABILITY TO INTERPRET CUSTOM COLOR MAPS
    """
    # get all key names, create a plotorder list in case maplayers is not an
    # ordered dict, making sure that anything called 'model' is first
    if keys is None:
        keys = list(maplayers.keys())
    plotorder = []

    try:
        limits = config[config.keys()[0]]['display_options']['lims']
        lims = []
    except:
        lims = None
        limits = None

    try:
        colors = config[config.keys()[0]]['display_options']['colors']
        colormaps = []
    except:
        colormaps = None
        colors = None

    try:
        logs = config[config.keys()[0]]['display_options']['logscale']
        logscale = []
    except:
        logscale = False
        logs = None

    try:
        masks = config[config.keys()[0]]['display_options']['maskthresholds']
        maskthreshes = []
    except:
        maskthreshes = None
        masks = None

    try:
        default = config[config.keys()[0]]['display_options']['colors']['default']
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


def modelMap(grids, shakefile=None, suptitle=None, inventory_shapefile=None,
             plotorder=None, maskthreshes=None, colormaps=None, boundaries=None,
             zthresh=0, scaletype='continuous', lims=None, logscale=False,
             ALPHA=0.7, maproads=True, mapcities=True, isScenario=False,
             roadfolder=None, topofile=None, cityfile=None, oceanfile=None,
             roadcolor='#6E6E6E', watercolor='#B8EEFF', countrycolor='#177F10',
             outputdir=None, outfilename=None, savepdf=True, savepng=True, showplots=False,
             roadref='unknown', cityref='unknown', oceanref='unknown',
             printparam=False, ds=True, dstype='mean', upsample=False):
    """
    This function creates maps of mapio grid layers (e.g. liquefaction or
    landslide models with their input layers)
    All grids must use the same bounds
    TO DO change so that all input layers do not have to have the same bounds,
    test plotting multiple probability layers, and add option so that if PDF and
    PNG aren't output, opens plot on screen using plt.show()

    :param grids: Dictionary of N layers and metadata formatted like:
        maplayers['layer name']={
        'grid': mapio grid2D object,
        'label': 'label for colorbar and top line of subtitle',
        'type': 'output or input to model',
        'description': 'detailed description of layer for subtitle'}.
      Layer names must be unique.
    :type name: Dictionary or Ordered dictionary - import collections;
      grids = collections.OrderedDict()
    :param shakefile: optional ShakeMap file (url or full file path) to extract information for labels and folder names
    :type shakefile: Shakemap Event Dictionary
    :param suptitle: This will be displayed at the top of the plots and in the
      figure names
    :type suptitle: string
    :param plotorder: List of keys describing the order to plot the grids, if
      None and grids is an ordered dictionary, it will use the order of the
      dictionary, otherwise it will choose order which may be somewhat random
      but it will always put a probability grid first
    :type plotorder: list
    :param maskthreshes: N x 1 array or list of lower thresholds for masking
      corresponding to order in plotorder or order of OrderedDict if plotorder
      is None. If grids is not an ordered dict and plotorder is not specified,
      this will not work right. If None (default), nothing will be masked
    :param colormaps: List of strings of matplotlib colormaps (e.g. cm.autumn_r)
      corresponding to plotorder or order of dictionary if plotorder is None.
      The list can contain both strings and None e.g. colormaps = ['cm.autumn',
      None, None, 'cm.jet'] and None's will default to default colormap
    :param boundaries: None to show entire study area, 'zoom' to zoom in on the
      area of action (only works if there is a probability layer) using zthresh
      as a threshold, or a dictionary defining lats and lons in the form of
      boundaries.xmin = minlon, boundaries.xmax = maxlon, boundaries.ymin =
      min lat, boundaries.ymax = max lat
    :param zthresh: threshold for computing zooming bounds, only used if
      boundaries = 'zoom'
    :type zthresh: float
    :param scaletype: Type of scale for plotting, 'continuous' or 'binned' -
      will be reflected in colorbar
    :type scaletype: string
    :param lims: None or Nx1 list of tuples or numpy arrays corresponding to
      plotorder defining the limits for saturating the colorbar (vmin, vmax) if
      scaletype is continuous or the bins to use (clev) if scaletype if binned.
      The list can contain tuples, arrays, and Nones, e.g. lims = [(0., 10.),
      None, (0.1, 1.5), np.linspace(0., 1.5, 15)]. When None is specified, the
      program will estimate the limits, when an array is specified but the scale
      type is continuous, vmin will be set to min(array) and vmax will be set
      to max(array)
    :param lims: None or Nx1 list of Trues and Falses corresponding to
      plotorder defining whether to use a linear or log scale (log10) for
      plotting the layer. This will be reflected in the labels
    :param ALPHA: Transparency for mapping, if there is a hillshade that will
      plot below each layer, it is recommended to set this to at least 0.7
    :type ALPHA: float
    :param maproads: Whether to show roads or not, default True, but requires
      that roadfile is specified and valid to work
    :type maproads: boolean
    :param mapcities: Whether to show cities or not, default True, but requires
      that cityfile is specified and valid to work
    :type mapcities: boolean
    :param isScenario: Whether this is a scenario (True) or a real event (False)
      (default False)
    :type isScenario: boolean
    :param roadfolder: Full file path to folder containing road shapefiles
    :type roadfolder: string
    :param topofile: Full file path to topography grid (GDAL compatible) - this
      is only needed to make a hillshade if a premade hillshade is not specified
    :type topofile: string
    :param cityfile: Full file path to Pager file containing city & population
      information
    :type cityfile: string
    :param roadcolor: Color to use for roads, if plotted, default #6E6E6E
    :type roadcolor: Hex color or other matplotlib compatible way of defining
      color
    :param watercolor: Color to use for oceans, lakes, and rivers, default
      #B8EEFF
    :type watercolor: Hex color or other matplotlib compatible way of defining
      color
    :param countrycolor: Color for country borders, default #177F10
    :type countrycolor: Hex color or other matplotlib compatible way of defining
      color
    :param outputdir: File path for outputting figures, if edict is defined, a
      subfolder based on the event id will be created in this folder. If None,
      will use current directory
    :param savepdf: True to save pdf figure, False to not
    :param savepng: True to save png figure, False to not
    :param ds: True to allow downsampling for display (necessary when arrays
      are quite large, False to not allow)
    :param dstype: What function to use in downsampling, options are 'min',
      'max', 'median', or 'mean'
    :param upsample: True to upsample the layer to the DEM resolution for better
      looking hillshades

    :returns:
        * PDF and/or PNG of map
        * Downsampled and trimmed version of input grids. If no
        modification was needed for plotting, this will be identical to grids but
        without the metadata
        * list of output filenames

    """

    if suptitle is None:
        suptitle = ' '

    if not showplots or (not savepdf and not savepng):  # display fig only if static fig not output
        plt.ioff()

    defaultcolormap = cm.CMRmap_r

    if shakefile is not None:
        edict = ShakeGrid.load(shakefile, adjust='res').getEventDict()
        temp = ShakeGrid.load(shakefile, adjust='res').getShakeDict()
        edict['eventid'] = temp['shakemap_id']
        edict['version'] = temp['shakemap_version']
    else:
        edict = None

    # Get output file location
    if outputdir is None:
        print('No output location given, using current directory for outputs\n')
        outputdir = os.getcwd()
    if edict is not None:
        outfolder = os.path.join(outputdir, edict['event_id'])
    else:
        outfolder = outputdir
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)

    # Get plotting order, if not specified
    if plotorder is None:
        plotorder = list(grids.keys())

    # Get boundaries to use for all plots
    cut = True
    if boundaries is None:
        cut = False
        keytemp = list(plotorder)
        boundaries = grids[keytemp[0]]['grid'].getGeoDict()
    elif boundaries == 'zoom':
        # Find probability layer (will just take the maximum bounds if there is
        # more than one)
        keytemp = list(plotorder)
        key1 = [key for key in keytemp if 'model' in key.lower()]
        if len(key1) == 0:
            print('Could not find model layer to use for zoom, using default boundaries')
            keytemp = list(plotorder)
            boundaries = grids[keytemp[0]]['grid'].getGeoDict()
        else:
            lonmax = -1.e10
            lonmin = 1.e10
            latmax = -1.e10
            latmin = 1.e10
            for key in key1:
                # get lat lons of areas affected and add, if no areas affected,
                # switch to shakemap boundaries
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
        keytemp = list(grids.keys())
        tempgdict = grids[keytemp[0]]['grid'].getGeoDict()
        if np.abs(tempgdict.xmin-boundaries['xmin']) < 0.05 and \
           np.abs(tempgdict.ymin-boundaries['ymin']) < 0.05 and \
           np.abs(tempgdict.xmax-boundaries['xmax']) < 0.05 and \
           np.abs(tempgdict.ymax - boundaries['ymax']) < 0.05:
            print('Input boundaries are almost the same as specified boundaries, no cutting needed')
            boundaries = tempgdict
            cut = False
        else:
            try:
                if boundaries['xmin'] > boundaries['xmax'] or \
                   boundaries['ymin'] > boundaries['ymax']:
                    print('Input boundaries are not usable, using default boundaries')
                    keytemp = list(grids.keys())
                    boundaries = grids[keytemp[0]]['grid'].getGeoDict()
                    cut = False
                else:
                    # Build dummy GeoDict
                    boundaries = GeoDict({'xmin': boundaries['xmin'],
                                          'xmax': boundaries['xmax'],
                                          'ymin': boundaries['ymin'],
                                          'ymax': boundaries['ymax'],
                                          'dx': 100.,
                                          'dy': 100.,
                                          'ny': 100.,
                                          'nx': 100.},
                                         adjust='res')
            except:
                print('Input boundaries are not usable, using default boundaries')
                keytemp = list(grids.keys())
                boundaries = grids[keytemp[0]]['grid'].getGeoDict()
                cut = False

    # Pull out bounds for various uses
    bxmin, bxmax, bymin, bymax = boundaries.xmin, boundaries.xmax, boundaries.ymin, boundaries.ymax

    # Determine if need a single panel or multi-panel plot and if multi-panel,
    # how many and how it will be arranged
    fig = plt.figure()
    numpanels = len(plotorder)
    if numpanels == 1:
        rowpan = 1
        colpan = 1
        # create the figure and axes instances.
        fig.set_figwidth(5)
    elif numpanels == 2 or numpanels == 4:
        rowpan = np.ceil(numpanels/2.)
        colpan = 2
        fig.set_figwidth(13)
    else:
        rowpan = np.ceil(numpanels/3.)
        colpan = 3
        fig.set_figwidth(15)
    if rowpan == 1:
        fig.set_figheight(rowpan*6.0)
    else:
        fig.set_figheight(rowpan*5.3)

    # Need to update naming to reflect the shakemap version once can get
    # getHeaderData to work, add edict['version'] back into title, maybe
    # shakemap id also?
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
                print(('Cutting failed, %s, continuing with full layers' % e))
                newgrids = grids
                continue
        del templayer
        gc.collect()
    else:
        newgrids = grids
    tempgdict = newgrids[list(grids.keys())[0]]['grid'].getGeoDict()

    # Upsample layers to same as topofile if desired for better looking hillshades
    if upsample is True and topofile is not None:
        try:
            topodict = GDALGrid.getFileGeoDict(topofile)
            if topodict.dx >= tempgdict.dx or topodict.dy >= tempgdict.dy:
                print('Upsampling not possible, resolution of results already smaller than DEM')
                pass
            else:
                tempgdict1 = GeoDict({'xmin': tempgdict.xmin-xbuff,
                                      'ymin': tempgdict.ymin-ybuff,
                                      'xmax': tempgdict.xmax+xbuff,
                                      'ymax': tempgdict.ymax+ybuff,
                                      'dx': topodict.dx,
                                      'dy': topodict.dy,
                                      'nx': topodict.nx,
                                      'ny': topodict.ny},
                                     adjust='res')
                tempgdict2 = tempgdict1.getBoundsWithin(tempgdict)
                for k, layer in enumerate(plotorder):
                    newgrids[layer]['grid'] = newgrids[layer]['grid'].subdivide(tempgdict2)
        except:
            print('Upsampling failed, continuing')

    # Downsample all of them for plotting, if needed, and replace them in
    # grids (to save memory)
    tempgrid = newgrids[list(grids.keys())[0]]['grid']
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
            dat = block_reduce(layergrid.getData().copy(),
                               block_size=(divy, divx),
                               cval=float('nan'),
                               func=func)
            if k == 0:
                lons = block_reduce(np.linspace(xmin, xmax, layergrid.getGeoDict().nx),
                                    block_size=(divx,),
                                    func=np.mean,
                                    cval=float('nan'))
                if math.isnan(lons[-1]):
                    lons[-1] = lons[-2] + (lons[1]-lons[0])
                lats = block_reduce(np.linspace(ymax, ymin, layergrid.getGeoDict().ny),
                                    block_size=(divy,),
                                    func=np.mean,
                                    cval=float('nan'))
                if math.isnan(lats[-1]):
                    lats[-1] = lats[-2] + (lats[1]-lats[0])
                gdict = GeoDict({'xmin': lons.min(),
                                 'xmax': lons.max(),
                                 'ymin': lats.min(),
                                 'ymax': lats.max(),
                                 'dx': np.abs(lons[1]-lons[0]),
                                 'dy': np.abs(lats[1]-lats[0]),
                                 'nx': len(lons),
                                 'ny': len(lats)},
                                adjust='res')
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
            oc = next(f)
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
    if mapcities is True and cityfile is not None:
        try:
            mycity = BasemapCities.loadFromGeoNames(cityfile=cityfile)
            bcities = mycity.limitByBounds((bxmin, bxmax, bymin, bymax))
            #bcities = bcities.limitByPopulation(40000)
            bcities = bcities.limitByGrid(nx=4, ny=4, cities_per_grid=2)
        except:
            print('Could not read in cityfile, not plotting cities')
            mapcities = False
            cityfile = None

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
    roadslist = None
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

    val = 1
    for k, layer in enumerate(plotorder):
        layergrid = newgrids[layer]['grid']
        if 'label' in list(grids[layer].keys()):
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
        # setup of basemap ('lcc' = lambert conformal conic, or tmerc is transverse mercator).
        # use major and minor sphere radii from WGS84 ellipsoid.
        m = Basemap(llcrnrlon=bxmin, llcrnrlat=bymin, urcrnrlon=bxmax, urcrnrlat=bymax,
                    rsphere=(6378137.00, 6356752.3142),
                    resolution='l', area_thresh=1000., projection='tmerc',
                    lat_0=clat, lon_0=clon, ax=ax)

        x1, y1 = m(llons1, llats1)  # get projection coordinates
        axsize = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        if k == 0:
            wid, ht = axsize.width, axsize.height
        default1 = True
        if len(colormaps) == 1 and len(newgrids) == 1:
            palette = colormaps
            default1 = False
        if colormaps is not None and len(colormaps) == len(newgrids) and colormaps[k] is not None:
            palette = colormaps[k]
            default1 = False
        if default1:  # Find preferred default color map for each type of layer if no colormaps found
            if 'prob' in layer.lower() or 'pga' in layer.lower() or \
               'pgv' in layer.lower() or 'cohesion' in layer.lower() or \
               'friction' in layer.lower() or 'fs' in layer.lower():
                palette = cm.CMRmap_r
            elif 'slope' in layer.lower():
                palette = cm.gnuplot2
            elif 'precip' in layer.lower():
                palette = cm2.s3pcpn
            else:
                palette = defaultcolormap

        if topodata is not None:
            if k == 0:
                #lons[lons > 0] = lons[lons > 0] - 360
                ptopo = m.transform_scalar(
                    np.flipud(topodata), lons+0.5*gdict.dx,
                    lats[::-1]-0.5*gdict.dy, int(np.round(300.*wid)),
                    int(np.round(300.*ht)), returnxy=False, checkbounds=False,
                    order=1, masked=False)
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

        # if logscale is not False and len(logscale) == len(newgrids):
        #     if logscale[k] is True:
        #         dat = np.log10(dat)
        #         label1 = r'$log_{10}$(' + label1 + ')'

        if scaletype.lower() == 'binned':
            if logscale is not False and len(logscale) == len(newgrids):
                if logscale[k] is True:
                    clev = 10.**(np.arange(np.floor(np.log10(np.nanmin(dat))), np.ceil(np.log10(np.nanmax(dat))), 0.25))
                else:
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
            else:
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
                dat[(dat >= clev[j]) & (dat < clev[j+1])] = (clev[j] + clev[j+1])/2.
            # So colorbar saturates at top
            dat[dat > clev[-1]] = clev[-1]
            #panelhandle = m.contourf(x1, y1, datm, clev, cmap=palette, linewidth=0., alpha=ALPHA, rasterized=True)
            vmin = clev[0]
            vmax = clev[-1]

        else:
            if isinstance(logscale, (bool)):
                # Put it in a list so it won't break things later
                logscale = [logscale]
            if lims is not None and len(lims) == len(newgrids):
                if lims[k] is None:
                    vmin = np.nanmin(dat)
                    vmax = np.nanmax(dat)
                else:
                    vmin = lims[k][0]
                    vmax = lims[k][-1]
            else:
                vmin = np.nanmin(dat)
                vmax = np.nanmax(dat)

        # Might need to move this up to before downsampling...might give illusion of no hazard in places where there is some that just got averaged out
        if maskthreshes is not None and len(maskthreshes) == len(newgrids):
            if maskthreshes[k] is not None:
                dat[dat <= maskthreshes[k]] = float('NaN')
                dat = np.ma.array(dat, mask=np.isnan(dat))

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
                if 'point' in str(type(in1)):
                    x, y = in1.xy
                    x = x[0]
                    y = y[0]
                    m.scatter(x, y, c='m', s=50, latlon=True, marker='^',
                              zorder=100001)
                else:
                    x, y = m(in1.exterior.xy[0], in1.exterior.xy[1])
                    xy = list(zip(x, y))
                    patch = Polygon(xy, facecolor='none', edgecolor='k', lw=0.5, zorder=10.)
                    #patches.append(Polygon(xy, facecolor=watercolor, edgecolor=watercolor, zorder=500.))
                    ax.add_patch(patch)
        palette.set_bad(clear_color, alpha=0.0)
        # Plot it up
        dat_im = m.transform_scalar(np.flipud(dat), lons+0.5*gdict.dx, lats[::-1]-0.5*gdict.dy,
                                    int(np.round(300.*wid)), int(np.round(300.*ht)), returnxy=False,
                                    checkbounds=False, order=0, masked=True)
        if logscale[k] is True:
            logsc = LogNorm(vmin=vmin, vmax=vmax)
        else:
            logsc = None
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
            m.imshow(draped_hsv, zorder=3., interpolation='none', norm=logsc)
            # This is just a dummy layer that will be deleted to make the
            # colorbar look right
            panelhandle = m.imshow(dat_im, cmap=palette, zorder=0.,
                                   vmin=vmin, vmax=vmax, norm=logsc, interpolation='none')
        else:
            panelhandle = m.imshow(dat_im, cmap=palette, zorder=3., norm=logsc,
                                   vmin=vmin, vmax=vmax, interpolation='none')
        #panelhandle = m.pcolormesh(x1, y1, dat, linewidth=0., cmap=palette, vmin=vmin, vmax=vmax, alpha=ALPHA, rasterized=True, zorder=2.);
        #panelhandle.set_edgecolors('face')
        # add colorbar
        cbfmt = '%1.2f'
        if vmax is not None and vmin is not None:
            if logscale is not False and len(logscale) == len(newgrids):
                if logscale[k] is True:
                    cbfmt = None
            elif (vmax - vmin) < 1.:
                cbfmt = '%1.2f'
            elif vmax > 5.:  # (vmax - vmin) > len(clev):
                cbfmt = '%1.0f'

        #norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        if scaletype.lower() == 'binned':
            cbar = fig.colorbar(panelhandle, spacing='proportional',
                                ticks=clev, boundaries=clev, fraction=0.036,
                                pad=0.04, format=cbfmt, extend='both', norm=logsc)  # extend='both'
            #cbar1 = ColorbarBase(cbar.ax, cmap=palette, norm=norm, spacing='proportional', ticks=clev, boundaries=clev, fraction=0.036, pad=0.04, format=cbfmt, extend='both', extendfrac='auto')

        else:
            cbar = fig.colorbar(panelhandle, fraction=0.036, pad=0.04,
                                extend='both', format=cbfmt, norm=logsc)  # extend='both'
            #cbar1 = ColorbarBase(cbar.ax, cmap=palette, norm=norm, fraction=0.036, pad=0.04, extend='both', extendfrac='auto', format=cbfmt)

        if topodata is not None:
            panelhandle.remove()

        cbar.set_label(label1, fontsize=10)
        cbar.ax.tick_params(labelsize=8)

        parallels = m.drawparallels(getMapLines(bymin, bymax, 3),
                                    labels=[1, 0, 0, 0], linewidth=0.5,
                                    labelstyle='+/-', fontsize=9, xoffset=-0.8,
                                    color='gray', zorder=100.)
        m.drawmeridians(getMapLines(bxmin, bxmax, 3), labels=[0, 0, 0, 1],
                        linewidth=0.5, labelstyle='+/-', fontsize=9,
                        color='gray', zorder=100.)
        for par in parallels:
            try:
                parallels[par][1][0].set_rotation(90)
            except:
                pass

        #draw roads on the map, if they were provided to us
        if maproads is True and roadslist is not None:
            try:
                for road in roadslist:
                    try:
                        xy = list(road['geometry']['coordinates'])
                        roadx, roady = list(zip(*xy))
                        mapx, mapy = m(roadx, roady)
                        m.plot(mapx, mapy, roadcolor, lw=0.5, zorder=9)
                    except:
                        continue
            except Exception as e:
                print(('Failed to plot roads, %s' % e))

        #add city names to map
        if mapcities is True and cityfile is not None:
            try:
                #fontname = 'Helvetica'
                fontsize = 8
                if k == 0:  # Only need to choose cities first time and then apply to rest
                    fcities = bcities.limitByMapCollision(
                        m, fontsize=fontsize)
                    ctlats, ctlons, names = fcities.getCities()
                    cxis, cyis = m(ctlons, ctlats)
                for ctlat, ctlon, cxi, cyi, name in zip(ctlats, ctlons, cxis, cyis, names):
                    m.scatter(ctlon, ctlat, c='k', latlon=True, marker='.',
                              zorder=100000)
                    ax.text(cxi, cyi, name,
                            fontsize=fontsize, zorder=100000)
            except Exception as e:
                print('Failed to plot cities, %s' % e)

        #draw star at epicenter
        plt.sca(ax)
        if edict is not None:
            elat, elon = edict['lat'], edict['lon']
            ex, ey = m(elon, elat)
            plt.plot(ex, ey, '*', markeredgecolor='k', mfc='None', mew=1.0,
                     ms=15, zorder=10000.)

        m.drawmapboundary(fill_color=watercolor)

        m.fillcontinents(color=clear_color, lake_color=watercolor)
        m.drawrivers(color=watercolor)
        ##m.drawcoastlines()

        #draw country boundaries
        m.drawcountries(color=countrycolor, linewidth=1.0)

        #add map scale
        m.drawmapscale((bxmax+bxmin)/2., (bymin+0.1*(bymax-bymin)), clon, clat,
                       np.round((((bxmax-bxmin)*111)/5)/10.)*10,
                       barstyle='simple', zorder=200)

        # Add border
        autoAxis = ax.axis()
        rec = Rectangle((autoAxis[0]-0.7, autoAxis[2]-0.2),
                        (autoAxis[1]-autoAxis[0])+1, (autoAxis[3]-autoAxis[2])+0.4,
                        fill=False, lw=1, zorder=1e8)
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
            print(paramstring)
            fig.text(0.01, 0.015, paramstring, fontsize=fontsizesmallest)
            plt.draw()
        except:
            print('Could not display model parameters')

    if edict is not None:
        eventid = edict['eventid']
    else:
        eventid = ''

    if outfilename is None:
        time1 = datetime.datetime.utcnow().strftime('%d%b%Y_%H%M')
        outfile = os.path.join(outfolder, '%s_%s_%s.pdf' % (eventid, suptitle, time1))
        pngfile = os.path.join(outfolder, '%s_%s_%s.png' % (eventid, suptitle, time1))
    else:
        outfile = os.path.join(outfolder, outfilename + '.pdf')
        pngfile = os.path.join(outfolder, outfilename + '.png')

    filenames = []
    if savepdf is True:
        print('Saving map output to %s' % outfile)
        plt.savefig(outfile, dpi=300)
        filenames.append(outfile)
    if savepng is True:
        print('Saving map output to %s' % pngfile)
        plt.savefig(pngfile)
        filenames.append(pngfile)
    if showplots is True:
        plt.show()
    else:
        plt.close(fig)

    return newgrids, filenames


def interactiveMap(grids, shakefile=None, plotorder=None, inventory_shapefile=None,
                   maskthreshes=None, colormaps=None, scaletype='continuous', lims=None, logscale=False,
                   ALPHA=0.7, isScenario=False, outputdir=None, outfilename=None, tiletype='Stamen Terrain',
                   printparam=False, ds=True, dstype='mean', smcontourfile=None):
    """
    This function creates interactive html plots of mapio grid layers (e.g. liquefaction or
    landslide models with their input layers)

    :param grids: Dictionary of N layers and metadata formatted like:
        maplayers['layer name']={
        'grid': mapio grid2D object,
        'label': 'label for colorbar and top line of subtitle',
        'type': 'output or input to model',
        'description': 'detailed description of layer for subtitle'}.
      Layer names must be unique.
    :type name: Dictionary or Ordered dictionary - import collections;
      grids = collections.OrderedDict()
    :param shakefile: optional ShakeMap file (url or full file path) to extract information for labels and folder names
    :type shakefile: Shakemap Event Dictionary
    :param plotorder: List of keys describing the order to plot the grids, if
      None and grids is an ordered dictionary, it will use the order of the
      dictionary, otherwise it will choose order which may be somewhat random
      but it will always put a probability grid first
    :type plotorder: list
    :param maskthreshes: N x 1 array or list of lower thresholds for masking
      corresponding to order in plotorder or order of OrderedDict if plotorder
      is None. If grids is not an ordered dict and plotorder is not specified,
      this will not work right. If None (default), nothing will be masked
    :param colormaps: List of strings of matplotlib colormaps (e.g. cm.autumn_r)
      corresponding to plotorder or order of dictionary if plotorder is None.
      The list can contain both strings and None e.g. colormaps = ['cm.autumn',
      None, None, 'cm.jet'] and None's will default to default colormap
    :param scaletype: Type of scale for plotting, 'continuous' or 'binned' -
      will be reflected in colorbar
    :type scaletype: string
    :param lims: None or Nx1 list of tuples or numpy arrays corresponding to
      plotorder defining the limits for saturating the colorbar (vmin, vmax) if
      scaletype is continuous or the bins to use (clev) if scaletype if binned.
      The list can contain tuples, arrays, and Nones, e.g. lims = [(0., 10.),
      None, (0.1, 1.5), np.linspace(0., 1.5, 15)]. When None is specified, the
      program will estimate the limits, when an array is specified but the scale
      type is continuous, vmin will be set to min(array) and vmax will be set
      to max(array)
    :param lims: None or Nx1 list of Trues and Falses corresponding to
      plotorder defining whether to use a linear or log scale (log10) for
      plotting the layer. This will be reflected in the labels
    :param logscale: None or Nx1 list of booleans defining whether to use a log scale or not for each layer
    :param ALPHA: Transparency for mapping, if there is a hillshade that will
      plot below each layer, it is recommended to set this to at least 0.7
    :type ALPHA: float
    :param isScenario: Whether this is a scenario (True) or a real event (False)
      (default False)
    :type isScenario: boolean
    :param outputdir: File path for outputting figures, if edict is defined, a
      subfolder based on the event id will be created in this folder. If None,
      will use current directory
    :param outfilename: File name for output without any file extensions
    :param tiletype: folium tile type:
        - "OpenStreetMap"
        - "Mapbox Bright" (Limited levels of zoom for free tiles)
        - "Mapbox Control Room" (Limited levels of zoom for free tiles)
        - "Stamen" (Terrain, Toner, and Watercolor)
        - "Cloudmade" (Must pass API key)
        - "Mapbox" (Must pass API key)
        - "CartoDB" (positron and dark_matter)
    :param printparam: True prints model parameters on figure - NOT IMPLEMENTED YET
    :param ds: True to allow downsampling for display (necessary when arrays
      are quite large, False to not allow) NOT IMPLEMENTED YET
    :param dstype: What function to use in downsampling, options are 'min',
      'max', 'median', or 'mean' NOT IMPLEMENTED YET
    :param smcontourfile: file extension to shakemap contour file to plot
     NOT FUNCTIONAL YET

    :returns:
        * Interactive plot (html file) of all grids listed in plotorder

    """
    plt.ioff()
    clear_color = [0, 0, 0, 0.0]
    if plotorder is None:
        plotorder = grids.keys()

    defaultcolormap = cm.CMRmap_r

    if shakefile is not None:
        edict = ShakeGrid.load(shakefile, adjust='res').getEventDict()
        temp = ShakeGrid.load(shakefile, adjust='res').getShakeDict()
        edict['eventid'] = temp['shakemap_id']
        edict['version'] = temp['shakemap_version']
    else:
        edict = None

    # Get output file location
    if outputdir is None:
        print('No output location given, using current directory for outputs\n')
        outputdir = os.getcwd()
    if edict is not None:
        outfolder = os.path.join(outputdir, edict['event_id'])
    else:
        outfolder = outputdir
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)
    # cbfolder = os.path.join(outfolder, 'cbfolder')
    # if not os.path.isdir(cbfolder):
    #     os.makedirs(cbfolder)

    # ADD IN DOWNSAMPLING CODE FROM MODELMAP HERE
    filenames = []
    for k, key in enumerate(plotorder):

        # Get simplified name of key for file naming
        RIDOF = '[+-]?(?=\d*[.eE])(?=\.?\d)\d*\.?\d*(?:[eE][+-]?\d+)?'
        OPERATORPAT = '[\+\-\*\/]*'
        keyS = re.sub(OPERATORPAT, '', key)
        #remove floating point numbers
        keyS = re.sub(RIDOF, '', keyS)
        #remove parentheses
        keyS = re.sub('[()]*', '', keyS)
        #remove any blank spaces
        keyS = keyS.replace(' ', '')

        grid = grids[key]['grid']

        # get labels and metadata info
        if 'label' in list(grids[key].keys()):
            label1 = grids[key]['label']
        else:
            label1 = key
        try:
            sref = grids[key]['description']['name']
        except:
            sref = ''

        if colormaps is not None and len(colormaps) == len(plotorder):
            if colormaps[k] is not None:
                palette = colormaps[k]
            else:
                palette = defaultcolormap
        else:  # Find preferred default color map for each type of layer
            if 'prob' in key.lower() or 'pga' in key.lower() or \
               'pgv' in key.lower() or 'cohesion' in key.lower() or \
               'friction' in key.lower() or 'fs' in key.lower():
                palette = cm.inferno
            elif 'slope' in key.lower():
                palette = cm.gnuplot2
            elif 'precip' in key.lower():
                palette = cm2.s3pcpn
            else:
                palette = defaultcolormap

        palette.set_bad(clear_color, alpha=0.0)

        dat = grid.getData().copy()

        if scaletype.lower() == 'binned':
            if logscale is not False and len(logscale) == len(plotorder):
                if logscale[k] is True:
                    clev = 10.**(np.arange(np.floor(np.log10(np.nanmin(dat))), np.ceil(np.log10(np.nanmax(dat))), 0.25))
                else:
                    # Find order of range to know how to scale
                    order = np.round(np.log(np.nanmax(dat) - np.nanmin(dat)))
                    if order < 1.:
                        scal = 10**-order
                    else:
                        scal = 1.

                    if lims is None or len(lims) != len(plotorder):
                        clev = (np.linspace(np.floor(scal*np.nanmin(dat)), np.ceil(scal*np.nanmax(dat)), 10))/scal
                    else:
                        if lims[k] is None:
                            clev = (np.linspace(np.floor(scal*np.nanmin(dat)), np.ceil(scal*np.nanmax(dat)), 10))/scal
                        else:
                            clev = lims[k]
            else:
                # Find order of range to know how to scale
                order = np.round(np.log(np.nanmax(dat) - np.nanmin(dat)))
                if order < 1.:
                    scal = 10**-order
                else:
                    scal = 1.

                if lims is None or len(lims) != len(plotorder):
                    clev = (np.linspace(np.floor(scal*np.nanmin(dat)), np.ceil(scal*np.nanmax(dat)), 10))/scal
                else:
                    if lims[k] is None:
                        clev = (np.linspace(np.floor(scal*np.nanmin(dat)), np.ceil(scal*np.nanmax(dat)), 10))/scal
                    else:
                        clev = lims[k]

            # Adjust to colorbar levels
            dat[dat < clev[0]] = clev[0]
            for j, level in enumerate(clev[:-1]):
                dat[(dat >= clev[j]) & (dat < clev[j+1])] = (clev[j] + clev[j+1])/2.
            # So colorbar saturates at top
            dat[dat > clev[-1]] = clev[-1]
            vmin = clev[0]
            vmax = clev[-1]

        else:
            if lims is not None and len(lims) == len(plotorder):
                if lims[k] is None:
                    vmin = np.nanmin(dat)
                    vmax = np.nanmax(dat)
                else:
                    vmin = lims[k][0]
                    vmax = lims[k][-1]
            else:
                vmin = np.nanmin(dat)
                vmax = np.nanmax(dat)
            clev = np.linspace(vmin, vmax, 10)

        if maskthreshes is not None and len(maskthreshes) == len(plotorder):
            if maskthreshes[k] is not None:
                dat[dat <= maskthreshes[k]] = float('NaN')
                dat = np.ma.array(dat, mask=np.isnan(dat))

        #turn data into an RGBA image
        #adjust data so scaled between vmin and vmax and between 0 and 1
        cmap = palette
        dat1 = dat.copy()
        dat1[dat1 < vmin] = vmin  # saturate at ends
        dat1[dat1 > vmax] = vmax
        dat1 = (dat1 - vmin)/(vmax-vmin)
        rgba_img = cmap(dat1)

        gd = grid.getGeoDict()
        minlat = gd.ymin - gd.dy/2.
        minlon = gd.xmin - gd.dx/2.
        maxlat = gd.ymax + gd.dy/2.
        maxlon = gd.xmax + gd.dx/2.

        # if logscale is not False and len(logscale) == len(plotorder):
        #     logsc = None
        #     if logscale[k] is True:
        #         logsc = LogNorm(vmin=vmin, vmax=vmax)
        # else:
        #     logsc = None

        # # Make colorbar figure

        # # This is just a dummy layer that will be deleted to make the colorbar look right
        # panelhandle = plt.imshow(dat1, cmap=palette, vmin=vmin, vmax=vmax, norm=logsc)

        # cbfmt = '%1.1f'
        # if vmax is not None and vmin is not None:
        #     if logscale is not False and len(logscale) == len(plotorder):
        #         if logscale[k] is True:
        #             cbfmt = '%1.0e'
        #     elif (vmax - vmin) < 1.:
        #         cbfmt = '%1.2f'
        #     elif vmax > 5.:  # (vmax - vmin) > len(clev):
        #         cbfmt = '%1.0f'

        # fig = plt.figure(figsize=(8., 2.5))

        # if scaletype.lower() == 'binned':
        #     cbar = fig.colorbar(panelhandle, fraction=0.8, pad=0., orientation='horizontal',
        #                         extend='both', format=cbfmt, spacing='proportional',
        #                         ticks=clev, boundaries=clev)
        # else:
        #     cbar = fig.colorbar(panelhandle, fraction=0.8, pad=0., orientation='horizontal',
        #                         extend='both', format=cbfmt)
        # cbar.set_label(label1, fontsize=14)
        # cbar.ax.tick_params(labelsize=14)
        # plt.axis('off')
        # panelhandle.remove()

        # if edict is not None:
        #     if isScenario:
        #         title = edict['event_description']
        #     else:
        #         timestr = edict['event_timestamp'].strftime('%b %d %Y')
        #         title = 'M%.1f %s v%i - %s' % (edict['magnitude'], timestr, edict['version'], edict['event_description'])
        #     plt.suptitle(title+'\n'+sref, fontsize=16)
        # else:
        #     plt.suptitle(sref, fontsize=16)

        sref_fix = sref
        if sref != '':
            sref_fix = sref_fix.replace(' (', '_')
            sref_fix = sref_fix.replace(')', '')
            sref_fix = sref_fix.replace(' ', '_')

        if outfilename is None:
            outfilename = '%s_%s' % (edict['event_id'], sref_fix)

        # plt.tight_layout()
        # fig.savefig(os.path.join(cbfolder, ('%s_%s_colorbar.png' % (outfilename, keyS))), transparent=True)  # This file has to move with the html files

        if inventory_shapefile is not None:
            reader = shapefile.Reader(inventory_shapefile)
            fields = reader.fields[1:]
            field_names = [field[0] for field in fields]
            buffer1 = []
            for sr in reader.shapeRecords():
                atr = dict(zip(field_names, sr.record))
                geom = sr.shape.__geo_interface__
                style_function = lambda x: {'fillColor': 'none', 'color': 'black', 'weight': 0.7}
                buffer1.append(dict(type="Feature", geometry=geom, properties=atr))

            # create geojson object
            invt = GeoJson({"type": "FeatureCollection", "features": buffer1}, style_function=style_function)

        map1 = folium.Map(location=[(maxlat+minlat)/2., (maxlon+minlon)/2.],
                          tiles=tiletype, zoom_start=9, max_zoom=14, min_lat=minlat, max_lat=maxlat,
                          min_lon=minlon, max_lon=maxlon, prefer_canvas=True)

        #map1.add_children(plugins.HeatMap(zip(lats, lons, dat1), radius=gd.dx))
        img = plugins.ImageOverlay(rgba_img, opacity=ALPHA, bounds=[[minlat, minlon],
                                   [maxlat, maxlon]], mercator_project=True)
        img.layer_name = label1
        map1.add_children(img)

        # plugins.FloatImage(os.path.join(cbfolder, ('%s_%s_colorbar.png' % (outfilename, keyS))), bottom=0, left=1).add_to(map1)
        if scaletype.lower() == 'binned':
            color1 = palette(clev/clev.max())
            cbar = cmb.StepColormap(color1, vmin=vmin, vmax=vmax, index=clev, caption=label1)
        else:
            color1 = [tuple(p) for p in palette._lut]
            cbar = cmb.LinearColormap(color1, vmin=vmin, vmax=vmax, caption=label1)
        map1.add_child(cbar)

        # map1.choropleth(threshold_scale=clev,
        #                 fill_color=color1, fill_opacity=0.7, line_opacity=0.5,
        #                 legend_name=label1,
        #                 reset=True)
        map1.add_child(folium.LatLngPopup())
        map1.add_child(RectangleMarker(bounds=[[minlat, minlon], [maxlat, maxlon]], fill_opacity=0.5,
                       weight=2, fill_color='none'))

        if inventory_shapefile is not None:
            map1.add_child(invt)
            invt.layer_name = 'Inventory'

        folium.LayerControl().add_to(map1)

        if smcontourfile is not None:
            style_function = lambda x: {'fillColor': 'none', 'color': 'white', 'weight': 0.7}
            smc = GeoJson(open(smcontourfile))
            #smc.layer_name = 'ShakeMap Contours'
            map1.add_child(smc)
            #for feature in smc.data['features']:
            #    label = ('%s (%s)') % (feature['properties']['value'], feature['properties']['units'].replace('pct', '%'))
            #    plugins.PolyLineTextPath(feature['geometry']['coordinates'], label,
            #                             center=True, attributes={'fill': 'white', 'font-size': '14'}).add_to(map1)

        #draw star at epicenter
        if edict is not None:
            folium.RegularPolygonMarker(location=[edict['lat'], edict['lon']], popup='Epicenter',
                                        fill_color='#769d96', number_of_sides=4, radius=6).add_to(map1)

        filen = os.path.join(outfolder, '%s_%s.html' % (outfilename, keyS))
        map1.save(filen)
        filenames.append(filen)
        plt.close('all')

    return map1, filenames


def InteractivePage(grids, keys=None, shakefile=None, inventory_shapefile=None,
                    maskthreshes=None, colormaps=None, isScenario=False,
                    zthresh=0, scaletype='continuous', lims=None, logscale=False,
                    ALPHA=0.7, outputdir=None, outfilename=None, tiletype='Stamen Terrain',
                    printparam=False, ds=True, dstype='mean'):
    """embed multiple interactive plots in one page
    """
    pass


def make_hillshade(topogrid, azimuth=315., angle_altitude=50.):
    """
    Computes a hillshade from a digital elevation model. Most of this script
    borrowed from
    <http://geoexamples.blogspot.com/2014/03/shaded-relief-images-using-gdal-python.html>
    last accessed 9/2/2015

    :param topogrid: Digital elevation model
    :type topogrid: Grid2D object
    :param azimuth: azimuth of illumination in degrees
    :type azimuth: float
    :param angle_altitude: altitude angle of illumination
    :type angle_altitude: float

    :returns: Hillshade map layer (Grid2D object)

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
        lon, lat = list(zip(*sequence))
        x, y = m(lon, lat)
        tlist.append(tuple(zip(x, y)))
    polyjson['coordinates'] = tuple(tlist)
    ppolygon = shape(polyjson)
    patch = PolygonPatch(ppolygon, facecolor=facecolor, edgecolor=edgecolor,
                         zorder=zorder, linewidth=lw, fill=True, visible=True)
    return patch


def roundToNearest(value, roundValue=1000):
    """Return the value, rounded to nearest roundValue (defaults to 1000).

    :param value: Value to be rounded.
    :param roundValue: Number to which the value should be rounded.

    :returns: rounded value
    """
    if roundValue < 1:
        ds = str(roundValue)
        nd = len(ds) - (ds.find('.')+1)
        value = value * 10**nd
        roundValue = roundValue * 10**nd
        value = int(round(float(value)/roundValue)*roundValue)
        value = float(value) / 10**nd
    else:
        value = int(round(float(value)/roundValue)*roundValue)
    return value


def floorToNearest(value, floorValue=1000):
    """Return the value, floored to nearest floorValue (defaults to 1000).

    :param value: Value to be floored.
    :param floorValue: Number to which the value should be floored.

    :returns: floored value
    """
    if floorValue < 1:
        ds = str(floorValue)
        nd = len(ds) - (ds.find('.')+1)
        value = value * 10**nd
        floorValue = floorValue * 10**nd
        value = int(np.floor(float(value)/floorValue)*floorValue)
        value = float(value) / 10**nd
    else:
        value = int(np.floor(float(value)/floorValue)*floorValue)
    return value


def ceilToNearest(value, ceilValue=1000):
    """Return the value, ceiled to nearest ceilValue (defaults to 1000).

    :param value: Value to be ceiled.
    :param ceilValue: Number to which the value should be ceiled.

    :returns: ceiling-ed value
    """
    if ceilValue < 1:
        ds = str(ceilValue)
        nd = len(ds) - (ds.find('.')+1)
        value = value * 10**nd
        ceilValue = ceilValue * 10**nd
        value = int(np.ceil(float(value)/ceilValue)*ceilValue)
        value = float(value) / 10**nd
    else:
        value = int(np.ceil(float(value)/ceilValue)*ceilValue)
    return value


if __name__ == '__main__':
    pass
