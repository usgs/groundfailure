#!/usr/bin/env python

# stdlib imports
import os
import gc
import math
import glob
import matplotlib as mpl

import collections
from datetime import datetime
#from configobj import ConfigObj

# third party imports
import matplotlib.cm as cm
import numpy as np
import matplotlib.pyplot as plt
import fiona
from shapely.geometry import mapping, shape
from shapely.geometry import Polygon as PolygonSH
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import cm as cm2
from mpl_toolkits.basemap import maskoceans
from matplotlib.patches import Polygon, Rectangle
from skimage.measure import block_reduce
from descartes import PolygonPatch
from mapio.shake import ShakeGrid


# local imports
from mapio.gmt import GMTGrid
from mapio.gdal import GDALGrid
from mapio.geodict import GeoDict
from mapio.grid2d import Grid2D
from mapio.basemapcity import BasemapCities


# Make fonts readable and recognizable by illustrator
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = \
    ['Helvetica', 'Arial', 'Bitstream Vera Serif', 'sans-serif']
plt.switch_backend('agg')


def modelMap(grids, shakefile=None,
             suptitle=None, inventory_shapefile=None,
             plotorder=None, maskthreshes=None, colormaps=None,
             boundaries=None, zthresh=0, scaletype='continuous', lims=None,
             logscale=False, ALPHA=0.7, maproads=True, mapcities=True,
             isScenario=False, roadfolder=None, topofile=None, cityfile=None,
             oceanfile=None, roadcolor='#6E6E6E', watercolor='#B8EEFF',
             countrycolor='#177F10', outputdir=None, outfilename=None,
             savepdf=True, savepng=True, showplots=False, printparam=False, ds=True,
             dstype='mean', upsample=False):
    """
    Create static maps of mapio grid layers (e.g. liquefaction or
    landslide models with their input layers).

    Args:
        grids (dict): Dictionary of N layers and metadata formatted like:

            .. code-block:: python

                {
                    'grid': mapio grid2D object,
                    'label': 'label for colorbar and top line of subtitle',
                    'type': 'output or input to model',
                    'description': 'description for subtitle'
                }

          Layer names must be unique. All grids must use the same bounds.
        shakefile (str): Optional ShakeMap file (url or full file path) to
            extract information for labels and folder names.
        suptitle (str): This will be displayed at the top of the plots and in
            the figure names.
        inventory_shapefile (str): Path to inventory file.
        plotorder (list): List of keys describing the order to plot the grids,
            if None and grids is an ordered dictionary, it will use the order
            of the dictionary, otherwise it will choose order which may be
            somewhat random but it will always put a probability grid first.
        maskthreshes (array): N x 1 array or list of lower thresholds for
            masking corresponding to order in plotorder or order of OrderedDict
            if plotorder is None. If grids is not an ordered dict and plotorder
            is not specified, this will not work right. If None (default),
            nothing will be masked.
        colormaps (list): List of strings of matplotlib colormaps (e.g.
            cm.autumn_r) corresponding to plotorder or order of dictionary if
            plotorder is None. The list can contain both strings and None e.g.
            colormaps = ['cm.autumn', None, None, 'cm.jet'] and None's will
            default to default colormap.
        boundaries (*): None to show entire study area, 'zoom' to zoom in on
            the area of action (only works if there is a probability layer)
            using ``zthresh`` as a threshold, or a dictionary with keys 'xmin',
            'xmax', 'ymin', and 'ymax'.
        zthresh (float): Threshold for computing zooming bounds, only used if
            boundaries = 'zoom'.
        scaletype (str): Type of scale for plotting, 'continuous' or 'binned'.
            Will be reflected in colorbar.
        lims (*): None or Nx1 list of tuples or numpy arrays corresponding to
            plotorder defining the limits for saturating the colorbar
            (vmin, vmax) if scaletype is continuous or the bins to use (clev)
            if scaletype if binned. The list can contain tuples, arrays, and
            Nones, e.g.

            .. code-block:: python

                [(0., 10.), None, (0.1, 1.5), np.linspace(0., 1.5, 15)]

            When None is specified, the program will
            estimate the limits, when an array is specified but the scale
            type is continuous, vmin will be set to min(array) and vmax will
            be set to max(array).
        logscale (*): boolean to apply log colorbar to all plotted layers or list
            of booleans corresponding to each layer in grid.
        ALPHA (float): Transparency for mapping, if there is a hillshade that
            will plot below each layer, it is recommended to set this to at
            least 0.7.
        maproads (bool): Whether to show roads or not, default True, but
            requires that roadfile is specified and valid to work.
        mapcities (bool): Whether to show cities or not, default True, but
            requires that cityfile is specified and valid to work.
        isScenario (bool): Whether this is a scenario (True) or a real event
            (False) (default False).
        roadfolder (str): Full file path to folder containing road shapefiles.
        topofile (str): Path to topography grid (GDAL compatible). This is
            only needed to make a hillshade if a premade hillshade is not
            specified.
        cityfile (str): Path to Pager file containing city & population
            information.
        oceanfile (str): Path to file ocean information.
        roadcolor (str): Color to use for roads, if plotted, default is
            '#6E6E6E'.
        watercolor (str): Color to use for oceans, lakes, and rivers, default
            is '#B8EEFF'.
        countrycolor (str): Color for country borders, default is '#177F10'.
        outputdir (str): Path for output figures, if edict is defined, a
            subfolder based on the event id will be created in this folder.
            If None, will use current directory.
        outfilename (str): Output file name.
        savepdf (bool): True to save pdf figure.
        savepng (bool): True to save png figure.
        showplots (bool): True to display plots
        printparam (bool): Show parameter values used to run model on plots.
        ds (bool): True to allow downsampling for display (necessary when
            arrays are quite large, False to not allow).
        dstype (str): What function to use in downsampling? Options are 'min',
            'max', 'median', or 'mean'.
        upsample (bool): True to upsample the layer to the DEM resolution for
            better looking hillshades.

    Returns:
        tuple: (newgrids, filenames), where:
            * newgrids: list of downsampled and trimmed version of input grids. If
                no modification was needed for plotting, this will be identical to grids
                but without the metadata.
            * filenames: a list of filenames that were created

    """
    # TODO:
    #     - Change so that all input layers do not have to have the same bounds,
    #       test plotting multiple probability layers, and add option so that if
    #       PDF and PNG aren't output, opens plot on screen using plt.show().

    if suptitle is None:
        suptitle = ' '

    # display fig only if static fig not output
    if not showplots or (not savepdf and not savepng):
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
        print('No output location given, using current directory '
              'for outputs\n')
        outputdir = os.getcwd()
        if edict is not None:
            outfolder = os.path.join(outputdir, edict['event_id'])
        else:
            outfolder = outputdir
    else:
        outfolder = outputdir

    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)

    # Get plotting order, if not specified
    if plotorder is None:
        plotorder = list(grids.keys())

    if colormaps is None:
        colormaps = [None] * len(plotorder)

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
            print('Could not find model layer to use for zoom, using '
                  'default boundaries')
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
                # backwards so it plots right
                lats = np.linspace(ymax, ymin, temp.getGeoDict().ny)
                row, col = np.where(temp.getData() > float(zthresh))
                lonmin = lons[col].min()
                lonmax = lons[col].max()
                latmin = lats[row].min()
                latmax = lats[row].max()

            # dummy fillers, only really care about bounds
            boundaries1 = {'dx': 100, 'dy': 100., 'nx': 100., 'ny': 100}
            # Add padding around zoom limits
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
        if (np.abs(tempgdict.xmin-boundaries['xmin']) < 0.05 and
                np.abs(tempgdict.ymin-boundaries['ymin']) < 0.05 and
                np.abs(tempgdict.xmax-boundaries['xmax']) < 0.05 and
                np.abs(tempgdict.ymax - boundaries['ymax']) < 0.05):
            print('Input boundaries are almost the same as specified '
                  'boundaries, no cutting needed')
            boundaries = tempgdict
            cut = False
        else:
            try:
                if (boundaries['xmin'] > boundaries['xmax'] or
                        boundaries['ymin'] > boundaries['ymax']):
                    print('Input boundaries are not usable, using '
                          'default boundaries')
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
                print('Input boundaries are not usable, using default '
                      'boundaries')
                keytemp = list(grids.keys())
                boundaries = grids[keytemp[0]]['grid'].getGeoDict()
                cut = False

    # Pull out bounds for various uses
    bxmin, bxmax, bymin, bymax = (boundaries.xmin, boundaries.xmax,
                                  boundaries.ymin, boundaries.ymax)

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
            title = 'M%.1f %s v%i - %s' % (edict['magnitude'],
                                           timestr, edict['version'],
                                           edict['event_description'])
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
                newgrids[layer] = {
                    'grid': templayer.cut(cutxmin, cutxmax,
                                          cutymin, cutymax, align=True)}
            except Exception as e:
                print(('Cutting failed, %s, continuing with full layers' % e))
                newgrids = grids
                continue
        del templayer
        gc.collect()
    else:
        newgrids = grids
    tempgdict = newgrids[list(grids.keys())[0]]['grid'].getGeoDict()

    # Upsample layers to same as topofile if desired for better looking
    # hillshades
    if upsample is True and topofile is not None:
        try:
            topodict = GDALGrid.getFileGeoDict(topofile)
            if topodict.dx >= tempgdict.dx or topodict.dy >= tempgdict.dy:
                print('Upsampling not possible, resolution of results '
                      'already smaller than DEM')
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
                    newgrids[layer]['grid'] = \
                        newgrids[layer]['grid'].subdivide(tempgdict2)
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
                lons = block_reduce(np.linspace(xmin, xmax,
                                                layergrid.getGeoDict().nx),
                                    block_size=(divx,),
                                    func=np.mean,
                                    cval=float('nan'))
                if math.isnan(lons[-1]):
                    lons[-1] = lons[-2] + (lons[1]-lons[0])
                lats = block_reduce(np.linspace(ymax, ymin,
                                                layergrid.getGeoDict().ny),
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
        # backwards so it plots right side up
        lats = np.linspace(ymax, ymin, ysize)

    # make meshgrid
    llons1, llats1 = np.meshgrid(lons, lats)

    # See if there is an oceanfile for masking
    bbox = PolygonSH(((cutxmin, cutymin),
                      (cutxmin, cutymax),
                      (cutxmax, cutymax),
                      (cutxmax, cutymin)))
    if oceanfile is not None:
        try:
            f = fiona.open(oceanfile)
            oc = next(f)
            f.close
            shapes = shape(oc['geometry'])
            # make boundaries into a shape
            ocean = shapes.intersection(bbox)
        except:
            print('Not able to read specified ocean file, will use default '
                  'ocean masking')
            oceanfile = None
    if inventory_shapefile is not None:
        try:
            f = fiona.open(inventory_shapefile)
            invshp = list(f.items(bbox=(bxmin, bymin, bxmax, bymax)))
            f.close()
            inventory = [shape(inv[1]['geometry']) for inv in invshp]
        except:
            print('unable to read inventory shapefile specified, will not '
                  'plot inventory')
            inventory_shapefile = None

    # Find cities that will be plotted
    if mapcities is True and cityfile is not None:
        try:
            mycity = BasemapCities.loadFromGeoNames(cityfile=cityfile)
            bcities = mycity.limitByBounds((bxmin, bxmax, bymin, bymax))
            bcities = bcities.limitByGrid(nx=3, ny=3, cities_per_grid=1)
        except:
            print('Could not read in cityfile, not plotting cities')
            mapcities = False
            cityfile = None

    # Load in topofile
    if topofile is not None:
        try:
            topomap = GDALGrid.load(topofile, resample=True, method='linear',
                                    samplegeodict=gdict)
        except:
            topomap = GMTGrid.load(topofile, resample=True, method='linear',
                                   samplegeodict=gdict)
        topodata = topomap.getData().copy()
        # mask oceans if don't have ocean shapefile
        if oceanfile is None:
            topodata = maskoceans(llons1, llats1, topodata, resolution='h',
                                  grid=1.25, inlands=True)
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
        # setup of basemap ('lcc' = lambert conformal conic, or tmerc
        # is transverse mercator).
        # use major and minor sphere radii from WGS84 ellipsoid.
        m = Basemap(llcrnrlon=bxmin, llcrnrlat=bymin,
                    urcrnrlon=bxmax, urcrnrlat=bymax,
                    rsphere=(6378137.00, 6356752.3142),
                    resolution='l', area_thresh=1000., projection='tmerc',
                    lat_0=clat, lon_0=clon, ax=ax)

        x1, y1 = m(llons1, llats1)  # get projection coordinates
        axsize = ax.get_window_extent().transformed(
            fig.dpi_scale_trans.inverted())
        if k == 0:
            wid, ht = axsize.width, axsize.height
        default1 = True
        if len(colormaps) == 1 and len(plotorder) == 1:
            if colormaps[k] is None:
                default1 = True
            else:
                default1 = False
            palette = colormaps[k]
        if (colormaps is not None and
                len(colormaps) == len(plotorder) and
                colormaps[k] is not None):
            palette = colormaps[k]
            default1 = False
        # Find preferred default color map for each type of layer if no
        # colormaps found
        if default1:
            if ('prob' in layer.lower() or 'pga' in layer.lower() or
                    'pgv' in layer.lower() or 'cohesion' in layer.lower() or
                    'friction' in layer.lower() or 'fs' in layer.lower()):
                palette = cm.CMRmap_r
            elif 'slope' in layer.lower():
                palette = cm.gnuplot2
            elif 'precip' in layer.lower():
                palette = cm2.s3pcpn
            else:
                palette = defaultcolormap

        if topodata is not None:
            if k == 0:
                ptopo = m.transform_scalar(
                    np.flipud(topodata), lons+0.5*gdict.dx,
                    lats[::-1]-0.5*gdict.dy, int(np.round(300.*wid)),
                    int(np.round(300.*ht)), returnxy=False, checkbounds=False,
                    order=1, masked=False)
                # use lightsource class to make our shaded topography
                ls = LightSource(azdeg=135, altdeg=45)
                ls1 = LightSource(azdeg=120, altdeg=45)
                ls2 = LightSource(azdeg=225, altdeg=45)
                intensity1 = ls1.hillshade(ptopo, fraction=0.25, vert_exag=1.)
                intensity2 = ls2.hillshade(ptopo, fraction=0.25, vert_exag=1.)
                intensity = intensity1*0.5 + intensity2*0.5

        # Get the data
        dat = layergrid.getData().copy()

        # mask out anything below any specified thresholds

        # if logscale is not False and len(logscale) == len(plotorder):
        #     if logscale[k] is True:
        #         dat = np.log10(dat)
        #         label1 = r'$log_{10}$(' + label1 + ')'

        if scaletype.lower() == 'binned':
            if logscale is not False and len(logscale) == len(plotorder):
                if logscale[k] is True:
                    clev = 10.**(np.arange(np.floor(np.log10(np.nanmin(dat))),
                                           np.ceil(np.log10(np.nanmax(dat))),
                                           0.25))
                else:
                    # Find order of range to know how to scale
                    order = np.round(np.log(np.nanmax(dat) - np.nanmin(dat)))
                    if order < 1.:
                        scal = 10**-order
                    else:
                        scal = 1.

                    if lims is None or len(lims) != len(plotorder):
                        clev = (np.linspace(np.floor(scal*np.nanmin(dat)),
                                            np.ceil(scal*np.nanmax(dat)),
                                            10))/scal
                    else:
                        if lims[k] is None:
                            clev = (np.linspace(np.floor(scal*np.nanmin(dat)),
                                                np.ceil(scal*np.nanmax(dat)),
                                                10))/scal
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
                    clev = (np.linspace(np.floor(scal*np.nanmin(dat)),
                                        np.ceil(scal*np.nanmax(dat)),
                                        10))/scal
                else:
                    if lims[k] is None:
                        clev = (np.linspace(np.floor(scal*np.nanmin(dat)),
                                            np.ceil(scal*np.nanmax(dat)),
                                            10))/scal
                    else:
                        clev = lims[k]

            # Adjust to colorbar levels
            dat[dat < clev[0]] = clev[0]
            for j, level in enumerate(clev[:-1]):
                dat[(dat >= clev[j]) & (dat < clev[j+1])] = \
                    (clev[j] + clev[j+1])/2.
            # So colorbar saturates at top
            dat[dat > clev[-1]] = clev[-1]
            vmin = clev[0]
            vmax = clev[-1]

        else:
            if isinstance(logscale, (bool)):
                # Put it in a list so it won't break things later
                logscale = [logscale]

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

        # Might need to move this up to before downsampling...
        # might give illusion of no hazard in places where there is some that
        # just got averaged out.
        if maskthreshes is not None and len(maskthreshes) == len(plotorder):
            if maskthreshes[k] is not None:
                dat[dat <= maskthreshes[k]] = float('NaN')
                dat = np.ma.array(dat, mask=np.isnan(dat))

        # Mask out cells overlying oceans or block with a shapefile if
        # available
        if oceanfile is None:
            dat = maskoceans(llons1, llats1, dat, resolution='h',
                             grid=1.25, inlands=True)
        else:
            if type(ocean) is PolygonSH:
                ocean = [ocean]
            for oc in ocean:
                patch = getProjectedPatch(oc, m, edgecolor="#006280",
                                          facecolor=watercolor, lw=0.5,
                                          zorder=4.)
                ax.add_patch(patch)

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
                    patch = Polygon(xy, facecolor='none', edgecolor='k',
                                    lw=0.5, zorder=10.)
                    ax.add_patch(patch)
        palette.set_bad(clear_color, alpha=0.0)
        # Plot it up
        dat_im = m.transform_scalar(np.flipud(dat),
                                    lons+0.5*gdict.dx,
                                    lats[::-1]-0.5*gdict.dy,
                                    int(np.round(300.*wid)),
                                    int(np.round(300.*ht)),
                                    returnxy=False,
                                    checkbounds=False,
                                    order=0,
                                    masked=True)
        if isinstance(logscale, (bool)):
            if logscale is True:
                logsc = LogNorm(vmin=vmin, vmax=vmax)
            else:
                logsc = None
        else:
            if len(logscale) > 1:
                if logscale[k] is True:
                    logsc = LogNorm(vmin=vmin, vmax=vmax)
                else:
                    logsc = None
            else:
                if logscale[0] is True:
                    logsc = LogNorm(vmin=vmin, vmax=vmax)
                else:
                    logsc = None

        # Drape over hillshade
        if topodata is not None:
            # turn data into an RGBA image
            cmap = palette
            # adjust data so scaled between vmin and vmax and between 0 and 1
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
                                   vmin=vmin, vmax=vmax, norm=logsc,
                                   interpolation='none')
        else:
            panelhandle = m.imshow(dat_im, cmap=palette, zorder=3., norm=logsc,
                                   vmin=vmin, vmax=vmax, interpolation='none')
        cbfmt = '%1.2f'
        if vmax is not None and vmin is not None:
            if logscale is not False and len(logscale) == len(plotorder):
                if logscale[k] is True:
                    cbfmt = None
            elif (vmax - vmin) <= 1.:
                cbfmt = '%1.2f'
            elif vmax > 5.:  # (vmax - vmin) > len(clev):
                cbfmt = '%1.0f'

        if scaletype.lower() == 'binned':
            cbar = fig.colorbar(panelhandle, spacing='proportional',
                                ticks=clev, boundaries=clev, fraction=0.036,
                                pad=0.04, format=cbfmt, extend='neither',
                                norm=logsc)

        else:
            cbar = fig.colorbar(panelhandle, fraction=0.036, pad=0.04,
                                extend='both', format=cbfmt, norm=logsc)

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

        # draw roads on the map, if they were provided to us
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

        # add city names to map
        if mapcities is True and cityfile is not None:
            try:
                fontsize = 8
                # Only need to choose cities first time and then apply to rest
                if k == 0:
                    fcities = bcities.limitByMapCollision(
                        m, fontsize=fontsize)
                    ctlats, ctlons, names = fcities.getCities()
                    cxis, cyis = m(ctlons, ctlats)
                for ctlat, ctlon, cxi, cyi, name in \
                        zip(ctlats, ctlons, cxis, cyis, names):
                    m.scatter(ctlon, ctlat, c='k', latlon=True, marker='.',
                              zorder=100000)
                    ax.text(cxi, cyi, name,
                            fontsize=fontsize, zorder=100000)
            except Exception as e:
                print('Failed to plot cities, %s' % e)

        # draw star at epicenter
        plt.sca(ax)
        if edict is not None:
            elat, elon = edict['lat'], edict['lon']
            ex, ey = m(elon, elat)
            plt.plot(ex, ey, '*', markeredgecolor='k', mfc='None', mew=1.0,
                     ms=15, zorder=10000.)

        m.drawmapboundary(fill_color=watercolor)

        m.fillcontinents(color=clear_color, lake_color=watercolor)
        m.drawrivers(color=watercolor)

        # draw country boundaries
        m.drawcountries(color=countrycolor, linewidth=1.0)

        # add map scale
        m.drawmapscale((bxmax+bxmin)/2., (bymin+0.1*(bymax-bymin)), clon, clat,
                       np.round((((bxmax-bxmin)*111)/5)/10.)*10,
                       barstyle='simple', zorder=200)

        # Add border
        autoAxis = ax.axis()
        rec = Rectangle((autoAxis[0]-0.7,
                         autoAxis[2]-0.2),
                        (autoAxis[1]-autoAxis[0])+1,
                        (autoAxis[3]-autoAxis[2])+0.4,
                        fill=False,
                        lw=1,
                        zorder=1e8)
        rec = ax.add_patch(rec)
        rec.set_clip_on(False)

        plt.draw()

        if sref is not None:
            label2 = '%s\nsource: %s' % (label1, sref)
        else:
            label2 = label1
        plt.title(label2, axes=ax, fontsize=fontsizesub)

        # draw scenario watermark, if scenario
        if isScenario:
            plt.sca(ax)
            cx, cy = m(clon, clat)
            plt.text(cx, cy, 'SCENARIO', rotation=45, alpha=0.10, size=72,
                     ha='center', va='center', color='red')

        # if ds: # Could add this to print "downsampled" on map
        #     plt.text()

        if k == 1 and rowpan == 1:
            # adjust single level plot
            axsize = ax.get_window_extent().transformed(
                fig.dpi_scale_trans.inverted())
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
        time1 = datetime.utcnow().strftime('%d%b%Y_%H%M')
        outfile = os.path.join(outfolder,
                               '%s_%s_%s.pdf'
                               % (eventid, suptitle, time1))
        pngfile = os.path.join(outfolder,
                               '%s_%s_%s.png'
                               % (eventid, suptitle, time1))
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


def getMapLines(dmin, dmax, nlines):
    """Get equally spaced lat or lon lines for mapping

    Args:
        dmin (float): minimum lat or lon
        dmax (float): maximum lat or lon
        nlines (int): number of lines

    Returns:
        array: location of lines

    """
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
        # make the increment the closest power of 10
        near = np.power(10, round(math.log10(drange)))
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
    """Project a polygon to coordinate system used in Basemap m

    Args:
        polygon: list of shapely polygons to project
        m: Basemap map to project polygon to
        edgecolor: color to use for polygon edge
        facecolor: color to use to fill polygon
        lw: line width for edge of polygon
        zorder: plotting order of polygons

    Returns:
        patch: patch of projected polygons that can be added to map figure
    """
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
    """
    Return the value, rounded to nearest roundValue (defaults to 1000).

    Args:
        value (float): Value to be rounded.
        roundValue (float): Number to which the value should be rounded.

    Returns:
        float: Rounded value.
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
    """
    Return the value, floored to nearest floorValue (defaults to 1000).

    Args:
        value (float): Value to be floored.
        floorValue (float): Number to which the value should be floored.

    Returns:
        float: Floored value.
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
    """
    Return the value, ceiled to nearest ceilValue (defaults to 1000).

    Args:
        value (float): Value to be ceiled.
        ceilValue (float): Number to which the value should be ceiled.

    Returns:
        float: Ceiling-ed value.
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
