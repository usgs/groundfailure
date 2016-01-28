#!/usr/bin/env python

#stdlib imports
import os.path
import math
import glob
from matplotlib.path import Path
from mpl_toolkits.basemap import maskoceans

#third party imports
import matplotlib.cm as cm
import numpy as np
import matplotlib.pyplot as plt
import fiona
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

#local imports
from mapio.gmt import GMTGrid
from neicmap.city import PagerCity
from neicutil.text import ceilToNearest, floorToNearest, roundToNearest


def parseMapConfig(config):
    # Parse config object
    # ADD PLOTORDER TO CONFIG? OTHER THINGS LIKE COLORMAPS?
    topofile = None
    hillshade = None
    roadfolder = None
    cityfile = None
    roadcolor = '6E6E6E'
    countrycolor = '177F10'
    watercolor = 'B8EEFF'
    ALPHA = 0.7
    outputdir = None

    try:
        config1 = config['MAPDATA']
        if 'topofile' in config1.keys():
            topofile = config1['topofile']
        if 'hillshade' in config1.keys():
            hillshade = config1['hillshade']
        if 'roadfolder' in config1.keys():
            roadfolder = config1['roadfolder']
            if os.path.exists(roadfolder) is False:
                print('roadfolder not valid - roads will not be displayed\n')
                roadfolder = None
        if 'cityfile' in config1.keys():
            cityfile = config1['cityfile']
            if os.path.exists(cityfile):
                try:
                    PagerCity(cityfile)
                except Exception as e:
                    print e
                    print('cities file not valid - cities will not be displayed\n')
                    cityfile = None
            else:
                print('cities file not valid - cities will not be displayed\n')
                cityfile = False
        if 'roadcolor' in config1.keys():
            roadcolor = config1['roadcolor']
        if 'countrycolor' in config1.keys():
            countrycolor = config1['countrycolor']
        if 'watercolor' in config1.keys():
            watercolor = config1['watercolor']
        if 'alpha' in config1.keys():
            ALPHA = float(config1['alpha'])
        outputdir = config['OUTPUT']['folder']
    except Exception as e:
        print e
        print('MAPDATA missing from or misformatted in config')

    countrycolor = '#'+countrycolor
    watercolor = '#'+watercolor
    roadcolor = '#'+roadcolor

    return topofile, hillshade, roadfolder, cityfile, roadcolor, countrycolor, watercolor, ALPHA, outputdir


def modelMap(grids, edict=None, modelname=None, plotorder=None, maskthreshes=None, colormaps=None, boundaries=None, zthresh=0, scaletype='continuous', lims=None, logscale=False, ALPHA=0.7, maproads=True, mapcities=True, isScenario=False, roadfolder=None, topofile=None, hillshade=None, cityfile=None, roadcolor='#6E6E6E', watercolor='#B8EEFF', countrycolor='#177F10', outputdir=None, savepdf=True, savepng=True):

    """
    This function creates maps of mapio grid layers (e.g. liquefaction or landslide models with their input layers)
    All grids must use the same bounds
    TO DO change so that all input layers do not have to have the same bounds, test plotting multiple probability layers, and add option so that if PDF and PNG aren't output, opens plot on screen using plt.show()

    :param grids: Dictionary of N layers and metadata formatted like maplayers['layer name']={'grid': mapio grid2D object, 'label': 'label for colorbar and top line of subtitle', 'type': 'output or input to model', 'description': 'detailed description of layer for subtitle'}. Layer names must be unique.
    :type name: Dictionary or Ordered dictionary - import collections; grids = collections.OrderedDict()
    :param edict: Optional event dictionary for ShakeMap, used for title and naming output folder
    :type edit: Shakemap Event Dictionary
    :param modelname: Name of model being run, this will be displayed at the top of the plots and in the figure names
    :type modelname: string
    :param plotorder: List of keys describing the order to plot the grids, if None and grids is an ordered dictionary, it will use the order of the dictionary, otherwise it will choose order which may be somewhat random but it will always put a probability grid first
    :type plotorder: list
    :param maskthreshes: N x 1 array or list of lower thresholds for masking corresponding to order in plotorder or order of OrderedDict if plotorder is None. If grids is not an ordered dict and plotorder is not specified, this will not work right. If None (default), nothing will be masked
    :param colormaps: List of strings of matplotlib colormaps (e.g. cm.autumn_r) corresponding to plotorder or order of dictionary if plotorder is None. The list can contain both strings and None e.g. colormaps = ['cm.autumn', None, None, 'cm.jet'] and None's will default to default colormap
    :param boundaries: None to show entire study area, 'zoom' to zoom in on the area of action (only works if there is a probability layer) using zthresh as a threshold, or a dictionary defining lats and lons in the form of boundaries['xmin'] = minlon, boundaries['xmax'] = maxlon, boundaries['ymin'] = min lat, boundaries['ymax'] = max lat
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
    :param hillshade: Full file path to hillshade grid (GDAL compatible)
    :type hillshade: string
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


    :returns:  PDF's or PNG's
    """
    if modelname is None:
        modelname = ' '

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

    # Get boundaries to use for all plots
    if boundaries is None:
        keytemp = grids.keys()
        boundaries = grids[keytemp[0]]['grid'].getGeoDict()
    elif boundaries == 'zoom':
        # Find probability layer (will just take the maximum bounds if there is more than one)
        keytemp = grids.keys()
        key1 = [key for key in keytemp if 'prob' in key.lower()]
        if len(key1) == 0:
            print('Could not find probability layer to use for zoom, using default boundaries')
            keytemp = grids.keys()
            boundaries = grids[keytemp[0]].getGeoDict()
        else:
            lonmax = -1.e10
            lonmin = 1.e10
            latmax = -1.e10
            latmin = 1.e10
            for key in key1:
                # get lat lons of areas affected and add, if no areas affected, switch to shakemap boundaries
                temp = grids[key]['grid']
                xmin, xmax, ymin, ymax = temp.getBounds()
                lons = np.arange(xmin, xmax+temp.getGeoDict()['xdim'], temp.getGeoDict()['xdim'])
                lons = lons[:temp.getGeoDict()['ncols']]  # make sure right length
                lats = np.arange(ymax, ymin-temp.getGeoDict()['ydim'], -temp.getGeoDict()['ydim'])  # backwards so it plots right
                lats = lats[:temp.getGeoDict()['nrows']]
                llons, llats = np.meshgrid(lons, lats)  # make meshgrid
                llons1 = llons[temp.getData() > float(zthresh)]
                llats1 = llats[temp.getData() > float(zthresh)]
                if llons1.min() < lonmin:
                    lonmin = llons1.min()
                if llons1.max() > lonmax:
                    lonmax = llons1.max()
                if llats1.min() < latmin:
                    latmin = llats1.min()
                if llats1.max() > latmax:
                    latmax = llats1.max()
            boundaries1 = {}
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
            boundaries = boundaries1
    else:
        try:
            if boundaries['xmin'] > boundaries['xmax'] or boundaries['ymin'] > boundaries['ymax']:
                print('Input boundaries are not usable, using default boundaries')
                keytemp = grids.keys()
                boundaries = grids[keytemp[0]].getGeoDict()
        except:
            print('Input boundaries are not usable, using default boundaries')
            keytemp = grids.keys()
            boundaries = grids[keytemp[0]].getGeoDict()

    # Load in topofile and load in or compute hillshade, resample to same grid as first layer NEED ERROR HANDLING TO MAKE SURE ALL ARE SAME BOUNDARIES
    gdict = grids[grids.keys()[0]]['grid'].getGeoDict()
    if hillshade is not None:
        hillsmap = GMTGrid.load(hillshade, resample=True, method='linear', preserve='shape', samplegeodict=gdict)
    elif topofile is not None and hillshade is None:
        topomap = GMTGrid.load(topofile, resample=True, method='linear', preserve='shape', samplegeodict=gdict)
        hillsmap = make_hillshade(topomap, 315, 50)
    else:
        print('no hillshade is possible\n')
        hillsmap = None
        ALPHA = 1.

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
    fig.set_figheight(rowpan*5.1)

    # Need to update naming to reflect the shakemap version once can get getHeaderData to work, add edict['version'] back into title, maybe shakemap id also?
    if isScenario:
        title = edict['event_description']
    else:
        timestr = edict['event_timestamp'].strftime('%b %d %Y')
        title = 'M%.1f %s v%i - %s' % (edict['magnitude'], timestr, edict['version'], edict['event_description'])
    plt.suptitle(title+'\n'+'model:'+modelname)

    clear_color = [0, 0, 0, 0.0]

    if plotorder is None:
        plotorder = grids.keys()
    val = 1
    for k, layer in enumerate(plotorder):
        layergrid = grids[layer]['grid']
        if 'label' in grids[layer].keys():
            label1 = grids[layer]['label']
        else:
            label1 = layer
        ax = fig.add_subplot(rowpan, colpan, val)
        val += 1
        xmin, xmax, ymin, ymax = boundaries['xmin'], boundaries['xmax'], boundaries['ymin'], boundaries['ymax']
        clat = ymin + (ymax-ymin)/2.0
        clon = xmin + (xmax-xmin)/2.0
        # setup of basemap ('lcc' = lambert conformal conic).
        # use major and minor sphere radii from WGS84 ellipsoid.
        m = Basemap(llcrnrlon=xmin, llcrnrlat=ymin, urcrnrlon=xmax, urcrnrlat=ymax,
                    rsphere=(6378137.00, 6356752.3142),
                    resolution='i', area_thresh=1000., projection='lcc',
                    lat_1=clat, lon_0=clon, ax=ax)

        dat = layergrid.getData().copy()
        if colormaps is not None and len(colormaps) == len(grids) and colormaps[k] is not None:
            palette = eval(colormaps[k])
        else:  # Find preferred default color map for each type of layer
            if 'prob' in layer.lower() or 'pga' in layer.lower() or 'pgv' in layer.lower() or 'cohesion' in layer.lower() or 'friction' in layer.lower() or 'fs' in layer.lower():
                palette = cm.jet
            elif 'slope' in layer.lower():
                palette = cm.gnuplot2
            elif 'precip' in layer.lower():
                palette = cm.s3pcpn
            else:
                palette = defaultcolormap

        xmin, xmax, ymin, ymax = layergrid.getBounds()
        lons = np.arange(xmin, xmax+layergrid.getGeoDict()['xdim'], layergrid.getGeoDict()['xdim'])
        lons = lons[:layergrid.getGeoDict()['ncols']]  # make sure it's the right length, sometimes it is one too long, sometimes it isn't because of GMTGrid issue
        lats = np.arange(ymax, ymin-layergrid.getGeoDict()['ydim'], -layergrid.getGeoDict()['ydim'])  # backwards so it plots right side up
        lats = lats[:layergrid.getGeoDict()['nrows']]  # make sure it's the right length, sometimes it is one too long, sometimes it isn't because of GMTGrid issue
        #make meshgrid
        llons1, llats1 = np.meshgrid(lons, lats)
        x1, y1 = m(llons1, llats1)  # get projection coordinates

        if hillsmap is not None:
            hillshm = hillsmap.getData()
            hillshm = maskoceans(llons1, llats1, hillshm, resolution='h', grid=1.25, inlands=True)
            m.pcolormesh(x1, y1, hillshm/np.abs(hillshm).max(), cmap='Greys', linewidth=0., rasterized=True, vmin=0., vmax=3., edgecolors='none', zorder=1)
            plt.draw()

        # mask out anything below any specified thresholds
        if maskthreshes is not None and len(maskthreshes) == len(grids):
            if maskthreshes[k] is not None:
                dat[dat <= maskthreshes[k]] = float('NaN')
                dat = np.ma.array(dat, mask=np.isnan(dat))

        if logscale is not None and len(logscale) == len(grids):
            if logscale[k] is True:
                dat = np.log10(dat)
                label1 = r'$log_{10}$(' + label1 + ')'
        if scaletype.lower() == 'binned':
            if lims is None or len(lims) != len(grids):
                clev = np.linspace(np.floor(dat.min()), np.ceil(dat.max()), 10)
            else:
                if lims[k] is None:
                    clev = np.linspace(np.floor(dat.min()), np.ceil(dat.max()), 10)
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
            if lims is not None and len(lims) == len(grids):
                if lims[k] is None:
                    vmin = None
                    vmax = None
                else:
                    vmin = lims[k][0]
                    vmax = lims[k][-1]
            else:
                vmin = None
                vmax = None
            # Tried to use imshow, but doesn't quite work right yet, not lining up properly so still using pcolormesh
            # interpolate to a rectangular map projection grid.
            # dat1 = m.transform_scalar(datm, lons, np.flipud(lats), 500, 500, returnxy=False, order=1, checkbounds=True, masked=True)
            # panelhandle = m.imshow(np.flipud(dat1), cmap=palette, vmin=vmin, vmax=vmax, alpha=ALPHA, rasterized=True, zorder=2)

        # Mask out cells overlying oceans
        datm = maskoceans(llons1, llats1, dat, resolution='h', grid=1.25, inlands=True)
        palette.set_bad(clear_color, alpha=0.0)
        # Plot it up
        panelhandle = m.pcolormesh(x1, y1, datm, linewidth=0., cmap=palette, vmin=vmin, vmax=vmax, alpha=ALPHA, rasterized=True)
        panelhandle.set_edgecolors('face')
        # add colorbar
        if scaletype.lower() == 'binned':
            cbar = fig.colorbar(panelhandle, spacing='proportional', ticks=clev, boundaries=clev, fraction=0.036, pad=0.04, format='%1.1f', extend='both')
        else:
            cbar = fig.colorbar(panelhandle, fraction=0.036, pad=0.04, extend='both', format='%1.1f')

        cbar.set_label(label1, fontsize=10)
        cbar.ax.tick_params(labelsize=8)

        parallels = m.drawparallels(getMapLines(ymin, ymax, 3), labels=[1, 0, 0, 0], linewidth=0.5, labelstyle='+/-', fontsize=6, xoffset=-0.8, color='gray')
        m.drawmeridians(getMapLines(xmin, xmax, 3), labels=[0, 0, 0, 1], linewidth=0.5, labelstyle='+/-', fontsize=6, color='gray')
        for par in parallels:
            try:
                parallels[par][1][0].set_rotation(90)
            except:
                pass

        #draw roads on the map, if they were provided to us
        if maproads is True and roadfolder is not None:
            try:
                roadslist = []
                for folder in os.listdir(roadfolder):
                    road1 = os.path.join(roadfolder, folder)
                    shpfiles = glob.glob(os.path.join(road1, '*.shp'))
                    if len(shpfiles):
                        shpfile = shpfiles[0]
                        f = fiona.open(shpfile)
                        shapes = list(f.items(bbox=(boundaries['xmin'], boundaries['ymin'], boundaries['xmax'], boundaries['ymax'])))
                        for shapeid, shapedict in shapes:
                            roadslist.append(shapedict)
                        f.close()
                for road in roadslist:
                    xy = list(road['geometry']['coordinates'])
                    roadx, roady = zip(*xy)
                    mapx, mapy = m(roadx, roady)
                    m.plot(mapx, mapy, roadcolor, lw=0.5, zorder=9)
            except Exception as e:
                print('Failed to plot roads, %s' % e)

        #add city names to map
        if mapcities is True and cityfile is not None:
            try:
                dmin = 0.04*(m.ymax-m.ymin)
                xyplotted = []
                cities = PagerCity(cityfile)
                #Find cities within bounding box
                boundcity = cities.findCitiesByRectangle(bounds=(boundaries['xmin'], boundaries['xmax'], boundaries['ymin'], boundaries['ymax']))
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
                    if not dist or np.min(dist) > dmin:
                        m.scatter(cit['lon'], cit['lat'], c='k', latlon=True, marker='.', zorder=100000)
                        ax.text(xi, yi, cit['name'], ha='right', va='top', fontsize=8, zorder=100000)
                        xyplotted.append((xi, yi))
            except Exception as e:
                print('Failed to plot cities, %s' % e)

        #draw star at epicenter
        plt.sca(ax)
        elat, elon = edict['lat'], edict['lon']
        ex, ey = m(elon, elat)
        plt.plot(ex, ey, '*', markeredgecolor='k', mfc='None', mew=1.0, ms=12)

        m.drawmapboundary(fill_color=watercolor)

        m.fillcontinents(color=clear_color, lake_color=watercolor)
        m.drawrivers(color=watercolor)

        #draw country boundaries
        m.drawcountries(color=countrycolor, linewidth=1.0)

        #add map scale
        # xmax, xmin, ymax, ymin = boundaries['xmin'], boundaries['xmax'], boundaries['ymin'], boundaries['ymax']
        # clat = ymin + (ymax-ymin)/2.0
        # clon = xmin + (xmax-xmin)/2.0
        # #m.drawmapscale((xmax+xmin)/2., (ymin+(ymax-ymin)/5.), clon, clat, np.round((((xmax-xmin)*110)/5)/10.)*10, barstyle='simple')
        # m.drawmapscale((xmax+xmin)/2., ymin, clon, clat, np.round((((xmax-xmin)*110)/5)/10.)*10, barstyle='simple')

        plt.title(label1, axes=ax)

        #draw scenario watermark, if scenario
        if isScenario:
            plt.sca(ax)
            cx, cy = m(clon, clat)
            plt.text(cx, cy, 'SCENARIO', rotation=45, alpha=0.10, size=72, ha='center', va='center', color='red')

    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    outfile = os.path.join(outfolder, '%s_%s.pdf' % (edict['eventid'], modelname))
    pngfile = os.path.join(outfolder, '%s_%s.png' % (edict['eventid'], modelname))
    if savepdf is True:
        print 'Saving map output to %s' % outfile
        plt.savefig(outfile, dpi=300)
    if savepng is True:
        print 'Saving map output to %s' % pngfile
        plt.savefig(pngfile)


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
    hillshade = topogrid.copy()
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


def makeMap_OSM(grids, edict, modelname=None):
    """
    Plot results over Open Street Map base
    """
    pass


def comparisonMap():
    """
    Compare output probability maps from different models, compile statistics
    """
    pass


def compareModelInventory():
    pass


def saveGrid():
    pass


if __name__ == '__main__':
    pass
