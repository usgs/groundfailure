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


def parseMapConfig():
    pass


def modelMap(grids, config, edict=None, modelname=None, maproads=True, mapcities=True, isScenario=False, lowthreshes=None, plotorder=None, colormaps=None, boundaries=None, zthresh=0, scaletype='continuous', lims=None):
    """This function creates maps of mapio grid layers (e.g. liquefaction or landslide models with their input layers) with the same boundaries with many plotting options, some in config file, others in

    :param grids: Grid of layers and metadata formatted like maplayers['layer name']={'grid': mapio grid2D object, 'label': 'label for colorbar and top line of subtitle', 'type': 'output or input to model', 'description': 'detailed description of layer for subtitle'}
    :type name: Ordered dictionary - import collections; grids = collections.OrderedDict()
    :param edict: Optional event dictionary for ShakeMap, used only for title
    :param state: Current state to be in.
    :type state: bool.
    :returns:  int -- the return code.
    :raises: AttributeError, KeyError

    """
    """
    Creates maps of either liquefaction or landslide model with options for showing roads, choosing colors, and showing cities
    config object
    grids must be dictionary of mapio grids, label of dictionary will be label on colorbar key and will appear in title. Easy to build, example -> grids['probability'] = lsgrid - using an ordered dictionary instead of a regular one will mean you don't need plotorder and the order of vlims, colormaps, etc. will be correct
    zthresh = theshold for computing zooming bounds, only used if boundaries = 'zoom'
    All grids must use same exact bounds (FUTURE ALTERATION, ALLOW THEM TO BE DIFFERENT?)
    boundaries = 'zoom' will automatically search for bounds around a probability file using zthresh, else a dictionary with xmin, xmax, ymin, ymax in it to use as boundaries for plotting, when boundaries=None, boundaries are taken from first grid in grids
    lowthreshes = list of lower threshold values for each grid (below this value will be transparent)
    plotorder = list of keys in order they should be plotted, if None, will search for probability grid to plot first and rest will be random
    list of colormaps corresponding to each grid in plotorder, defaults set to search for certain types of outputs using their key and label
    scaletype - binned or continuous scale to apply to all
    lims = if scaletype is continuous, this is a list of tuples of vmin and vmax corresponding to each grid, if None, will search for this info in config file, otherwise will just use automatic bounds
        if scaletype is binned, this is a list of arrays of bin boundaries corresponding to each layer, any layers that you don't want to specify, insert None, or None for all

    """
    # Set some defaults that may be overwritten by config file
    roadcolor = '6E6E6E'
    countrycolor = '177F10'
    defaultcolormap = cm.jet
    watercolor = 'B8EEFF'
    cityfile = None
    roadfolder = None
    hillshade = None
    topofile = None
    ALPHA = 0.7

    if modelname is None:
        modelname = ' '

    # Parse config object
    try:
        config1 = config['MAPDATA']
        if 'topofile' in config1.keys():
            topofile = config1['topofile']
        if 'hillshade' in config1.keys():
            hillshade = config1['hillshade']
        if 'topofile' not in config1.keys() and 'hillshade' not in config1.keys():
            print('No hillshade or topomap provided, you will have a sad flat map\n')
        if 'roadfolder' in config1.keys() and maproads is True:
            roadfolder = config1['roadfolder']
            if os.path.exists(roadfolder) is False:
                print('roadfolder not valid - roads will not be displayed\n')
                maproads = False
        if 'cityfile' in config1.keys() and mapcities is True:
            cityfile = config1['cityfile']
            if os.path.exists(cityfile):
                try:
                    cities = PagerCity(cityfile)
                except Exception as e:
                    print e
                    print('cities file not valid - cities will not be displayed\n')
                    mapcities = False
            else:
                print('cities file not valid - cities will not be displayed\n')
                mapcities = False
        if 'roadcolor' in config1.keys():
            roadcolor = config1['roadcolor']
        if 'countrycolor' in config1.keys():
            countrycolor = config1['countrycolor']
        if 'watercolor' in config1.keys():
            watercolor = config1['watercolor']
        if 'defaultcolormap' in config1.keys():
            defaultcolormap = eval(config1['defaultcolormap'])
        if 'alpha' in config1.keys():
            ALPHA = float(config1['alpha'])
    except Exception as e:
        print e
        print('MAPDATA missing from or misformatted in config, map will be missing components\n')

    countrycolor = '#'+countrycolor
    watercolor = '#'+watercolor
    roadcolor = '#'+roadcolor

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
                llons1 = llons[temp.getData() > zthresh]
                llats1 = llats[temp.getData() > zthresh]
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
        hillsmap = hillshade(topomap, 315, 50)
    else:
        print('no hillshade is possible\n')
        hillsmap = None
        ALPHA = 1.

    # Get output file location
    if 'OUTPUT'in config.keys():
        outfile = config['OUTPUT']['folder']
    else:  # Use current directory
        print('No output location given, using current directory for outputs\n')
        outfile = os.getcwd()
    outfolder = os.path.join(outfile, edict['eventid'])
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)

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
        if colormaps is not None and len(colormaps) == len(grids):
            palette = colormaps[k]
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
        if lowthreshes is not None and len(lowthreshes) == len(grids):
            if lowthreshes[k] is not None:
                dat[dat <= lowthreshes[k]] = float('NaN')
                dat = np.ma.array(dat, mask=np.isnan(dat))
        #datm = np.ma.masked_values(dat, -1.)
        #import pdb; pdb.set_trace()

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
            # extent_xy = (m.xmin, m.xmax, m.ymin, m.ymax)
            # panelhandle = m.imshow(np.flipud(dat1), cmap=palette, vmin=vmin, vmax=vmax, alpha=ALPHA, rasterized=True, zorder=2, extent=extent_xy)
        datm = maskoceans(llons1, llats1, dat, resolution='h', grid=1.25, inlands=True)
        palette.set_bad(clear_color, alpha=0.0)
        panelhandle = m.pcolormesh(x1, y1, datm, linewidth=0.0, cmap=palette, vmin=vmin, vmax=vmax, alpha=ALPHA, rasterized=True, edgecolor='None')
        # add colorbar
        if scaletype.lower() == 'binned':
            cbar = fig.colorbar(panelhandle, spacing='proportional', ticks=clev, boundaries=clev, format='%1.1f', extend='both')
        else:
            cbar = fig.colorbar(panelhandle, fraction=0.036, pad=0.04, extend='both')

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
        if maproads is True:
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
            plt.draw()

        #add city names to map
        if mapcities is True:
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
        xmax, xmin, ymax, ymin = boundaries['xmin'], boundaries['xmax'], boundaries['ymin'], boundaries['ymax']
        clat = ymin + (ymax-ymin)/2.0
        clon = xmin + (xmax-xmin)/2.0
        #m.drawmapscale((xmax+xmin)/2., (ymin+(ymax-ymin)/5.), clon, clat, np.round((((xmax-xmin)*110)/5)/10.)*10, barstyle='simple')

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
    plt.savefig(outfile, dpi=300)
    plt.savefig(pngfile)
    print 'Saving map output to %s' % outfile
    print 'Saving map output to %s' % pngfile


def hillshade(topogrid, azimuth, angle_altitude):
    """
    Most of this script borrowed from http://geoexamples.blogspot.com/2014/03/shaded-relief-images-using-gdal-python.html last accessed 9/2/2015
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
    return 255*(shaded + 1)/2


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


def makeMap_KML(grids, edict, modelname=None, filename=None):
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
