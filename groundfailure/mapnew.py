#!/usr/bin/env python

#stdlib imports
import os.path
import math
from configobj import ConfigObj
import glob

#third party imports
import matplotlib.cm as cm
import numpy as np
import matplotlib.pyplot as plt
import fiona
from mpl_toolkits.basemap import Basemap

#local imports
from mapio.gmt import GMTGrid
from neicmap.city import PagerCity
from neicutil.text import ceilToNearest, floorToNearest, roundToNearest

SEA_LEVEL = 0


def hillshade(topogrid, azimuth, angle_altitude):
    """
    Most of this script borrowed from http://geoexamples.blogspot.com/2014/03/shaded-relief-images-using-gdal-python.html last accessed 9/2/2015
    """
    topotmp = topogrid.getData().copy()
    #make a masked array
    topotmp = np.ma.array(topotmp)
    topodat = np.ma.masked_where(np.isnan(topotmp), topotmp)
    topodat = np.ma.masked_where(topodat == SEA_LEVEL, topodat)
    x, y = np.gradient(topodat)
    slope = np.pi/2. - np.arctan(np.sqrt(x*x + y*y))
    aspect = np.arctan2(-x, y)
    azimuthrad = azimuth*np.pi / 180.
    altituderad = angle_altitude*np.pi / 180.
    shaded = np.sin(altituderad) * np.sin(slope) + np.cos(altituderad) * np.cos(slope) * np.cos(azimuthrad - aspect)
    return 255*(shaded + 1)/2


def getMapLines(dmin, dmax):
    NLINES = 4
    drange = dmax-dmin
    if drange > 4:
        near = 1
    else:
        if drange >= 0.5:
            near = 0.25
        else:
            near = 0.125
    inc = roundToNearest(drange/NLINES, near)
    if inc == 0:
        near = np.power(10, round(math.log10(drange)))  # make the increment the closest power of 10
        inc = ceilToNearest(drange/NLINES, near)
        newdmin = floorToNearest(dmin, near)
        newdmax = ceilToNearest(dmax, near)
    else:
        newdmin = ceilToNearest(dmin, near)
        newdmax = floorToNearest(dmax, near)
    darray = np.arange(newdmin, newdmax+inc, inc)
    if darray[-1] > dmax:
        darray = darray[0:-1]
    return darray


def latstr(parallel):
    if parallel < 0:
        parstr = '%.2f' % (-1*parallel) + '$\degree$ S'
    else:
        parstr = '%.2f' % (parallel) + '$\degree$ N'
    return parstr


def lonstr(meridian):
    if meridian < 0:
        merstr = '%.2f' % (-1*meridian) + '$\degree$ W'
    else:
        merstr = '%.2f' % (meridian) + '$\degree$ E'
    return merstr


def getMapTicks(m, xmin, xmax, ymin, ymax):
    meridians = getMapLines(xmin, xmax)
    parallels = getMapLines(ymin, ymax)
    #do tick stuff
    xlabels = [lonstr(mer) for mer in meridians]
    ylabels = [latstr(par) for par in parallels]
    xticks = []
    yticks = []
    for i in range(0, len(meridians)):
        lat = ymin
        lon = meridians[i]
        x, y = m(lon, lat)
        xticks.append(x)
    for i in range(0, len(parallels)):
        lon = xmin
        lat = parallels[i]
        x, y = m(lon, lat)
        yticks.append(y)

    return (xticks, xlabels, yticks, ylabels)


def makeMap_OSM(grids, edict, modelname=None):
    """
    Plot results over Open Street Map base
    """
    pass


def makeMap_KML(grids, edict, modelname=None, filename=None):
    pass


def makeMap(grids, edict, configfile, modelname=None, maproads=True, mapcities=True, isScenario=False, vlims=None, lowthreshes=None, plotorder=None, colormaps=None, boundaries=None):
    """
    Creates single map of either liquefaction or landslide model with options for showing roads, choosing colors, and showing cities
    grids must be dictionary of mapio grids, label of dictionary will be label on colorbar key and will appear in title. Easy to build, example -> grids['probability'] = lsgrid
    All grids must use same exact bounds (FUTURE ALTERATION, ALLOW THEM TO BE DIFFERENT?)
    boundaries = dictionary with xmin, xmax, ymin, ymax in it to use as boundaries for plotting, when boundaries=None, boundaries are taken from first grid in grids
    vlims = list of tuples of vmin and vmax corresponding to each grid
    lowthreshes = list of lower threshold values for each grid (below this value won't be plotted)
    plotorder = list of keys in order they should be plotted
    list of colormaps corresponding to each grid, None uses default jet

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
        modelname = 'unknown'

    # Parse config file
    try:
        config = ConfigObj(configfile)['MAPDATA']
        if 'topofile' in config.keys():
            topofile = config['topofile']
        if 'hillshade' in config.keys():
            hillshade = config['hillshade']
        if 'topofile' not in config.keys() and 'hillshade' not in config.keys():
            print('No hillshade or topomap provided, you will have a sad flat map\n')
        if 'roadfolder' in config.keys() and maproads is True:
            roadfolder = config['roadfolder']
            if os.path.exists(roadfolder) is False:
                print('roadfolder not valid - roads will not be displayed\n')
                maproads = False
        if 'cityfile' in config.keys() and mapcities is True:
            cityfile = config['cityfile']
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
        if 'roadcolor' in config.keys():
            roadcolor = config['roadcolor']
        if 'countrycolor' in config.keys():
            countrycolor = config['countrycolor']
        if 'watercolor' in config.keys():
            watercolor = config['watercolor']
        if 'defaultcolormap' in config.keys():
            defaultcolormap = config['defaultcolormap']
        if 'alpha' in config.keys():
            ALPHA = config['alpha']
    except Exception as e:
        print e
        print('MAPDATA missing from or misformatted in configfile, map will be missing most components\n')

    countrycolor = '#'+countrycolor
    watercolor = '#'+watercolor
    roadcolor = '#'+roadcolor

    # Get boundaries to use for all plots
    if boundaries is None:
        keytemp = grids.keys()
        boundaries = grids[keytemp[0]].getGeoDict()

    # Load in topofile and load in or compute hillshade
    if topofile is not None and hillshade is None:
        gdict = GMTGrid.getBoundsWithin(topofile, boundaries)
        # if no hillshade expand a little to avoid edge cutoff
        lonrange = gdict['xmax'] - gdict['xmin']
        latrange = gdict['ymax'] - gdict['ymin']
        gdict['xmin'] = gdict['xmin'] - lonrange*0.2
        gdict['xmax'] = gdict['xmax'] + lonrange*0.2
        gdict['ymin'] = gdict['ymin'] - latrange*0.2
        gdict['ymax'] = gdict['ymax'] + latrange*0.2
        topomap = GMTGrid.load(topofile, samplegeodict=gdict)
    if hillshade is not None:
        gdict = GMTGrid.getBoundsWithin(hillshade, boundaries)
        # expand a little to avoid edge cutoff
        lonrange = gdict['xmax'] - gdict['xmin']
        latrange = gdict['ymax'] - gdict['ymin']
        gdict['xmin'] = gdict['xmin'] - lonrange*0.2
        gdict['xmax'] = gdict['xmax'] + lonrange*0.2
        gdict['ymin'] = gdict['ymin'] - latrange*0.2
        gdict['ymax'] = gdict['ymax'] + latrange*0.2
        hillsmap = GMTGrid.load(hillshade, samplegeodict=gdict)
    elif topofile is not None:
            hillsmap = hillshade(topomap, 315, 50)
    else:
        print('no hillshade is possible\n')
        hillsmap = None
        ALPHA = 1.
    if hillsmap is not None:  # Make mesh for later
        xmin, xmax, ymin, ymax = hillsmap.getBounds()
        lons = np.arange(xmin, xmax, hillsmap.getGeoDict()['xdim'])
        lons = lons[:hillsmap.getGeoDict()['ncols']]  # make sure right length
        lats = np.arange(ymax, ymin, -hillsmap.getGeoDict()['ydim'])  # backwards so it plots right
        lats = lats[:hillsmap.getGeoDict()['nrows']]
        llons, llats = np.meshgrid(lons, lats)  # make meshgrid

    # Get output file location
    if 'OUTPUT'in ConfigObj(configfile).keys():
        outfile = ConfigObj(configfile)['OUTPUT']['folder']
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
    elif numpanels == 4:
        rowpan = np.ceil(numpanels/2.)
        colpan = 2
        fig.set_figwidth(10)
    else:
        rowpan = np.ceil(numpanels/3.)
        colpan = 3
        fig.set_figwidth(15)
    fig.set_figheight(rowpan*5.1)
    if isScenario:
        title = edict['loc']
    else:
        timestr = edict['time'].strftime('%b %d %Y')
        title = 'M%.1f %s v%i - %s' % (edict['mag'], timestr, edict['version'], edict['loc'])
    plt.suptitle(title+'\n'+'model:'+modelname)

    clear_color = [0, 0, 0, 0.0]

    if plotorder is None:
        plotorder = grids.keys()
    val = 1
    for k, layer in enumerate(plotorder):
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
        if hillsmap is not None:
            x, y = m(llons, llats)  # get projection coordinates
            if topofile is not None:
                topomap5 = GMTGrid.load(topofile, resample=True, samplegeodict=hillsmap.getGeoDict())
                hillsh = hillsmap.getData()
                hillsh[topomap5.getData() == SEA_LEVEL] = 0
                hillshm = np.ma.masked_equal(hillsh, 0)
            else:
                hillshm = hillsmap.getData()
            m.pcolormesh(x, y, hillshm/np.abs(hillshm).max(), cmap='Greys', lw=0, rasterized=True, vmin=0., vmax=3.)
            plt.draw()

        dat = grids[layer].getData().copy()
        if colormaps is not None and len(colormaps) == len(grids):
            palette = colormaps[k]
        else:
            palette = defaultcolormap
        if topofile is not None:
            # resample topomap to the same as current layer
            topomap2 = GMTGrid.load(topofile, resample=True, method='linear', samplegeodict=grids[layer].getGeoDict())
            iwater = np.where(topomap2.getData() == SEA_LEVEL)
            dat[iwater] = 0
        if lowthreshes is not None and len(lowthreshes) == len(grids):
            dat[dat <= lowthreshes[k]] = 0
        else:
            dat[dat <= 0.] = 0
        datm = np.ma.masked_equal(dat, 0)
        palette.set_bad(clear_color, alpha=0.0)

        xmin, xmax, ymin, ymax = grids[layer].getBounds()
        lons = np.arange(xmin, xmax, grids[layer].getGeoDict()['xdim'])
        lons = lons[:grids[layer].getGeoDict()['ncols']]  # make sure it's the right length, sometimes it is one too long, sometimes it isn't because of GMTGrid issue
        lats = np.arange(ymax, ymin, -grids[layer].getGeoDict()['ydim'])  # backwards so it plots right side up
        lats = lats[:grids[layer].getGeoDict()['nrows']]  # make sure it's the right length, sometimes it is one too long, sometimes it isn't because of GMTGrid issue
        #make meshgrid
        llons1, llats1 = np.meshgrid(lons, lats)
        x1, y1 = m(llons1, llats1)  # get projection coordinates
        if vlims is not None and len(vlims) == len(grids):
            panelhandle = m.pcolormesh(x1, y1, datm, lw=0, cmap=palette, vmin=vlims[k][0], vmax=vlims[k][1], alpha=ALPHA, rasterized=True)
        else:
            panelhandle = m.pcolormesh(x1, y1, datm, lw=0, cmap=palette, alpha=ALPHA, rasterized=True)

        # add colorbar
        cbar = fig.colorbar(panelhandle, fraction=0.046, pad=0.04)
        cbar.set_label(layer, fontsize=10)
        cbar.ax.tick_params(labelsize=8)

        #this stuff has to be added after something has been rendered on the map
        #draw the map ticks on outside of all edges
        fig.canvas.draw()  # have to do this for tick stuff to show
        xticks, xlabels, yticks, ylabels = getMapTicks(m, boundaries['xmin'], boundaries['xmax'], boundaries['ymin'], boundaries['ymax'])
        plt.sca(ax)
        plt.tick_params(axis='both', direction='in', right='on', colors='gray')
        plt.xticks(xticks, xlabels, size=6)
        plt.yticks(yticks, ylabels, size=6)
        for tick in ax.axes.yaxis.get_major_ticks():
            tick.set_pad(-38)
            tick.label2.set_horizontalalignment('right')
        for tick in ax.axes.xaxis.get_major_ticks():
            tick.set_pad(-15)
            tick.label2.set_verticalalignment('top')
        [j.set_color("gray") for j in plt.gca().get_xticklabels()]
        [j.set_color("gray") for j in plt.gca().get_yticklabels()]

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
                m.plot(mapx, mapy, roadcolor, lw=0.5)
            plt.draw()

        #add city names to map with population >50,000 (add option later)
        if mapcities is True:
            dmin = 0.04*(m.ymax-m.ymin)
            xyplotted = []
            cities = PagerCity(cityfile)
            #Find cities within bounding box
            boundcity = cities.findCitiesByRectangle(bounds=(boundaries['xmin'], boundaries['xmax'], boundaries['ymin'], boundaries['ymax']))
            #Just keep 5 biggest cities
            thresh = sorted([cit['pop'] for cit in boundcity])[-5]
            plotcity = [cit for cit in boundcity if cit['pop'] >= thresh]
            #For cities that are more than one xth of the xwidth apart, keep only the larger one
            pass  # do later
            #Plot cities
            for cit in plotcity:  # should sort so it plots them in order of population so larger cities are preferentially plotted - do later
                xi, yi = m(cit['lon'], cit['lat'])
                dist = [np.sqrt((xi-x0)**2+(yi-y0)**2) for x0, y0 in xyplotted]
                if not dist or np.min(dist) > dmin:
                    m.scatter(cit['lon'], cit['lat'], c='k', latlon=True, marker='.')
                    ax.text(xi, yi, cit['name'], ha='right', va='top', fontsize=8)
                    xyplotted.append((xi, yi))

            #draw a map boundary, fill in oceans with water
            m.drawmapboundary(fill_color=watercolor)

            # add the subtitle
            #ax.set_title(layer) # DONT NEED BECAUSE COLORBAR IS LABELED

            #draw star at epicenter
            plt.sca(ax)
            elat, elon = edict['epicenter']
            ex, ey = m(elon, elat)
            plt.plot(ex, ey, '*', markeredgecolor='k', mfc='None', mew=1.5, ms=24)

        #fill in the lakes and rivers
        m.fillcontinents(color=clear_color, lake_color=watercolor)
        m.drawrivers(color=watercolor)

        #draw country boundaries
        m.drawcountries(color=countrycolor, linewidth=1.0)

        #add map scale
        xmax, xmin, ymax, ymin = boundaries['xmin'], boundaries['xmax'], boundaries['ymin'], boundaries['ymax']
        clat = ymin + (ymax-ymin)/2.0
        clon = xmin + (xmax-xmin)/2.0
        m.drawmapscale((xmax+xmin)/2., (ymin+(ymax-ymin)/10.), clon, clat, np.round((((xmax-xmin)*110)/5)/10.)*10, barstyle='simple')
        #draw coastlines
        m.drawcoastlines(color='#476C91', linewidth=0.5)

        #draw scenario watermark, if scenario
        if isScenario:
            plt.sca(ax)
            cx, cy = m(clon, clat)
            plt.text(cx, cy, 'SCENARIO', rotation=45, alpha=0.10, size=72, ha='center', va='center', color='red')

    plt.tight_layout()
    #plt.title(ptitle,axes=ax)
    outfile = os.path.join(outfolder, '%s_%s.pdf' % (edict['eventid'], modelname))
    pngfile = os.path.join(outfolder, '%s_%s.png' % (edict['eventid'], modelname))
    plt.savefig(outfile, dpi=300)
    plt.savefig(pngfile)
    print 'Saving map output to %s' % outfile
    print 'Saving map output to %s' % pngfile


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


def saveMap():
    pass


if __name__ == '__main__':
    pass
