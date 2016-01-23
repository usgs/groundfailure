#!/usr/bin/env python

#stdlib imports
import os.path
import math
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


def modelMap(grids, edict, config, modelname=None, maproads=True, mapcities=True, isScenario=False, lowthreshes=None, plotorder=None, colormaps=None, boundaries=None, zthresh=0, scaletype='continuous', lims=None):
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
            defaultcolormap = config1['defaultcolormap']
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
                lons = np.arange(xmin, xmax, temp.getGeoDict()['xdim'])
                lons = lons[:temp.getGeoDict()['ncols']]  # make sure right length
                lats = np.arange(ymax, ymin, -temp.getGeoDict()['ydim'])  # backwards so it plots right
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
            if xmin < lonmin-0.1*(lonmax-lonmin):
                boundaries1['xmin'] = lonmin-0.1*(lonmax-lonmin)
            else:
                boundaries1['xmin'] = xmin
            if xmax > lonmax+0.1*(lonmax-lonmin):
                boundaries1['xmax'] = lonmax+0.1*(lonmax-lonmin)
            else:
                boundaries1['xmax'] = xmax
            if ymin < latmin-0.1*(latmax-latmin):
                boundaries1['ymin'] = latmin-0.1*(latmax-latmin)
            else:
                boundaries1['ymin'] = ymin
            if ymax > latmax+0.1*(latmax-latmin):
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

    # Load in topofile and load in or compute hillshade
    if topofile is not None and hillshade is None:
        gdict = GMTGrid.getBoundsWithin(topofile, boundaries)
        # Only load in if no hillshade is available, expand a little to avoid edge cutoff
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

        dat = layergrid.getData().copy()
        if colormaps is not None and len(colormaps) == len(grids):
            palette = colormaps[k]
        else:  # Find preferred default color map for each type of layer
            if 'prob' in layer.lower() or 'pga' in layer.lower() or 'pgv' in layer.lower() or 'cohesion' in layer.lower() or 'friction' in layer.lower():
                palette = cm.jet
            elif 'slope' in layer.lower():
                palette = cm.gnuplot2
            elif 'precip' in layer.lower():
                palette = cm.s3pcpn
            else:
                palette = defaultcolormap
        if topofile is not None:
            # resample topomap to the same as current layer
            topomap2 = GMTGrid.load(topofile, resample=True, method='linear', samplegeodict=layergrid.getGeoDict())
            iwater = np.where(topomap2.getData() == SEA_LEVEL)
            dat[iwater] = 0
        if lowthreshes is not None and len(lowthreshes) == len(grids):
            dat[dat <= lowthreshes[k]] = 0
        else:
            dat[dat <= 0.] = 0
        datm = np.ma.masked_equal(dat, 0)
        palette.set_bad(clear_color, alpha=0.0)

        xmin, xmax, ymin, ymax = layergrid.getBounds()
        lons = np.arange(xmin, xmax+layergrid.getGeoDict()['xdim'], layergrid.getGeoDict()['xdim'])
        lons = lons[:layergrid.getGeoDict()['ncols']]  # make sure it's the right length, sometimes it is one too long, sometimes it isn't because of GMTGrid issue
        lats = np.arange(ymax, ymin-layergrid.getGeoDict()['ydim'], -layergrid.getGeoDict()['ydim'])  # backwards so it plots right side up
        lats = lats[:layergrid.getGeoDict()['nrows']]  # make sure it's the right length, sometimes it is one too long, sometimes it isn't because of GMTGrid issue
        #make meshgrid
        llons1, llats1 = np.meshgrid(lons, lats)
        x1, y1 = m(llons1, llats1)  # get projection coordinates
        #x2, y2 = m(lons, lats)
        if scaletype == 'binned':
            if lims is None or len(lims) != len(grids):
                clev = np.linspace(np.floor(datm.min()), np.ceil(datm.max()), 6)
            else:
                if lims[k] is None:
                    clev = np.linspace(np.floor(datm.min()), np.ceil(datm.max()), 6)
                else:
                    clev = lims[k]
            panelhandle = m.contourf(x1, y1, datm, clev, cmap=palette, lw=0, alpha=ALPHA, rasterized=True)
        else:
            if lims is not None and len(lims) == len(grids):
                vmin = lims[k][0]
                vmax = lims[k][1]
            else:
                vmin = None
                vmax = None
            panelhandle = m.pcolormesh(x1, y1, datm, lw=0, cmap=palette, vmin=vmin, vmax=vmax, alpha=ALPHA, rasterized=True)
        # add colorbar
        cbar = fig.colorbar(panelhandle, fraction=0.046, pad=0.04)
        cbar.set_label(label1, fontsize=10)
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
            elat, elon = edict['lat'], edict['lon']
            ex, ey = m(elon, elat)
            plt.plot(ex, ey, '*', markeredgecolor='k', mfc='None', mew=1.5, ms=12)

        #fill in the lakes and rivers
        m.fillcontinents(color=clear_color, lake_color=watercolor)
        m.drawrivers(color=watercolor)

        #draw country boundaries
        m.drawcountries(color=countrycolor, linewidth=1.0)

        #add map scale
        xmax, xmin, ymax, ymin = boundaries['xmin'], boundaries['xmax'], boundaries['ymin'], boundaries['ymax']
        clat = ymin + (ymax-ymin)/2.0
        clon = xmin + (xmax-xmin)/2.0
        #m.drawmapscale((xmax+xmin)/2., (ymin+(ymax-ymin)/5.), clon, clat, np.round((((xmax-xmin)*110)/5)/10.)*10, barstyle='simple')
        #draw coastlines
        m.drawcoastlines(color='#476C91', linewidth=0.5)

        #draw scenario watermark, if scenario
        if isScenario:
            plt.sca(ax)
            cx, cy = m(clon, clat)
            plt.text(cx, cy, 'SCENARIO', rotation=45, alpha=0.10, size=72, ha='center', va='center', color='red')

    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    #plt.title(ptitle,axes=ax)
    outfile = os.path.join(outfolder, '%s_%s.pdf' % (edict['eventid'], modelname))
    pngfile = os.path.join(outfolder, '%s_%s.png' % (edict['eventid'], modelname))
    plt.savefig(outfile, dpi=300)
    plt.savefig(pngfile)
    print 'Saving map output to %s' % outfile
    print 'Saving map output to %s' % pngfile


def comparisonMap():
    """
    Compare maps from different models
    """
    pass


def compareModelInventory():
    pass


def saveGrid():
    pass


if __name__ == '__main__':
    pass
