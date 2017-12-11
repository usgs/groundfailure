#!/usr/bin/env python

#stdlib imports
import os
import fiona
from shapely.geometry import shape
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from scipy import interpolate
from sklearn.metrics import roc_curve, roc_auc_score, auc
from skimage.feature import match_template
import copy
import collections
import urllib
import json

#local imports
from groundfailure.sample import pointsFromShapes
from mapio.gdal import GDALGrid
from mapio.geodict import GeoDict
from mapio.grid2d import Grid2D
from gfail.stats import computeHagg

# Make fonts readable and recognizable by illustrator
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = ['Arial', 'Bitstream Vera Serif', 'sans-serif']


def concatenateModels(modellist, astitle='id', includeunc=False):
    """
    Put several models together into dictionary in format for modelSummary
    :param astitle: 'id' to use shakemap id or 'model' to use model name
    :param includeunc: include modelmin and modelmax if present, will include in same dictionary as model
    :returns: dictionary containing all models combined into a single dictionary
    """
    newdict = collections.OrderedDict()
    for model in modellist:
        if astitle == 'id':
            title = model['model']['description']['shakemap']
        else:
            title = model['model']['description']['name']
        newdict[title] = model['model']
        if len(model) == 3 and includeunc:
            newdict[title + '_min'] = model['modelmin']
            newdict[title + '_max'] = model['modelmax']
    return newdict


def modelSummary(models, titles=None, eids=None, outputtype='unknown', cumulative=False, histtype='bar', bounds=None, bins=25,
                 semilogy=False, normed=True, thresh=0., showplots=True, csvfile=None, saveplots=False, filepath=None,
                 xlims=[0., 1.], ylims=None, summary_figure=True, individual_plots=True, fileprefix=''):
    """
    Function for creating a summary histogram of a model output - can only do cumulative step plot if combine_events=True

    :param models: model results (ordered dict ) or list of model results (list of ordered dicts) of multiple models
                   (Model results should be in the original output format from the model codes) - only the first key
                    will be used, the rest are ignored
    :param titles: List of titles to use for each model, in same order as model. If none, will use key of each model (may be non-unique)
    :param outputtype: Type of model output, just used for label (e.g.'probability', 'coverage', 'index')
    :param cumulative: True for cumulative histogram, false for non-cumulative
    :param histtype: ‘bar’, ‘barstacked’, ‘step’, ‘stepfilled’ (same as for plt.hist) - only applies to individual plots, all summary plots are step
                    # Add scatter plot?
    :param bounds: Bounding box of area to include in summary. Input as dictionary e.g. {'xmin': -119.2, 'xmax': -118.,
                   'ymin': 34., 'ymax': 34.7}. If None, will use entire area in model
    :param bins: bins to use for histogram. if a single integer, will create that number of bins using min
                 and max of data, if a numpy array, will use those as bin edges
    :param semilogy: uses log scale instead of linear on y axis if True
    :param thresh: threshold for a nonoccurrence, default is zero but for models that never have nonzero values, can set
                   to what you decide is insignificant
    :param csvfile: name of csvfile of summary of outputs, if None, not output
    :param showplots: if True, will display the plots
    :param saveplots: if True, will save the plots
    :param filepath: Filepath for saved plots, if None, will save in current directory. Files are named with test name
                     and time stamp
    :param getquakenames: Get earthquake names from comcat using event id (deleted at least temporarily)
    :param includeunc: Include uncertainty in summary, if files are present?
    :param combine_events: if True, put all events on same plot
    :returns: means, medians, totareas, titles, n, bins
    """
    if showplots is False:
        plt.ioff()

    means = []
    medians = []
    totareas = []
    vallist = []
    means_min = []
    medians_min = []
    totareas_min = []
    vallist_min = []
    means_max = []
    medians_max = []
    totareas_max = []
    vallist_max = []
    if type(models) != list:
        models = [models]
    keylist = [list(mod.keys())[0] for mod in models]  # will only get first one, not min and max
    if titles is None:
        titles = keylist
    else:
        if len(titles) != len(models):
            raise Exception('Length of titles provided are not equal to length of models provided')
    if eids is None:
        eids = keylist
    else:
        if len(eids) != len(models):
            raise Exception('Length of eids provided not equal to length of models provided')

    # for i in np.arange(len(titles)):
    #     names = []
    #     times = []
    #     magnitudes = []
    #     if getquakenames:
    #         #try:
    #         id1 = titles[i].split('_')[0]
    #         name, time, magnitude = getQuakeInfo(id1)
    #         names.append(name)
    #         times.append(time)
    #         magnitudes.append(magnitude)
    #         # except Exception as e:
    #         #     print(e)
    #         #     print('setting quake info to unknown')
    #         #     names.append('unknown')
    #         #     times.append('unknown')
    #         #     magnitudes.append('unknown')
    #     else:
    #         names.append('unknown')
    #         times.append('unknown')
    #         magnitudes.append('unknown')

    for k, mod in enumerate(models):
        # Trim model, if needed
        if bounds is not None and len(bounds) == 4:
            model1 = copy.deepcopy(mod)  # to avoid changing original file
            model1 = model1.cut(bounds['xmin'], bounds['xmax'], bounds['ymin'], bounds['ymax'], align=True)
        else:
            model1 = mod
        grid = model1[keylist[k]]['grid'].getData()
        allvals = grid[~np.isnan(grid)]
        # compute area
        totareas.append(computeHagg(model1[keylist[k]]['grid'], probthresh=thresh))
        total = len(allvals)
        totalnonz = len(allvals[allvals > float(thresh)])
        allvals = allvals[allvals > float(thresh)]
        means.append(np.mean(allvals))
        medians.append(np.median(allvals))
        vallist.append(allvals)
        try:
            # first - 1std
            gridmin = model1[keylist[k] + '_min']['grid'].getData()
            allvalsmin = gridmin[~np.isnan(gridmin)]
            # compute area
            totareas_min.append(computeHagg(model1[keylist[k] + '_min']['grid'], probthresh=thresh))
            #totalmin = len(allvalsmin)
            #totalnonzmin = len(allvalsmin[allvalsmin > float(thresh)])
            allvalsmin = allvalsmin[allvalsmin > float(thresh)]
            means_min.append(np.mean(allvalsmin))
            medians_min.append(np.median(allvalsmin))
            vallist_min.append(allvalsmin)
            # and +1 std
            gridmax = model1[keylist[k] + '_max']['grid'].getData()
            allvalsmax = gridmax[~np.isnan(gridmax)]
            # compute area
            totareas_max.append(computeHagg(model1[keylist[k] + '_max']['grid'], probthresh=thresh))
            #totalmax = len(allvalsmax)
            #totalnonzmax = len(allvalsmax[allvalsmax > float(thresh)])
            allvalsmax = allvalsmax[allvalsmax > float(thresh)]
            means_max.append(np.mean(allvalsmax))
            medians_max.append(np.median(allvalsmax))
            vallist_max.append(allvalsmax)
        except Exception:  # as e:
            #print(e)
            print('Unable to include uncertainty for %s' % keylist[k])
            vallist_max.append(0)  # (np.nan*(np.zeros(np.shape(allvals))))
            vallist_min.append(0)  # (np.nan*(np.zeros(np.shape(allvals))))
            totareas_min.append(float('nan'))
            means_min.append(float('nan'))
            medians_min.append(float('nan'))
            # and +1 std
            totareas_max.append(float('nan'))
            means_max.append(float('nan'))
            medians_max.append(float('nan'))

        # for n, m, t, ti, tot in zip(names, magnitudes, times, titles, totareas):
        #     if 'unknown' in n:
        #         labels = ['%s - %1.1e km2' % (t, m) for t, m in zip(titles, totareas)]
        #     else:
        #         labels = ['M%s %s - %s - %1.1e km2\nid: %s' % (m, n, t, tot, ti)]
        labels = ['%s - %1.1e km2' % (t, m) for t, m in zip(titles, totareas)]

    if summary_figure:  # make summary figure
        fig = plt.figure(figsize=(16, 10))
        ax = fig.add_subplot(111)
        #labels1 = [lab.split('\n')[0] for lab in labels]
        n, bins, rects = ax.hist(tuple(vallist), bins=bins, normed=normed, cumulative=cumulative, histtype='step', range=(0., 1.),
                                 label=labels)
        ## Shows uncertainty, but made plot too busy
        # figjunk = plt.figure(frameon=False)
        # axjunk = figjunk.add_subplot(111)
        # if type(n) != list:
        #     n = [n]
        # for vmin, vmax, nt, r in zip(vallist_max, vallist_min, n, rects[::-1]):  # for some reason, need to flip rectangle order to get colors right
        #     if type(vmin) is int:  # skip if there are not uncertainties
        #         continue
        #     ymin, bins, rects2 = axjunk.hist(vmin, bins=bins, normed=normed, cumulative=cumulative, histtype='step')
        #     ymax, bins, rects2 = axjunk.hist(vmax, bins=bins, normed=normed, cumulative=cumulative, histtype='step')
        #     mid = 0.5*(bins[1:] + bins[:-1])
        #     try:  # draw transparent boxes the same color around line instead of error bars
        #         color = r.get_facecolor()
        #     except:
        #         color = r[0].get_facecolor()
        #     x = np.concatenate(([0.], mid, mid[::-1], [0.]))
        #     y = np.concatenate(([0.], ymin, ymax[::-1], [0.]))
        #     #import pdb; pdb.set_trace()
        #     poly = Polygon(np.vstack([x, y]).T, facecolor=color, edgecolor='none', alpha=0.2)
        #     ax.add_patch(poly)

        #plt.close(figjunk)
        #plt.ion()
        if cumulative and histtype == 'step':  # remove ugly vertical line at end
            for r in rects:
                try:
                    r.set_xy(r.get_xy()[:-1])
                except:
                    for rect in r:
                        rect.set_xy(rect.get_xy()[:-1])
        if normed:
            ax.set_ylabel('Proportion of cells')
        else:
            ax.set_ylabel('Total cells')
        ax.set_xlabel(outputtype)
        if semilogy:
            ax.set_yscale('log')

        if cumulative:
            cumul = 'Cumulative'
        else:
            cumul = 'Non-cumulative'
        if thresh > 0:
            plt.suptitle('%s Summary Plot\nExcluded that were < threshold of %0.2f' % (cumul, thresh))

        # Put a legend below plot
        if len(models) > 10:
            ax.legend(loc='center right', bbox_to_anchor=(1.3, 0.5), prop={'size': 6})
            plt.subplots_adjust(left=0.05, right=0.8)
        else:
            ax.legend(loc='upper right', prop={'size': 8})

        #plt.tight_layout()
        #plt.subplots_adjust(top=0.95)

        ax.set_xlim(xlims)
        ax.set_ylim(ylims)

        if saveplots is True:
            if filepath is None:
                filepath = os.getcwd()
            import datetime
            time1 = datetime.datetime.utcnow().strftime('%d%b%Y_%H%M')
            fig.savefig(os.path.join(filepath, '%sSummary_%s.png' % (fileprefix, time1,)))

    if individual_plots:  # Make individual summary plots
        for val, vmin, vmax, label, id1 in zip(vallist, vallist_max, vallist_min, labels, titles):
            if len(val) == 0:
                print('All values zero for %s, going to next event' % (label))
                continue
            fig = plt.figure()
            ax = fig.add_subplot(111)
            n, bins, rects = ax.hist(val, bins=bins, range=(0., 1.), normed=normed, cumulative=cumulative, histtype=histtype, label=label)
            if type(vmin) is not int:  # skip this if no uncertainties
                ymin, bins, rects2 = plt.hist(vmin, bins=bins, normed=normed, cumulative=cumulative, histtype=histtype)
                ymax, bins, rects2 = plt.hist(vmax, bins=bins, normed=normed, cumulative=cumulative, histtype=histtype)
                mid = 0.5*(bins[1:] + bins[:-1])
                if histtype == 'step':  # draw transparent boxes the same color around line instead of error bars
                    try:
                        color = r.get_facecolor()
                    except:
                        color = r[0].get_facecolor()
                    x = np.concatenate(([0.], mid, mid[::-1], [0.]))
                    y = np.concatenate(([0.], ymin, ymax[::-1], [0.]))
                    poly = Polygon(np.vstack([x, y]).T, facecolor=color, edgecolor='none', alpha=0.2)
                    ax.add_patch(poly)
                else:
                    yerr = np.vstack((np.abs(n-ymin), np.abs(ymax - n)))
                    ax.errorbar(mid, n, yerr=yerr, fmt='none', ecolor='k')

            if cumulative and histtype == 'step':  # remove ugly vertical line at end
                for r in rects:
                    try:
                        r.set_xy(r.get_xy()[:-1])
                    except:
                        for rect in r:
                            rect.set_xy(rect.get_xy()[:-1])
            if normed:
                ax.set_ylabel('Proportion of cells')
            else:
                ax.set_ylabel('Total cells')
            ax.set_xlabel(outputtype)
            ax.set_xlim(xlims)
            ax.set_ylim(ylims)
            if semilogy:
                ax.set_yscale('log')

            if cumulative:
                cumul = 'Cumulative'
            else:
                cumul = 'Non-cumulative'
            if thresh > 0:
                plt.suptitle('%s\n %s Summary Plot - Excluded %d out of %d cells (%0.1f%%) that were < threshold of %0.2f' % (label, cumul,
                             totalnonz, total, totalnonz/float(total) * 100., thresh))
            else:
                plt.suptitle('%s\n %s Summary Plot' % (label, cumul))
            #plt.tight_layout()
            #plt.subplots_adjust(top=0.95)
            plt.subplots_adjust(bottom=0.25)

            if saveplots is True:
                if filepath is None:
                    filepath = os.getcwd()
                import datetime
                time1 = datetime.datetime.utcnow().strftime('%d%b%Y_%H%M')
                fig.savefig(os.path.join(filepath, '%sSummary_%s_%s.png' % (fileprefix, id1, time1)))

    if showplots is True:
        plt.show()
    else:
        plt.close('all')
    plt.ion()

    if csvfile is not None:
        #binrange = (bins[:-1] + bins[1:])/2.
        # if combine_events:
        #     binsS = ['%0.2f - %0.2f' % (b0, b1) for b0, b1 in zip(bins[:-1], bins[1:])]
        # else:
        #     binsS = []
        import csv
        with open(csvfile, 'w') as csvfile1:
            writer = csv.writer(csvfile1)
            writer.writerow(['Id', 'Description', 'Mean', 'Meanmin', 'Meanmax', 'Median', 'Medianmin',
                             'Medianmax', 'Area affected', 'Area affected_min', 'Area affected_max',
                             'Threshold'])  # + binsS)
            for i, ti in enumerate(titles):
                # if combine_events:
                #     nvals = ['%0.2f' % nval for nval in n[i]]
                # else:
                #     nvals = []
                writer.writerow([eids[i], titles[i], means[i], means_min[i], means_max[i], medians[i],
                                 medians_min[i], medians_max[i], totareas[i], totareas_min[i], totareas_max[i],
                                 thresh])  # + nvals)

    return means, medians, totareas, titles, means_min, means_max, medians_min, medians_max, totareas_min, totareas_max


def getQuakeInfo(id):
    BASEURL = 'http://earthquake.usgs.gov/fdsnws/event/1/query?'
    #here we're using the request submodule to open a url, just like we would use open() to open a file.
    contribs = ['', 'us', 'nc', 'usp', 'atlas', 'ci', 'ak', 'at', 'cgs', 'hv', 'ismp', 'ld', 'mb', 'nc', 'nm', 'nn', 'np', 'pr', 'pt', 'se', 'us', 'uu', 'uw']
    for con in contribs:
        indict = {'format': 'geojson',
                  'eventid': con + id,
                  'producttype': 'origin',
                  'includesuperseded': 'false'}
        #the urllib module (in Python 3) encapsulates everything the standard library knows about urls
        params = urllib.parse.urlencode(indict)
        try:
            #assemble complete url
            url = '%s%s' % (BASEURL, params)
            f = urllib.request.urlopen(url)
            break
        except:
            continue
    #urlopen() returns a file-like object, which means that it behaves just like a file object does, including
    #allowing you to read in all the data from the thing that has been opened.
    #note the decode() method, which is a new necessity in Python 3, in order to convert a string of bytes
    #into ASCII (utf-8).
    data = f.read().decode('utf-8')
    #Always close your file!
    f.close()
    jsondict = json.loads(data)
    title = jsondict['properties']['title']
    time = jsondict['properties']['products']['origin'][0]['properties']['eventtime']
    magnitude = jsondict['properties']['mag']
    return title, time, magnitude


def convert2Coverage(gdict, inventory, numdiv=30., method='nearest', proj='moll'):
    """Fast method to produce grid of area actually affected by landsliding in each cell defined by geodict

    :param gdict: geodict, likely taken from model to compare inventory against
    :param inventory: full file path to shapefile of inventory, must be in geographic coordinates, WGS84
    :type inventory: string
    :param numdiv: Approximate amount to subdivide each cell of geodict by to compute areas (higher number slower but more accurate)
    :return inventorygrid: Grid2D object reporting proportional area of landsliding inside each cell defined by geodict
    :param method: method for resampling when projecting back to geographic coordinates, nearest recommended but not perfect. Cubic not recommended.
    :param proj: projection to use to obtain equal area,  'moll'  mollweide, or 'laea' lambert equal area
    :returns: Grid2D object reporting approximate areal coverage of input inventory corresponding to geodict
    """

    lat0 = np.mean((gdict.ymin, gdict.ymax))
    lon0 = np.mean((gdict.xmin, gdict.xmax))
    gdsubdiv = {'xmin': gdict.xmin - 3*gdict.dx, 'xmax': gdict.xmax + 3*gdict.dx, 'ymin': gdict.ymin - 3*gdict.dy, 'ymax': gdict.ymax + 3*gdict.dy, 'dx': gdict.dx/numdiv,
                'dy': gdict.dy/numdiv, 'ny': gdict.ny*numdiv, 'nx': gdict.nx*numdiv}
    subgd = GeoDict(gdsubdiv, adjust='res')
    #subgd = Grid2D.bufferBounds(gdict, subgd, doPadding=True, buffer_pixels=2)

    with fiona.open(inventory) as f:
        invshp = list(f.items())
    shapes = [shape(inv[1]['geometry']) for inv in invshp]

    # Rasterize with oversampled area
    rast = Grid2D.rasterizeFromGeometry(shapes, subgd, fillValue=0., burnValue=1.0, mustContainCenter=True)
    #rast = GDALGrid.copyFromGrid(rast)
    #rast.save('test7.tif')

    # Transform to equal area projection
    #projs = '+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs'
    #projs = '+proj=%s +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs' % (proj)
    #projs = '+proj=%s +lat_0=%f +lon_0=%f +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs' % (proj, lat0, lon0)
    #projs = '+proj=%s +datum=WGS84 +lat_0=%0.5f +lon_0=%0.5F +units=m +x_0=0 +y_0=0' % (proj, lat0, lon0)
    projs = '+proj=%s +lat_0=%f +lon_0=%f +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs' % (proj, lat0, lon0)
    equal_area = rast.project(projection=projs, method='nearest')
    egdict = equal_area.getGeoDict()

    # Use block mean to get %coverage of larger grid
    bm, xdiff, ydiff = block_mean(equal_area.getData(), factor=numdiv, replacenan=0.0)
    #dat = equal_area.getData()
    #dat[np.isnan(dat)] = 0.0
    #bm1 = block_reduce(dat, block_size=(numdiv, numdiv), func=np.mean, cval=0.0)

    # complicated adjustment to get the alignment right
    gdds = {'xmin': egdict.xmin + 0.5 * egdict.dx * (numdiv - 1), 'xmax': egdict.xmax - egdict.dx*(xdiff-0.5+numdiv) + (numdiv * egdict.dx)/2., 'ymin': egdict.ymin + egdict.dy*(ydiff-0.5+numdiv) - (egdict.dy*numdiv)/2., 'ymax': egdict.ymax + 0.5*egdict.dy - 0.5*numdiv*egdict.dy/2, 'dx': egdict.dx*numdiv,
            'dy': egdict.dy*numdiv, 'ny': np.floor(egdict.ny/numdiv), 'nx': np.floor(egdict.nx/numdiv), 'projection': projs}
    dsgd = GeoDict(gdds, adjust='res')
    dsgd.setProjection(projs)
    eabig = Grid2D(data=bm, geodict=dsgd)
    eabig = GDALGrid.copyFromGrid(eabig)

    #TODO add block-mean to mapio as an interpolation method so it can be called like this
    #eabig = equal_area.interpolateToGrid(dsgd, method='block_mean')

    # Project back
    eabigproj = eabig.project(projection=gdict.projection)

    # Resample to original grid
    inventorygrid = eabigproj.interpolateToGrid(gdict, method='linear')

    return inventorygrid


def convert2Prob(gdict, inventory, mustContainCenter=False):
    """Convert inventory shapefile to binary grid (with geodict of gdict) with 1 if any part of landslide was in a cell, 0 if not

    :param gdict: geodict, likely taken from model to compare inventory against
    :param inventory: full file path to shapefile of inventory, must be in geographic coordinates, WGS84
    :type inventory: string
    :returns rast: Grid2D object containing rasterized version of inventory where cell is 1. if any part of a landslide was in the cell
    """
    with fiona.open(inventory) as f:
        invshp = list(f.items())

    shapes = [shape(inv[1]['geometry']) for inv in invshp]

    # Rasterize with allTouch
    rast = Grid2D.rasterizeFromGeometry(shapes, gdict, fillValue=0., burnValue=1.0, mustContainCenter=mustContainCenter)

    return rast


def statsCoverage(modelgrid, inventorygrid, bins=None, showplots=True, saveplots=False, filepath=None):
    """TO DO - FIND MORE TESTS THAT ARE EASY TO COMPARE WITH EACH OTHER
    Compute stats and make comparison plots specific to models that output areal coverage like Godt et al 2008

    :param modelgrid: Grid2D object of model results
    :param inventorygrid: Grid2D object of areal coverage of inventory computed on same grid as modelgrid using, for example, convert2Coverage
    :param bins: bin edges to use for various binning and threshold statistical calculations. if None bins = np.linspace(0, np.max((inv.max(), model.max())), 11)
    :param showplots: if True, will display the plots
    :param saveplots: if True, will save the plots
    :param filepath: Filepath for saved plots, if None, will save in current directory. Files are named with test name and time stamp
    :returns:
        * results: dictionary of results of tests.
                      {'Compare coverage': dictionary,
                        'RMS': float,
                        'RMS_nonzero': float,
                        'Percent in bin': dictionary,
                        'Direct Comparison': dictionary]}
        * invminusmod: Grid2D object of difference between inventory and model (inventory - model)

    """

    inv = inventorygrid.getData()
    model = modelgrid.getData()
    invminusmod = GDALGrid(inv-model, inventorygrid.getGeoDict())

    plt.ioff()

    # Make statistical comparisons
    results = {}
    if bins is None:
        bins = np.linspace(0, np.max((inv.max(), model.max())), 11)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist((model, inv), bins=bins)
    ax.set_yscale('log')
    ax.legend(('Model', 'Inventory'))
    ax.set_ylabel('Total # of cells')
    ax.set_xlabel('Predicted coverage bin')
    results['Direct Comparison'] = {'bins': bins, 'model': model, 'inventory': inv}

    # Statistical comparison
    perc = []
    areatot = []
    for i, bin in enumerate(bins[:-1]):
        idx = (model > bin) & (model <= bins[i+1])
        areain = inv[idx]
        totalin = sum((areain > bin) & (areain <= bins[i+1]))
        perc.append(float(totalin)/len(idx)*100)
        areatot.append(np.mean(areain))
    areatot = np.nan_to_num(areatot)

    # get centerpoint of bins
    binvec = []
    for i in range(len(bins[:-1])):
        binvec.append(bins[i]+(bins[i+1]-bins[i])/2)
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(binvec, perc, 'o-')
    ax1.set_ylabel('% in right bin')
    ax1.set_xlabel('Predicted coverage bin')
    ax1.set_title('Actual coverage in predicted bin')
    results['Percent in bin'] = {'bincenters': binvec, 'perc': perc}

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(binvec, areatot, 'o-')
    ax2.plot([0., 1.], [0., 1.], '--', color='gray')
    ax2.set_ylabel('Actual coverage')
    ax2.set_xlabel('Predicted coverage')
    ax2.set_title('Actual vs. Predicted coverage')
    results['Compare coverage'] = {'bincenters': binvec, 'areatot': areatot}

    if showplots is True:
        plt.show()
    if saveplots is True:
        if filepath is None:
            filepath = os.getcwd()
        import datetime
        time1 = datetime.datetime.utcnow().strftime('%d%b%Y_%H%M')
        fig.savefig(os.path.join(filepath, 'Direct_compare_%s.pdf' % (time1,)))
        fig1.savefig(os.path.join(filepath, 'Percent_in_bin_%s.pdf' % (time1,)))
        fig2.savefig(os.path.join(filepath, 'Compare_coverage_%s.pdf' % (time1,)))

    # RMS
    idx = np.union1d(model.nonzero()[0], inv.nonzero()[0])
    results['RMS'] = np.sqrt((model - inv)**2).mean()
    print(('RMS: %0.3f' % results['RMS']))
    results['RMS_nonzero'] = np.sqrt((model[idx] - inv[idx])**2).mean()
    print(('RMS_nonzero: %0.3f' % results['RMS_nonzero']))

    return invminusmod, results


def stats(modelgrid, inventory, dx=100., Nsamp=None, method='nearest', extent='inventory', bins=None, runtests=True, showplots=True, saveplots=False, filepath=None):
    """Run through suite of tests for models that output probability or index that varies between 0 and 1

    :param modelgrid: Grid2D object of model results
    :param inventory: full file path to shapefile of inventory, must be in geographic coordinates, WGS84
    :type inventory: string
    :param dx: Approximate sample spacing in meters, overwritten if Nsamp is defined
    :type dx: float
    :param Nsamp: Total number of samples desired - will choose optimal dx to get slightly more than this number of samples delete samples outside bounds and randomly delete others until sample number is exactly Nsamp
    :type Nsamp: integer
    :param method: method used for interp2d when transforming sampled model values back from projected coordinates to geographic coordinates - 'nearest', 'linear', or 'cubic'
    :type method: string
    :param extent: extent to include in sampling - 'inventory' or 'model' or custom bounds as tuple (xmin, ymin, xmax, ymax) - in lats and lons
    :param bins: bin edges to use for various binning and threshold statistical calculations. if None bins = [0, 0.2, 0.4, 0.6, 0.8, 1.]
    :param runtests: if True, will run various statistical tests, if False will just output sampled values
    :param showplots: if True, will disply the plots
    :param saveplots: if True, will save the plots
    :param filepath: Filepath for saved plots, if None, will save in current directory. Files are named with test name and time stamp
    :returns:
        * yespoints: Nx2 array of geographic coordinates of positive sample locations
        * nopoints: Nx2 array of geographic coordinates of negative sample locations
        * modelvalyes: N model output values corresponding to yespoints
        * modelvalno: N model output values corresponding to nopoints
        * results: dictionary of results of statistical tests. Will be empty if runtests=False
        {'Occ_nonocc': dict,
        'SRC': dict,
        'ROC': dict,
        'AUC_ROC': float,
        'Log_loss': float,
        'GFC': dict,
        'Pred_vs_Obs': dict,
        'Brier': float,
        'Brier_no': float,
        'Brier_yes': float}
    """
    plt.close('all')
    f = fiona.collection(inventory, 'r')
    shapes = list(f)
    bxmin, bymin, bxmax, bymax = f.bounds
    gdict = modelgrid.getGeoDict()

    if extent == 'model':
        extent = gdict.xmin, gdict.ymin, gdict.xmax, gdict.ymax
    elif extent == 'inventory':
        extent = bxmin, bymin, bxmax, bymax
    #yespoints, junk, nopoints, junk2, xvar, yvar, pshapes, proj = sampleFromShapes(shapes, extent, dx=dx, Nsamp=Nsamp, testPercent=100.)

    yespoints, nopoints, xvar, yvar, pshapes, proj = pointsFromShapes(shapes, extent, dx=dx, Nsamp=Nsamp)
    yesptx = [pt[0] for pt in yespoints]
    yespty = [pt[1] for pt in yespoints]
    noptx = [pt[0] for pt in nopoints]
    nopty = [pt[1] for pt in nopoints]
    #import pdb; pdb.set_trace()
    # Get values of model at those points
    lons = np.linspace(gdict.xmin, gdict.xmax, gdict.nx)
    lats = np.linspace(gdict.ymax, gdict.ymin, gdict.ny)
    if method.lower() == 'nearest':
        modelvalyes = []
        modelvalno = []
        for XX, YY in zip(yesptx, yespty):
            row = (np.abs(lats - YY)).argmin()
            col = (np.abs(lons - XX)).argmin()
            modelvalyes.append(modelgrid.getData()[row, col])
        for XX, YY in zip(noptx, nopty):
            row = (np.abs(lats - YY)).argmin()
            col = (np.abs(lons - XX)).argmin()
            modelvalno.append(modelgrid.getData()[row, col])
    else:
        func = interpolate.interp2d(lons, lats, modelgrid.getData(), kind=method.lower())
        modelvalyes = np.array([float(func(XX, YY)) for XX, YY in zip(yesptx, yespty)])
        modelvalno = np.array([float(func(XX, YY)) for XX, YY in zip(noptx, nopty)])

    modelvalyes = np.nan_to_num(np.array(modelvalyes))  # replace nan with zeros
    modelvalno = np.nan_to_num(np.array(modelvalno))  # replace nan with zeros

    # Now run the desired tests and make the desired plots
    results = {}

    if runtests is True:
        # Brier score
        N = len(yespoints) + len(nopoints)
        yessum = np.sum([(val-1)**2 for val in modelvalyes])
        nosum = np.sum([(val)**2 for val in modelvalno])
        results['Brier_yes'] = yessum/len(modelvalyes)
        results['Brier_no'] = nosum/len(modelvalno)
        results['Brier'] = (yessum + nosum)/N
        print(('Brier scores: overall %0.3f\nBrier_yes score: %0.3f\nBrier_no score %0.3f' % (results['Brier'], results['Brier_yes'], results['Brier_no'])))

        # Logarithmic score
        tempno = np.array(modelvalno).copy()
        tempyes = np.array(modelvalyes).copy()
        tempno[tempno == 0] = 1.e-15
        tempyes[tempyes == 0] = 1.e-15
        results['Log_loss'] = -(np.sum(np.log(tempyes)) + np.sum(np.log(1.-tempno)))/N
        print(('Log loss score: %0.3f' % (results['Log_loss'],)))

        if bins is None:
            bins = [0, 0.2, 0.4, 0.6, 0.8, 1.]
        binvec = []
        observed = []
        percyes = []
        percno = []
        overall_tot = len(modelvalyes) + len(modelvalno)
        for i in range(len(bins[:-1])):
            binvec.append(bins[i]+(bins[i+1]-bins[i])/2)
            yestot = np.sum([(modelvalyes > bins[i]) & (modelvalyes < bins[i+1])])
            notot = np.sum([(modelvalno > bins[i]) & (modelvalno < bins[i+1])])
            if notot+yestot != 0:
                observed.append(float(yestot)/(yestot+notot))
            else:
                observed.append('nan')
            percyes.append((yestot/float(overall_tot))*100.)
            percno.append((notot/float(overall_tot))*100.)

        plt.ioff()

        # Predicted vs. Observed ratios
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(binvec, observed, '-o')
        ax.plot([0]+binvec, [0]+binvec, '--', color='gray')
        ax.set_xlabel('Expected ratio')
        ax.set_ylabel('Observed ratio')
        ax.set_xlim([bins[0], bins[-1]])
        ax.set_title('Predicted vs. Observed')
        results['Pred_vs_Obs'] = {'binvec': binvec, 'observed': observed}

        # Ground failure occurrence/nonoccurrence
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        wid = (bins[1]-bins[0])/2.5
        rects1 = ax1.bar(np.array(bins[:-1]), percyes, width=wid)
        rects2 = ax1.bar(np.array(bins[:-1])+wid, percno, width=wid, color='r')
        ax1.set_xlabel('Predicted susceptibility range')
        ax1.set_ylabel('% of samples')
        ax1.legend((rects1[0], rects2[0]), ('Occurrence', 'Nonoccurrence'))
        ax1.set_title('Occurrence vs. Nonoccurrence')
        results['Occ_nonocc'] = {'bins': bins, 'percyes': percyes, 'percno': percno}

        # Ground failure capture for various thresholds
        gfc = []
        for val in bins:
            gfc.append(np.sum([modelvalyes > val])/float(len(yespoints)))
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        ax2.plot(bins, gfc, 'o-')
        ax2.set_xlabel('Threshold')
        ax2.set_ylabel(r'%GFC')
        ax2.set_title('Ground Failure Capture')
        results['GFC'] = {'thresholds': bins, 'gfc': gfc}

        # ROC curves
        fpr, tpr, thresholds = roc_curve(np.concatenate((np.ones(len(yespoints)), np.zeros(len(nopoints)))), np.concatenate((modelvalyes, modelvalno)))
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        ax3.plot(fpr, tpr)
        ax3.set_xlabel('False positive rate')
        ax3.set_ylabel('True positive rate')
        ax3.set_xlim([0, 1.])
        ax3.set_ylim([0, 1.])
        ax3.plot(fpr, fpr, '--', color='gray')
        ax3.set_title('ROC curve')
        results['ROC'] = {'thresholds': bins, 'gfc': gfc}

        results['AUC_ROC'] = roc_auc_score(np.concatenate((np.ones(len(yespoints)), np.zeros(len(nopoints)))), np.concatenate((modelvalyes, modelvalno)))
        print(('AUC_ROC: %0.3f' % (results['AUC_ROC'],)))
        ax3.text(0.8, 0.2, 'AUC: %0.3f' % results['AUC_ROC'])

        # Success rate curves
        sucbin = np.linspace(0, 1., 100)
        prop = []
        realvals = np.concatenate((np.ones(len(yespoints)), np.zeros(len(nopoints))))
        predvals = np.concatenate((modelvalyes, modelvalno))
        indx = np.argsort(predvals)
        predvals = predvals[indx]
        realvals = realvals[indx]
        for val in sucbin:
            prop.append(np.sum(realvals[predvals < val])/len(yespoints))
        fig4 = plt.figure()
        ax4 = fig4.add_subplot(111)
        ax4.plot(sucbin, prop)
        ax4.set_xlabel('Success Rate Curve')
        ax4.set_ylabel('Proportion of actual occurrences')
        ax4.set_title('Proportion of Study Area')
        AUC = auc(sucbin, prop)
        print(('AUC_SRC: %0.3f' % AUC))
        ax4.text(0.8, 0.2, 'AUC: %0.3f' % AUC)
        ax4.set_xlim([0, 1.])
        ax4.set_ylim([0, 1.])
        results['SRC'] = {'xvals': sucbin, 'proportion': prop, 'auc': AUC}

        if showplots is True:
            plt.show()
        if saveplots is True:
            if filepath is None:
                filepath = os.getcwd()
            import datetime
            time1 = datetime.datetime.utcnow().strftime('%d%b%Y_%H%M')
            fig.savefig(os.path.join(filepath, 'Pred_vs_obs_%s.pdf' % (time1,)))
            fig1.savefig(os.path.join(filepath, 'Occ_nonocc_%s.pdf' % (time1,)))
            fig2.savefig(os.path.join(filepath, 'GFC_%s.pdf' % (time1,)))
            fig3.savefig(os.path.join(filepath, 'ROC_%s.pdf' % (time1,)))
            fig4.savefig(os.path.join(filepath, 'SRC_%s.pdf' % (time1,)))

    return yespoints, nopoints, modelvalyes, modelvalno, results


def normXcorr(model, inventory):
    """Perform normalized cross correlation (no shifts) to assess how well the spatial extent was predicted, regardless of overall values

    :param model: Grid2D file of model
    :param inventory: Grid 2D file of inventory processed to simulate what model is supposed to predict (using convert2Prob or convert2Coverage)
    :returns: normalized cross correlation coefficient (between 0 and 1)
    """
    if model.getGeoDict() != inventory.getGeoDict():
        raise Exception('model and inventory files are not identical')

    modat = model.getData()
    invdat = inventory.getData()
    modat[np.isnan(modat)] = 0.
    invdat[np.isnan(invdat)] = 0.
    result = match_template(modat, invdat)
    xcorrcoeff = result[0, 0]
    return xcorrcoeff


def block_mean(ar, factor, replacenan=None):
    """
    Block mean for downsampling 2d array
    From here http://stackoverflow.com/questions/18666014/downsample-array-in-python

    :param ar: 2D array to downsample using a block mean
    :param factor: factor by which to downsample
    :returns: res: downsampled array
    xdiff: number of cells added to right edge
    ydiff: number of cells added to bottom edge
    """
    from scipy import ndimage
    if replacenan is not None:
        ar[np.isnan(ar)] = replacenan
    sx, sy = ar.shape
    # get new dimensions (cutting off extra on edges)
    newx = int(np.floor(sx/float(factor)))
    newy = int(np.floor(sy/float(factor)))
    # build regions over which to average
    vec = np.reshape(np.arange(newx*newy), (newx, newy))
    regions = vec.repeat(factor, 0).repeat(factor, 1)
    # Patch on edges (just expand existing cells on the edges a little) - to account for rounding down earlier
    xdiff = ar.shape[0] - regions.shape[0]
    if xdiff != 0:
        row = regions[-1, :]
        regions = np.row_stack((regions, np.repeat(np.reshape(row, (1, len(row))), xdiff, axis=0)))  # add onto the bottom
    ydiff = ar.shape[1] - regions.shape[1]
    if ydiff != 0:
        col = regions[:, -1]
        regions = np.column_stack((regions, np.repeat(np.reshape(col, (len(col), 1)), ydiff, axis=1)))  # add onto the right side
    res = ndimage.mean(ar, labels=regions, index=np.arange(regions.max() + 1))
    res = res.reshape(newx, newy)
    return res, xdiff, ydiff

#
#
#def line_mean(ar, fact):
#    """
#    Same as block_mean but for 1d array
#    """
#    from scipy import ndimage
#    sx = len(ar)
#    newx = int(np.floor(sx/float(fact)))
#    vec = np.arange(newx)
#    regions = vec.repeat(fact)
#    # Patch on edges (just expand existing cells on the edges a little)
#    xdiff = ar.shape[0] - regions.shape[0]
#    if xdiff != 0:
#        row = regions[-1]
#        regions = np.concatenate((regions, np.repeat(row, xdiff)))
#    res = ndimage.mean(ar, labels=regions, index=np.arange(regions.max() + 1))
#    return res
#

# def computeCoverage_accurate(gdict, inventory, numdiv=10.):
#     """
#     VERY SLOW!!
#     Slow but more accurate method to produce grid of area actually affected by landsliding in each cell defined by geodict

#     :param gdict: geodict, likely taken from model to compare inventory against
#     :param inventory: full file path to shapefile of inventory, must be in geographic coordinates, WGS84
#     :type inventory: string
#     :param numdiv: Approximate amount to subdivide each cell of geodict by to compute areas (higher number slower but more accurate)

#     :returns: Grid2D object reporting areal coverage of landsliding inside each cell defined by geodict
#     """

#     f = fiona.collection(inventory, 'r')
#     shapes = list(f)
#     bxmin, bymin, bxmax, bymax = f.bounds

#     lons = np.linspace(gdict.xmin, gdict.xmax, gdict.nx)
#     lats = np.linspace(gdict.ymax, gdict.ymin, gdict.ny)
#     llons, llats = np.meshgrid(lons, lats)

#     spacing = np.round(np.abs(((lats[1]-lats[0])*111.12*1000.)/numdiv))  # in meters
#     yespoints, nopoints, xvar, yvar, pshapes, proj = pointsFromShapes(shapes, bounds=(gdict.xmin, gdict.ymin, gdict.xmax, gdict.ymax), dx=spacing)

#     # Loop over lat lon pairs that are within boundaries of yes and no points
#     ptlonmax = (np.max((yespoints[:, 0].max(), nopoints[:, 0].max())))
#     ptlonmin = (np.max((yespoints[:, 0].min(), nopoints[:, 0].min())))
#     ptlatmax = (np.max((yespoints[:, 1].max(), nopoints[:, 1].max())))
#     ptlatmin = (np.max((yespoints[:, 1].min(), nopoints[:, 1].min())))

#     subllons = llons[(llons >= ptlonmin) & (llons <= ptlonmax) & (llats >= ptlatmin) & (llats <= ptlatmax)]
#     subllats = llats[(llons >= ptlonmin) & (llons <= ptlonmax) & (llats >= ptlatmin) & (llats <= ptlatmax)]

#     import time
#     # Contains points method
#     t1 = time.clock()
#     dx = gdict.dx
#     area = np.zeros(np.shape(llons))
#     numpts = area.copy()
#     numyes = area.copy()
#     for lat1, lon1 in zip(subllats, subllons):
#         # Find ratio of yes points to no points
#         bbPath = mplPath.Path(np.array([[lon1-0.5*dx, lat1-0.5*dx], [lon1-0.5*dx, lat1+0.5*dx], [lon1+0.5*dx, lat1+0.5*dx], [lon1+0.5*dx, lat1-0.5*dx]]))
#         yesin = sum(bbPath.contains_points(yespoints))  # sum([(yes0 > lon1-0.5*dx) & (yes0 <= lon1+0.5*dx) & (yes1 > lat1-0.5*dx) & (yes1 <= lat1+0.5*dx)])
#         noin = sum(bbPath.contains_points(nopoints))  # sum([(no0 > lon1-0.5*dx) & (no0 <= lon1+0.5*dx) & (no1 > lat1-0.5*dx) & (no1 <= lat1+0.5*dx)])
#         total = yesin + noin
#         if total == 0.:
#             continue
#         # get indices
#         row = np.where(lats == lat1)
#         col = np.where(lons == lon1)
#         # Store total number of points in matrix
#         numpts[row, col] = total
#         # Store area
#         numyes[row, col] = yesin
#     t2 = time.clock()
#     print(('Time elapsed %0.2f seconds' % (t2-t1)))

#     # Correct for incompletely sampled squared (all unsampled points would be no points)
#     numpts[numpts < (numpts[numpts != 0].mean() - numpts[numpts != 0].std())] = np.median(numpts[numpts != 0])  # Will change zeros to nonzeros, but yeses will be 0 in those cells so it doesn't matter
#     area = numyes/numpts

#     inventorygrid = GDALGrid(area, gdict)

#     return inventorygrid

