#!/usr/bin/env python
# -*- coding: utf-8 -*-

# stdlib imports
import os
import fiona
from shapely.geometry import shape
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from scipy import interpolate
from sklearn.metrics import roc_curve, roc_auc_score, auc
from skimage.feature import match_template
from scipy import ndimage
import copy
import collections
import urllib
import json
import csv
import datetime

# local imports
from gfail.sample import pointsFromShapes
from mapio.gdal import GDALGrid
from mapio.geodict import GeoDict
from mapio.grid2d import Grid2D
from mapio.shake import ShakeGrid

# Make fonts readable and recognizable by illustrator
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = ['Arial',
                                   'Bitstream Vera Serif',
                                   'sans-serif']


def concatenateModels(modellist, astitle='id', includeunc=False):
    """
    Put several models together into dictionary in format for modelSummary.

    Args:
        modellist: List of model results.
        astitle: 'id' to use shakemap id or 'model' to use model name.
        includeunc: include modelmin and modelmax if present, will include in
            same dictionary as model.

    Returns:
        Dictionary containing all models combined into a single dictionary.
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


def modelSummary(models, titles=None, eids=None, outputtype='unknown',
                 cumulative=False, histtype='bar', bounds=None, bins=25,
                 semilogy=False, normed=True, thresh=0., showplots=True,
                 csvfile=None, saveplots=False, filepath=None,
                 xlims=[0., 1.], ylims=None, summary_figure=True,
                 individual_plots=True, fileprefix=''):
    """
    Function for creating a summary histogram of a model output. Can only do
    cumulative step plot if `combine_events=True`.

    Args:
        models (ordered dict): Model results or list of model results (list of
            ordered dicts) of multiple models (Model results should be in the
            original output format from the model codes). Only the first key
            will be used, the rest are ignored.
        titles: List of titles to use for each model, in same order as model.
            If none, will use key of each model (may be non-unique).
        outputtype: Type of model output, just used for label (e.g.,
            'probability', 'coverage', 'index').
        cumulative: True for cumulative histogram, false for non-cumulative.
        histtype: ‘bar’, ‘barstacked’, ‘step’, ‘stepfilled’ (same as for
            plt.hist). Only applies to individual plots, all summary plots are
            step.
        bounds: None, will use entire area, or a dictionary, e.g.

            .. code-block:: python

                {
                    'xmin': -119.2,
                    'xmax': -118.,
                    'ymin': 34.0,
                    'ymax': 34.7
                }

        bins: Bins to use for histogram. if a single integer, will create that
            number of bins using min and max of data, if a numpy array, will
            use those as bin edges.
        semilogy: Uses log scale instead of linear on y axis if True.
        thresh: Threshold for a nonoccurrence, default is zero but for models
            that never have nonzero values, can set to what you decide is
            insignificant.
        csvfile: Name of csvfile of summary of outputs, if None, no output is
            created.
        showplots: Display plots?
        saveplots: Save the plots?
        filepath: Filepath for saved plots, if None, will save in current
            directory. Files are named with test name and time stamp.
        getquakenames: Get earthquake names from comcat using event id
            (deleted at least temporarily).
        includeunc: Include uncertainty in summary, if files are present?
        combine_events: If True, put all events on same plot.

    Returns:
        A tuple with the following elements:
            - means
            - medians
            - totareas
            - titles
            - means_min
            - means_max,
            - medians_min
            - medians_max
            - totareas_min
            - totareas_max
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

    # will only get first one, not min and max
    keylist = [list(mod.keys())[0] for mod in models]
    if titles is None:
        titles = keylist
    else:
        if len(titles) != len(models):
            raise Exception('Length of titles provided are not equal to '
                            'length of models provided')
    if eids is None:
        eids = keylist
    else:
        if len(eids) != len(models):
            raise Exception('Length of eids provided not equal to length of '
                            'models provided')

    for k, mod in enumerate(models):
        # Trim model, if needed
        if bounds is not None and len(bounds) == 4:
            model1 = copy.deepcopy(mod)  # to avoid changing original file
            model1 = model1.cut(bounds['xmin'], bounds['xmax'],
                                bounds['ymin'], bounds['ymax'], align=True)
        else:
            model1 = mod
        grid = model1[keylist[k]]['grid'].getData()
        allvals = grid[~np.isnan(grid)]

        # compute area
        totareas.append(computeHagg(
            model1[keylist[k]]['grid'], probthresh=thresh))
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
            totareas_min.append(computeHagg(
                model1[keylist[k] + '_min']['grid'], probthresh=thresh))
            allvalsmin = allvalsmin[allvalsmin > float(thresh)]
            means_min.append(np.mean(allvalsmin))
            medians_min.append(np.median(allvalsmin))
            vallist_min.append(allvalsmin)

            # and +1 std
            gridmax = model1[keylist[k] + '_max']['grid'].getData()
            allvalsmax = gridmax[~np.isnan(gridmax)]

            # compute area
            totareas_max.append(computeHagg(
                model1[keylist[k] + '_max']['grid'], probthresh=thresh))
            allvalsmax = allvalsmax[allvalsmax > float(thresh)]
            means_max.append(np.mean(allvalsmax))
            medians_max.append(np.median(allvalsmax))
            vallist_max.append(allvalsmax)
        except Exception:  # as e:
            # print(e)
            print('Unable to include uncertainty for %s' % keylist[k])
            vallist_max.append(0)
            vallist_min.append(0)
            totareas_min.append(float('nan'))
            means_min.append(float('nan'))
            medians_min.append(float('nan'))

            # and +1 std
            totareas_max.append(float('nan'))
            means_max.append(float('nan'))
            medians_max.append(float('nan'))

        labels = ['%s - %1.1e km2' % (t, m) for t, m in zip(titles, totareas)]

    if summary_figure:  # make summary figure
        fig = plt.figure(figsize=(16, 10))
        ax = fig.add_subplot(111)
        n, bins, rects = ax.hist(tuple(vallist), bins=bins, normed=normed,
                                 cumulative=cumulative, histtype='step',
                                 range=(0., 1.),
                                 label=labels)

        # remove ugly vertical line at end
        if cumulative and histtype == 'step':
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
            plt.suptitle('%s Summary Plot\nExcluded that were < '
                         'threshold of %0.2f' % (cumul, thresh))

        # Put a legend below plot
        if len(models) > 10:
            ax.legend(loc='center right', bbox_to_anchor=(
                1.3, 0.5), prop={'size': 6})
            plt.subplots_adjust(left=0.05, right=0.8)
        else:
            ax.legend(loc='upper right', prop={'size': 8})

        ax.set_xlim(xlims)
        ax.set_ylim(ylims)

        if saveplots is True:
            if filepath is None:
                filepath = os.getcwd()
            time1 = datetime.datetime.utcnow().strftime('%d%b%Y_%H%M')
            fig.savefig(os.path.join(filepath, '%sSummary_%s.png' %
                                     (fileprefix, time1,)))

    # Make individual summary plots
    if individual_plots:
        # NOTE: it looks like min and max are slipped here...
        # Probably better to just do a loop witn an index rather than zip
        for val, vmin, vmax, label, id1 in zip(vallist, vallist_max, vallist_min, labels, titles):
            if len(val) == 0:
                print('All values zero for %s, going to next event' % (label))
                continue
            fig = plt.figure()
            ax = fig.add_subplot(111)
            n, bins, rects = ax.hist(
                    val, bins=bins, range=(0., 1.), normed=normed,
                    cumulative=cumulative, histtype=histtype, label=label)
            if type(vmin) is not int:  # skip this if no uncertainties
                ymin, bins, rects2 = plt.hist(
                    vmin, bins=bins, normed=normed, cumulative=cumulative,
                    histtype=histtype)
                ymax, bins, rects2 = plt.hist(
                    vmax, bins=bins, normed=normed, cumulative=cumulative,
                    histtype=histtype)
                mid = 0.5 * (bins[1:] + bins[:-1])
                # draw transparent boxes the same color around line instead
                # of error bars
                if histtype == 'step':
                    try:
                        color = r.get_facecolor()
                    except:
                        color = r[0].get_facecolor()
                    x = np.concatenate(([0.], mid, mid[::-1], [0.]))
                    y = np.concatenate(([0.], ymin, ymax[::-1], [0.]))
                    poly = Polygon(np.vstack([x, y]).T, facecolor=color,
                                   edgecolor='none', alpha=0.2)
                    ax.add_patch(poly)
                else:
                    yerr = np.vstack((np.abs(n - ymin), np.abs(ymax - n)))
                    ax.errorbar(mid, n, yerr=yerr, fmt='none', ecolor='k')

            # remove ugly vertical line at end
            if cumulative and histtype == 'step':
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
                plt.suptitle('%s\n %s Summary Plot - Excluded %d out of %d '
                             'cells (%0.1f%%) that were < threshold of %0.2f'
                             % (label, cumul, totalnonz, total,
                                totalnonz / float(total) * 100., thresh))
            else:
                plt.suptitle('%s\n %s Summary Plot' % (label, cumul))
            plt.subplots_adjust(bottom=0.25)

            if saveplots is True:
                if filepath is None:
                    filepath = os.getcwd()
                time1 = datetime.datetime.utcnow().strftime('%d%b%Y_%H%M')
                file_name = '%sSummary_%s_%s.png' % (fileprefix, id1, time1)
                fig.savefig(os.path.join(filepath, file_name))

    if showplots is True:
        plt.show()
    else:
        plt.close('all')
    plt.ion()

    if csvfile is not None:
        with open(csvfile, 'w') as csvfile1:
            writer = csv.writer(csvfile1)
            writer.writerow(['Id', 'Description', 'Mean', 'Meanmin', 'Meanmax',
                             'Median', 'Medianmin', 'Medianmax',
                             'Area affected', 'Area affected_min',
                             'Area affected_max', 'Threshold'])
            for i, ti in enumerate(titles):
                writer.writerow([eids[i], titles[i], means[i], means_min[i],
                                 means_max[i], medians[i], medians_min[i],
                                 medians_max[i], totareas[i], totareas_min[i],
                                 totareas_max[i], thresh])

    return (means, medians, totareas, titles, means_min, means_max,
            medians_min, medians_max, totareas_min, totareas_max)


def computeHagg(grid2D, proj='moll', probthresh=0.0, shakefile=None,
                shakethreshtype='pga', shakethresh=0.0):
    """
    Computes the Aggregate Hazard (Hagg) which is equal to the
    probability * area of grid cell For models that compute areal coverage,
    this is equivalant to the total predicted area affected in km2.

    Args:
        grid2D: grid2D object of model output.
        proj: projection to use to obtain equal area, 'moll'  mollweide, or
            'laea' lambert equal area.
        probthresh: Probability threshold, any values less than this will not
            be included in aggregate hazard estimation.
        shakefile: Optional, path to shakemap file to use for ground motion
            threshold.
        shakethreshtype: Optional, Type of ground motion to use for
            shakethresh, 'pga', 'pgv', or 'mmi'.
        shakethresh: Optional, Float or list of shaking thresholds in %g for
            pga, cm/s for pgv, float for mmi.

    Returns:
        Float if no shakethresh defined or only one shakethresh defined,
        otherwise, a list of aggregate hazard for all shakethresh values.
    """
    Hagg = []
    bounds = grid2D.getBounds()
    lat0 = np.mean((bounds[2], bounds[3]))
    lon0 = np.mean((bounds[0], bounds[1]))
    projs = ('+proj=%s +lat_0=%f +lon_0=%f +x_0=0 +y_0=0 +ellps=WGS84 '
             '+units=km +no_defs' % (proj, lat0, lon0))
    geodict = grid2D.getGeoDict()

    if shakefile is not None:
        if type(shakethresh) is float or type(shakethresh) is int:
            shakethresh = [shakethresh]
        for shaket in shakethresh:
            if shaket < 0.:
                raise Exception('shaking threshold must be equal or greater '
                                'than zero')
        # resample shakemap to grid2D
        temp = ShakeGrid.load(shakefile, samplegeodict=geodict, resample=True,
                              doPadding=True, adjust='res')
        shk = temp.getLayer(shakethreshtype)
        if shk.getGeoDict() != geodict:
            raise Exception('shakemap was not resampled to exactly the same '
                            'geodict as the model')

    if probthresh < 0.:
        raise Exception('probability threshold must be equal or greater '
                        'than zero')

    grid = grid2D.project(projection=projs)
    geodictRS = grid.getGeoDict()
    cell_area_km2 = geodictRS.dx * geodictRS.dy
    model = grid.getData()
    model[np.isnan(model)] = -1.
    if shakefile is not None:
        for shaket in shakethresh:
            modcop = model.copy()
            shkgrid = shk.project(projection=projs)
            shkdat = shkgrid.getData()
            # use -1 to avoid nan errors and warnings, will always be thrown
            # out because default is 0.
            shkdat[np.isnan(shkdat)] = -1.
            modcop[shkdat < shaket] = -1.
            Hagg.append(np.sum(modcop[modcop >= probthresh] * cell_area_km2))
    else:
        Hagg.append(np.sum(model[model >= probthresh] * cell_area_km2))
    if len(Hagg) == 1:
        Hagg = Hagg[0]
    return Hagg


def computeHagg2(grid2D, proj='moll', probthresh=0.0, shakefile=None,
                 shakethreshtype='pga', shakethresh=0.0):
    """
    Alternative Aggregate Hazard (Hagg), which is equal to the
    the sum of the area of grid cell that exceeds a given probability.

    Args:
        grid2D: grid2D object of model output.
        proj: projection to use to obtain equal area, 'moll'  mollweide, or
            'laea' lambert equal area.
        probthresh: Probability threshold, any values less than this will not
            be included in aggregate hazard estimation.
        shakefile: Optional, path to shakemap file to use for ground motion
            threshold.
        shakethreshtype: Optional, Type of ground motion to use for
            shakethresh, 'pga', 'pgv', or 'mmi'.
        shakethresh: Optional, Float or list of shaking thresholds in %g for
            pga, cm/s for pgv, float for mmi.

    Returns:
        Float if no shakethresh defined or only one shakethresh defined,
        otherwise, a list of aggregate hazard for all shakethresh values.
    """
    Hagg = []
    bounds = grid2D.getBounds()
    lat0 = np.mean((bounds[2], bounds[3]))
    lon0 = np.mean((bounds[0], bounds[1]))
    projs = ('+proj=%s +lat_0=%f +lon_0=%f +x_0=0 +y_0=0 +ellps=WGS84 '
             '+units=km +no_defs' % (proj, lat0, lon0))
    geodict = grid2D.getGeoDict()

    if shakefile is not None:
        if type(shakethresh) is float or type(shakethresh) is int:
            shakethresh = [shakethresh]
        for shaket in shakethresh:
            if shaket < 0.:
                raise Exception('shaking threshold must be equal or greater '
                                'than zero')
        # resample shakemap to grid2D
        temp = ShakeGrid.load(shakefile, samplegeodict=geodict, resample=True,
                              doPadding=True, adjust='res')
        shk = temp.getLayer(shakethreshtype)
        if shk.getGeoDict() != geodict:
            raise Exception('shakemap was not resampled to exactly the same '
                            'geodict as the model')

    if probthresh < 0.:
        raise Exception('probability threshold must be equal or greater '
                        'than zero')

    grid = grid2D.project(projection=projs)
    geodictRS = grid.getGeoDict()
    cell_area_km2 = geodictRS.dx * geodictRS.dy
    model = grid.getData()
    model[np.isnan(model)] = -1.
    if shakefile is not None:
        for shaket in shakethresh:
            modcop = model.copy()
            shkgrid = shk.project(projection=projs)
            shkdat = shkgrid.getData()
            # use -1 to avoid nan errors and warnings, will always be thrown
            # out because default probthresh is 0 and must be positive.
            shkdat[np.isnan(shkdat)] = -1.
            modcop[shkdat < shaket] = -1.
            one_mat = np.ones_like(modcop)
            Hagg.append(np.sum(one_mat[modcop >= probthresh] * cell_area_km2))
    else:
        one_mat = np.ones_like(model)
        Hagg.append(np.sum(one_mat[model >= probthresh] * cell_area_km2))
    if len(Hagg) == 1:
        Hagg = Hagg[0]
    return Hagg


def getQuakeInfo(id):
    BASEURL = 'http://earthquake.usgs.gov/fdsnws/event/1/query?'
    # here we're using the request submodule to open a url, just like we
    # would use open() to open a file.
    contribs = ['', 'us', 'nc', 'usp', 'atlas', 'ci', 'ak', 'at', 'cgs', 'hv',
                'ismp', 'ld', 'mb', 'nc', 'nm', 'nn', 'np', 'pr', 'pt', 'se',
                'us', 'uu', 'uw']
    for con in contribs:
        indict = {'format': 'geojson',
                  'eventid': con + id,
                  'producttype': 'origin',
                  'includesuperseded': 'false'}
        # the urllib module (in Python 3) encapsulates everything the standard
        # library knows about urls
        params = urllib.parse.urlencode(indict)
        try:
            # assemble complete url
            url = '%s%s' % (BASEURL, params)
            f = urllib.request.urlopen(url)
            break
        except:
            continue
    # urlopen() returns a file-like object, which means that it behaves just
    # like a file object does, including allowing you to read in all the data
    # from the thing that has been opened. Note the decode() method, which is
    # a new necessity in Python 3, in order to convert a string of bytes
    # into ASCII (utf-8).
    data = f.read().decode('utf-8')
    # Always close your file!
    f.close()
    jsondict = json.loads(data)
    title = jsondict['properties']['title']
    origin = jsondict['properties']['products']['origin'][0]
    time = origin['properties']['eventtime']
    magnitude = jsondict['properties']['mag']
    return title, time, magnitude


def convert2Coverage(gdict, inventory, numdiv=30.0, method='nearest',
                     proj='moll'):
    """
    Fast method to produce grid of area actually affected by landsliding in
    each cell defined by geodict

    Args:
        gdict: Geodict, likely taken from model to compare inventory against.
        inventory: Path to shapefile of inventory, must be in geographic
            coordinates, WGS84.
        numdiv: Approximate amount to subdivide each cell of geodict by to
            compute areas (higher number slower but more accurate).
        method: Method for resampling when projecting back to geographic
            coordinates, nearest recommended but not perfect. Cubic is not
            recommended.
        proj: Projection to use to obtain equal area,  'moll'  mollweide, or
            'laea' lambert equal area.

    Returns:
        Grid2D object reporting approximate areal coverage of input inventory
        corresponding to geodict
    """

    lat0 = np.mean((gdict.ymin, gdict.ymax))
    lon0 = np.mean((gdict.xmin, gdict.xmax))
    gdsubdiv = {
        'xmin': gdict.xmin - 3 * gdict.dx,
        'xmax': gdict.xmax + 3 * gdict.dx,
        'ymin': gdict.ymin - 3 * gdict.dy,
        'ymax': gdict.ymax + 3 * gdict.dy,
        'dx': gdict.dx / numdiv,
        'dy': gdict.dy / numdiv,
        'ny': gdict.ny * numdiv,
        'nx': gdict.nx * numdiv
    }
    subgd = GeoDict(gdsubdiv, adjust='res')

    with fiona.open(inventory) as f:
        invshp = list(f.items())
    shapes = [shape(inv[1]['geometry']) for inv in invshp]

    # Rasterize with oversampled area
    rast = Grid2D.rasterizeFromGeometry(shapes, subgd, fillValue=0.,
                                        burnValue=1.0, mustContainCenter=True)
    projs = ('+proj=%s +lat_0=%f +lon_0=%f +x_0=0 +y_0=0 +ellps=WGS84 '
             '+units=m +no_defs' % (proj, lat0, lon0))
    equal_area = rast.project(projection=projs, method='nearest')
    egdict = equal_area.getGeoDict()

    # Use block mean to get %coverage of larger grid
    bm, xdiff, ydiff = block_mean(equal_area.getData(), factor=numdiv,
                                  replacenan=0.0)

    # complicated adjustment to get the alignment right
    gdds = {
        'xmin': egdict.xmin + 0.5 * egdict.dx * (numdiv - 1),
        'xmax': egdict.xmax - egdict.dx * (xdiff - 0.5 + numdiv)
                + (numdiv * egdict.dx) / 2.,
        'ymin': egdict.ymin + egdict.dy * (ydiff - 0.5 + numdiv)
                - (egdict.dy * numdiv) / 2.,
        'ymax': egdict.ymax + 0.5 * egdict.dy - 0.5 * numdiv * egdict.dy / 2,
        'dx': egdict.dx * numdiv,
        'dy': egdict.dy * numdiv,
        'ny': np.floor(egdict.ny / numdiv),
        'nx': np.floor(egdict.nx / numdiv),
        'projection': projs
    }
    dsgd = GeoDict(gdds, adjust='res')
    dsgd.setProjection(projs)
    eabig = Grid2D(data=bm, geodict=dsgd)
    eabig = GDALGrid.copyFromGrid(eabig)

    # Project back
    eabigproj = eabig.project(projection=gdict.projection)

    # Resample to original grid
    inventorygrid = eabigproj.interpolateToGrid(gdict, method='linear')

    return inventorygrid


def convert2Prob(gdict, inventory, mustContainCenter=False):
    """
    Convert inventory shapefile to binary grid (with geodict of gdict) with 1
    if any part of landslide was in a cell, 0 if not.


    Args:
        gdict: Geodict, likely taken from model to compare inventory against.
        inventory: Path to shapefile of inventory, must be in geographic
            coordinates, WGS84.
        mustContainCenter: Bool indicating whether the geometry must touch
            the center of the cell versus be inside the cell extent in order
            to set the value. Note that when false the class balance is not
            retained, but when true, the class balance of the reuslting
            grid is unbiased with respect to teh inventory.

    Returns:
        Grid2D object containing rasterized version of inventory where cell is
        1. if any part of a landslide was in the cell.
    """
    with fiona.open(inventory) as f:
        invshp = list(f.items())

    shapes = [shape(inv[1]['geometry']) for inv in invshp]

    # Rasterize with allTouch
    rast = Grid2D.rasterizeFromGeometry(shapes, gdict, fillValue=0.,
                                        burnValue=1.0,
                                        mustContainCenter=mustContainCenter)

    return rast


def statsCoverage(modelgrid, inventorygrid, bins=None, showplots=True,
                  saveplots=False, filepath=None):
    """
    Compute stats and make comparison plots specific to models that output
    areal coverage like Godt et al 2008.

    Args:
        modelgrid: Grid2D object of model results.
        inventorygrid: Grid2D object of areal coverage of inventory computed on
            same grid as modelgrid using, for example, convert2Coverage.
        bins: Bin edges to use for various binning and threshold statistical
            calculations. If None

            .. code-block:: python

                bins = np.linspace(0, np.max((inv.max(), model.max())), 11)

        showplots: Display the plots?
        saveplots: Save the plots?
        filepath: Filepath for saved plots. If None, will save in current
            directory. Files are named with test name and time stamp.

    Returns:
        tuple: Two dictionaries:
            * results: dictionary of results of tests.

                .. code-block:: python

                    {
                        'Compare coverage': dictionary,
                        'RMS': float,
                        'RMS_nonzero': float,
                        'Percent in bin': dictionary,
                        'Direct Comparison': dictionary
                    }

            * invminusmod: Grid2D object of difference between inventory and
              model (inventory - model)

    """

    inv = inventorygrid.getData()
    model = modelgrid.getData()
    invminusmod = GDALGrid(inv - model, inventorygrid.getGeoDict())

    plt.ioff()

    # Make statistical comparisons
    results = {}
    if bins is None:
        bins = np.linspace(0, np.max((np.nanmax(inv), np.nanmax(model))), 11)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    notnans = ~np.isnan(model) & ~np.isnan(inv)
    ax.hist((model[notnans], inv[notnans]), bins=bins)
    ax.set_yscale('log')
    ax.legend(('Model', 'Inventory'))
    ax.set_ylabel('Total # of cells')
    ax.set_xlabel('Predicted coverage bin')
    results['Direct Comparison'] = {
        'bins': bins,
        'model': model,
        'inventory': inv
    }

    # Statistical comparison
    perc = []
    areatot = []
    for i, bin in enumerate(bins[:-1]):
        idx = (model > bin) & (model <= bins[i + 1])
        areain = inv[idx]
        totalin = sum((areain > bin) & (areain <= bins[i + 1]))
        perc.append(float(totalin) / len(idx) * 100)
        areatot.append(np.mean(areain))
    areatot = np.nan_to_num(areatot)

    # get centerpoint of bins
    binvec = []
    for i in range(len(bins[:-1])):
        binvec.append(bins[i] + (bins[i + 1] - bins[i]) / 2)
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(binvec, perc, 'o-')
    ax1.set_ylabel('% in right bin')
    ax1.set_xlabel('Predicted coverage bin')
    ax1.set_title('Actual coverage in predicted bin')
    results['Percent in bin'] = {
        'bincenters': binvec,
        'perc': perc
    }

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(binvec, areatot, 'o-')
    ax2.plot([0., 1.], [0., 1.], '--', color='gray')
    ax2.set_ylabel('Actual coverage')
    ax2.set_xlabel('Predicted coverage')
    ax2.set_title('Actual vs. Predicted coverage')
    results['Compare coverage'] = {
        'bincenters': binvec,
        'areatot': areatot
    }

    if showplots is True:
        plt.show()
    if saveplots is True:
        if filepath is None:
            filepath = os.getcwd()

        time1 = datetime.datetime.utcnow().strftime('%d%b%Y_%H%M')
        fig.savefig(os.path.join(filepath, 'Direct_compare_%s.pdf' % (time1,)))
        fig1.savefig(os.path.join(
            filepath, 'Percent_in_bin_%s.pdf' % (time1,)))
        fig2.savefig(os.path.join(
            filepath, 'Compare_coverage_%s.pdf' % (time1,)))

    # RMS
    idx = np.union1d(model.nonzero()[0], inv.nonzero()[0])
    results['RMS'] = np.sqrt((model - inv)**2).mean()
    print(('RMS: %0.3f' % results['RMS']))
    results['RMS_nonzero'] = np.sqrt((model[idx] - inv[idx])**2).mean()
    print(('RMS_nonzero: %0.3f' % results['RMS_nonzero']))

    return invminusmod, results


def stats(modelgrid, inventory, dx=100.0, Nsamp=100, method='nearest',
          extent='inventory', bins=None, runtests=True, showplots=True,
          saveplots=False, filepath=None):
    """
    Run through suite of tests for models that output probability or index
    that varies between 0 and 1.

    Args:
        modelgrid: Grid2D object of model results.
        inventory: Path to shapefile of inventory, must be in geographic
            coordinates, WGS84.
        inventory: Path to inventory file.
        dx: Approximate sample spacing in meters, overwritten if Nsamp is
            defined.
        Nsamp: Total number of samples desired. Will choose optimal dx to get
            slightly more than this number of samples delete samples outside
            bounds and randomly delete others until sample number is exactly
            Nsamp.
        method: Method used for interp2d when transforming sampled model
            values back from projected coordinates to geographic coordinates.
            Can be 'nearest', 'linear', or 'cubic'.
        extent: Extent to include in sampling. Can be 'inventory', 'model',
            or a tuple like `(xmin, ymin, xmax, ymax)`.
        bins: Bin edges to use for various binning and threshold statistical
            calculations. If None `bins = [0, 0.2, 0.4, 0.6, 0.8, 1.]`
        runtests: If True, will run various statistical tests, if False will
            just output sampled values.
        showplots: Disply the plots?
        saveplots: Save the plots?
        filepath: Filepath for saved plots. If None, will save in current
            directory. Files are named with test name and time stamp.

    Returns:
        A tuple including:
            * yespoints: Nx2 array of geographic coordinates of positive
              sample locations.
            * nopoints: Nx2 array of geographic coordinates of negative sample
              locations.
            * modelvalyes: N model output values corresponding to yespoints.
            * modelvalno: N model output values corresponding to nopoints.
            * results: Dictionary of results of statistical tests. Will be
              empty if runtests=False

            .. code-block:: python

                {
                    'Occ_nonocc': dict,
                    'SRC': dict,
                    'ROC': dict,
                    'AUC_ROC': float,
                    'Log_loss': float,
                    'GFC': dict,
                    'Pred_vs_Obs': dict,
                    'Brier': float,
                    'Brier_no': float,
                    'Brier_yes': float
                }
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

    yespoints, nopoints, xvar, yvar, pshapes, proj = pointsFromShapes(
        shapes, extent, dx=dx, Nsamp=Nsamp)
    yesptx = [pt[0] for pt in yespoints]
    yespty = [pt[1] for pt in yespoints]
    noptx = [pt[0] for pt in nopoints]
    nopty = [pt[1] for pt in nopoints]
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
        func = interpolate.interp2d(
            lons, lats, modelgrid.getData(), kind=method.lower())
        modelvalyes = np.array([float(func(XX, YY))
                                for XX, YY in zip(yesptx, yespty)])
        modelvalno = np.array([float(func(XX, YY))
                               for XX, YY in zip(noptx, nopty)])

    modelvalyes = np.nan_to_num(
        np.array(modelvalyes))  # replace nan with zeros
    modelvalno = np.nan_to_num(np.array(modelvalno))  # replace nan with zeros

    # Now run the desired tests and make the desired plots
    results = {}

    if runtests is True:
        # Brier score
        N = len(yespoints) + len(nopoints)
        yessum = np.sum([(val - 1)**2 for val in modelvalyes])
        nosum = np.sum([(val)**2 for val in modelvalno])
        results['Brier_yes'] = yessum / len(modelvalyes)
        results['Brier_no'] = nosum / len(modelvalno)
        results['Brier'] = (yessum + nosum) / N
        print(('Brier scores: overall %0.3f\n'
               'Brier_yes score: %0.3f\n'
               'Brier_no score %0.3f'
               % (results['Brier'],
                  results['Brier_yes'],
                  results['Brier_no'])))

        # Logarithmic score
        tempno = np.array(modelvalno).copy()
        tempyes = np.array(modelvalyes).copy()
        tempno[tempno == 0] = 1.e-15
        tempyes[tempyes == 0] = 1.e-15
        results['Log_loss'] = - \
            (np.sum(np.log(tempyes)) + np.sum(np.log(1. - tempno))) / N
        print(('Log loss score: %0.3f' % (results['Log_loss'],)))

        if bins is None:
            bins = [0, 0.2, 0.4, 0.6, 0.8, 1.]
        binvec = []
        observed = []
        percyes = []
        percno = []
        overall_tot = len(modelvalyes) + len(modelvalno)
        for i in range(len(bins[:-1])):
            binvec.append(bins[i] + (bins[i + 1] - bins[i]) / 2)
            yestot = np.sum([(modelvalyes > bins[i]) &
                             (modelvalyes < bins[i + 1])])
            notot = np.sum(
                [(modelvalno > bins[i]) & (modelvalno < bins[i + 1])])
            if notot + yestot != 0:
                observed.append(float(yestot) / (yestot + notot))
            else:
                observed.append('nan')
            percyes.append((yestot / float(overall_tot)) * 100.)
            percno.append((notot / float(overall_tot)) * 100.)

        plt.ioff()

        # Predicted vs. Observed ratios
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(binvec, observed, '-o')
        ax.plot([0] + binvec, [0] + binvec, '--', color='gray')
        ax.set_xlabel('Expected ratio')
        ax.set_ylabel('Observed ratio')
        ax.set_xlim([bins[0], bins[-1]])
        ax.set_title('Predicted vs. Observed')
        results['Pred_vs_Obs'] = {'binvec': binvec, 'observed': observed}

        # Ground failure occurrence/nonoccurrence
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        wid = (bins[1] - bins[0]) / 2.5
        rects1 = ax1.bar(np.array(bins[:-1]), percyes, width=wid)
        rects2 = ax1.bar(np.array(bins[:-1]) +
                         wid, percno, width=wid, color='r')
        ax1.set_xlabel('Predicted susceptibility range')
        ax1.set_ylabel('% of samples')
        ax1.legend((rects1[0], rects2[0]), ('Occurrence', 'Nonoccurrence'))
        ax1.set_title('Occurrence vs. Nonoccurrence')
        results['Occ_nonocc'] = {'bins': bins,
                                 'percyes': percyes, 'percno': percno}

        # Ground failure capture for various thresholds
        gfc = []
        for val in bins:
            gfc.append(np.sum([modelvalyes > val]) / float(len(yespoints)))
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        ax2.plot(bins, gfc, 'o-')
        ax2.set_xlabel('Threshold')
        ax2.set_ylabel(r'%GFC')
        ax2.set_title('Ground Failure Capture')
        results['GFC'] = {'thresholds': bins, 'gfc': gfc}

        # ROC curves
        y_true = np.concatenate((np.ones(len(yespoints)),
                                 np.zeros(len(nopoints))))
        y_score = np.concatenate((modelvalyes, modelvalno))
        fpr, tpr, thresholds = roc_curve(y_true, y_score)
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

        results['AUC_ROC'] = roc_auc_score(y_true, y_score)
        print(('AUC_ROC: %0.3f' % (results['AUC_ROC'],)))
        ax3.text(0.8, 0.2, 'AUC: %0.3f' % results['AUC_ROC'])

        # Success rate curves
        sucbin = np.linspace(0, 1., 100)
        prop = []
        realvals = np.concatenate(
            (np.ones(len(yespoints)), np.zeros(len(nopoints))))
        predvals = np.concatenate((modelvalyes, modelvalno))
        indx = np.argsort(predvals)
        predvals = predvals[indx]
        realvals = realvals[indx]
        for val in sucbin:
            prop.append(np.sum(realvals[predvals < val]) / len(yespoints))
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
            time1 = datetime.datetime.utcnow().strftime('%d%b%Y_%H%M')
            fig.savefig(os.path.join(
                filepath, 'Pred_vs_obs_%s.pdf' % (time1,)))
            fig1.savefig(os.path.join(
                filepath, 'Occ_nonocc_%s.pdf' % (time1,)))
            fig2.savefig(os.path.join(filepath, 'GFC_%s.pdf' % (time1,)))
            fig3.savefig(os.path.join(filepath, 'ROC_%s.pdf' % (time1,)))
            fig4.savefig(os.path.join(filepath, 'SRC_%s.pdf' % (time1,)))

    return yespoints, nopoints, modelvalyes, modelvalno, results


def normXcorr(model, inventory):
    """
    Perform normalized cross correlation (no shifts) to assess how well the
    spatial extent was predicted, regardless of overall values.

    Args:
        model: Grid2D file of model.
        inventory: Grid 2D file of inventory processed to simulate what model
            is supposed to predict (using convert2Prob or convert2Coverage).

    Returns: Normalized cross correlation coefficient (between 0 and 1).
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
    From here:
    http://stackoverflow.com/questions/18666014/downsample-array-in-python

    Args:
        ar: 2D array to downsample using a block mean
        factor: factor by which to downsample

    Returns:
        tuple of:
            - downsampled array
            - number of cells added to right edge
            - number of cells added to bottom edge
    """

    if replacenan is not None:
        ar[np.isnan(ar)] = replacenan
    sx, sy = ar.shape
    # get new dimensions (cutting off extra on edges)
    newx = int(np.floor(sx / float(factor)))
    newy = int(np.floor(sy / float(factor)))
    # build regions over which to average
    vec = np.reshape(np.arange(newx * newy), (newx, newy))
    regions = vec.repeat(factor, 0).repeat(factor, 1)
    # Patch on edges (just expand existing cells on the edges a little)
    # to account for rounding down earlier
    xdiff = ar.shape[0] - regions.shape[0]
    if xdiff != 0:
        row = regions[-1, :]
        regions = np.row_stack((regions, np.repeat(np.reshape(
            row, (1, len(row))), xdiff, axis=0)))  # add onto the bottom
    ydiff = ar.shape[1] - regions.shape[1]
    if ydiff != 0:
        col = regions[:, -1]
        regions = np.column_stack((regions, np.repeat(np.reshape(
            col, (len(col), 1)), ydiff, axis=1)))  # add onto the right side
    res = ndimage.mean(ar, labels=regions, index=np.arange(regions.max() + 1))
    res = res.reshape(newx, newy)
    return res, xdiff, ydiff
