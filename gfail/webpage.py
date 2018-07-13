#!/usr/bin/env python


import matplotlib.pyplot as plt
import matplotlib as mpl
import json
import numpy as np
import os
from configobj import ConfigObj
from gfail.makemaps import setupsync
from gfail.utilities import parseConfigLayers
from gfail.stats import computeStats
from folium.utilities import mercator_transform
from mapio.shake import ShakeGrid
import matplotlib.cm as cm

from impactutils.textformat.text import set_num_precision
from impactutils.time.ancient_time import HistoricTime as ShakeDateTime
import pytz

from gfail.utilities import get_event_comcat, loadlayers
from gfail.utilities import is_grid_point_source


# temporary until mapio is updated
import warnings
warnings.filterwarnings('ignore')

plt.switch_backend('agg')

# # hex versions:
# DFCOLORS = [
#     '#efefb34D',  # 30% opaque 4D
#     '#e5c72f66',  # 40% opaque 66
#     '#ea720780',  # 50% opaque 80
#     '#c0375c99',  # 60% opaque 99
#     '#5b28b299',  # 60% opaque 99
#     '#1e1e6499'   # 60% opaque 99
# ]

DFCOLORS = [
    [0.94, 0.94, 0.70, 0.3],
    [0.90, 0.78, 0.18, 0.4],
    [0.92, 0.45, 0.03, 0.5],
    [0.75, 0.22, 0.36, 0.6],
    [0.36, 0.16, 0.70, 0.6],
    [0.12, 0.12, 0.39, 0.6]
]

DFBINS = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]


def hazdev(maplayerlist, configs, shakemap, outfolder=None, alpha=0.7,
           shakethreshtype='pga', probthresh=None, shakethresh=10.,
           prefLS='Nowicki Jessee and others (2017)', prefLQ='Zhu and others (2017)',
           pop_file=None, defaultcolors=True):
    """Create all files needed for product page creation
    Assumes gfail has been run already with -w flag

    Args:
        maplayerlist (list): List of model outputs from gfail.
        configs (list): List of dictionaries of config files corresponding to
            each model in maplayerlist and in the same order.
        shakemap (str): path to shakemap .xml file.
        outfolder (str): Location in which to save outputs. If None, will use
            current directory.
        alpha (float): Transparency to use for overlay pngs, value from 0 to 1.
        shakethreshtype (str): Type of ground motion to use for shakethresh,
            'pga', 'pgv', or 'mmi'.
        probthresh: Optional. Float or list of probability thresholds to apply
            before computing stats.
        shakethresh: Float or list of shaking thresholds in %g for pga, cm/s
            for pgv, float for mmi. Used for Hagg and Exposure computation.
        prefLS (str): shortref of "preferred" landslide model.
        prefLQ (str): shortref of "preferred" liquefaction model.
        pop_filt (str): file path to population file used to compute
            population-based alert levels.
        defaultcolors (bool): If True, will use DFCOLORS for all layers instead
            of determining new ones. This will crash if any of the layers have
            a different number of bins than the number of DFCOLORS

    Returns:
        Files that need to be sent to comcat for hazdev to create the product
            webpage including:
                - info.json
                - transparent png overlays of all models
    """
    event_id = maplayerlist[0]['model']['description']['event_id']

    if pop_file is None:
        # Read in default paths to get location of the population grid
        default_file = os.path.join(os.path.expanduser('~'), '.gfail_defaults')
        defaults = ConfigObj(default_file)
        pop_file = defaults['popfile']

    if outfolder is None:
        outfolder = os.path.join(os.getcwd(), event_id)

    filenames = []

    # Separate the LS and LQ models

    concLS = []
    concLQ = []

    lsmodels = []
    lqmodels = []
    logLS = []
    limLS = []
    colLS = []
    logLQ = []
    limLQ = []
    colLQ = []

    for conf, maplayer in zip(configs, maplayerlist):
        mdict = maplayer['model']['description']
        # config = ConfigObj(conf)

        if 'landslide' in mdict['parameters']['modeltype'].lower():
            title = maplayer['model']['description']['name']
            plotorder, logscale, lims, colormaps, maskthreshes = \
                parseConfigLayers(maplayer, conf, keys=['model'])
            logLS.append(logscale[0])
            limLS.append(lims[0])
            colLS.append(colormaps[0])
            concLS.append(title)

            if 'godt' in maplayer['model']['description']['name'].lower():
                statprobthresh = None
                id1 = 'godt_2008'
            else:
                # Since logistic models can't equal one, need to eliminate
                # placeholder zeros before computing stats
                statprobthresh = 0.0
                if 'jessee' in maplayer['model']['description']['name'].lower():
                    id1 = 'jessee_2017'
                else:
                    id1 = 'nowicki_2014_global'

            stats = computeStats(maplayer['model']['grid'],
                                 probthresh=probthresh,
                                 shakefile=shakemap,
                                 shakethresh=shakethresh,
                                 statprobthresh=statprobthresh,
                                 pop_file=pop_file,
                                 shakethreshtype=shakethreshtype)

            metadata = maplayer['model']['description']
            if len(maplayer) > 1:
                inputs = {}
                inkeys = list(maplayer.keys())
                for key in inkeys:
                    if key != 'model':
                        newkey = maplayer[key]['label']
                        inputs[newkey] = maplayer[key]['description']
                metadata['inputs'] = inputs

            if title == prefLS:
                on = True
                ls_haz_alert, ls_pop_alert, _, _, ls_alert, _ = get_alert(
                    stats['Hagg_0.10g'], 0.,
                    stats['exp_pop_0.10g'], 0.)
                lsext = get_extent(maplayer['model']['grid'])
                ls_haz_value = set_num_precision(stats['Hagg_0.10g'], 2, 'float')
                ls_pop_value = set_num_precision(stats['exp_pop_0.10g'], 2, 'int')
            else:
                on = False
                ls_haz_alert = None
                ls_pop_alert = None
                ls_haz_value = None
                ls_pop_value = None
                ls_alert = None
                lsext = None

            edict = {
                'id': id1,
                'title': metadata['name'],
                'overlay': '%s.png' % id1,
                'extent': '%s_extent.json' % id1,
                'units': metadata['units'],
                'preferred': on,
                'alert': ls_alert,
                'hazard_alert': {
                    'color': ls_haz_alert,
                    'value': ls_haz_value,
                    'parameter': 'Aggregate Hazard',
                    'units': 'km^2'
                },
                'population_alert': {
                    'color': ls_pop_alert,
                    'value': ls_pop_value,
                    'parameter': 'Population exposure',
                    'units': 'people'
                },
                'bin_edges': list(lims[0]),
                'probability': {
                    'max': float("%.2f" % stats['Max']),
                    'std': float("%.2f" % stats['Std']),
                    'hagg0.1g': float("%.2f" % stats['Hagg_0.10g']),
                    'popexp0.1g': float("%.2f" % stats['exp_pop_0.10g'])
                },
                'longref': metadata['longref'],
                'parameters': metadata['parameters'],
                'zoomext': lsext
            }

            lsmodels.append(edict)

        elif 'liquefaction' in mdict['parameters']['modeltype'].lower():
            title = maplayer['model']['description']['name']
            plotorder, logscale, lims, colormaps, maskthreshes = \
                parseConfigLayers(maplayer, conf, keys=['model'])
            logLQ.append(logscale[0])
            limLQ.append(lims[0])
            colLQ.append(colormaps[0])
            concLQ.append(title)

            if '2015' in maplayer['model']['description']['name'].lower():
                id1 = 'zhu_2015'
            elif '2017' in maplayer['model']['description']['name'].lower():
                id1 = 'zhu_2017_general'

            stats = computeStats(maplayer['model']['grid'],
                                 probthresh=probthresh,
                                 shakefile=shakemap,
                                 shakethresh=shakethresh,
                                 pop_file=pop_file,
                                 shakethreshtype=shakethreshtype)

            metadata = maplayer['model']['description']
            if len(maplayer) > 1:
                inputs = {}
                inkeys = list(maplayer.keys())
                for key in inkeys:
                    if key != 'model':
                        newkey = maplayer[key]['label']
                        inputs[newkey] = maplayer[key]['description']
                metadata['inputs'] = inputs

            if title == prefLQ:
                on = True
                _, _, lq_haz_alert, lq_pop_alert, _, lq_alert = get_alert(
                    0., stats['Hagg_0.10g'],
                    0., stats['exp_pop_0.10g'])
                lqext = get_extent(maplayer['model']['grid'])
                lq_haz_value = set_num_precision(stats['Hagg_0.10g'], 2, 'float')
                lq_pop_value = set_num_precision(stats['exp_pop_0.10g'], 2, 'int')

            else:
                on = False
                lq_haz_alert = None
                lq_pop_alert = None
                lq_haz_value = None
                lq_pop_value = None
                lq_alert = None
                lqext = None

            edict = {
                'id': id1,
                'title': metadata['name'],
                'overlay': '%s.png' % id1,
                'extent': '%s_extent.json' % id1,
                'units': metadata['units'],
                'preferred': on,
                'alert': lq_alert,
                'hazard_alert': {
                    'color': lq_haz_alert,
                    'value': lq_haz_value,
                    'parameter': 'Aggregate Hazard',
                    'units': 'km^2'
                },
                'population_alert': {
                    'color': lq_pop_alert,
                    'value': lq_pop_value,
                    'parameter': 'Population exposure',
                    'units': 'people'
                },
                'bin_edges': list(lims[0]),
                'probability': {
                    'max': float("%.2f" % stats['Max']),
                    'std': float("%.2f" % stats['Std']),
                    'hagg0.1g': float("%.2f" % stats['Hagg_0.10g']),
                    'popexp0.1g': float("%.2f" % stats['exp_pop_0.10g'])
                },
                'longref': metadata['longref'],
                'parameters': metadata['parameters'],
                'zoomext': lqext
            }

            lqmodels.append(edict)

        else:
            raise Exception("model type is undefined, check "
                            "maplayer['model']['parameters']"
                            "['modeltype'] to ensure it is defined")

    if defaultcolors:
        for ls in lsmodels:
            ls['bin_colors'] = DFCOLORS
        for lq in lqmodels:
            lq['bin_colors'] = DFCOLORS

    else:
        defaultcolormap = cm.CMRmap_r

        # Get colors and stuff into dictionaries
        sync, colorlistLS, reflims = setupsync(prefLS, concLS, limLS, colLS,
                                               defaultcolormap, logscale=logLS,
                                               alpha=alpha)

        if reflims is None:
            raise Exception('Check input config files, they must all have the '
                            'same number of bin edges')
        else:
            # Stuff colors into dictionary
            for ls in lsmodels:
                ls['bin_colors'] = list(colorlistLS)

        sync, colorlistLQ, reflims = setupsync(prefLQ, concLQ, limLQ, colLQ,
                                               defaultcolormap, logscale=logLQ,
                                               alpha=alpha)
        if reflims is None:
            raise Exception('Check input config files, they must all have the '
                            'same number of bin edges')
        else:
            # Stuff colors into dictionary
            for lq in lqmodels:
                lq['bin_colors'] = list(colorlistLQ)

    # Create pngs
    pngfiles = create_png(outfolder, lsmodels, lqmodels)
    filenames.append(pngfiles)

    # Create info.json
    infojson = create_info(outfolder, lsmodels, lqmodels)
    filenames.append(infojson)

    return filenames


def create_png(event_dir, lsmodels=None, lqmodels=None, mercator=True):
    """
    Creates transparent PNG file for website.

    Args:
        event_dir (srt): Directory containing ground failure results.
        lsmodels (list): List of dictionaries of model summary info compiled
            by the hazdev function. If not specified, code will search for
            the hdf5 files for the preferred model and will create this dictionary
            and will apply default colorbars and bins.
        lqmodels (list): Same as above for liquefaction.
        mercator (bool): Project raster to web mercator

    Returns:
        .png map overlays and .json files specifying their mapped extents
    """
    filenames = []
    files = os.listdir(event_dir)
    if lsmodels is None:
        # Read in the "preferred" model for landslides
        ls_mod_file = [f for f in files if 'jessee_2017.hdf5' in f]
        if len(ls_mod_file) == 1:
            ls_file = os.path.join(event_dir, ls_mod_file[0])
            ls_mod = loadlayers(ls_file)
            levels = DFBINS
            colors1 = DFCOLORS
            filesnippet = 'jessee_2017'
            ls_grid = ls_mod['model']['grid']
            ls_data = ls_grid.getData()
            ls_geodict = ls_grid.getGeoDict()
            ls_extent = [
                ls_geodict.xmin - 0.5*ls_geodict.dx,
                ls_geodict.xmax + 0.5*ls_geodict.dx,
                ls_geodict.ymin - 0.5*ls_geodict.dy,
                ls_geodict.ymax + 0.5*ls_geodict.dy,
            ]
            filen = os.path.join(event_dir, '%s_extent.json' % filesnippet)
            filenames.append(filen)
            with open(filen, 'w') as f:
                json.dump(ls_extent, f)

            lmin = levels[0]
            lmax = levels[-1]
            ls_data2 = np.clip(ls_data, lmin, lmax)
            cmap = mpl.colors.ListedColormap(colors1)
            norm = mpl.colors.BoundaryNorm(levels, cmap.N)
            ls_data2 = np.ma.array(ls_data2, mask=np.isnan(ls_data))
            #scalarMap = cm.ScalarMappable(norm=norm, cmap=cmap)
            #rgba_img = scalarMap.to_rgba_array(ls_data2, alpha=alpha)
            rgba_img = cmap(norm(ls_data2))
            if mercator:
                rgba_img = mercator_transform(rgba_img, (ls_extent[2], ls_extent[3]),
                                              origin='upper')

            filen = os.path.join(event_dir, '%s.png' % filesnippet)
            plt.imsave(filen,
                       rgba_img,
                       vmin=lmin,
                       vmax=lmax,
                       cmap=cmap
                       )
        else:
            raise OSError("Preferred landslide model result (%s) not found." % ls_mod_file)
    else:
        for lsm in lsmodels:
            #if lsm['preferred']:
            levels = lsm['bin_edges']
            colors1 = lsm['bin_colors']
            filesnippet = lsm['id']

            fsh = '%s.hdf5' % filesnippet
            ls_mod_file = [f for f in files if fsh in f]
            if len(ls_mod_file) == 1:
                ls_file = os.path.join(event_dir, ls_mod_file[0])
                ls_mod = loadlayers(ls_file)
            else:
                raise OSError("Specified landslide model result (%s) not found." % fsh)

            ls_grid = ls_mod['model']['grid']
            ls_data = ls_grid.getData()
            ls_geodict = ls_grid.getGeoDict()
            ls_extent = [
                ls_geodict.xmin - 0.5*ls_geodict.dx,
                ls_geodict.xmax + 0.5*ls_geodict.dx,
                ls_geodict.ymin - 0.5*ls_geodict.dy,
                ls_geodict.ymax + 0.5*ls_geodict.dy,
            ]
            filen = os.path.join(event_dir, '%s_extent.json' % filesnippet)
            filenames.append(filen)
            with open(filen, 'w') as f:
                json.dump(ls_extent, f)

            lmin = levels[0]
            lmax = levels[-1]
            ls_data2 = np.clip(ls_data, lmin, lmax)
            cmap = mpl.colors.ListedColormap(colors1)
            norm = mpl.colors.BoundaryNorm(levels, cmap.N)
            ls_data2 = np.ma.array(ls_data2, mask=np.isnan(ls_data))
            #scalarMap = cm.ScalarMappable(norm=norm, cmap=cmap)
            #rgba_img = scalarMap.to_rgba_array(ls_data2, alpha=alpha)
            rgba_img = cmap(norm(ls_data2))
            if mercator:
                rgba_img = mercator_transform(rgba_img, (ls_extent[2], ls_extent[3]),
                                              origin='upper')

            filen = os.path.join(event_dir, '%s.png' % filesnippet)
            plt.imsave(filen,
                       rgba_img,
                       vmin=lmin,
                       vmax=lmax,
                       cmap=cmap
                       )

    if lqmodels is None:
        # read in preferred model for liquefaction if none specified
        lq_mod_file = [f2 for f2 in files if 'zhu_2017_general.hdf5' in f2]
        if len(lq_mod_file) == 1:
            lq_file = os.path.join(event_dir, lq_mod_file[0])
            lq_mod = loadlayers(lq_file)
            levels = DFBINS
            colors1 = DFCOLORS
            filesnippet = 'zhu_2017_general'
            lq_grid = lq_mod['model']['grid']
            lq_data = lq_grid.getData()
            lq_geodict = lq_grid.getGeoDict()
            lq_extent = [
                lq_geodict.xmin - 0.5*lq_geodict.dx,
                lq_geodict.xmax + 0.5*lq_geodict.dx,
                lq_geodict.ymin - 0.5*lq_geodict.dy,
                lq_geodict.ymax + 0.5*lq_geodict.dy,
            ]
            filen = os.path.join(event_dir, '%s_extent.json' % filesnippet)
            filenames.append(filen)
            with open(filen, 'w') as f:
                json.dump(lq_extent, f)

            lmin = levels[0]
            lmax = levels[-1]
            lq_data2 = np.clip(lq_data, lmin, lmax)
            cmap = mpl.colors.ListedColormap(colors1)
            norm = mpl.colors.BoundaryNorm(levels, cmap.N)
            lq_data2 = np.ma.array(lq_data2, mask=np.isnan(lq_data))
            rgba_img = cmap(norm(lq_data2))
            if mercator:
                rgba_img = mercator_transform(rgba_img, (lq_extent[2], lq_extent[3]),
                                              origin='upper')
            filen = os.path.join(event_dir, '%s.png' % filesnippet)
            plt.imsave(filen,
                       rgba_img,
                       vmin=lmin,
                       vmax=lmax,
                       cmap=cmap
                       )
            filenames.append(filen)
        else:
            raise OSError("Preferred liquefaction model result (%s) not found." % lq_mod_file)
    else:
        for lqm in lqmodels:
            if lqm['preferred']:
                levels = lqm['bin_edges']
                colors1 = lqm['bin_colors']
                filesnippet = lqm['id']

            fsh = '%s.hdf5' % filesnippet
            lq_mod_file = [f2 for f2 in files if fsh in f2]
            if len(lq_mod_file) == 1:
                lq_file = os.path.join(event_dir, lq_mod_file[0])
                lq_mod = loadlayers(lq_file)
            else:
                raise OSError("Specified liquefaction model result (%s) not found." % fsh)

            lq_grid = lq_mod['model']['grid']
            lq_data = lq_grid.getData()
            lq_geodict = lq_grid.getGeoDict()
            lq_extent = [
                lq_geodict.xmin - 0.5*lq_geodict.dx,
                lq_geodict.xmax + 0.5*lq_geodict.dx,
                lq_geodict.ymin - 0.5*lq_geodict.dy,
                lq_geodict.ymax + 0.5*lq_geodict.dy,
            ]
            filen = os.path.join(event_dir, '%s_extent.json' % filesnippet)
            filenames.append(filen)
            with open(filen, 'w') as f:
                json.dump(lq_extent, f)

            lmin = levels[0]
            lmax = levels[-1]
            lq_data2 = np.clip(lq_data, lmin, lmax)
            cmap = mpl.colors.ListedColormap(colors1)
            norm = mpl.colors.BoundaryNorm(levels, cmap.N)
            lq_data2 = np.ma.array(lq_data2, mask=np.isnan(lq_data))
            rgba_img = cmap(norm(lq_data2))
            if mercator:
                rgba_img = mercator_transform(rgba_img, (lq_extent[2], lq_extent[3]),
                                              origin='upper')
            filen = os.path.join(event_dir, '%s.png' % filesnippet)
            plt.imsave(filen,
                       rgba_img,
                       vmin=lmin,
                       vmax=lmax,
                       cmap=cmap
                       )
            filenames.append(filen)

    return filenames


def create_info(event_dir, lsmodels=None, lqmodels=None):
    """Create info.json for ground failure product.

    Args:
        event_dir (srt): Directory containing ground failure results.
        lsmodels (list): List of dictionaries of model summary info compiled
            by the hazdev function. If not specified, code will search for
            the hdf5 files for the preferred model and will create this
            dictionary and will apply default colorbars and bins.
        lqmodels (list): Same as above for liquefaction.

    Returns:
        creates info.json for this event
    """
    filenames = []
    # Find the shakemap grid.xml file
    with open(os.path.join(event_dir, 'shakefile.txt'), 'r') as f:
        shakefile = f.read()

    files = os.listdir(event_dir)

    if lsmodels is None and lqmodels is None:

        # Read in the "preferred" model for landslides and liquefaction
        ls_mod_file = [f2 for f2 in files if 'jessee_2017.hdf5' in f2]
        if len(ls_mod_file) == 1:
            ls_file = os.path.join(event_dir, ls_mod_file[0])
            ls_mod = loadlayers(ls_file)
            # get extents
            lsext = get_extent(ls_mod['model']['grid'])
        else:
            raise OSError("Preferred landslide model result not found.")
        lq_mod_file = [f2 for f2 in files if 'zhu_2017_general.hdf5' in f2]
        if len(lq_mod_file) == 1:
            lq_file = os.path.join(event_dir, lq_mod_file[0])
            lq_mod = loadlayers(lq_file)
            # get extents
            lqext = get_extent(lq_mod['model']['grid'])
        else:
            raise OSError("Preferred liquefaction model result not found.")

        # Read in extents
        ls_extent_file = [
            f2 for f2 in files if 'jessee_2017_extent.json' in f2]
        if len(ls_extent_file) == 1:
            ls_file = os.path.join(event_dir, ls_extent_file[0])
            with open(ls_file) as f:
                jessee_extent = json.load(f)
        else:
            raise OSError("Landslide extent not found.")
        lq_extent_file = [
            f2 for f2 in files if 'zhu_2017_general_extent.json' in f2]
        if len(lq_extent_file) == 1:
            lq_file = os.path.join(event_dir, lq_extent_file[0])
            with open(lq_file) as f:
                zhu_extent = json.load(f)
        else:
            raise OSError("Liquefaction extent not found.")

        # Read in default paths to get location of the population grid
        default_file = os.path.join(os.path.expanduser('~'), '.gfail_defaults')
        defaults = ConfigObj(default_file)
        pop_file = defaults['popfile']

        # Landslide alert statistics
        ls_stats = computeStats(
            ls_mod['model']['grid'],
            probthresh=None,
            shakefile=shakefile,
            shakethresh=10.0,
            shakethreshtype='pga',
            statprobthresh=None,
            pop_file=pop_file)

        # Liquefaction alert statistics
        lq_stats = computeStats(
            lq_mod['model']['grid'],
            probthresh=None,
            shakefile=shakefile,
            shakethresh=10.0,
            shakethreshtype='pga',
            statprobthresh=None,
            pop_file=pop_file)

        # Get alert levels
        ls_haz_level = ls_stats['Hagg_0.10g']
        lq_haz_level = lq_stats['Hagg_0.10g']
        ls_pop_level = ls_stats['exp_pop_0.10g']
        lq_pop_level = lq_stats['exp_pop_0.10g']

        # If hazard alert level is less than 0.1, zero it out
        # (due to rounding to 2 sig digits later, this can give
        #  overly precise results, e.g., 0.000012 if we don't clip,
        #  but this doesn't happen with pop alerts because they are
        #  integers)
        if ls_haz_level < 0.1:
            ls_haz_level = 0.0
        if lq_haz_level < 0.1:
            lq_haz_level = 0.0

        # Convert levels into categories
        alert_info = get_alert(ls_haz_level, lq_haz_level,
                               ls_pop_level, lq_pop_level)
        # Unpack info (I think we are now assuming that the statements will be
        # constructed on the website and so we don't need them here)
        ls_haz_alert, ls_pop_alert, lq_haz_alert, lq_pop_alert, \
            ls_alert, lq_alert = alert_info

        if lsmodels is None:
            lsmodels = [{
                'id': 'nowicki_jessee_2017',
                'title': 'Nowicki Jessee and others (2017)',
                'overlay': 'jessee_2017.png',
                'extent': jessee_extent,
                'units': "Proportion of area affected",
                'preferred': True,
                'alert': ls_alert,
                'hazard_alert': {
                    'color': ls_haz_alert,
                    'value': set_num_precision(ls_haz_level, 2, 'float'),
                    'parameter': 'Aggregate Hazard',
                    'units': 'km^2'
                },
                'population_alert': {
                    'color': ls_pop_alert,
                    'value': set_num_precision(ls_pop_level, 2, 'int'),
                    'parameter': 'Population exposure',
                    'units': 'people'
                },
                'probability': {
                    'max': float("%.2f" % ls_stats['Max']),
                    'std': float("%.2f" % ls_stats['Std']),
                    'hagg0.1g': float("%.2f" % ls_stats['Hagg_0.10g']),
                    'popexp0.1g': float("%.2f" % ls_stats['exp_pop_0.10g'])
                }
            }]
        if lqmodels is None:
            lqmodels = [{
                'id': 'zhu_2017',
                'title': 'Zhu and others (2017)',
                'overlay': 'zhu_2017.png',
                'extent': zhu_extent,
                'units': "Proportion of area affected",
                'preferred': True,
                'alert': lq_alert,
                'hazard_alert': {
                    'color': lq_haz_alert,
                    'value': set_num_precision(lq_haz_level, 2, 'float'),
                    'parameter': 'Aggregate Hazard',
                    'units': 'km^2'
                },
                'population_alert': {
                    'color': lq_pop_alert,
                    'value': set_num_precision(lq_pop_level, 2, 'int'),
                    'parameter': 'Population exposure',
                    'units': 'people'
                },
                'probability': {
                    'max': float("%.2f" % lq_stats['Max']),
                    'std': float("%.2f" % lq_stats['Std']),
                    'hagg0.1g': float("%.2f" % ls_stats['Hagg_0.10g']),
                    'popexp0.1g': float("%.2f" % ls_stats['exp_pop_0.10g'])
                }
            }]
    else:
        # Get all info from dictionaries of preferred events, add in extent
        # and filename
        for lsm in lsmodels:
            # Add extent and filename for preferred model
            if lsm['preferred']:
                filesnippet = lsm['id']
                # Read in extents
                flnm = '%s_extent.json' % filesnippet
                ls_extent_file = [f for f in files if flnm in f]
                if len(ls_extent_file) == 1:
                    ls_file = os.path.join(event_dir, ls_extent_file[0])
                    with open(ls_file) as f:
                        ls_extent = json.load(f)
                else:
                    raise OSError("Landslide extent not found.")
                lsm['extent'] = ls_extent
                #lsm['filename'] = flnm
                lsext = lsm['zoomext']  # Get zoom extent
                ls_alert = lsm['alert']
                rmkeys = ['bin_edges', 'bin_colors', 'zoomext']
            else:
                # Remove any alert keys
                rmkeys = ['bin_edges', 'bin_colors', 'zoomext', 'population_alert',
                          'alert', 'hazard_alert']
            for key in rmkeys:
                if key in lsm:
                    lsm.pop(key)

        for lqm in lqmodels:
            if lqm['preferred']:
                filesnippet = lqm['id']
                # Read in extents
                flnm = '%s_extent.json' % filesnippet
                lq_extent_file = [f2 for f2 in files if flnm in f2]
                if len(lq_extent_file) == 1:
                    lq_file = os.path.join(event_dir, lq_extent_file[0])
                    with open(lq_file) as f:
                        lq_extent = json.load(f)
                else:
                    raise OSError("Liquefaction extent not found.")
                lqm['extent'] = lq_extent
                #lqm['filename'] = flnm
                lqext = lqm['zoomext']  # Get zoom extent
                lq_alert = lqm['alert']
                rmkeys = ['bin_edges', 'bin_colors', 'zoomext']
            else:
                # Remove any alert keys
                rmkeys = ['bin_edges', 'bin_colors', 'zoomext', 'population_alert',
                          'alert', 'hazard_alert']
            for key in rmkeys:
                if key in lqm:
                    lqm.pop(key)

    # Try to get event info
    shake_grid = ShakeGrid.load(shakefile, adjust='res')
    event_dict = shake_grid.getEventDict()
    sm_dict = shake_grid.getShakeDict()
    base_url = 'https://earthquake.usgs.gov/earthquakes/eventpage/'

    # Is this a point source?
    point = is_grid_point_source(shake_grid)

    try:
        #Hopefully this will eventually be more reliable once we get the
        #comcat info directly from the shakemap grid, rather than rely on
        #magnitude/location/time association.
        shakemap_info, detail, temp = get_event_comcat(shakefile)
        event_url = detail.url
        code = detail['code']
        net = detail['net']
        utc = pytz.utc
        detail_time = ShakeDateTime.fromtimestamp(detail['time']/1000.0, utc).strftime('%Y-%m-%dT%H:%M:%SZ')

    except:
        # Hopefully we can eventually remove this....
        event_url = '%s%s#executive' % (base_url, event_dict['event_id'])
        code = 'unknown'
        net = 'unknown'
        detail_time = -999

    # Get extents that work for both unless one is green and the other isn't
    if lq_alert == 'green' and ls_alert != 'green' and ls_alert is not None:
        xmin = lsext['xmin']
        xmax = lsext['xmax']
        ymin = lsext['ymin']
        ymax = lsext['ymax']
    elif lq_alert != 'green' and ls_alert == 'green' and lq_alert is not None:
        xmin = lqext['xmin']
        xmax = lqext['xmax']
        ymin = lqext['ymin']
        ymax = lqext['ymax']
    else:
        xmin = np.min((lqext['xmin'], lsext['xmin']))
        xmax = np.max((lqext['xmax'], lsext['xmax']))
        ymin = np.min((lqext['ymin'], lsext['ymin']))
        ymax = np.max((lqext['ymax'], lsext['ymax']))

    # Should we display the warning about point source?
    rupture_warning = False
    if point and event_dict['magnitude'] > 6.5:
        rupture_warning = True

    # Create info.json for website rendering and metadata purposes
    info_dict = {
        'Summary': {
            'code': code,
            'net': net,
            'magnitude': event_dict['magnitude'],
            'depth': event_dict['depth'],
            'time': detail_time,
            'lat': event_dict['lat'],
            'lon': event_dict['lon'],
            'event_url': event_url,
            'shakemap_version': sm_dict['shakemap_version'],
            'rupture_warning': rupture_warning,
            'point_source': point,
            'zoom_extent': [xmin, xmax, ymin, ymax]
        },
        'Landslides': lsmodels,
        'Liquefaction': lqmodels

    }

    info_file = os.path.join(event_dir, 'info.json')
    with open(info_file, 'w') as f:
        json.dump(info_dict, f)
    filenames.append(info_file)
    return filenames


def get_alert(paramalertLS, paramalertLQ, parampopLS, parampopLQ,
              hazbinLS=[1., 10., 100.], popbinLS=[100, 1000, 10000],
              hazbinLQ=[10., 100., 1000.], popbinLQ=[1000, 10000, 100000]):
    """
    Get alert levels

    Args:
        paramalertLS (float): Hazard statistic of preferred landslide model
        paramalertLQ (float): Hazard statistic of preferred liquefaction model
        parampopLS (float): Exposure statistic of preferred landslide model
        parampopLQ (float): Exposure statistic of preferred liquefaction model
        hazbinLS (list): 3 element list of bin edges for landslide
            hazard alert between Green and Yellow, Yellow and Orange, and
            Orange and Red.
        popbinLS (list): same as above but for population exposure
        hazbinLQ (list): 3 element list of bin edges for liquefaction hazard
            alert between Green and Yellow, Yellow and Orange, and Orange
            and Red.
        popbinLQ (list): same as above but for population exposure

    Returns:
        Returns:
            tuple: (hazLS, popLS, hazLQ, popLQ, LS, LQ) where:
                * hazLS: the landslide hazard alert level (str)
                * popLS: the landslide population alert level (str)
                * hazLQ: the liquefaction hazard alert level (str)
                * popLQ: the liquefaction population alert level (str)
                * LS: the overall landslide alert level (str)
                * LQ: the overall liquefaction alert level (str)

    """
    if paramalertLS is None:
        hazLS = None
    elif paramalertLS < hazbinLS[0]:
        hazLS = 'green'
    elif paramalertLS >= hazbinLS[0] and paramalertLS < hazbinLS[1]:
        hazLS = 'yellow'
    elif paramalertLS >= hazbinLS[1] and paramalertLS < hazbinLS[2]:
        hazLS = 'orange'
    elif paramalertLS > hazbinLS[2]:
        hazLS = 'red'
    else:
        hazLS = None

    if parampopLS is None:
        popLS = None
    elif parampopLS < popbinLS[0]:
        popLS = 'green'
    elif parampopLS >= popbinLS[0] and parampopLS < popbinLS[1]:
        popLS = 'yellow'
    elif parampopLS >= popbinLS[1] and parampopLS < popbinLS[2]:
        popLS = 'orange'
    elif parampopLS > popbinLS[2]:
        popLS = 'red'
    else:
        popLS = None

    if paramalertLQ is None:
        hazLQ = None
    elif paramalertLQ < hazbinLQ[0]:
        hazLQ = 'green'
    elif paramalertLQ >= hazbinLQ[0] and paramalertLQ < hazbinLQ[1]:
        hazLQ = 'yellow'
    elif paramalertLQ >= hazbinLQ[1] and paramalertLQ < hazbinLQ[2]:
        hazLQ = 'orange'
    elif paramalertLQ > hazbinLQ[2]:
        hazLQ = 'red'
    else:
        hazLQ = None

    if parampopLQ is None:
        popLQ = None
    elif parampopLQ < popbinLQ[0]:
        popLQ = 'green'
    elif parampopLQ >= popbinLQ[0] and parampopLQ < popbinLQ[1]:
        popLQ = 'yellow'
    elif parampopLQ >= popbinLQ[1] and parampopLQ < popbinLQ[2]:
        popLQ = 'orange'
    elif parampopLQ > popbinLQ[2]:
        popLQ = 'red'
    else:
        popLQ = None

    num2color = {
        '1': 'green',
        '2': 'yellow',
        '3': 'orange',
        '4': 'red'
    }
    col2num = dict((v, k) for k, v in num2color.items())

    LSnum1 = col2num[hazLS]
    LSnum2 = col2num[popLS]
    LSnum = str(np.max([int(LSnum1), int(LSnum2)]))
    LS = num2color[LSnum]

    LQnum1 = col2num[hazLQ]
    LQnum2 = col2num[popLQ]
    LQnum = str(np.max([int(LQnum1), int(LQnum2)]))
    LQ = num2color[LQnum]

    return hazLS, popLS, hazLQ, popLQ, LS, LQ


def get_extent(grid, propofmax=0.3):
    """
    Get the extent that contains all values with probabilities exceeding a threshold
    in order to determine ideal zoom level for interactive map
    If nothing is above the threshold, uses the full extent

    Args:
        grid: grid2d of model output
        propofmax (float): Proportion of maximum that should be fully included
            within the bounds.

    Returns:
        tuple: (boundaries, zoomed) where,
            * boundaries: a dictionary with keys 'xmin', 'xmax', 'ymin', and 'ymax' that
                defines the boundaries in geographic coordinates.

    """
    maximum = np.nanmax(grid.getData())

    xmin, xmax, ymin, ymax = grid.getBounds()
    lons = np.linspace(xmin, xmax, grid.getGeoDict().nx)
    lats = np.linspace(ymax, ymin, grid.getGeoDict().ny)

    if maximum <= 0.:
        boundaries1 = dict(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
        # If nothing is above the threshold, use full extent
        return boundaries1

    threshold = propofmax * maximum

    row, col = np.where(grid.getData() > float(threshold))
    lonmin = lons[col].min()
    lonmax = lons[col].max()
    latmin = lats[row].min()
    latmax = lats[row].max()

    boundaries1 = {}

    if xmin < lonmin:
        boundaries1['xmin'] = lonmin
    else:
        boundaries1['xmin'] = xmin
    if xmax > lonmax:
        boundaries1['xmax'] = lonmax
    else:
        boundaries1['xmax'] = xmax
    if ymin < latmin:
        boundaries1['ymin'] = latmin
    else:
        boundaries1['ymin'] = ymin
    if ymax > latmax:
        boundaries1['ymax'] = latmax
    else:
        boundaries1['ymax'] = ymax

    return boundaries1
