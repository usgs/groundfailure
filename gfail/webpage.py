#!/usr/bin/env python


import matplotlib.pyplot as plt
import json
import numpy as np
import os
from configobj import ConfigObj
from gfail.makemaps import setupsync, get_zoomextent, make_rgba, DFCOLORS,\
    DFBINS, make_legend
from gfail.utilities import parseConfigLayers, get_alert
from gfail.stats import computeStats
from mapio.shake import ShakeGrid
import matplotlib.cm as cm

from impactutils.textformat.text import set_num_precision
# from impactutils.time.ancient_time import HistoricTime as ShakeDateTime
# import pytz

from gfail.utilities import loadlayers
# from gfail.utilities import is_grid_point_source


# temporary until mapio is updated
import warnings
warnings.filterwarnings('ignore')

plt.switch_backend('agg')


def hazdev(maplayerlist, configs, shakemap, outfolder=None, alpha=0.7,
           shakethreshtype='pga', probthresh=None, shakethresh=10.,
           prefLS='Nowicki Jessee and others (2017)',
           prefLQ='Zhu and others (2017)',
           pop_file=None, defaultcolors=True, point=True,
           pager_alert='', eventsource='', eventsourcecode=''):
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
        point (bool): if True, event is a point source and warning should be
            displayed
        pager_alert (str): PAGER alert level, e.g., 'green'. 'pending', ...

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
                if 'jessee' in maplayer['model']['description']['name'].lower():
                    id1 = 'jessee_2017'
                    statprobthresh = 0.002
                else:
                    id1 = 'nowicki_2014_global'
                    statprobthresh = 0.0

            stats = computeStats(
                maplayer['model']['grid'],
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
                lsext = get_zoomextent(maplayer['model']['grid'])
                ls_haz_value = set_num_precision(
                    stats['Hagg_0.10g'], 2, 'float')
                ls_pop_value = set_num_precision(
                    stats['exp_pop_0.10g'], 2, 'int')
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
                statprobthresh = 0.0
            elif '2017' in maplayer['model']['description']['name'].lower():
                id1 = 'zhu_2017_general'
                statprobthresh = 0.005

            stats = computeStats(
                maplayer['model']['grid'],
                probthresh=probthresh,
                shakefile=shakemap,
                shakethresh=shakethresh,
                pop_file=pop_file,
                shakethreshtype=shakethreshtype,
                statprobthresh=statprobthresh)

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
                lqext = get_zoomextent(maplayer['model']['grid'])
                lq_haz_value = set_num_precision(
                    stats['Hagg_0.10g'], 2, 'float')
                lq_pop_value = set_num_precision(
                    stats['exp_pop_0.10g'], 2, 'int')

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
        sync, colorlistLS, reflims = setupsync(
            prefLS, concLS, limLS, colLS,
            defaultcolormap, logscale=logLS,
            alpha=alpha)

        if reflims is None:
            raise Exception('Check input config files, they must all have the '
                            'same number of bin edges')
        else:
            # Stuff colors into dictionary
            for ls in lsmodels:
                ls['bin_colors'] = list(colorlistLS)

        sync, colorlistLQ, reflims = setupsync(
            prefLQ, concLQ, limLQ, colLQ,
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

    # If PAGER alert is pending, overwrite our alerts
    if pager_alert == 'pending':
        for ls in lsmodels:
            ls['alert'] = 'pending'
            ls['hazard_alert']['color'] = 'pending'
            ls['population_alert']['color'] = 'pending'
        for lq in lqmodels:
            lq['alert'] = 'pending'
            lq['hazard_alert']['color'] = 'pending'
            lq['population_alert']['color'] = 'pending'

    # Create info.json
    infojson = create_info(outfolder, lsmodels, lqmodels, eventsource,
                           eventsourcecode, point)
    filenames.append(infojson)

    return filenames


def create_png(event_dir, lsmodels=None, lqmodels=None, mercator=True,
               lsmask=0.002, lqmask=0.005, legends=False,
               eventsource='', eventsourcecode=''):
    """
    Creates transparent PNG file for website.

    Args:
        event_dir (srt): Directory containing ground failure results.
        lsmodels (list): List of dictionaries of model summary info compiled
            by the hazdev function. If not specified, code will search for
            the hdf5 files for the preferred model and will create this
            dictionary and will apply default colorbars and bins.
        lqmodels (list): Same as above for liquefaction.
        mercator (bool): Project raster to web mercator
        lsmask (float): Mask all landslide cells with probabilities below this
            threshold
        lqmask (float): Mask all liquefaction cells with probabilities below
            this threshold
        legends (bool): if True, will produce png files of legends for each
            preferred model

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
            filesnippet = 'jessee_2017'
            out = make_rgba(ls_mod['model']['grid'], mask=lsmask,
                            mercator=mercator, levels=DFBINS,
                            colorlist=DFCOLORS)
            rgba_img, ls_extent, lmin, lmax, cmap = out
            filen = os.path.join(event_dir, '%s_extent.json' % filesnippet)
            filenames.append(filen)
            with open(filen, 'w') as f:
                json.dump(ls_extent, f)

            filen = os.path.join(event_dir, '%s.png' % filesnippet)
            plt.imsave(filen,
                       rgba_img,
                       vmin=lmin,
                       vmax=lmax,
                       cmap=cmap
                       )
        else:
            raise OSError("Preferred landslide model result (%s) not found." %
                          ls_mod_file)
    else:
        for lsm in lsmodels:
            # if lsm['preferred']:
            filesnippet = lsm['id']

            fsh = '%s.hdf5' % filesnippet
            ls_mod_file = [f for f in files if fsh in f]
            if len(ls_mod_file) == 1:
                ls_file = os.path.join(event_dir, ls_mod_file[0])
                ls_mod = loadlayers(ls_file)
            else:
                raise OSError(
                    "Specified landslide model result (%s) not found." % fsh)

            out = make_rgba(ls_mod['model']['grid'], mask=lsmask,
                            mercator=mercator, levels=DFBINS,
                            colorlist=DFCOLORS)
            rgba_img, ls_extent, lmin, lmax, cmap = out

            filen = os.path.join(event_dir, '%s_extent.json' % filesnippet)
            filenames.append(filen)
            with open(filen, 'w') as f:
                json.dump(ls_extent, f)

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
            filesnippet = 'zhu_2017_general'
            tempbins = DFBINS
            tempbins[0] = 0.005
            out = make_rgba(lq_mod['model']['grid'], mask=lqmask,
                            mercator=mercator, levels=tempbins,
                            colorlist=DFCOLORS)
            rgba_img, lq_extent, lmin, lmax, cmap = out            

            filen = os.path.join(event_dir, '%s_extent.json' % filesnippet)
            filenames.append(filen)
            with open(filen, 'w') as f:
                json.dump(lq_extent, f)

            filen = os.path.join(event_dir, '%s.png' % filesnippet)
            plt.imsave(filen,
                       rgba_img,
                       vmin=lmin,
                       vmax=lmax,
                       cmap=cmap
                       )
            filenames.append(filen)
        else:
            raise OSError("Preferred liquefaction model result (%s) not found."
                          % lq_mod_file)
    else:
        for lqm in lqmodels:
            #if lqm['preferred']:
            filesnippet = lqm['id']
            tempbins = DFBINS
            tempbins[0] = 0.005
            fsh = '%s.hdf5' % filesnippet
            lq_mod_file = [f2 for f2 in files if fsh in f2]
            if len(lq_mod_file) == 1:
                lq_file = os.path.join(event_dir, lq_mod_file[0])
                lq_mod = loadlayers(lq_file)
            else:
                raise OSError(
                    "Specified liquefaction model result (%s) not found."
                    % fsh)

            out = make_rgba(lq_mod['model']['grid'], mask=lqmask,
                            mercator=mercator, levels=tempbins,
                            colorlist=DFCOLORS)
            rgba_img, lq_extent, lmin, lmax, cmap = out
            
            filen = os.path.join(event_dir, '%s_extent.json' % filesnippet)
            filenames.append(filen)
            with open(filen, 'w') as f:
                json.dump(lq_extent, f)

            filen = os.path.join(event_dir, '%s.png' % filesnippet)
            plt.imsave(filen,
                       rgba_img,
                       vmin=lmin,
                       vmax=lmax,
                       cmap=cmap
                       )
            filenames.append(filen)

    if legends:
        lsname, lqname = make_legends(
            lqmin=lqmask, lsmin=lsmask, outfolder=event_dir)
        filenames.append(lsname)
        filenames.append(lqname)

    return filenames


def create_info(event_dir, lsmodels=None, lqmodels=None,
                eventsource='', eventsourcecode='', point=True):
    """Create info.json for ground failure product.

    Args:
        event_dir (srt): Directory containing ground failure results.
        lsmodels (list): List of dictionaries of model summary info compiled
            by the hazdev function. If not specified, code will search for
            the hdf5 files for the preferred model and will create this
            dictionary and will apply default colorbars and bins.
        lqmodels (list): Same as above for liquefaction.
        point (bool): if True, event is a point source and warning should be
            displayed

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
            lsext = get_zoomextent(ls_mod['model']['grid'])
        else:
            raise OSError("Preferred landslide model result not found.")
        lq_mod_file = [f2 for f2 in files if 'zhu_2017_general.hdf5' in f2]
        if len(lq_mod_file) == 1:
            lq_file = os.path.join(event_dir, lq_mod_file[0])
            lq_mod = loadlayers(lq_file)
            # get extents
            lqext = get_zoomextent(lq_mod['model']['grid'])
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
                # lsm['filename'] = flnm
                lsext = lsm['zoomext']  # Get zoom extent
                ls_alert = lsm['alert']
                rmkeys = ['bin_edges', 'bin_colors', 'zoomext']
            else:
                # Remove any alert keys
                rmkeys = ['bin_edges', 'bin_colors', 'zoomext',
                          'population_alert', 'alert', 'hazard_alert']
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
                # lqm['filename'] = flnm
                lqext = lqm['zoomext']  # Get zoom extent
                lq_alert = lqm['alert']
                rmkeys = ['bin_edges', 'bin_colors', 'zoomext']
            else:
                # Remove any alert keys
                rmkeys = ['bin_edges', 'bin_colors', 'zoomext',
                          'population_alert', 'alert', 'hazard_alert']
            for key in rmkeys:
                if key in lqm:
                    lqm.pop(key)

    # Try to get event info
    shake_grid = ShakeGrid.load(shakefile, adjust='res')
    event_dict = shake_grid.getEventDict()
    sm_dict = shake_grid.getShakeDict()
    base_url = 'https://earthquake.usgs.gov/earthquakes/eventpage/'

    # Is this a point source?
    # point = is_grid_point_source(shake_grid)
    # Temporarily hard code this until we can get a better solution via
    # new grid.xml attributes.
    #point = True

    net = eventsource
    code = eventsourcecode
    time = event_dict['event_timestamp'].strftime('%Y-%m-%dT%H:%M:%SZ')

    event_url = '%s%s%s#executive' % (base_url, net, code)

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
            'time': time,
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
        json.dump(info_dict, f)  # allow_nan=False)
    filenames.append(info_file)
    return filenames


def make_legends(lqmin=0.005, lsmin=0.002, outfolder=None,
                 orientation='horizontal', transparent=True):
    """Make png file of legends to go with pngs using DFCOLORS and DFBINS

    Args:
        lqmin (float): minimum visible value of liquefaction probability
        lsmin (float): same as above for landslides
        outfolder (float): folder to place pngs of legends
        orientation (str): orientation of colorbar, 'horizontal' or 'vertical'
        transparent (bool): if True, background will be transparent


    Returns:
        tuple: (lsfilename, lqfilename), locations of files that were created.

    """

    # Make liquefaction colorbar
    binedges = DFBINS
    binedges[0] = lqmin
    if outfolder is None:
        lqfilename = 'legend_liquefaction.png'
    else:
        lqfilename = os.path.join(outfolder, 'legend_liquefaction.png')

    make_legend(binedges, DFCOLORS, filename=lqfilename,
                orientation=orientation, title='Liquefaction probability',
                mask=lqmin, transparent=transparent)

    # --------------------------
    # Make landslide legend

    if outfolder is None:
        lsfilename = 'legend_landslide.png'
    else:
        lsfilename = os.path.join(outfolder, 'legend_landslide.png')

    make_legend(DFBINS, DFCOLORS, filename=lqfilename,
                orientation=orientation, title='Landslide probability',
                mask=lsmin, transparent=transparent)

    return lsfilename, lqfilename
