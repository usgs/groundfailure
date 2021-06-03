#!/usr/bin/env python


# stdlib imports
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as colors
import json
import numpy as np
import os
from configobj import ConfigObj
import tempfile
from datetime import datetime, timedelta

# third party imports
import simplekml
from folium.utilities import mercator_transform
from gfail.utilities import parseConfigLayers, get_alert, get_rangebeta
from gfail.stats import computeStats
from mapio.shake import ShakeGrid


from impactutils.textformat.text import set_num_precision
# from impactutils.time.ancient_time import HistoricTime as ShakeDateTime
# import pytz

from gfail.utilities import loadlayers

# from gfail.utilities import is_grid_point_source


# temporary until mapio is updated
import warnings
warnings.filterwarnings('ignore')

# plt.switch_backend('agg')

DFCOLORS = [
    [0.94, 0.94, 0.70, 0.7],
    [0.90, 0.78, 0.18, 0.7],
    [0.92, 0.45, 0.03, 0.7],
    [0.75, 0.22, 0.36, 0.7],
    [0.36, 0.16, 0.70, 0.7],
    [0.12, 0.12, 0.39, 0.7]
]

DFBINS = [0.002, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]


def hazdev(maplayerlist, configs, shakemap, outfolder=None, alpha=0.7,
           shakethreshtype='pga', shakethresh=10.,
           prefLS='Nowicki Jessee and others (2018)',
           prefLQ='Zhu and others (2017)',
           pop_file=None, defaultcolors=True, point=True,
           pager_alert='', eventsource='', eventsourcecode='',
           createpngs=True, gf_version=1, pdlcall=False):
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
        shakethresh: Float or list of shaking thresholds in %g for pga, cm/s
            for pgv, float for mmi. Used for Hagg and Exposure computation.
        prefLS (str): shortref of "preferred" landslide model.
        prefLQ (str): shortref of "preferred" liquefaction model.
        pop_file (str): file path to population file used to compute
            population-based alert levels.
        defaultcolors (bool): If True, will use DFCOLORS for all layers instead
            of determining new ones. This will crash if any of the layers have
            a different number of bins than the number of DFCOLORS
        point (bool): if True, event is a point source and warning should be
            displayed
        pager_alert (str): PAGER alert level, e.g., 'green'. 'pending', ...
        eventsource (str): net id (e.g., 'us')
        eventsourcecode (str): event code (e.g. '123456pq')
        createpngs (bool): if True, create pngs for web map
        gf_version (int): ground failure version
        pdlcall (bool): True if callgf was called by pdl automatically,
            false if called manually (or otherwise).

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
                probthresh = None
                id1 = 'godt_2008'
                maxP = 1.
            else:
                # Since logistic models can't equal one, need to eliminate
                # placeholder zeros before computing stats
                if 'jessee' in maplayer['model']['description']['name'].lower(
                ):
                    id1 = 'jessee_2018'
                    probthresh = 0.002
                    maxP = 0.26
                else:
                    id1 = 'nowicki_2014_global'
                    probthresh = 0.0
                    maxP = 1.

            if 'std' in list(maplayer.keys()):
                stdgrid2D = maplayer['std']['grid']
            else:
                stdgrid2D = None

            stats = computeStats(
                maplayer['model']['grid'],
                stdgrid2D=stdgrid2D,
                shakefile=shakemap,
                shakethresh=shakethresh,
                probthresh=probthresh,
                pop_file=pop_file,
                shakethreshtype=shakethreshtype,
                maxP=maxP)

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
                    stats['hagg_0.10g'], 0.,
                    stats['exp_pop_0.10g'], 0.)
                lsext = get_zoomextent(maplayer['model']['grid'])
                ls_haz_value = set_num_precision(
                    stats['hagg_0.10g'], 2, 'float')
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

            ls_haz_std = None
            ls_pop_std = None
            ls_hlim = None
            ls_elim = None
            ls_hp = None
            ls_hq = None
            ls_ep = None
            ls_eq = None
            ls_haz_1std_range = None
            ls_haz_2std_range = None
            ls_pop_1std_range = None
            ls_pop_2std_range = None

            if stdgrid2D is not None and title == prefLS:
                ph = stats['p_hagg_0.10g']
                qh = stats['q_hagg_0.10g']
                pe = stats['p_exp_0.10g']
                qe = stats['q_exp_0.10g']
                hmax = stats['hlim_0.10g']
                emax = stats['elim_0.10g']

                ls_haz_std = float("%.4f" % stats['hagg_std_0.10g'])
                ls_pop_std = float("%.4f" % stats['exp_std_0.10g'])
                ls_hlim = float("%.4f" % hmax)
                ls_elim = float("%.4f" % emax)
                ls_hp = float("%.4f" % ph)
                ls_hq = float("%.4f" % qh)
                ls_ep = float("%.4f" % pe)
                ls_eq = float("%.4f" % qe)

                # Add bar uncertainty extents here using p and q if applicable
                # make sure not a non-event/placeholder
                if ph > 0. and qh > 0.:
                    h68 = get_rangebeta(ph, qh, prob=0.6827, minlim=0.,
                                        maxlim=hmax)
                    h95 = get_rangebeta(ph, qh, prob=0.9545, minlim=0.,
                                        maxlim=hmax)
                    ls_haz_1std_range = h68
                    ls_haz_2std_range = h95
                # make sure not a non-event/placeholder
                if pe > 0. and qe > 0.:
                    e68 = get_rangebeta(pe, qe, prob=0.6827, minlim=0.,
                                        maxlim=emax)
                    e95 = get_rangebeta(pe, qe, prob=0.9545, minlim=0.,
                                        maxlim=emax)
                    ls_pop_1std_range = e68
                    ls_pop_2std_range = e95

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
                    'std': ls_haz_std,
                    'parameter': 'Aggregate Hazard',
                    'units': 'km^2'
                },
                'population_alert': {
                    'color': ls_pop_alert,
                    'value': ls_pop_value,
                    'std': ls_pop_std,
                    'parameter': 'Population exposure',
                    'units': 'people'
                },
                'bin_edges': list(lims[0]),
                'probability': {
                    'max': float("%.2f" % stats['Max']),
                    'std': float("%.2f" % stats['Std']),
                    'hagg0.1g': float("%.2f" % stats['hagg_0.10g']),
                    'popexp0.1g': float("%.2f" % stats['exp_pop_0.10g'],),
                    'hagg0.1g_std': ls_haz_std,
                    'popexp0.1g_std': ls_pop_std,
                    'hlim0.1g': ls_hlim,
                    'elim0.1g': ls_elim,
                    'p_hagg': ls_hp,
                    'q_hagg': ls_hq,
                    'p_exp': ls_ep,
                    'q_exp': ls_eq,
                    'hagg_1std': ls_haz_1std_range,
                    'hagg_2std': ls_haz_2std_range,
                    'pop_1std': ls_pop_1std_range,
                    'pop_2std': ls_pop_2std_range,
                },
                'longref': metadata['longref'],
                'parameters': metadata['parameters'],
                'zoomext': lsext
            }
            # Replace any Nans with Nones so info.json is valid
            for key in edict['probability']:
                if np.isscalar(edict['probability'][key]):
                    if np.isnan(edict['probability'][key]):
                        edict['probability'][key] = None
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
                probthresh = 0.0
                maxP = 1.
            elif '2017' in maplayer['model']['description']['name'].lower():
                id1 = 'zhu_2017_general'
                probthresh = 0.005
                maxP = 0.487
            else:
                probthresh = None
                maxP = 1.
                id1 = None

            if 'std' in list(maplayer.keys()):
                stdgrid2D = maplayer['std']['grid']
            else:
                stdgrid2D = None

            stats = computeStats(
                maplayer['model']['grid'],
                stdgrid2D=stdgrid2D,
                shakefile=shakemap,
                shakethresh=shakethresh,
                pop_file=pop_file,
                shakethreshtype=shakethreshtype,
                probthresh=probthresh,
                maxP=maxP)
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
                    0., stats['hagg_0.10g'],
                    0., stats['exp_pop_0.10g'])
                lqext = get_zoomextent(maplayer['model']['grid'])
                lq_haz_value = set_num_precision(
                    stats['hagg_0.10g'], 2, 'float')
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

            lq_haz_std = None
            lq_pop_std = None
            lq_hlim = None
            lq_elim = None
            lq_hp = None
            lq_hq = None
            lq_ep = None
            lq_eq = None
            lq_haz_1std_range = None
            lq_haz_2std_range = None
            lq_pop_1std_range = None
            lq_pop_2std_range = None

            if stdgrid2D is not None and title == prefLQ:
                ph = stats['p_hagg_0.10g']
                qh = stats['q_hagg_0.10g']
                pe = stats['p_exp_0.10g']
                qe = stats['q_exp_0.10g']
                hmax = stats['hlim_0.10g']
                emax = stats['elim_0.10g']

                lq_haz_std = float("%.2f" % stats['hagg_std_0.10g'])
                lq_pop_std = float("%.2f" % stats['exp_std_0.10g'])
                lq_hlim = float("%.4f" % hmax)
                lq_elim = float("%.4f" % emax)
                lq_hp = float("%.4f" % ph)
                lq_hq = float("%.4f" % qh)
                lq_ep = float("%.4f" % pe)
                lq_eq = float("%.4f" % qe)

                # Add bar uncertainty extents here using p and q if applicable
                # make sure not a non-event/placeholder
                if ph > 0. and qh > 0.:
                    h68 = get_rangebeta(ph, qh, prob=0.6827, minlim=0.,
                                        maxlim=hmax)
                    h95 = get_rangebeta(ph, qh, prob=0.9545, minlim=0.,
                                        maxlim=hmax)
                    lq_haz_1std_range = h68
                    lq_haz_2std_range = h95
                # make sure not a non-event/placeholder
                if pe > 0. and qe > 0.:
                    e68 = get_rangebeta(pe, qe, prob=0.6827, minlim=0.,
                                        maxlim=emax)
                    e95 = get_rangebeta(pe, qe, prob=0.9545, minlim=0.,
                                        maxlim=emax)
                    lq_pop_1std_range = e68
                    lq_pop_2std_range = e95

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
                    'std': lq_haz_std,
                    'parameter': 'Aggregate Hazard',
                    'units': 'km^2'
                },
                'population_alert': {
                    'color': lq_pop_alert,
                    'value': lq_pop_value,
                    'std': lq_pop_std,
                    'parameter': 'Population exposure',
                    'units': 'people'
                },
                'bin_edges': list(lims[0]),
                'probability': {
                    'max': float("%.2f" % stats['Max']),
                    'std': float("%.2f" % stats['Std']),
                    'hagg0.1g': float("%.2f" % stats['hagg_0.10g']),
                    'popexp0.1g': float("%.2f" % stats['exp_pop_0.10g']),
                    'hagg0.1g_std': lq_haz_std,
                    'popexp0.1g_std': lq_pop_std,
                    'hlim0.1g': lq_hlim,
                    'elim0.1g': lq_elim,
                    'p_hagg': lq_hp,
                    'q_hagg': lq_hq,
                    'p_exp': lq_ep,
                    'q_exp': lq_eq,
                    'hagg_1std': lq_haz_1std_range,
                    'hagg_2std': lq_haz_2std_range,
                    'pop_1std': lq_pop_1std_range,
                    'pop_2std': lq_pop_2std_range,

                },
                'longref': metadata['longref'],
                'parameters': metadata['parameters'],
                'zoomext': lqext
            }
            # Replace any Nans with Nones so info.json is valid
            for key in edict['probability']:
                if np.isscalar(edict['probability'][key]):
                    if np.isnan(edict['probability'][key]):
                        edict['probability'][key] = None
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
        sync, colorlistLS, reflims = setupcolors(
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

        sync, colorlistLQ, reflims = setupcolors(
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
    if createpngs:
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
    infojson = create_info(outfolder, lsmodels, lqmodels, gf_version,
                           eventsource, eventsourcecode, point, pdlcall)
    filenames.append(infojson)

    return filenames


def create_png(event_dir, lsmodels=None, lqmodels=None, mercator=True,
               lsmask=0.002, lqmask=0.005, legends=False):
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
        ls_mod_file = [f for f in files if 'jessee_2018.hdf5' in f]
        if len(ls_mod_file) == 1:
            ls_file = os.path.join(event_dir, ls_mod_file[0])
            ls_mod = loadlayers(ls_file)
            filesnippet = 'jessee_2018'
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
            # if lqm['preferred']:
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


def create_info(event_dir, lsmodels, lqmodels, gf_version=1,
                eventsource='', eventsourcecode='', point=True,
                pdlcall=False):
    """Create info.json for ground failure product.

    Args:
        event_dir (srt): Directory containing ground failure results.
        lsmodels (list): List of dictionaries of model summary info compiled
            by the hazdev function. If not specified, code will search for
            the hdf5 files for the preferred model and will create this
            dictionary and will apply default colorbars and bins.
        lqmodels (list): Same as above for liquefaction.
        gf_version (int): ground failure version
        eventsource (str): net id (e.g., 'us')
        eventsourcecode (str): event code (e.g. '123456pq')
        point (bool): if True, event is a point source and warning should be
            displayed
        pdlcall (bool): True if callgf was called by pdl automatically

    Returns:
        creates info.json for this event
    """
    filenames = []
    # Find the shakemap grid.xml file
    with open(os.path.join(event_dir, 'shakefile.txt'), 'r') as f:
        shakefile = f.read()

    files = os.listdir(event_dir)

    # Get all info from dictionaries of preferred events, add in extent
    # and filename
    ls_alert = None
    lq_alert = None
    lsext = None
    lqext = None
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
            # Deal with extent
            try:
                with open(os.path.join(event_dir, lsm['extent'])) as ff:
                    extent1 = json.load(ff)
            except BaseException:
                extent1 = None
            lsm['extent'] = extent1
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
            # Deal with extent
            try:
                with open(os.path.join(event_dir, lqm['extent'])) as ff:
                    extent1 = json.load(ff)
            except BaseException:
                extent1 = None
            lqm['extent'] = extent1
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

    # Figure out if realtime (pdl run and same week as event)
    if (pdlcall and datetime.utcnow() - event_dict['event_timestamp'] <
            timedelta(days=7)):
        realtime = True
    else:
        realtime = False

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
            'location': event_dict['event_description'],
            'event_url': event_url,
            'gf_version': gf_version,
            'gf_time': datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ'),
            'shakemap_version': sm_dict['shakemap_version'],
            'shakemap_source': sm_dict['shakemap_originator'],
            'shakemap_time': sm_dict['process_timestamp'].strftime('%Y-%m-%dT%H:%M:%SZ'),
            'rupture_warning': rupture_warning,
            'point_source': point,
            'zoom_extent': [xmin, xmax, ymin, ymax],
            'realtime': realtime
        },
        'Landslides': lsmodels,
        'Liquefaction': lqmodels
    }

    info_file = os.path.join(event_dir, 'info.json')
    with open(info_file, 'w') as f:
        json.dump(info_dict, f)
    fixfile(info_file)
    filenames.append(info_file)
    return filenames


def fixfile(filename):
    """Replace any Nan or Infinity values with null in info.json

    Args:
        filename (str): full path to info.json file

    Returns: same file with invalid entries replaced with valid

    """
    if os.path.exists(filename):
        with open(filename, 'r') as file:
            filedata = file.read()
        if 'NaN' in filedata:
            filedata = filedata.replace('NaN', 'null')
        if 'Infinity' in filedata:
            filedata = filedata.replace('-Infinity', 'null')
            filedata = filedata.replace('Infinity', 'null')
        with open(filename, 'w') as file:
            file.write(filedata)


def make_legends(lqmin=0.005, lsmin=0.002, outfolder=None,
                 orientation='horizontal', transparent=True):
    """Make png file of legends to go with pngs using DFCOLORS and DFBINS

    Args:
        lqmin (float): minimum visible value of liquefaction probability
        lsmin (float): same as above for landslides
        outfolder (str): folder to place pngs of legends
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
                transparent=transparent)

    # --------------------------
    # Make landslide legend
    binedges2 = DFBINS
    binedges2[0] = lsmin
    if outfolder is None:
        lsfilename = 'legend_landslide.png'
    else:
        lsfilename = os.path.join(outfolder, 'legend_landslide.png')

    make_legend(DFBINS, DFCOLORS, filename=lqfilename,
                orientation=orientation, title='Landslide probability',
                transparent=transparent)

    return lsfilename, lqfilename


def create_kmz(maplayer, outfile, mask=None, levels=None,
               colorlist=None, qdict=None):
    """
    Create kmz files of models

    Args:
        maplayer (dict):
            Dictionary of one model result formatted like:

            .. code-block:: python

                {
                    'grid': mapio grid2D object,
                    'label': 'label for colorbar and top line of subtitle',
                    'type': 'output or input to model',
                    'description': 'description for subtitle'
                }
        outfile (str):
            File extension
        mask (float):
            make all cells below this value transparent
        levels (array):
            list of bin edges for each color, must be same length
        colorlist (array):
            list of colors for each bin, should be length one less than levels
        qdict (dict):
            dictionary of quantile grids

    Returns:
        kmz file
    """
    # Figure out lims
    if levels is None:
        levels = DFBINS
    if colorlist is None:
        colorlist = DFCOLORS

    if len(levels) - 1 != len(colorlist):
        raise Exception('len(levels) must be one longer than len(colorlist)')

    # Make place to put temporary files
    temploc = tempfile.TemporaryDirectory()

    # Figure out file names
    name, ext = os.path.splitext(outfile)
    basename = os.path.basename(name)
    if ext != '.kmz':
        ext = '.kmz'
    filename = '%s%s' % (name, ext)
    mapfile = os.path.join(temploc.name, '%s.tiff' % basename)
    legshort = '%s_legend.png' % basename
    legfile = os.path.join(temploc.name, legshort)

    # Make colored geotiff
    out = make_rgba(maplayer['grid'], mask=mask,
                    levels=levels, colorlist=colorlist)
    rgba_img, extent, lmin, lmax, cmap = out
    # Save as a tiff
    plt.imsave(mapfile, rgba_img, vmin=lmin, vmax=lmax, cmap=cmap)
    if qdict is not None:
        for quant, quantdict in qdict.items():
            out = make_rgba(quantdict['grid'], mask=mask,
                            levels=levels, colorlist=colorlist)
            rgba_img, extent, lmin, lmax, cmap = out
            qfile = os.path.join(
                temploc.name, '%s_%s.tiff' % (basename, quant))
            plt.imsave(qfile, rgba_img, vmin=lmin, vmax=lmax, cmap=cmap)

    # Start creating kmz
    L = simplekml.Kml()

    # Set zoom window
    doc = L.document  # have to put lookat in root document directory
    doc.altitudemode = simplekml.AltitudeMode.relativetoground
    boundaries1 = get_zoomextent(maplayer['grid'])
    doc.lookat.latitude = np.mean([boundaries1['ymin'], boundaries1['ymax']])
    doc.lookat.longitude = np.mean([boundaries1['xmax'], boundaries1['xmin']])
    doc.lookat.altitude = 0.
    # dist in m from point
    doc.lookat.range = (boundaries1['ymax'] -
                        boundaries1['ymin']) * 111. * 1000.
    doc.description = (
        'USGS near-real-time earthquake-triggered %s model for event id %s'
        % (maplayer['description']['parameters']['modeltype'],
           maplayer['description']['event_id'])
    )

    prob = L.newgroundoverlay(name=maplayer['label'])
    prob.icon.href = 'files/%s.tiff' % basename
    prob.latlonbox.north = extent[3]
    prob.latlonbox.south = extent[2]
    prob.latlonbox.east = extent[1]
    prob.latlonbox.west = extent[0]
    L.addfile(mapfile)
    if qdict is not None:
        for quant, quantdict in qdict.items():
            qfile = os.path.join(
                temploc.name, '%s_%s.tiff' % (basename, quant))
            qlayer = L.newgroundoverlay(
                name=quantdict['label'], visibility=0)
            qlayer.icon.href = 'files/%s_%s.tiff' % (basename, quant)
            qlayer.latlonbox.north = extent[3]
            qlayer.latlonbox.south = extent[2]
            qlayer.latlonbox.east = extent[1]
            qlayer.latlonbox.west = extent[0]
            L.addfile(qfile)

    # Add legend and USGS icon as screen overlays
    # Make legend
    make_legend(levels, colorlist, filename=legfile, title=maplayer['label'])

    size1 = simplekml.Size(x=0.3, xunits=simplekml.Units.fraction)
    leg = L.newscreenoverlay(name='Legend', size=size1)
    leg.icon.href = 'files/%s' % legshort
    leg.screenxy = simplekml.ScreenXY(x=0.2, y=0.05,
                                      xunits=simplekml.Units.fraction,
                                      yunits=simplekml.Units.fraction)
    L.addfile(legfile)

    size2 = simplekml.Size(x=0.15, xunits=simplekml.Units.fraction)
    icon = L.newscreenoverlay(name='USGS', size=size2)
    icon.icon.href = 'files/USGS_ID_white.png'
    icon.screenxy = simplekml.ScreenXY(x=0.8, y=0.95,
                                       xunits=simplekml.Units.fraction,
                                       yunits=simplekml.Units.fraction)
    L.addfile(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           os.pardir, 'content', 'USGS_ID_white.png'))

    L.savekmz(filename)
    return filename


# noinspection PyArgumentList
def get_zoomextent(grid, propofmax=0.3):
    """
    Get the extent that contains all values with probabilities exceeding
    a threshold in order to determine ideal zoom level for interactive map
    If nothing is above the threshold, uses the full extent

    Args:
        grid: grid2d of model output
        propofmax (float): Proportion of maximum that should be fully included
            within the bounds.

    Returns:
        * boundaries: a dictionary with keys 'xmin', 'xmax', 'ymin', and
         'ymax' that defines the zoomed boundaries in geographic coordinates.

    """
    maximum = np.nanmax(grid.getData())

    xmin, xmax, ymin, ymax = grid.getBounds()

    if np.isnan(maximum):
        # If no finite values, use entire extent for zoom
        return dict(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

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


def make_rgba(grid2D, levels, colorlist, mask=None,
              mercator=False):
    """
    Make an rgba (red, green, blue, alpha) grid out of raw data values and
    provide extent and limits needed to save as an image file

    Args:
        grid2D: Mapio Grid2D object of result to mape
        levels (list): list of bin edges for each color, must be same length
        colorlist (list): list of colors for each bin, should be length one
            less than levels
        mask (float): mask all values below this value
        mercator (bool): project to web mercator (needed for leaflet, not
                 for kmz)

    Returns:
        tuple: (rgba_img, extent, lmin, lmax, cmap), where:
            * rgba_img: rgba (red green blue alpha) image
            * extent: list of outside corners of image,
                [minlat, maxlat, minlon, maxlon]
            * lmin: lowest bin edge
            * lmax: highest bin edge
            * cmap: colormap corresponding to image
    """

    data1 = grid2D.getData()
    if mask is not None:
        data1[data1 < mask] = float('nan')
    geodict = grid2D.getGeoDict()
    extent = [
        geodict.xmin - 0.5 * geodict.dx,
        geodict.xmax + 0.5 * geodict.dx,
        geodict.ymin - 0.5 * geodict.dy,
        geodict.ymax + 0.5 * geodict.dy,
    ]

    lmin = levels[0]
    lmax = levels[-1]
    data2 = np.clip(data1, lmin, lmax)
    cmap = mpl.colors.ListedColormap(colorlist)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    data2 = np.ma.array(data2, mask=np.isnan(data1))
    rgba_img = cmap(norm(data2))
    if mercator:
        rgba_img = mercator_transform(
            rgba_img, (extent[2], extent[3]), origin='upper')

    return rgba_img, extent, lmin, lmax, cmap


def make_legend(levels, colorlist, filename=None, orientation='horizontal',
                title=None, transparent=False):
    """Make legend file

    Args:

        levels (list): list of bin edges for each color, must be same length
        colorlist (list): list of colors for each bin, should be length one
            less than levels
        filename (str): File extension of legend file
        orientation (str): orientation of colorbar, 'horizontal' or 'vertical'
        title (str): title of legend (usually units)
        transparent (bool): if True, background will be transparent

    Returns:
        figure of legend

    """
    fontsize = 16
    labels = ['< %1.1f%%' % (levels[0] * 100.,)]
    for db in levels:
        if db < 0.01:
            labels.append('%1.1f' % (db * 100,))
        else:
            labels.append('%1.0f' % (db * 100.,))

    if orientation == 'vertical':
        # Flip order to darker on top
        labels = labels[::-1]
        colors1 = colorlist[::-1]
        fig, axes = plt.subplots(len(colors1) + 1, 1,
                                 figsize=(3., len(colors1) - 1.7))
        clearind = len(axes) - 1
        maxind = 0
    else:
        colors1 = colorlist
        fig, axes = plt.subplots(1, len(colorlist) + 1,
                                 figsize=(len(colorlist) + 1.7, 0.8))
        # DPI = fig.get_dpi()
        # fig.set_size_inches(440/DPI, 83/DPI)
        clearind = 0
        maxind = len(axes) - 1

    for i, ax in enumerate(axes):
        ax.set_ylim((0., 1.))
        ax.set_xlim((0., 1.))
        # draw square
        if i == clearind:
            color1 = colors1[0]
            color1[-1] = 0.  # make completely transparent
            if orientation == 'vertical':
                label = labels[i + 1]
            else:
                label = labels[0]
        else:
            if orientation == 'vertical':
                label = '%s-%s%%' % (labels[i + 1], labels[i])
                color1 = colors1[i]
            else:
                label = '%s-%s%%' % (labels[i], labels[i + 1])
                color1 = colors1[i - 1]
            color1[-1] = 0.8  # make less transparent
            if i == maxind:
                label = '> %1.0f%%' % (levels[-2] * 100.)
        ax.set_facecolor(color1)
        if orientation == 'vertical':
            ax.text(1.1, 0.5, label, fontsize=fontsize,
                    rotation='horizontal', va='center')
        else:
            ax.set_xlabel(label, fontsize=fontsize,
                          rotation='horizontal')

        ax.set_yticks([])
        ax.set_xticks([])
        plt.setp(ax.get_yticklabels(), visible=False)
        plt.setp(ax.get_xticklabels(), visible=False)

    if orientation == 'vertical':
        fig.suptitle(title.title(), weight='bold', fontsize=fontsize + 2)
        plt.subplots_adjust(hspace=0.01, right=0.4, top=0.82)
    else:
        fig.suptitle(title.title(), weight='bold', fontsize=fontsize + 2)
        # , left=0.01, right=0.99, top=0.99, bottom=0.01)
        plt.subplots_adjust(wspace=0.1, top=0.6)
    # plt.tight_layout()
    if filename is not None:
        fig.savefig(filename, bbox_inches='tight', transparent=transparent)
    else:
        plt.show()


def setupcolors(sync, plotorder, lims, colormaps, defaultcolormap=cm.CMRmap_r,
                logscale=None, alpha=None):
    """Get colors that will be used for all colorbars from reference grid

    Args:
        sync(str): If False, will exit program, else corresponds to the
            shortref of the model which should serve as the template for
            the colorbars used by all other models. All other models must
            have the exact same number of bins
        plotorder (list): List of keys of shortrefs of the grids that will be
            plotted.
        lims (*): Nx1 list of tuples or numpy arrays corresponding to
            plotorder defining the bin edges to use for each model.
            Example:
            .. code-block:: python
                [(0., 0.1, 0.2, 0.3), np.linspace(0., 1.5, 15)]
        colormaps (list): List of strings of matplotlib colormaps (e.g.
            cm.autumn_r) corresponding to plotorder
        defaultcolormap (matplotlib colormap): Colormap to use if
            colormaps is not defined. default cm.CMRmap_r
        logscale (*): If not None, then a list of booleans corresponding to
            plotorder stating whether to use log scaling in determining colors
        alpha (*): list of transparencies

    Returns:
        tuple: (sync, colorlist, lim1) where:
            * sync (bool): whether or not colorbars are/can be synced
            * colorlist (list): list of rgba colors that will be applied to all
                models regardless of bin edge values
            * lim1 (array): bin edges of model to which others are synced

    """

    if not sync:
        return False, None, None

    elif sync in plotorder:
        k = [indx for indx, key in enumerate(plotorder) if key in sync][0]
        # Make sure lims exist and all have the same number of bins'

        if logscale is not None:
            logs = logscale[k]
        else:
            logs = False

        lim1 = np.array(lims[k])
        sum1 = 0
        for lim in lims:
            if lim is None:
                sum1 += 1
                continue
            if len(lim) != len(lim1):
                sum1 += 1
                continue
        if sum1 > 0:
            print('Cannot sync colorbars, different number of bins or lims '
                  'not specified')
            sync = False
            return sync, None, None

        if colormaps[k] is not None:
            palette1 = colormaps[k]
        else:
            palette1 = defaultcolormap
        #palette1.set_bad(clear_color, alpha=0.0)
        if logs:
            cNorm = colors.LogNorm(vmin=lim1[0], vmax=lim1[-1])
            # geometric mean for midpoints
            midpts = np.sqrt(lim1[1:] * lim1[:-1])
        else:
            cNorm = colors.Normalize(vmin=lim1[0], vmax=lim1[-1])
            midpts = (lim1[1:] - lim1[:-1]) / 2 + lim1[:-1]
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=palette1)
        colorlist = []
        for value in midpts:
            colorlist.append(scalarMap.to_rgba(value, alpha=alpha))
        sync = True

    else:
        print('Cannot sync colorbars, different number of bins or lims not '
              'specified')
        sync = False
        colorlist = None
        lim1 = None
    return sync, colorlist, lim1


if __name__ == '__main__':
    pass
