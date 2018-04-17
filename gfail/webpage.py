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
from mapio.shake import ShakeGrid
import matplotlib.cm as cm

from impactutils.textformat.text import set_num_precision
from impactutils.time.ancient_time import HistoricTime as ShakeDateTime
import pytz

from gfail.utilities import get_event_comcat, loadlayers


# temporary until mapio is updated
import warnings
warnings.filterwarnings('ignore')

plt.switch_backend('agg')

# DFCOLORS = [[0.9403921568627451, 0.9403921568627451, 0.7019607843137254, 0.7],
#            [0.9, 0.781764705882353, 0.18470588235294128, 0.7],
#            [0.92, 0.45, 0.03, 0.7],
#            [0.7552941176470588, 0.21941176470588236, 0.36411764705882355, 0.7],
#            [0.35882352941176465, 0.15980392156862744, 0.7009803921568627, 0.7],
#            [0.11764705882352941, 0.11764705882352941, 0.39215686274509803, 0.7]]

# hex versions:
DFCOLORS = [
    '#efefb3',
    '#e5c72f',
    '#ea7207',
    '#c0375c',
    '#5b28b2',
    '#1e1e64'
]

DFBINS = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]

# reformat for info.json
color_bins = []
for i in range(len(DFCOLORS)):
    color_bins.append({
        'min': DFBINS[i],
        'max': DFBINS[i+1],
        'color': DFCOLORS[i]
    })


def hazdev(maplayerlist, configs, shakemap, outfolder=None, alpha=0.7,
           shakethreshtype='pga', probthresh=None, shakethresh=10.,
           prefLS='Nowicki Jessee (2017)', prefLQ='Zhu and others (2017)',
           pop_file=None):
    """Create all files needed for product page creation
    Assumes gfail has been run already with -w flag

    Args:
        maplayerlist (list): List of model outputs from gfail
        configs (list): List of dictionaries of config files corresponding to
            each model in maplayerlist and in the same order
        shakemap (str): path to shakemap .xml file
        outfolder (str): Location in which to save outputs. If None, will use current directory 
        alpha (float): Transparency to use for overlay pngs, value from 0 to 1.
        shakethreshtype (str): Type of ground motion to use for shakethresh, 'pga', 'pgv', or 'mmi'.
        probthresh: Optional. Float or list of probability thresholds to apply before computing stats.
        shakethresh: Float or list of shaking thresholds in %g for pga, cm/s for pgv, float for mmi.
            Used for Hagg and Exposure computation
        prefLS (str): shortref of "preferred" landslide model
        prefLQ (str): shortref of "preferred" liquefaction model
        pop_filt (str): file path to population file used to compute population-based
            alert levels.

    Returns:
        Files that need to be sent to comcat for hazdev to create the product webpage including:
            info.json
            transparent png overlays of all models
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
            else:
                # Since logistic models can't equal one, need to eliminate
                # placeholder zeros before computing stats
                statprobthresh = 0.0

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
                ls_haz_alert, ls_pop_alert, _, _ = get_alert(
                    stats['Hagg_0.10g'], 0.,
                    stats['exp_pop_0.10g'], 0.)
            else:
                on = False
                ls_haz_alert = None
                ls_pop_alert = None

            edict = dict(model=metadata['name'],  preferred=on,
                         probability_max=float("%.2f" % stats['Max']),
                         probability_std=float("%.2f" % stats['Std']),
                         units=metadata['units'], bin_edges=list(lims[0]),
                         bin_colors=[], filename=[], extent=[],
                         hazard_alert_value=stats['Hagg_0.10g'],
                         populate_alert_value=stats['exp_pop_0.10g'],
                         hazard_alert_parameter='Hagg_0.10g',
                         population_alert_parameter='exp_pop_0.10g',
                         parameters=metadata['parameters'],
                         longref=metadata['longref'],
                         hazard_alert=ls_haz_alert,
                         population_alert=ls_pop_alert,)

            lsmodels.append(edict)

        elif 'liquefaction' in mdict['parameters']['modeltype'].lower():
            title = maplayer['model']['description']['name']
            plotorder, logscale, lims, colormaps, maskthreshes = \
                parseConfigLayers(maplayer, conf, keys=['model'])
            logLQ.append(logscale[0])
            limLQ.append(lims[0])
            colLQ.append(colormaps[0])
            concLQ.append(title)

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
                _, _, lq_haz_alert, lq_pop_alert = get_alert(
                    0., stats['Hagg_0.10g'],
                    0., stats['exp_pop_0.10g'])
            else:
                on = False
                lq_haz_alert = None
                lq_pop_alert = None

            edict = dict(model=metadata['name'],  preferred=on,
                         probability_max=float("%.2f" % stats['Max']),
                         probability_std=float("%.2f" % stats['Std']),
                         units=metadata['units'], bin_edges=list(lims[0]),
                         bin_colors=[], filename=[], extent=[],
                         hazard_alert_value=stats['Hagg_0.10g'],
                         populate_alert_value=stats['exp_pop_0.10g'],
                         hazard_alert_parameter='Hagg_0.10g',
                         population_alert_parameter='exp_pop_0.10g',
                         parameters=metadata['parameters'],
                         longref=metadata['longref'],
                         hazard_alert=lq_haz_alert,
                         population_alert=lq_pop_alert,)

            lqmodels.append(edict)

        else:
            raise Exception("model type is undefined, check "
                            "maplayer['model']['parameters']"
                            "['modeltype'] to ensure it is defined")
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

    # Create invisible pngs
    out = create_png(outfolder, lsmodels=lsmodels, lqmodels=lqmodels)
    filenames += out

    # Create info.json
    out = create_info(outfolder, lsmodels=lsmodels, lqmodels=lqmodels)
    filenames += out

    return filenames


def create_png(event_dir, lsmodels=None, lqmodels=None):
    """
    Creates transparent PNG file for website.

    Args:
        event_dir (srt): Directory containing ground failure results.
        lsmodels (list): List of dictionaries of model summary info compiled
            by the hazdev function. If not specified, code will search for
            the hdf5 files for the preferred model and will create this dictionary
            and will apply default colorbars and bins.
        lqmodels (list): Same as above for liquefaction.

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
        else:
            raise OSError("Preferred landslide model result not found.")
    else:
        for lsm in lsmodels:
            if lsm['preferred']:
                levels = lsm['bin_edges']
                colors1 = lsm['bin_colors']
                shortname = lsm['model']
                if 'Jessee (2017)' in shortname:
                    filesnippet = 'jessee_2017'
                elif 'Nowicki and others (2014)' in shortname:
                    filesnippet = 'nowicki_2014_global'
                elif 'godt' in shortname.lower():
                    filesnippet = 'godt_2008'
        fsh = '%s.hdf5' % filesnippet
        ls_mod_file = [f for f in files if fsh in f]
        if len(ls_mod_file) == 1:
            ls_file = os.path.join(event_dir, ls_mod_file[0])
            ls_mod = loadlayers(ls_file)
        else:
            raise OSError("Preferred landslide model result not found.")

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
    rgba_img = cmap(norm(ls_data2))
    filen = os.path.join(event_dir, '%s.png' % filesnippet)
    plt.imsave(filen,
               rgba_img,
               vmin=lmin,
               vmax=lmax,
               cmap=cmap
               )

    if lqmodels is None:
        # read in preferred model for liquefaction if none specified
        lq_mod_file = [f for f in files if 'zhu_2017_general.hdf5' in f]
        if len(lq_mod_file) == 1:
            lq_file = os.path.join(event_dir, lq_mod_file[0])
            lq_mod = loadlayers(lq_file)
            levels = DFBINS
            colors1 = DFCOLORS
            filesnippet = 'zhu_2017_general'
        else:
            raise OSError("Preferred liquefaction model result not found.")
    else:
        for lqm in lqmodels:
            if lqm['preferred']:
                levels = lqm['bin_edges']
                colors1 = lqm['bin_colors']
                shortname = lqm['model']
                if 'Zhu and others (2017)' in shortname:
                    filesnippet = 'zhu_2017_general'
                elif 'Zhu and others (2015)' in shortname:
                    filesnippet = 'zhu_2015'

        fsh = '%s.hdf5' % filesnippet
        lq_mod_file = [f for f in files if fsh in f]
        if len(lq_mod_file) == 1:
            lq_file = os.path.join(event_dir, lq_mod_file[0])
            lq_mod = loadlayers(lq_file)
        else:
            raise OSError("Preferred liquefaction model result not found.")

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
            the hdf5 files for the preferred model and will create this dictionary
            and will apply default colorbars and bins.
        lqmodels (list): Same as above for liquefaction.

    Returns:
        creates info.json for this event
    """
    filenames = []
    # Find the shakemap grid.xml file
    with open(os.path.join(event_dir, 'shakefile.txt'), 'r') as f:
        shakefile = f.read()

    files = os.listdir(event_dir)

    if lsmodels is None or lqmodels is None:

        # Read in the "preferred" model for landslides and liquefaction
        ls_mod_file = [f for f in files if 'jessee_2017.hdf5' in f]
        if len(ls_mod_file) == 1:
            ls_file = os.path.join(event_dir, ls_mod_file[0])
            ls_mod = loadlayers(ls_file)
        else:
            raise OSError("Preferred landslide model result not found.")
        lq_mod_file = [f for f in files if 'zhu_2017_general.hdf5' in f]
        if len(lq_mod_file) == 1:
            lq_file = os.path.join(event_dir, lq_mod_file[0])
            lq_mod = loadlayers(lq_file)
        else:
            raise OSError("Preferred liquefaction model result not found.")

        # Read in extents
        ls_extent_file = [f for f in files if 'jessee_2017_extent.json' in f]
        if len(ls_extent_file) == 1:
            ls_file = os.path.join(event_dir, ls_extent_file[0])
            with open(ls_file) as f:
                jessee_extent = json.load(f)
        else:
            raise OSError("Landslide extent not found.")
        lq_extent_file = [
            f for f in files if 'zhu_2017_general_extent.json' in f]
        if len(lq_extent_file) == 1:
            lq_file = os.path.join(event_dir, lq_extent_file[0])
            with open(lq_file) as f:
                zhu_extent = json.load(f)
        else:
            raise OSError("Landslide extent not found.")

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
        ls_haz_alert, ls_pop_alert, lq_haz_alert, lq_pop_alert = alert_info

        if lsmodels is None:
            lsmodels = [{
                'id': 'nowicki_jessee_2017',
                'title': 'Nowicki Jessee (2017)',
                'overlay': 'jessee_2017.png',
                'extent': jessee_extent,
                'units': "Proportiona of area affected",
                'preferred': True,
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
                'color_bins': color_bins,
                'probability': {
                    'max': float("%.2f" % ls_stats['Max']),
                    'std': float("%.2f" % ls_stats['Std'])
                }
            }]
        if lqmodels is None:
            lqmodels = [{
                'id': 'zhu_2017',
                'title': 'Zhu and others (2017)',
                'overlay': 'zhu_2017.png',
                'extent': zhu_extent,
                'units': "Proportiona of area affected",
                'preferred': True,
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
                'color_bins': color_bins,
                'probability': {
                    'max': float("%.2f" % lq_stats['Max']),
                    'std': float("%.2f" % lq_stats['Std'])
                }
            }]
    else:
        # Get all info from dictionaries of preferred events, add in extent
        # and filename
        for lsm in lsmodels:
            # Add extent and filename for preferred model
            if lsm['preferred']:
                shortname = lsm['model']
                if 'Jessee (2017)' in shortname:
                    filesnippet = 'jessee_2017'
                elif 'Nowicki and others (2014)' in shortname:
                    filesnippet = 'nowicki_2014_global'
                elif 'godt' in shortname.lower():
                    filesnippet = 'godt_2008'
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
                lsm['filename'] = flnm

        for lqm in lqmodels:
            if lqm['preferred']:
                shortname = lqm['model']
                if 'Zhu and others (2017)' in shortname:
                    filesnippet = 'zhu_2017_general'
                elif 'Zhu and others (2015)' in shortname:
                    filesnippet = 'zhu_2015'
                # Read in extents
                flnm = '%s_extent.json' % filesnippet
                lq_extent_file = [f for f in files if flnm in f]
                if len(lq_extent_file) == 1:
                    lq_file = os.path.join(event_dir, lq_extent_file[0])
                    with open(lq_file) as f:
                        lq_extent = json.load(f)
                else:
                    raise OSError("Liquefaction extent not found.")
                lqm['extent'] = lq_extent
                lqm['filename'] = flnm

    # Try to get event info
    event_dict = ShakeGrid.load(shakefile, adjust='res').getEventDict()
    sm_dict = ShakeGrid.load(shakefile, adjust='res').getShakeDict()
    base_url = 'https://earthquake.usgs.gov/earthquakes/eventpage/'
    try:
        # Hopefully this will eventually be more reliable once we get the
        # comcat info directly from the shakemap grid, rather than rely on
        # magnitude/location/time association.
        shakemap_info, detail, temp = get_event_comcat(shakefile)
        event_url = detail.url
        code = detail['code']
        net = detail['net']
        utc = pytz.utc
        detail_time = ShakeDateTime.fromtimestamp(detail['time']/1000.0, utc)

        # ---------------------------------------------------------------------
        # Finite fault stuff:
        #    Other sections of the code do some relatively complicated stuff to
        #    try to sort out the finite fault. Here, I'm just simplifying it so
        #    that it checks comcat for a finite fault file.
        # ---------------------------------------------------------------------
        fault_file = shakemap_info['input']['event_information']['faultfiles']
        if len(fault_file) > 0:
            point = False
        else:
            point = True

    except:
        # Hopefully we can eventually remove this....
        event_url = '%s%s#executive' % (base_url, event_dict['event_id'])
        code = 'unknown'
        net = 'unknown'
        detail_time = -999
        point = False

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
            'time': detail_time.strftime('%Y-%m-%dT%H:%M:%SZ'),
            'lat': event_dict['lat'],
            'lon': event_dict['lon'],
            'event_url': event_url,
            'shakemap_version': sm_dict['shakemap_version'],
            'rupture_warning': rupture_warning,
            'point_source': point
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
              hazbinLQ=[10., 100., 1000.], popbinLQ=[100, 1000, 10000]):
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
            tuple: alertLS, popalertLS, alertLQ, popalertLQ, where
                * alertLS is the landslide hazard alert level (str)
                * popalertLS is the landslide population alert level (str)
                * alertLQ is the liquefaction hazard alert level (str)
                * popalertLQ is the liquefaction population alert level (str)

    """
    if paramalertLS is None:
        alertLS = None
    elif paramalertLS < hazbinLS[0]:
        alertLS = 'green'
    elif paramalertLS >= hazbinLS[0] and paramalertLS < hazbinLS[1]:
        alertLS = 'yellow'
    elif paramalertLS >= hazbinLS[1] and paramalertLS < hazbinLS[2]:
        alertLS = 'orange'
    elif paramalertLS > hazbinLS[2]:
        alertLS = 'red'
    else:
        alertLS = None

    if parampopLS is None:
        popalertLS = None
    elif parampopLS < popbinLS[0]:
        popalertLS = 'green'
    elif parampopLS >= popbinLS[0] and parampopLS < popbinLS[1]:
        popalertLS = 'yellow'
    elif parampopLS >= popbinLS[1] and parampopLS < popbinLS[2]:
        popalertLS = 'orange'
    elif parampopLS > popbinLS[2]:
        popalertLS = 'red'
    else:
        popalertLS = None

    if paramalertLQ is None:
        alertLQ = None
    elif paramalertLQ < hazbinLQ[0]:
        alertLQ = 'green'
    elif paramalertLQ >= hazbinLQ[0] and paramalertLQ < hazbinLQ[1]:
        alertLQ = 'yellow'
    elif paramalertLQ >= hazbinLQ[1] and paramalertLQ < hazbinLQ[2]:
        alertLQ = 'orange'
    elif paramalertLQ > hazbinLQ[2]:
        alertLQ = 'red'
    else:
        alertLQ = None

    if parampopLQ is None:
        popalertLQ = None
    elif parampopLQ < popbinLQ[0]:
        popalertLQ = 'green'
    elif parampopLQ >= popbinLQ[0] and parampopLQ < popbinLQ[1]:
        popalertLQ = 'yellow'
    elif parampopLQ >= popbinLQ[1] and parampopLQ < popbinLQ[2]:
        popalertLQ = 'orange'
    elif parampopLQ > popbinLQ[2]:
        popalertLQ = 'red'
    else:
        popalertLQ = None

    return alertLS, popalertLS, alertLQ, popalertLQ
