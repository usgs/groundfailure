# stdlib imports

import os
import numpy as np
import urllib
import json
from datetime import timedelta
import collections

# local imports
from mapio.shake import getHeaderData
from libcomcat.search import get_event_by_id, search
from libcomcat.classes import VersionOption
from mapio.basemapcity import BasemapCities
from mapio.multihaz import MultiHazardGrid

import matplotlib.cm as cm  # Don't delete this, it's needed in an eval function


def get_event_comcat(shakefile, timewindow=60, degwindow=0.3, magwindow=0.2):
    """
    """
    header_dicts = getHeaderData(shakefile)
    grid_dict = header_dicts[0]
    event_dict = header_dicts[1]
    version = grid_dict['shakemap_version']
    try:
        eid = event_dict['event_id']
        net = 'us'
        if 'event_network' in event_dict:
            net = event_dict['event_network']
        if not eid.startswith(net):
            eid = net + eid
        detail = get_event_by_id(eid, includesuperseded=True)
    except:
        lat = event_dict['lat']
        lon = event_dict['lon']
        mag = event_dict['magnitude']
        time = event_dict['event_timestamp']
        starttime = time - timedelta(seconds=timewindow)
        endtime = time + timedelta(seconds=timewindow)
        minlat = lat - degwindow
        minlon = lon - degwindow
        maxlat = lat + degwindow
        maxlon = lon + degwindow
        minmag = max(0, mag - magwindow)
        maxmag = min(10, mag + magwindow)
        events = search(starttime=starttime,
                        endtime=endtime,
                        minmagnitude=minmag,
                        maxmagnitude=maxmag,
                        minlatitude=minlat,
                        minlongitude=minlon,
                        maxlatitude=maxlat,
                        maxlongitude=maxlon)
        if not len(events):
            return None
        detail = events[0].getDetailEvent()
    allversions = detail.getProducts('shakemap', version=VersionOption.ALL)
    # Find the right version
    vers = [allv.version for allv in allversions]
    idx = np.where(np.array(vers) == version)[0][0]
    shakemap = allversions[idx]
    infobytes, url = shakemap.getContentBytes('info.json')
    info = json.loads(infobytes.decode('utf-8'))
    return info, detail, shakemap

def parseMapConfig(config, fileext=None):
    """
    Parse config for mapping options.

    Args:
        config (ConfigObj): ConfigObj object.
        fileext (str): File extension to add to relative filepaths, will be
            prepended to any file paths in config.

    Returns:
        dict: Dictionary of map options pulled from config file.
    """
    topofile = None
    roadfolder = None
    cityfile = None
    roadcolor = '6E6E6E'
    countrycolor = '177F10'
    watercolor = 'B8EEFF'
    ALPHA = 0.7
    oceanfile = None
    oceanref = None
    roadref = None
    cityref = None

    if fileext is None:
        fileext = '.'
    if 'dem' in config:
        topofile = os.path.join(fileext, config['dem']['file'])
        if os.path.exists(topofile) is False:
            print('DEM not valid - hillshade will not be possible\n')
    if 'ocean' in config:
        oceanfile = os.path.join(fileext, config['ocean']['file'])
        try:
            oceanref = config['ocean']['shortref']
        except:
            oceanref = 'unknown'
    if 'roads' in config:
        roadfolder = os.path.join(fileext, config['roads']['file'])
        if os.path.exists(roadfolder) is False:
            print('roadfolder not valid - roads will not be displayed\n')
            roadfolder = None
        try:
            roadref = config['roads']['shortref']
        except:
            roadref = 'unknown'
    if 'cities' in config:
        cityfile = os.path.join(fileext, config['cities']['file'])
        try:
            cityref = config['cities']['shortref']
        except:
            cityref = 'unknown'
        if os.path.exists(cityfile):
            try:
                BasemapCities.loadFromGeoNames(cityfile=cityfile)
            except Exception as e:
                print(e)
                print('cities file not valid - cities will not be displayed\n')
                cityfile = None
        else:
            print('cities file not valid - cities will not be displayed\n')
            cityfile = None
    if 'roadcolor' in config['colors']:
        roadcolor = config['colors']['roadcolor']
    if 'countrycolor' in config['colors']:
        countrycolor = config['colors']['countrycolor']
    if 'watercolor' in config['colors']:
        watercolor = config['colors']['watercolor']
    if 'alpha' in config['colors']:
        ALPHA = float(config['colors']['alpha'])

    countrycolor = '#'+countrycolor
    watercolor = '#'+watercolor
    roadcolor = '#'+roadcolor

    mapin = {'topofile': topofile, 'roadfolder': roadfolder,
             'cityfile': cityfile, 'roadcolor': roadcolor,
             'countrycolor': countrycolor, 'watercolor': watercolor,
             'ALPHA': ALPHA, 'roadref': roadref,
             'cityref': cityref, 'oceanfile': oceanfile, 'oceanref': oceanref}

    return mapin


def parseConfigLayers(maplayers, config, keys=None):
    """
    TODO:
        - Add ability to interpret custom color maps.

    Parse things that need to coodinate with each layer (like lims, logscale,
    colormaps etc.) from config file, in right order, where the order is from
    maplayers.

    Args:
        maplayers (dict): Dictionary containing model output.
        config (ConfigObj): Config object describing options for specific
            model.
        keys (list): List of keys of maplayers to process, e.g. ``['model']``.

    Returns:
        list: List of the following:
            * plotorder: maplayers keys in order of plotting.
            * logscale: list of logscale options from config corresponding to
              keys in plotorder (same order).
            * lims: list of colorbar limits from config corresponding to keys
              in plotorder (same order).
            * colormaps: list of colormaps from config corresponding to keys
              in plotorder (same order),
            * maskthreshes: list of mask thresholds from config corresponding
              to keys in plotorder (same order).

    """
    # get all key names, create a plotorder list in case maplayers is not an
    # ordered dict, making sure that anything called 'model' is first
    if keys is None:
        keys = list(maplayers.keys())
    plotorder = []

    try:
        limits = config[config.keys()[0]]['display_options']['lims']
        lims = []
    except:
        lims = None
        limits = None

    try:
        colors = config[config.keys()[0]]['display_options']['colors']
        colormaps = []
    except:
        colormaps = None
        colors = None

    try:
        logs = config[config.keys()[0]]['display_options']['logscale']
        logscale = []
    except:
        logscale = False
        logs = None

    try:
        masks = config[config.keys()[0]]['display_options']['maskthresholds']
        maskthreshes = []
    except:
        maskthreshes = None
        masks = None

    try:
        default = \
            config[config.keys()[0]]['display_options']['colors']['default']
        default = eval(default)
    except:
        default = None

    for i, key in enumerate(keys):
        plotorder += [key]
        if limits is not None:
            found = False
            for l in limits:
                getlim = None
                if l in key:
                    if type(limits[l]) is list:
                        getlim = np.array(limits[l]).astype(np.float)
                    else:
                        try:
                            getlim = eval(limits[l])
                        except:
                            getlim = None
                    lims.append(getlim)
                    found = True
            if not found:
                lims.append(None)

        if colors is not None:
            found = False
            for c in colors:
                if c in key:
                    getcol = colors[c]
                    colorobject = eval(getcol)
                    if colorobject is None:
                        colorobject = default
                    colormaps.append(colorobject)
                    found = True
            if not found:
                colormaps.append(default)

        if logs is not None:
            found = False
            for g in logs:
                getlog = False
                if g in key:
                    if logs[g].lower() == 'true':
                        getlog = True
                    logscale.append(getlog)
                    found = True
            if not found:
                logscale.append(False)

        if masks is not None:
            found = False
            for m in masks:
                if m in key:
                    getmask = eval(masks[m])
                    maskthreshes.append(getmask)
                    found = True
            if not found:
                maskthreshes.append(None)

    # Reorder everything so model is first, if it's not already
    if plotorder[0] != 'model':
        indx = [idx for idx, key in enumerate(plotorder) if key == 'model']
        if len(indx) == 1:
            indx = indx[0]
            firstpo = plotorder.pop(indx)
            plotorder = [firstpo] + plotorder
            firstlog = logscale.pop(indx)
            logscale = [firstlog] + logscale
            firstlim = lims.pop(indx)
            lims = [firstlim] + lims
            firstcol = colormaps.pop(indx)
            colormaps = [firstcol] + colormaps

    return plotorder, logscale, lims, colormaps, maskthreshes


def text_to_json(input1):
    """
    Simplification of text_to_json from shakelib.rupture.factory
    
    Args:
        input1 (str): url or filepath to text file
    """
    if os.path.exists(input1):
        with open(input1, 'r') as f:
            lines = f.readlines()
    else:
        with urllib.request.urlopen(input1) as f:
            lines = f.readlines()

    x = []
    y = []
    z = []
    reference = ''
    # convert to geojson
    for line in lines:
        sline = line.strip()
        if sline.startswith('#'):
            reference += sline.strip('#').strip('Source: ')
            continue
        if sline.startswith('>'):
            if len(x):  # start of new line segment
                x.append(np.nan)
                y.append(np.nan)
                z.append(np.nan)
                continue
            else:  # start of file
                continue
        if not len(sline.strip()):
            continue
        parts = sline.split()

        y.append(float(parts[0]))
        x.append(float(parts[1]))
        if len(parts) >= 3:
            z.append(float(parts[2]))
        else:
            print('Fault file has no depths, assuming zero depth')
            z.append(0.0)
        coords = []
        poly = []
        for lon, lat, dep in zip(x, y, z):
            if np.isnan(lon):
                coords.append(poly)
                poly = []
            else:
                poly.append([lon, lat, dep])
        if poly != []:
            coords.append(poly)

    d = {
        "type": "FeatureCollection",
        "metadata": {
            'reference': reference
        },
        "features": [
            {
                "type": "Feature",
                "properties": {
                    "rupture type": "rupture extent"
                },
                "geometry": {
                    "type": "MultiPolygon",
                    "coordinates": [coords]
                }
            }
        ]
    }
    return json.dumps(d)


def write_floats(filename, grid2d):
    """Create a binary (with acc. header file) version of a Grid2D object.

    Given a filename input of "probability.flt", this function will
    create that file, plus a text file called "probability.hdr".

    Args:
        filename (str): String filename to write (i.e., 'probability.flt')
        grid2d (Grid2D): MapIO Grid2D object.
    """
    geodict = grid2d.getGeoDict().asDict()
    array = grid2d.getData().astype('float32')
    np.save(filename, array)
    npyfilename = filename + '.npy'
    os.rename(npyfilename, filename)
    fpath, fname = os.path.split(filename)
    fbase, _ = os.path.splitext(fname)
    hdrfile = os.path.join(fpath, fbase + '.hdr')
    f = open(hdrfile, 'wt')
    for key, value in geodict.items():
        if isinstance(value, int):
            fmt = '%s = %i\n'
        elif isinstance(value, float):
            fmt = '%s = %.4f\n'
        else:
            fmt = '%s = %s\n'
        f.write(fmt % (key, value))
    f.close()


def savelayers(grids, filename):
    """
    Save ground failure layers object as a MultiHazard HDF file, preserving
    metadata structures. Must all have the same geodict.
    Args:
        grids: Ground failure layers object.
        filename (str): Path to where you want to save this file.
    """
    layers = collections.OrderedDict()
    metadata = collections.OrderedDict()
    for key in list(grids.keys()):
        layers[key] = grids[key]['grid'].getData()
        metadata[key] = {
            'description': grids[key]['description'],
            'type': grids[key]['type'],
            'label': grids[key]['label']
        }
    origin = {}
    header = {}
    mgrid = MultiHazardGrid(layers, grids[key]['grid'].getGeoDict(),
                            origin,
                            header,
                            metadata=metadata)
    mgrid.save(filename)


def loadlayers(filename):
    """
    Load a MultiHazard HDF file back in as a ground failure layers object in
    active memory (must have been saved for this purpose).
    Args:
        filename (str): Path to layers file.
    """
    mgrid = MultiHazardGrid.load(filename)
    grids = collections.OrderedDict()
    for key in mgrid.getLayerNames():
        grids[key] = {
            'grid': mgrid.getData()[key],
            'description': mgrid.getMetadata()[key]['description'],
            'type': mgrid.getMetadata()[key]['type'],
            'label': mgrid.getMetadata()[key]['label']
        }

    return grids