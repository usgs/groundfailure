#!/usr/bin/env python

import os
import matplotlib.pyplot as plt
from mapio.shake import ShakeGrid
from gfail import makemaps
from gfail.stats import computeStats
from impactutils.io.cmd import get_command_output
from shutil import copy
from datetime import datetime
from bs4 import BeautifulSoup
import shutil
import glob
import json
import collections
from configobj import ConfigObj
import numpy as np

from datetime import timedelta
from libcomcat.search import get_event_by_id, search
from libcomcat.classes import VersionOption
from mapio.shake import getHeaderData


# temporary until mapio is updated
import warnings
warnings.filterwarnings('ignore')

plt.switch_backend('agg')


def makeWebpage(maplayerlist, configs, web_template, shakemap, outfolder=None,
                includeunc=False, cleanup=True, includeAlert=False,
                alertkeyHAZ='Hagg_0.10g', alertkeyPOP='exp_pop_0.10g',
                faultfile=None,
                shakethreshtype='pga', point=False, pop_file=None,
                statlist=['Max', 'Std', 'Hagg_0.10g', 'exp_pop_0.10g'],
                probthresh=None, shakethresh=[5., 10.], statement=None):
    """
    Create a webpage that summarizes ground failure results (both landslides
        and liquefaction)

    Args:
        maplayerlist (list): list of model output structures to include.
        configs (list): list of paths to config files corresponding to each
            of the models in maplayerlist in the same order.
        web_template (str): Path to location of pelican template
            (final folder should be "theme").
        shakemap (str): path to shakemap .xml file for the current event.
        outfolder (str, optional): path to folder where output should be
            placed.
        includeunc (bool, optional): include uncertainty, NOT IMPLEMENTED.
        cleanup (bool, optional): cleanup all unneeded intermediate files that
            pelican creates, default True.
        includeAlert (bool, optional): if True, computes and reports alert
            level, default False.
        alertkey (str): stat key used for alert calculation
        faultfile (str, optional): GeoJson file of finite fault to display on
            interactive maps
        shakethreshtype (str, optional): Type of ground motion to use for stat
            thresholds, 'pga', 'pgv', or 'mmi'
        point (bool): True if it is known that the ShakeMap used a point
            source, False does not assume anything about source type.
        pop_file (str): Path to population file used for statistics
        statlist (list): list of strings indicating which stats to show on
            webpage.
        probthresh (float, optional): List of probability thresholds for which
            to compute Parea.
        shakethresh (list, optional): List of ground motion thresholds for
            which to compute Hagg, units corresponding to shakethreshtype.
        statement (str): Text to include in the summary section of the web
            page. Alert statements will be appended prior to this statement,
            Point source warnings for >M7 will be appended after this
            statement. If None, will use a generic explanatory statement.

    Returns:
        Folder where webpage files are located
    """
    print('Creating webpages')
    # get event id
    event_id = maplayerlist[0]['model']['description']['event_id']

    if outfolder is None:
        outfolder = os.path.join(os.getcwd(), event_id)
    fullout = os.path.join(outfolder, 'temp')
    finalout = os.path.join(outfolder, 'webpage')
    content = os.path.join(fullout, 'content')
    articles = os.path.join(content, 'articles')
    pages = os.path.join(content, 'pages')
    images = os.path.join(content, 'images')
    theme = web_template
    static = os.path.join(theme, 'static')
    if not os.path.exists(outfolder):
        os.mkdir(outfolder)
    if not os.path.exists(fullout):
        os.mkdir(fullout)
    if not os.path.exists(content):
        os.mkdir(content)
    if not os.path.exists(images):
        os.mkdir(images)
    if not os.path.exists(pages):
        os.mkdir(pages)
    if not os.path.exists(articles):
        os.mkdir(articles)
    if os.path.exists(finalout):
        shutil.rmtree(finalout)

    peliconf = os.path.join(fullout, 'pelicanconf.py')
    copy(os.path.join(os.path.dirname(web_template),
                      'pelicanconf.py'),
         peliconf)
    outjsfileLS = os.path.join(static, 'js', 'mapLS.js')
    with open(outjsfileLS, 'w') as f:
        f.write('\n')
    outjsfileLQ = os.path.join(static, 'js', 'mapLQ.js')
    with open(outjsfileLQ, 'w') as f:
        f.write('\n')

    if statement is None:
        statement = (
            "This product provides an early understanding of the "
            "landslides and liquefaction that may have been triggered "
            "by this earthquake until first responders and experts "
            "are able to survey the actual damage that has occurred. "
            "Results provide regional estimates and "
            "do not predict specific occurrences.")

    # Separate the LS and LQ models

    concLS = collections.OrderedDict()
    concLQ = collections.OrderedDict()

    lsmodels = {}
    lqmodels = {}
    logLS = []
    limLS = []
    colLS = []
    logLQ = []
    limLQ = []
    colLQ = []

    il = 0
    iq = 0

    filenames = []

    for conf, maplayer in zip(configs, maplayerlist):
        mdict = maplayer['model']['description']
        config = ConfigObj(conf)
        filename = '%s_%s' % (event_id, config.keys()[0])
        outfilebase = os.path.join(outfolder, filename)

        if 'landslide' in mdict['parameters']['modeltype'].lower():
            title = maplayer['model']['description']['name']
            plotorder, logscale, lims, colormaps, maskthreshes = \
                makemaps.parseConfigLayers(maplayer, conf, keys=['model'])
            logLS.append(logscale[0])
            limLS.append(lims[0])
            colLS.append(colormaps[0])
            concLS[title] = maplayer['model']

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
                                 pop_file=pop_file)

            if il == 0:
                on = True
            else:
                on = False
            metadata = maplayer['model']['description']
            if len(maplayer) > 1:
                inputs = {}
                inkeys = list(maplayer.keys())
                for key in inkeys:
                    if key != 'model':
                        newkey = maplayer[key]['label']
                        inputs[newkey] = maplayer[key]['description']
                metadata['inputs'] = inputs
            metad2 = json.dumps(metadata)
            filenames.append(outfilebase + '.json')
            with open(outfilebase + '.json', mode='w') as f3:
                f3.write(metad2)

            lsmodels[maplayer['model']['description']['name']] = {
                'geotiff_file': os.path.basename(outfilebase) + '.tif',
                'bin_edges': list(lims[0]),
                'metadata': metadata,
                'stats': dict(stats),
                'layer_on': on
            }
            il += 1
            filenames.append(outfilebase + '.tif')
        elif 'liquefaction' in mdict['parameters']['modeltype'].lower():
            title = maplayer['model']['description']['name']
            plotorder, logscale, lims, colormaps, maskthreshes = \
                makemaps.parseConfigLayers(maplayer, conf, keys=['model'])
            logLQ.append(logscale[0])
            limLQ.append(lims[0])
            colLQ.append(colormaps[0])
            concLQ[title] = maplayer['model']

            stats = computeStats(maplayer['model']['grid'],
                                 probthresh=probthresh,
                                 shakefile=shakemap,
                                 shakethresh=shakethresh)

            if iq == 0:
                on = True
            else:
                on = False
            metadata = maplayer['model']['description']
            if len(maplayer) > 1:
                inputs = {}
                inkeys = list(maplayer.keys())
                for key in inkeys:
                    if key != 'model':
                        newkey = maplayer[key]['label']
                        inputs[newkey] = maplayer[key]['description']
                metadata['inputs'] = inputs
            # Save metadata separately for each model
            metad2 = json.dumps(metadata)
            filenames.append(outfilebase + '.json')
            with open(outfilebase + '.json', mode='w') as f3:
                f3.write(metad2)

            lqmodels[maplayer['model']['description']['name']] = {
                'geotiff_file': os.path.basename(outfilebase) + '.tif',
                'bin_edges': list(lims[0]),
                'metadata': metadata,
                'stats': dict(stats),
                'layer_on': on
            }
            iq += 1
            filenames.append(outfilebase + '.tif')
        else:
            raise Exception("model type is undefined, check "
                            "maplayer['model']['parameters']"
                            "['modeltype'] to ensure it is defined")

    # Make interactive maps for each
    if il > 0:
        mapLS, filenameLS = makemaps.interactiveMap(
            concLS, shakefile=shakemap, scaletype='binned',
            colormaps=colLS, lims=limLS, clear_zero=False,
            logscale=logLS, separate=False, outfilename='LS_%s' % event_id,
            mapid='LS', savefiles=True, outputdir=images,
            sepcolorbar=True, floatcb=False, faultfile=faultfile,
            sync='Nowicki Jessee (2017)')
        makemaps.createColorsJson(concLS, colormaps=colLS,
                     lims=limLS, logscale=logLS, alpha=0.7,
                     outputdir=outfolder, outfilename='colorsLS.json',
                     sync='Nowicki Jessee (2017)')
        filenameLS = filenameLS[0]
    else:
        filenameLS = None

    if iq > 0:
        mapLQ, filenameLQ = makemaps.interactiveMap(
            concLQ, shakefile=shakemap, scaletype='binned',
            colormaps=colLQ, lims=limLQ, clear_zero=False,
            logscale=logLQ, separate=False, outfilename='LQ_%s' % event_id,
            savefiles=True, mapid='LQ', outputdir=images,
            sepcolorbar=True, floatcb=False, faultfile=faultfile,
            sync='Zhu and others (2017)')
        makemaps.createColorsJson(concLQ, colormaps=colLQ,
                     lims=limLQ, logscale=logLQ, alpha=0.7,
                     outputdir=outfolder, outfilename='colorsLQ.json',
                     sync='Zhu and others (2017)')
        filenameLQ = filenameLQ[0]
    else:
        filenameLQ = None

    # Get alert levels
    # TODO update to exact name of Hagg to use
    if includeAlert:
        try:
            paramalertLS = lsmodels['Nowicki Jessee (2017)']['stats'][alertkeyHAZ]
        except:
            paramalertLS = None
        try:
            parampopLS = lsmodels['Nowicki Jessee (2017)']['stats'][alertkeyPOP]
        except:
            parampopLS = None

        try:
            paramalertLQ = lqmodels['Zhu and others (2017)']['stats'][alertkeyHAZ]
        except:
            paramalertLQ = None

        try:
            parampopLQ = lqmodels['Zhu and others (2017)']['stats'][alertkeyPOP]
        except:
            parampopLQ = None

        alertLS, popalertLS, alertLQ, popalertLQ, alertstatementLS, alertstatementLQ = get_alert(
            paramalertLS, paramalertLQ, parampopLS, parampopLQ)
        topfileLQ = make_alert_img(alertLQ, popalertLQ, 'liquefaction', images)
        topfileLS = make_alert_img(alertLS, popalertLS, 'landslide', images)
        #statement = '%s<br>%s' % (alertstatementLQ, statement)
    else:
        alertLS = None
        alertLQ = None
        popalertLS = None
        popalertLQ = None
        topfileLQ = None
        topfileLS = None
        paramalertLS = None
        paramalertLQ = None
        parampopLS = None
        parampopLQ = None
        alertstatementLS = None
        alertstatementLQ = None

    if faultfile is not None:
        finitefault = True
    else:
        finitefault = False

    # Try to get urls
    try:
        info1, detail, temp = get_event_comcat(shakemap)
        eventurl = detail.url
        shakemapurl = detail.url + '#shakemap'
    except:
        shakemapurl = 'https://earthquake.usgs.gov/earthquakes/eventpage/%s#shakemap' % event_id
        eventurl = 'https://earthquake.usgs.gov/earthquakes/eventpage/%s#executive' % event_id

    sks = write_summary(shakemap, pages, images, point=point, event_url=eventurl,
                        statement=statement, finitefault=finitefault, shake_url=shakemapurl)

    # Create webpages for each type of ground motion
    write_individual(lsmodels, articles, 'Landslides',
                     interactivehtml=filenameLS, outjsfile=outjsfileLS,
                     topimage=topfileLS, statlist=statlist,
                     statement=alertstatementLS)
    write_individual(lqmodels, articles, 'Liquefaction',
                     interactivehtml=filenameLQ, outjsfile=outjsfileLQ,
                     topimage=topfileLQ, statlist=statlist,
                     statement=alertstatementLQ)

    # Create info.json for website rendering and metadata purposes
    web_dict = {
        'Summary': {
            'magnitude': sks['magnitude'],
            'depth': sks['depth'],
            'lat': sks['lat'],
            'lon': sks['lon'],
            'name': sks['name'],
            'date': sks['date'],
            'event_id': sks['event_id'],
            'event_url': eventurl,
            'shakemap_url': shakemapurl,
            'shakemap_version': sks['shakemap_version'],
            'statement': sks['statement'],
        },
        'Landslides': {
            'models': lsmodels,
            'alert': alertLS,
            'popalert': popalertLS,
            'alertkeyHAZ': alertkeyHAZ,
            'alertvalueHAZ': paramalertLS,
            'alertkeyPOP': alertkeyPOP,
            'alertvaluePOP': parampopLS,

        },
        'Liquefaction': {
            'models': lqmodels,
            'alert': alertLQ,
            'popalert': popalertLQ,
            'alertkeyHAZ': alertkeyHAZ,
            'alertvalueHAZ': paramalertLQ,
            'alertkeyPOP': alertkeyPOP,
            'alertvaluePOP': parampopLQ
        }
    }

    web_file = os.path.join(outfolder, 'info.json')
    filenames.append(web_file)

    with open(web_file, 'w') as f:
        json.dump(web_dict, f)

    # run website
    retcode, stdout, stderr = get_command_output(
        ('pelican -s %s -o %s -t %s')
        % (peliconf, os.path.join(fullout, 'output'), theme))
    print(stderr)

    if cleanup:  # delete everything except what is needed to make website
        files = glob.glob(os.path.join(fullout, '*'))
        for filen in files:
            if os.path.basename(filen) not in 'output':
                if os.path.isdir(filen):
                    shutil.rmtree(filen)
                else:
                    os.remove(filen)
            else:
                ofiles = glob.glob(os.path.join(fullout, 'output', '*'))
                for ofilen in ofiles:
                    if os.path.basename(ofilen) not in "index.htmlimagestheme":
                        if os.path.isdir(ofilen):
                            shutil.rmtree(ofilen)
                        else:
                            os.remove(ofilen)
        shutil.copytree(os.path.join(fullout, 'output'), finalout)
        shutil.rmtree(fullout)
        # Get rid of mapLS.js and mapLQ.js files from theme directory
        try:
            os.remove(os.path.join(theme, 'static', 'js', 'mapLQ.js'))
        except Exception as e:
            print(e)
        try:
            os.remove(os.path.join(theme, 'static', 'js', 'mapLS.js'))
        except Exception as e:
            print(e)
    filenames.append(finalout)
    return filenames


def write_individual(concatmods, outputdir, modeltype, topimage=None,
                     staticmap=None, interactivehtml=None,
                     outjsfile=None, statlist=None, statement=None):
    """
    Write markdown file for landslides or liquefaction.

    Args:
        concatmods (float or list): Ordered dictionary of models with
            fields required for write_individual populated (stats in
            particular)
        outputdir (str): Path to output directory.
        modeltype (str): 'Landslides' for landslide model, otherwise it is
            a liquefaction model.
        topimage (str, optional): Path to image for top of page.
        staticmap (str, optional): Path to static map.
        interactivehtml (str, optional): Path to interactive map file.
        outjsfile (str, optional): Path for output javascript file.
        stats (list): List of stats keys to include in the table, if None,
            it will include all of them
    """

    if modeltype == 'Landslides':
        id1 = 'LS'
    else:
        id1 = 'LQ'

    if len(concatmods) > 0:
        # If single model and not in list form, turn into lists
        modelnames = [key for key, value in concatmods.items()]
        # TODO Extract stats
        stattable = collections.OrderedDict()
        stattable['Model'] = modelnames
        if statlist is None:
            statlist = list(concatmods[modelnames[0]]['stats'].keys())

        # initialize empty lists for each
        for st in statlist:
            stattable[st] = []

        # put stats in
        for i, mod in enumerate(modelnames):
            for st in statlist:
                try:
                    stattable[st].append(concatmods[mod]['stats'][st])
                except:
                    stattable[st].append(float('nan'))

    if outjsfile is None:
        outjsfile = 'map.js'

    with open(os.path.join(outputdir, modeltype + '.md'), 'w') as file1:
        file1.write('title: %s\n' % modeltype.title())
        file1.write('date: %s\n'
                    % datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        file1.write('<h2>%s</h2>' % modeltype.title())

        if len(concatmods) > 0:

            if statement is not None:
                file1.write('<p>%s</p>\n' % statement)

            if topimage is not None:
                file1.write('<img src="images%s" width="250" '
                            'href="images%s"/>\n'
                            % (topimage.split('images')[-1],
                               topimage.split('images')[-1]))
            if interactivehtml is not None:
                # Extract js and move to map.js
                with open(interactivehtml) as f:
                    soup = BeautifulSoup(f, 'html.parser')
                    soup.prettify(encoding='utf-8')
                    scs = soup.find_all('script')
                    if scs is not None:
                        for sc in scs:
                            if 'var' in str(sc):
                                temp = str(sc)
                                temp = temp.strip(
                                    '<script>').strip('</script>')
                                with open(outjsfile, 'a') as f2:
                                    f2.write(temp)
                # Embed and link to fullscreen
                fileloc = interactivehtml.split('images')[-1]
                file1.write('<br><a href="images%s">Full '
                            'interactive map</a>\n'
                            % fileloc)

                file1.write('<center><div class="folium-map" id="map_%s">'
                            '</div></center>\n' % id1)

                cbname = fileloc.split('.html')[0] + '_colorbar' + '.png'
                file1.write('<center><img src="images%s" width="410" '
                            'href="images%s"/></center>\n'
                            % (cbname, cbname))
                if staticmap is not None:
                    file1.write('<center><img src="images%s" width="450" '
                                'href="images%s"/></center>\n'
                                % (staticmap.split('images')[-1],
                                   staticmap.split('images')[-1]))

                file1.write('<hr>\n')
                file1.write('<h3>%s Model Summary Statistics</h3>' %
                            modeltype.title())
                file1.write('<table style="width:100%">')
                file1.write('<tr><th>Model</th>')
                file1.write('<th>Output Type</th>')

                for st in statlist:
                    # if 'Hagg_' in st:
                    #    thresh = float(st.split('_')[-1].replace('g',''))
                    #    titl = 'Aggregate Hazard (km^2)' % (100.*thresh,)
                    if 'Hagg' in st:
                        titl = 'Aggregate Hazard (km^2)'
                    # elif 'exp_pop_' in st:
                    #    thresh = float(st.split('_')[-1].replace('g',''))
                    #    titl = 'People at risk (>%2.0f g)' % (100.*thresh,)
                    elif 'exp_pop' in st:
                        titl = 'People at risk'
                    elif 'Max' in st:
                        titl = 'Maximum probability'
                    elif 'Std' in st:
                        titl = 'Standard deviation'

                    file1.write('<th>%s</th>' % titl)
                file1.write('\n')

                # Write each row
                for i, mod in enumerate(modelnames):
                    file1.write('<tr>')
                    file1.write('<td>%s</td>' % mod.title())
                    # Get output type
                    if '2014' in mod:
                        type1 = 'Probability of any occurrence in cell'
                    else:
                        type1 = 'Proportion of area affected'
                    file1.write('<td>%s</td>' % type1)
                    for st in statlist:
                        if 'hagg' in st.lower() or 'exp' in st.lower() or 'parea' in st.lower():
                            file1.write('<td>%1.0f</td>' % stattable[st][i])
                        else:
                            file1.write('<td>%1.2f</td>' % stattable[st][i])

                    file1.write('</tr>\n')

                file1.write('</table>')
        else:
            file1.write('<h3>No results</h3>')


def write_summary(shakemap, outputdir, imgoutputdir, statement=None,
                  finitefault=False, point=False, event_url=None,
                  shake_url=None):
    """
    Write markdown file summarizing event

    Args:
        shakemap (str): path to shakemap .xml file for the current event
        outputdir (str): path to folder where output should be placed
        imgoutputdir (str): path to folder where images should be placed
            and linked
        HaggLS (float, optional): Aggregate hazard of preferred landslide model
        HaggLQ (float, optional): Aggregate hazard of preferred liquefaction
            model

    Returns:
        Markdown file
    """
    edict = ShakeGrid.load(shakemap, adjust='res').getEventDict()
    smdict = ShakeGrid.load(shakemap, adjust='res').getShakeDict()

    if event_url is None:
        event_url = 'https://earthquake.usgs.gov/earthquakes/eventpage/%s#executive' % edict['event_id']
    if shake_url is None:
        shake_url = event_url + '#shakemap'

    # NEED TO ADD FIX HERE IN CASE NO FINITE FAULT FILE MODEL IS AVAILABLE
    # BUT WAS USED IN SHAKEMAP
    if finitefault and not point:
        faulttype = '(finite fault model)'
    elif point and not finitefault:
        faulttype = '(point source model)'
        if edict['magnitude'] > 6.5:
            statement = ('%s<br>ShakeMap is currently approximating '
                         'this earthquake as a point source. '
                         'This may underestimate the extent of strong '
                         'shaking for larger earthquakes. '
                         'Please interpret Ground Failure results with '
                         'caution until ShakeMap has been updated '
                         'with a fault model.') % statement
    else:
        faulttype = ''

    with open(os.path.join(outputdir, 'Summary.md'), 'w') as file1:
        file1.write('title: summary\n')
        file1.write('date: 2017-06-09\n')
        file1.write('modified: 2017-06-09\n')
        file1.write('<h1 align="left">Ground Failure</h1>\n')
        if 'scenario' in smdict['shakemap_event_type'].lower():
            file1.write('<h2 align="left"><a href=%s>Magnitude %1.1f Scenario Earthquake - %s</a></h2>\n'
                        % (event_url, edict['magnitude'],
                           edict['event_description']))
        else:
            file1.write('<h2 align="left"><a href=%s>Magnitude %1.1f - %s</a></h2>\n'
                        % (event_url, edict['magnitude'],
                           edict['event_description']))

        writeline = '<h3 align="left"> %s (UTC) | %1.4f&#176,  %1.4f&#176 | %1.1f km depth</h3>\n' \
                    % (edict['event_timestamp'].strftime('%Y-%m-%dT%H:%M:%S'),
                        edict['lat'], edict['lon'], edict['depth'])
        file1.write(writeline)

        file1.write('<p align="left">Last updated at: %s (UTC)<br>'
                    % datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S'))
        file1.write('Based on ground motion estimates from '
                    '<a href=%s>ShakeMap</a> version %1.1f %s<br></p>'
                    % (shake_url, smdict['shakemap_version'], faulttype))
        statement = ('%s<br>Refer to the <a href=https://dev-earthquake.cr.usgs.gov/data/'
                     'grdfailure/background.php>Ground Failure Background</a>'
                     ' page for more details.' % statement)

        if statement is not None:
            file1.write('<h2 align="left">Summary</h2>\n')
            file1.write('<p align="left">%s</p>' % statement)

        file1.write('<hr>')

    shakesummary = {'magnitude': edict['magnitude'],
                    'shakemap_version': smdict['shakemap_version'],
                    'date': edict['event_timestamp'].strftime('%Y-%m-%dT%H:%M:%S'),
                    'lat': edict['lat'],
                    'lon': edict['lon'],
                    'depth': edict['depth'],
                    'name': 'Magnitude %1.1f - %s' % (edict['magnitude'],
                                                      edict['event_description']),
                    'statement': statement,
                    'event_id': edict['event_id'],
                    'shakemap_id': smdict['shakemap_id'],
                    'event_url': event_url
                    }

    return shakesummary


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
            hazard alert between Green and Yellow, Yellow and Orange, and Orange and Red.
        popbinLS (list): same as above but for population exposure
        hazbinLQ (list): 3 element list of bin edges for liquefaction hazard
            alert between Green and Yellow, Yellow and Orange, and Orange
            and Red.
        popbinLQ (list): same as above but for population exposure

    Returns:
        Returns:
            tuple: alertLS, alertLQ, statement, where
                * alertLS is the landslide alert level (str)
                * alertLQ is the liquefaction alert level (str)
                * statement is a sentence describing the ground failure hazard
                    based on the alert levels (str)
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

    if alertLS is not None:
        statementLS = ('Landslide hazard for this event is expected to be '
                       '%s with %s population at risk.'
                       % (get_word(alertLS), get_word(popalertLS)))
    else:
        statementLS = None

    if alertLQ is not None:
        statementLQ = ('Liquefaction hazard for this event is expected to be '
                       '%s with %s population at risk.'
                       % (get_word(alertLQ), get_word(popalertLQ)))
    else:
        statementLQ = None

    return alertLS, popalertLS, alertLQ, popalertLQ, statementLS, statementLQ


def get_word(color):
    """
    Get the alert-based word describing the hazard level.

    Args:
        color (str): Alert level; either 'green', 'yellow', 'orange', or 'red'.

    Returns:
        str: Phrase or word desribing hazard level.
    """
    if color is None:
        word = 'unknown'
    elif color in 'green':
        word = 'little to no'
    elif color in 'yellow':
        word = 'limited'
    elif color in 'orange':
        word = 'significant'
    elif color in 'red':
        word = 'extensive'
    else:
        word = 'unknown'
    return word


def make_alert_img(colorHAZ, colorPOP, type1, outfolder):
    """
    Construct alert image.

    Args:
        color (str): Alert color.
        type1 (str): Alert type, indicating landslide vs liquefaction.
        outfolder (str): Path for output file.

    Returns:
        str: Output file name.
    """
    if colorHAZ is None:
        outfilename = None
    else:
        fig = plt.figure(figsize=(6, 1.8))
        ax = fig.add_subplot(111)
        ax.add_artist(plt.Rectangle((0.1, 0.5), 0.2, 0.2, facecolor=colorHAZ,
                                    edgecolor='black', lw=1))
        #ax.text(0.1, 0.8, type1.title(), fontsize=25)
        ax.text(0.35, 0.52, 'Hazard %s' % colorHAZ, fontsize=25)
        ax.add_artist(plt.Rectangle((0.1, 0.2), 0.2, 0.2, facecolor=colorPOP,
                                    edgecolor='black', lw=1))
        ax.text(0.35, 0.22, 'Population exposure %s' % colorPOP, fontsize=25)

        ax.axis('off')
        ax.set_xlim([0, 2.0])
        ax.set_ylim([0., 0.8])
        plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)
        outfilename = os.path.join(outfolder, '%s_alert.png' % type1)
        fig.savefig(outfilename, transparent=True)
        plt.close()
    return outfilename


def get_event_comcat(shakefile, timewindow=60, degwindow=0.3, magwindow=0.2):
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
