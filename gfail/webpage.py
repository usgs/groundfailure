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

# temporary until mapio is updated
import warnings
warnings.filterwarnings('ignore')

plt.switch_backend('agg')


def makeWebpage(maplayerlist, configs, web_template, shakemap, outfolder=None,
                includeunc=False, cleanup=True, includeAlert=False, alertkey='Hagg_0.05g',
                faultfile=None, shakethreshtype='pga',
                statlist=['Max', 'Std', 'Hagg_0.05g', 'Hagg_0.10g', 'Parea_0.10', 'Parea_0.30'],
                probthresh=[0.1, 0.3], shakethresh=[5., 10.]):
    """
    Create a webpage that summarizes ground failure results (both landslides
        and liquefaction)

    Args:
        maplayerlist (list): list of model output structures to include
        configs (list): list of paths to config files corresponding to each
            of the models in maplayerlist in the same order
        web_template (str): Path to location of pelican template
            (final folder should be "theme")
        shakemap (str): path to shakemap .xml file for the current event
        outfolder (str, optional): path to folder where output should be placed
        includeunc (bool, optional): include uncertainty, NOT IMPLEMENTED
        cleanup (bool, optional): cleanup all unneeded intermediate files that
            pelican creates, default True.
        includeAlert (bool, optional): if True, computes and reports alert level, default
            False
        alertkey (str): stat key used for alert calculation
        faultfile (str, optional): GeoJson file of finite fault to display on
            interactive maps
        shakethreshtype (str, optional): Type of ground motion to use for stat thresholds,
            'pga', 'pgv', or 'mmi'
        statlist (list): list of strings indicating which stats to show on webpage
        probthresh (float, optional): List of probability thresholds for which to compute
            Parea.
        shakethresh (list, optional): List of ground motion thresholds for which
            to compute Hagg, units corresponding to shakethreshtype.

    Returns:
        Folder where webpage files are located
     """
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
                statprobthresh = 0.0  # Since logistic models can't equal one, need to eliminate placeholder zeros before computing stats

            stats = computeStats(maplayer['model']['grid'],
                                 probthresh=probthresh,
                                 shakefile=shakemap,
                                 shakethresh=shakethresh,
                                 statprobthresh=statprobthresh)

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

            lsmodels[maplayer['model']['description']['name']] = {'geotiff_file': os.path.basename(outfilebase) + '.tif',
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

            lqmodels[maplayer['model']['description']['name']] = {'geotiff_file': os.path.basename(outfilebase) + '.tif',
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
                sepcolorbar=True, floatcb=False, faultfile=faultfile)
            filenameLS = filenameLS[0]
        else:
            filenameLS = None

        if iq > 0:
            mapLQ, filenameLQ = makemaps.interactiveMap(
                concLQ, shakefile=shakemap, scaletype='binned',
                colormaps=colLQ, lims=limLQ, clear_zero=False,
                logscale=logLQ, separate=False, outfilename='LQ_%s' % event_id,
                savefiles=True, mapid='LQ', outputdir=images,
                sepcolorbar=True, floatcb=False, faultfile=faultfile)
            filenameLQ = filenameLQ[0]
        else:
            filenameLQ = None

        # Get alert levels
        #TODO update to exact name of Hagg to use
        if includeAlert:
            try:
                paramalertLS = lsmodels['Nowicki and others (2014)']['stats'][alertkey]
            except:
                paramalertLS = None

            try:
                paramalertLQ = lqmodels['Zhu and others (2017)']['stats'][alertkey]
            except:
                paramalertLQ = None

            alertLS, alertLQ, statement = get_alert(paramalertLS, paramalertLQ)
            topfileLQ = make_alert_img(alertLQ, 'liquefaction', images)
            topfileLS = make_alert_img(alertLS, 'landslide', images)
        else:
            alertLS = None
            alertLQ = None
            statement = None
            topfileLQ = None
            topfileLS = None
            paramalertLS = None
            paramalertLQ = None

    if faultfile is not None:
        finitefault = True
    else:
        finitefault = False
    sks = write_summary(shakemap, pages, images,
                        alert=includeAlert,
                        alertLS=alertLS, alertLQ=alertLQ,
                        statement=statement, finitefault=finitefault)

    # Create webpages for each type of ground motion
    write_individual(lsmodels, articles, 'Landslides',
                     interactivehtml=filenameLS, outjsfile=outjsfileLS,
                     topimage=topfileLS, statlist=statlist)
    write_individual(lqmodels, articles, 'Liquefaction',
                     interactivehtml=filenameLQ, outjsfile=outjsfileLQ,
                     topimage=topfileLQ, statlist=statlist)

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
            'event_url': sks['event_url'],
            'shakemap_url': 'https://earthquake.usgs.gov/earthquakes/eventpage/%s#shakemap' % sks['shakemap_id'],
            'shakemap_version': sks['shakemap_version'],
            'statement': sks['statement'],
            'scibackground': 'url placeholder'
            },
        'Landslides': {
            'models': lsmodels,
            'alert': sks['alertLS'],
            'alertkey': alertkey,
            'alertvalue': paramalertLS
            },
        'Liquefaction': {
            'models': lqmodels,
            'alert': sks['alertLQ'],
            'alertkey': alertkey,
            'alertvalue': paramalertLQ
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
                     outjsfile=None, statlist=None):
    """
    Write markdown file for landslides or liquefaction.

    Args:
        concatmods (float or list): Ordered dictionary of models with
            fields required for write_individual populated (stats in particular)
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
        #TODO Extract stats
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
                stattable[st].append(concatmods[mod]['stats'][st])

    if outjsfile is None:
        outjsfile = 'map.js'

    with open(os.path.join(outputdir, modeltype + '.md'), 'w') as file1:
        file1.write('title: %s\n' % modeltype.title())
        file1.write('date: %s\n'
                    % datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        file1.write('<center><h2>%s</h2></center>' % modeltype.title())

        if len(concatmods) > 0:
            if topimage is not None:
                file1.write('<center><img src="images%s" width="250" '
                            'href="images%s"/></center>\n'
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
                                temp = temp.strip('<script>').strip('</script>')
                                with open(outjsfile, 'a') as f2:
                                    f2.write(temp)
                # Embed and link to fullscreen
                fileloc = interactivehtml.split('images')[-1]
                file1.write('<center><a href="images%s">Click here for full '
                            'interactive map</a></center>\n'
                            % fileloc)

                file1.write('<center><div class="folium-map" id="map_%s">'
                            '</div></center>\n' % id1)

                cbname = fileloc.split('.html')[0] + '_colorbar' + '.png'
                file1.write('<center><img src="images%s" width="300" '
                            'href="images%s"/></center>\n'
                            % (cbname, cbname))
                if staticmap is not None:
                    file1.write('<center><img src="images%s" width="450" '
                                'href="images%s"/></center>\n'
                                % (staticmap.split('images')[-1],
                                   staticmap.split('images')[-1]))

                file1.write('<hr>\n')
                file1.write('<center><h3>%s Summary</h3></center>' % modeltype.title())
                file1.write('<table style="width:100%">')
                file1.write('<tr><th>Model</th>')
                for st in statlist:
                    file1.write('<th>%s</th>' % st)
                file1.write('\n')

                # Write each row
                for i, mod in enumerate(modelnames):
                    file1.write('<tr>')
                    file1.write('<td>%s</td>' % mod.title())
                    for st in statlist:
                        if 'hagg' in st.lower() or 'parea' in st.lower():
                            file1.write('<td>%1.0f</td>' % stattable[st][i])
                        else:
                            file1.write('<td>%1.2f</td>' % stattable[st][i])

                    file1.write('</tr>\n')

                file1.write('</table>')
        else:
            file1.write('<center><h3>No results</h3></center>')


def write_summary(shakemap, outputdir, imgoutputdir, alert=False,
                  alertLS=None, alertLQ=None, statement=None,
                  finitefault=False):
    """
    Write markdown file summarizing event

    Args:
        shakemap (str): path to shakemap .xml file for the current event
        outputdir (str): path to folder where output should be placed
        imgoutputdir (str): path to folder where images should be placed
            and linked
        HaggLS (float, optional): Aggregate hazard of preferred landslide model
        HaggLQ (float, optional): Aggregate hazard of preferred liquefaction model

    Returns:
        Markdown file
    """
    edict = ShakeGrid.load(shakemap, adjust='res').getEventDict()
    smdict = ShakeGrid.load(shakemap, adjust='res').getShakeDict()
    event_url = 'https://earthquake.usgs.gov/earthquakes/eventpage/%s#executive' % edict['event_id']

    if finitefault:
        faulttype = '(finite fault model)'
    else:
        faulttype = '(point source model)'
    with open(os.path.join(outputdir, 'Summary.md'), 'w') as file1:
        file1.write('title: summary\n')
        file1.write('date: 2017-06-09\n')
        file1.write('modified: 2017-06-09\n')
        file1.write('<h1>Ground Failure</h1>\n')
        if 'scenario' in smdict['shakemap_event_type'].lower():
            file1.write('<h2><a href=%s>Magnitude %1.1f Scenario Earthquake - %s</a></h2>\n'
                        % (event_url, edict['magnitude'],
                           edict['event_description']))
        else:
            file1.write('<h2><a href=%s>Magnitude %1.1f - %s</a></h2>\n'
                        % (event_url, edict['magnitude'],
                           edict['event_description']))

        writeline = '<h3> %s (UTC) | %1.4f&#176,  %1.4f&#176 | %1.1f km depth</h3>\n' \
                    % (edict['event_timestamp'].strftime('%Y-%m-%dT%H:%M:%S'),
                        edict['lat'], edict['lon'], edict['depth'])
        file1.write(writeline)

        file1.write('<p>Last updated at: %s (UTC)</p>\n'
                    % datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        file1.write('<p>Based on ground motion estimates from '
                    'ShakeMap version %1.1f %s</p>\n'
                    % (smdict['shakemap_version'], faulttype))
        if alert:
            file1.write('<h2>Summary</h2>\n')
            file1.write('<p>%s</p>' % statement)
            
        else:
            statement = None
            alertLS = None
            alertLQ = None

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
                    'alertLS': alertLS,
                    'alertLQ': alertLQ,
                    'event_id': edict['event_id'],
                    'shakemap_id': smdict['shakemap_id'],
                    'event_url': event_url
                    }

    return shakesummary


def get_alert(HaggLS, HaggLQ, binLS=[100., 850., 4000.],
              binLQ=[70., 500., 1200.]):
    """
    Get alert levels

    Args:
        HaggLS (float): Aggregate hazard (km2) of preferred landslide model
        HaggLQ (float): Aggregate hazard (km2) of preferred liquefaction model
        binLS (list, optional): 3 element list of bin edges for landslide alert
            between Green and Yellow, Yellow and Orange, and Orange and Red.
        binLQ (list, optional): 3 element list of bin edges for liquefaction
            alert between Green and Yellow, Yellow and Orange, and Orange
            and Red.

    Returns:
        Returns:
            tuple: alertLS, alertLQ, statement, where
                * alertLS is the landslide alert level (str)
                * alertLQ is the liquefaction alert level (str)
                * statement is a sentence describing the ground failure hazard
                    based on the alert levels (str)
    """
    if HaggLS is None:
        alertLS = None
    elif HaggLS < binLS[0]:
        alertLS = 'green'
    elif HaggLS >= binLS[0] and HaggLS < binLS[1]:
        alertLS = 'yellow'
    elif HaggLS >= binLS[1] and HaggLS < binLS[2]:
        alertLS = 'orange'
    elif HaggLS > binLS[2]:
        alertLS = 'red'
    else:
        alertLS = None

    if HaggLQ is None:
        alertLQ = None
    elif HaggLQ < binLQ[0]:
        alertLQ = 'green'
    elif HaggLQ >= binLQ[0] and HaggLQ < binLQ[1]:
        alertLQ = 'yellow'
    elif HaggLQ >= binLQ[1] and HaggLQ < binLQ[2]:
        alertLQ = 'orange'
    elif HaggLQ > binLQ[2]:
        alertLQ = 'red'
    else:
        alertLQ = None

    statement = ('Global ground failure models estimate that this earthquake '
                 'likely triggered %s liquefaction and %s landsliding'
                 % (get_word(alertLQ), get_word(alertLS)))

    return alertLS, alertLQ, statement


def get_word(color):
    """
    Get the alert-based word describing the hazard level.

    Args:
        color (str): Alert level; either 'green', 'yellow', 'orange', or 'red'.

    Returns:
        str: Phrase or word desribing hazard level.
    """
    if color is None:
        word = 'unknown levels of'
    elif color in 'green':
        word = 'little to no'
    elif color in 'yellow':
        word = 'localized'
    elif color in 'orange':
        word = 'substantial'
    elif color in 'red':
        word = 'extensive'
    else:
        word = 'unknown levels of'
    return word


def make_alert_img(color, type1, outfolder):
    """
    Construct alert image.

    Args:
        color (str): Alert color.
        type1 (str): Alert type, indicating landslide vs liquefaction.
        outfolder (str): Path for output file.

    Returns:
        str: Output file name.
    """
    if color is None:
        outfilename = None
    else:
        fig = plt.figure(figsize=(6, 1.8))
        ax = fig.add_subplot(111)
        ax.add_artist(plt.Circle((0.4, 0.3), 0.15, facecolor=color,
                                 edgecolor='black', lw=1))
        ax.text(0.7, 0.25, '%s alert %s' % (type1.title(), color), fontsize=25)
        ax.axis('off')
        ax.set_xlim([0, 2.0])
        ax.set_ylim([0., 0.6])
        plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)
        outfilename = os.path.join(outfolder, '%s_alert.png' % type1)
        fig.savefig(outfilename, transparent=True)
        plt.close()
    return outfilename
