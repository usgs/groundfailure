#!/usr/bin/env python

import os
import matplotlib.pyplot as plt
from mapio.shake import ShakeGrid
from gfail import makemaps
from gfail.assess_models import concatenateModels as concM
from gfail.assess_models import computeHagg
import numpy as np
from impactutils.io.cmd import get_command_output
from shutil import copy
from datetime import datetime
from bs4 import BeautifulSoup
import shutil
import glob
import json

# temporary until mapio is updated
import warnings
warnings.filterwarnings('ignore')


def makeWebpage(maplayerlist, configs, web_template, shakemap, outfolder=None,
                includeunc=False, cleanup=True, includeAlert=False,
                shakethreshtype='pga', shakethresh=0.0, faultfile=None):
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
        shakethreshtype (str, optional): Type of ground motion to use for Hagg threshold,
            'pga', 'pgv', or 'mmi'
        shakethresh (float, optional): Ground motion threshold corresponding to
            gmthreshtype. If None (default), no threshold will be used
        faultfile (str, optional): GeoJson file of finite fault to display on
            interactive maps

    Returns:
        Folder where webpage files are located

     TODO add in logic to deal with when one of the model types is missing
     """
    # get ShakeMap id
    sm_id = maplayerlist[0]['model']['description']['shakemap']
    if outfolder is None:
        outfolder = os.path.join(os.getcwd(), sm_id)
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
    LS = []
    LQ = []
    confLS = []
    confLQ = []
    for conf, maplayer in zip(configs, maplayerlist):
        mdict = maplayer['model']['description']
        if 'landslide' in mdict['parameters']['modeltype'].lower():
            LS.append(maplayer)
            confLS.append(conf)
        elif 'liquefaction' in mdict['parameters']['modeltype'].lower():
            LQ.append(maplayer)
            confLQ.append(conf)
        else:
            raise Exception("model type is undefined, check "
                            "maplayer['model']['parameters']"
                            "['modeltype'] to ensure it is defined")

    if len(LS) > 0:
        HaggLS = []
        maxLS = []
        logLS = []
        limLS = []
        colLS = []
        namesLS = []

        for conf, L in zip(confLS, LS):
            HaggLS.append(computeHagg(L['model']['grid'], probthresh=0.0,
                          shakefile=shakemap, shakethreshtype=shakethreshtype,
                          shakethresh=shakethresh))
            maxLS.append(np.nanmax(L['model']['grid'].getData()))
            plotorder, logscale, lims, colormaps, maskthreshes = \
                makemaps.parseConfigLayers(L, conf, keys=['model'])
            logLS.append(logscale[0])
            limLS.append(lims[0])
            colLS.append(colormaps[0])
            namesLS.append(L['model']['description']['name'])

        mapLS, filenameLS = makemaps.interactiveMap(
            concM(LS, astitle='model', includeunc=includeunc),
            colormaps=colLS, lims=limLS, clear_zero=False,
            logscale=logLS, separate=False, outfilename='LS_%s' % sm_id,
            mapid='LS', savefiles=True, outputdir=images,
            sepcolorbar=True, floatcb=False, faultfile=faultfile)

    if len(LQ) > 0:
        HaggLQ = []
        maxLQ = []
        logLQ = []
        limLQ = []
        colLQ = []
        namesLQ = []

        for conf, L in zip(confLQ, LQ):
            HaggLQ.append(computeHagg(L['model']['grid'], probthresh=0.0,
                          shakefile=shakemap, shakethreshtype=shakethreshtype,
                          shakethresh=shakethresh))
            maxLQ.append(np.nanmax(L['model']['grid'].getData()))
            plotorder, logscale, lims, colormaps, maskthreshes = \
                makemaps.parseConfigLayers(L, conf, keys=['model'])
            logLQ.append(logscale[0])
            limLQ.append(lims[0])
            colLQ.append(colormaps[0])
            namesLQ.append(L['model']['description']['name'])
        mapLQ, filenameLQ = makemaps.interactiveMap(
            concM(LQ, astitle='model', includeunc=includeunc),
            colormaps=colLQ, lims=limLQ, clear_zero=False,
            logscale=logLQ, separate=False, outfilename='LQ_%s' % sm_id,
            savefiles=True, mapid='LQ', outputdir=images,
            sepcolorbar=True, floatcb=False, faultfile=faultfile)

    write_summary(shakemap, pages, images,
                  HaggLS=HaggLS[namesLS.index('Nowicki and others (2014)')],
                  HaggLQ=HaggLQ[namesLQ.index('Zhu and others (2017)')])
    if includeAlert:
        alertLS, alertLQ, statement = get_alert(
            HaggLS=HaggLS[namesLS.index('Nowicki and others (2014)')],
            HaggLQ=HaggLQ[namesLQ.index('Zhu and others (2017)')])
        topfileLQ = make_alert_img(alertLQ, 'liquefaction', images)
        topfileLS = make_alert_img(alertLS, 'landslide', images)
    else:
        topfileLQ = None
        topfileLS = None
        alertLS = None
        alertLQ = None

    write_individual(HaggLQ, maxLQ, namesLQ, articles, 'Liquefaction',
                     interactivehtml=filenameLQ[0], outjsfile=outjsfileLQ,
                     topimage=topfileLQ)
    write_individual(HaggLS, maxLS, namesLS, articles, 'Landslides',
                     interactivehtml=filenameLS[0], outjsfile=outjsfileLS,
                     topimage=topfileLS)

    # Save some alert info to output directory for use later
    alert_dict = {
        'alertLS': alertLS,
        'alertLQ': alertLQ,
        'HaggLS': HaggLS[namesLS.index('Nowicki and others (2014)')],
        'HaggLQ': HaggLQ[namesLQ.index('Zhu and others (2017)')]}
    alert_file = os.path.join(outfolder, 'alert.json')
    with open(alert_file, 'w') as f:
        json.dump(alert_dict, f)

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

    return finalout


def write_individual(Hagg, maxprobs, modelnames, outputdir, modeltype,
                     topimage=None, staticmap=None, interactivehtml=None,
                     outjsfile=None):
    """
    Write markdown file for landslides or liquefaction.

    Args:
        Hagg (float or list): Aggregate hazard.
        maxprobs (float or list): Maximum probability.
        modelnames (str or list): Model name.
        outputdir (str): Path to output directory.
        modeltype (str): 'Landslides' for landslide model, otherwise it is
            a liquefaction model.
        topimage (str, optional): Path to image for top of page.
        staticmap (str, optional): Path to static map.
        interactivehtml (str, optional): Path to interactive map file.
        outjsfile (str, optional): Path for output javascript file.

    """
    if modeltype == 'Landslides':
        id1 = 'LS'
    else:
        id1 = 'LQ'

    # If single model and not in list form, turn into lists
    if type(Hagg) is float:
        Hagg = [Hagg]
        maxprobs = [maxprobs]
        modelnames = [modelnames]

    if outjsfile is None:
        outjsfile = 'map.js'

    with open(os.path.join(outputdir, modeltype + '.md'), 'w') as file1:
        file1.write('title: %s\n' % modeltype.title())
        file1.write('date: %s\n'
                    % datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        file1.write('<center><h2>%s</h2></center>' % modeltype.title())
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
        file1.write('<tr><th>Model</th><th>Aggregate Hazard</th><th>Max. '
                    'Probability</th></tr>\n')
        for H, m, n in zip(Hagg, maxprobs, modelnames):
            file1.write('<tr><td>%s</td><td>%0.2f km<sup>2</sup>'
                        '</td><td>%0.2f</td></tr>\n'
                        % (n.title(), H, m))
        file1.write('</table>')


def write_summary(shakemap, outputdir, imgoutputdir, HaggLS=None, HaggLQ=None):
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
    temp = ShakeGrid.load(shakemap, adjust='res').getShakeDict()
    edict['eventid'] = temp['shakemap_id']
    edict['version'] = temp['shakemap_version']
    alertLS, alertLQ, statement = get_alert(HaggLS, HaggLQ)

    with open(os.path.join(outputdir, 'Summary.md'), 'w') as file1:
        file1.write('title: summary\n')
        file1.write('date: 2017-06-09\n')
        file1.write('modified: 2017-06-09\n')
        file1.write('# Ground Failure\n')

        file1.write('## Magnitude %1.1f - %s\n'
                    % (edict['magnitude'],
                       edict['event_description']))

        file1.write('### %s (UTC) | %1.4f°,  %1.4f° | %1.1f km depth\n'
                    % (edict['event_timestamp'].strftime('%Y-%m-%dT%H:%M:%S'),
                       edict['lat'],
                       edict['lon'], edict['depth']))

        file1.write('Last updated at: %s (UTC)\n\n'
                    % datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        file1.write('Based on ground motion estimates from '
                    'ShakeMap version %1.1f\n'
                    % edict['version'])

        file1.write('## Summary\n')
        file1.write(statement)
        file1.write('<hr>')


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
    if color in 'green':
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
