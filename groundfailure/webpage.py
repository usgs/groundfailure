import os
from mapio.shake import ShakeGrid
from configobj import ConfigObj
from groundfailure.logisticmodel import LogisticModel
from groundfailure.conf import correct_config_filepaths
from groundfailure.makemaps import parseConfigLayers, parseMapConfig
from groundfailure import makemaps
from groundfailure.assessmodels import concatenateModels as concM
from groundfailure.assessmodels import computeHagg
#from groundfailure.newmark import godt2008
import numpy as np
from impactutils.io.cmd import get_command_output
from shutil import copy
from datetime import datetime


def makeWebpage(maplayerlist, configs, web_template, shakemap, outfolder=None, includeunc=False,
                cleanup=False):
    """
    :param maplayers: list of maplayer outputs from multiple models
    """
    # get ShakeMap id
    sm_id = maplayerlist[0]['model']['description']['shakemap']
    if outfolder is None:
        outfolder = os.getcwd()

    fullout = os.path.join(outfolder, sm_id)
    content = os.path.join(fullout, 'content')
    articles = os.path.join(content, 'articles')
    hidden_pages = os.path.join(content, 'hidden_pages')
    pages = os.path.join(content, 'pages')
    images = os.path.join(content, 'images')
    #images1 = os.path.join(images, sm_id)
    theme = os.path.join(web_template, 'theme')
    static = os.path.join(fullout, 'output', 'static')
    try:
        os.mkdir(fullout)
    except Exception as e:
        print(e)
    try:
        os.mkdir(content)
    except Exception as e:
        print(e)
    try:
        os.mkdir(images)
    except Exception as e:
        print(e)
    try:
        os.mkdir(pages)
    except Exception as e:
        print(e)
    try:
        os.mkdir(articles)
    except Exception as e:
        print(e)

    peliconf = os.path.join(fullout, 'pelicanconf.py')
    copy(os.path.join(web_template, 'pelicanconf.py'), peliconf)

    # Separate the LS and LQ models
    LS = []
    LQ = []
    confLS = []
    confLQ = []
    for conf, maplayer in zip(configs, maplayerlist):
        if 'landslide' in maplayer['model']['parameters']['modeltype'].lower():
            LS.append(maplayer)
            confLS.append(conf)
        elif 'liquefaction' in maplayer['model']['parameters']['modeltype'].lower():
            LQ.append(maplayer)
            confLQ.append(conf)
        else:
            raise Exception("model type is undefined, check maplayer['model']['parameters']['modeltype'] to ensure it is defined")

    if len(LS) > 0:
        HaggLS = []
        maxLS = []
        logLS = []
        limLS = []
        colLS = []
        maskLS = []
        namesLS = []

        for conf, L in zip(confLS, LS):
            #TODO add threshold option for Hagg
            HaggLS.append(computeHagg(L['model']['grid']))
            maxLS.append(np.nanmax(L['model']['grid'].getData()))
            plotorder, logscale, lims, colormaps, maskthreshes = parseConfigLayers(L, conf, keys=['model'])
            logLS.append(logscale)
            limLS.append(lims)
            colLS.append(colormaps)
            maskLS.append(maskthreshes)
            namesLS = L['model']['description']['name']

        # TURN INTO JSON AND INJECT THAT INTO TEMPLATE .JS INSTEAD OF MAKING HTML FILE
        map1, filenameLS = makemaps.interactiveMap(concM(LS, astitle='model', includeunc=includeunc),
                                                   maskthreshes=maskLS, colormaps=colLS, lims=limLS,
                                                   logscale=logLS)
        write_individual(HaggLS, maxLS, namesLS, articles, 'Landslides', interactivehtml=filenameLS)

    if len(LQ) > 0:
        HaggLQ = []
        maxLQ = []
        logLQ = []
        limLQ = []
        colLQ = []
        maskLQ = []
        namesLQ = []

        for conf, L in zip(confLQ, LQ):
            HaggLQ.append(computeHagg(L['model']['grid']))
            maxLQ.append(np.nanmax(L['model']['grid'].getData()))
            plotorder, logscale, lims, colormaps, maskthreshes = parseConfigLayers(L, conf, keys=['model'])
            logLQ.append(logscale)
            limLQ.append(lims)
            colLQ.append(colormaps)
            maskLQ.append(maskthreshes)
            namesLQ = L['model']['description']['name']
        map2, filenameLQ = makemaps.interactiveMap(concM(LQ, astitle='model', includeunc=includeunc),
                                                   maskthreshes=maskLQ, colormaps=colLQ, lims=limLQ,
                                                   logscale=logLQ)

        write_individual(HaggLQ, maxLQ, namesLQ, articles, 'Liquefaction', interactivehtml=filenameLQ)

    write_scibackground(LS, LQ)
    #statement = get_statement(HaggLS, HaggLQ)
    write_summary(shakemap, pages, statement=None)

    # run website
    retcode, stdout, stderr = get_command_output(('pelican -s %s -o %s -t %s') %
                                                (peliconf, os.path.join(fullout, 'output'), theme))
    print(stderr)
    write_static_map(filenameLS, filenameLQ, static)

    if cleanup:
        #delete the content folder
        pass


def write_individual(Hagg, maxprobs, modelnames, outputdir, modeltype,
                     topimage=None, staticmap=None, map1=None, interactivehtml=None):
    """
    write markdown file for landslides or liquefaction
    """
    # If single model and not in list form, turn into lists
    if type(Hagg) is float:
        Hagg = [Hagg]
        maxprobs = [maxprobs]
        modelnames = [modelnames]

    with open(os.path.join(outputdir, modeltype + '.md'), 'w') as file1:
        file1.write('title: %s\n' % modeltype.title())
        file1.write('date: %s\n' % datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        file1.write('''
            <center>
            <h2 class="uk-panel-title">%s</h2>
            </center>''' % modeltype.title())

        if topimage is not None:
            file1.write('  <img src="/images/%s" width="300" />\n' % topimage)
        file1.write('''<center>
            Model | Aggregate Hazard | Maximum Probability
            :---: | :---: | :---:\n''')
        for H, m, n in zip(Hagg, maxprobs, modelnames):
            file1.write('%s | %0.2f km<sup>2</sup> | %0.2f\n' % (n.title(), H, m))
        file1.write('</center>\n<hr>\n')

        if interactivehtml is not None:
            file1.write('<center><div class="folium-map" id="map%s"></div></center>' % modeltype)
                #file1.write('    <center><object type="text/html" data=images%s height=450 width=450></object></center>\n'
                #            % interactivehtml.split('images')[-1])
            file1.write('    <center><a href="images%s">Click here for full interactive map</a></center>'
                        % interactivehtml.split('images')[-1])
        if staticmap is not None:
            #file1.write('<center> <h2>Static Map<h2> </center>\n')
            file1.write('    <center><img src="images%s" width="450" href="images%s"/></center>\n' %
                        (staticmap.split('images')[-1], staticmap.split('images')[-1]))


def write_static_map(filenameLS, filenameLQ, static):
    # assign 'map%s' % modeltype to id
    # css between <style> markers
    # js between <script> markers (exclude if src)
    pass



def write_scibackground(configLS, configLQ):
    """
    write markdown file describing model background and references
    """
    if modelname is None:
        modelname = 'unspecified'
    with open(os.path.join(outputdir, modeltype + '.md'), 'w') as file1:
        file1.write('title: %s\n' % modeltype)
        file1.write('date: %s\n' % datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        file1.write('<center> <h2>%s</h2> </center>\n' % modeltype.title())
        file1.write('<center> <h3>Model: %s<h3> </center>\n' % modelname)
        if topimage is not None:
            file1.write('  <img src="/images/%s" width="300" />\n' % topimage)
        file1.write('<center> Aggregate Hazard: %0.2f km<sup>2</sup> </center>\n' % Hagg)
        file1.write('<center> Maximum probability: %0.2f </center>\n' % max1)
        file1.write('<p float="center">\n')

        if interactivehtml is not None:
            file1.write('<center><object type="text/html" data=images%s height=450 width=450></object></center>\n' % interactivehtml.split('images')[-1])
            file1.write('<center><a href="images%s">Click here for full interactive map</a></center>' % interactivehtml.split('images')[-1])
            #cbar = interactivehtml.split('images')[-1].split('.html')[0] + '_colorbar.png'
            #file1.write('  <center><img src="images%s" width="400" href="images%s"/></center>\n' % (cbar, cbar))
        if staticmap is not None:
            #file1.write('<center> <h2>Static Map<h2> </center>\n')
            file1.write('  <center><img src="images%s" width="450" href="images%s"/></center>\n' % (staticmap.split('images')[-1], staticmap.split('images')[-1]))
        file1.write('</p>')


def write_summary(shakemap, outputdir, statement):
    edict = ShakeGrid.load(shakemap, adjust='res').getEventDict()
    temp = ShakeGrid.load(shakemap, adjust='res').getShakeDict()
    edict['eventid'] = temp['shakemap_id']
    edict['version'] = temp['shakemap_version']

    with open(os.path.join(outputdir, 'Summary.md'), 'w') as file1:
        file1.write('title: summary\n')
        file1.write('date: 2017-06-09\n')
        file1.write('modified: 2017-06-09\n')
        file1.write('# Ground failure\n')
        file1.write('### Last updated at: %s (UTC)\n' % datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        file1.write('### Based on ground motion estimates from ShakeMap version %1.1f\n' % edict['version'])
        file1.write('## Magnitude %1.1f - %s\n' % (edict['magnitude'], edict['event_description']))
        file1.write('### %s (UTC) | %1.4f°,  %1.4f° | %1.1f km\n' %
                    (edict['event_timestamp'].strftime('%Y-%m-%dT%H:%M:%S'), edict['lat'],
                     edict['lon'], edict['depth']))

        file1.write('### Summary\n')
        file1.write(statement)


def get_statement(HaggLS=None, HaggLQ=None):
    """
    get standardized statement based on Hagg of landslides and liquefaction
    """
    return 'Some automatically produced statement about the levels of landslide and liquefaction hazard'


def runmodels(shakemap, configLS, configLQ, fileloc, mapconfig=None, mapdatadir=None, pgalim=None):
    """
    """
    plotorder, logscale, limsLS, colormapsLS, maskthreshesLS = parseConfigLayers(LS, configLS)
    HaggLS = computeHagg(LS['model']['grid'])
    maxLS = np.nanmax(LS['model']['grid'].getData())
    map1 = makemaps.interactiveMap(LS, shakefile=shakemap, maskthreshes=maskthreshesLS,
                                   colormaps=colormapsLS, lims=limsLS,
                                   outputdir=fileloc, outfilename='Landslides', sepcolorbar=True)
    if pgalim is not None:
        bounds = get_bounds(shakemap, parameter='pga', threshold=pgalim)
    else:
        bounds = None
    if mapconfig is not None and mapdatadir is not None:
        makemaps.modelMap(LS, maskthreshes=maskthreshesLS, colormaps=colormapsLS, lims=limsLS,
                          outputdir=fileloc, outfilename='Landslides', savepng=True,
                          savepdf=False, boundaries=bounds, maproads=False, **kwargs)
    configLQ = ConfigObj(configLQ)
    configLQ = correct_config_filepaths(datadir, configLQ)
    lmn = LogisticModel(shakemap, configLQ, saveinputs=False)
    LQ = lmn.calculate()
    plotorder, logscale, limsLQ, colormapsLQ, maskthreshesLQ = parseConfigLayers(LQ, configLQ)
    HaggLQ = computeHagg(LQ['model']['grid'])
    maxLQ = np.nanmax(LQ['model']['grid'].getData())
    map2 = makemaps.interactiveMap(LQ, shakefile=shakemap, maskthreshes=maskthreshesLQ,
                                   colormaps=colormapsLQ, lims=limsLQ,
                                   outputdir=fileloc, outfilename='Liquefaction', sepcolorbar=True)
    if mapconfig is not None and mapdatadir is not None:
        makemaps.modelMap(LQ, maskthreshes=maskthreshesLQ, colormaps=colormapsLQ, lims=limsLQ,
                          outputdir=fileloc, outfilename='Liquefaction', savepng=True,
                          savepdf=False, boundaries=bounds, maproads=False, **kwargs)
    return HaggLS, HaggLQ, maxLS, maxLQ



