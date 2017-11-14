import os
from mapio.shake import ShakeGrid
from groundfailure.makemaps import parseConfigLayers
from groundfailure import makemaps
from groundfailure.assessmodels import concatenateModels as concM
from groundfailure.assessmodels import computeHagg
#from groundfailure.newmark import godt2008
import numpy as np
from impactutils.io.cmd import get_command_output
from shutil import copy
from datetime import datetime
from bs4 import BeautifulSoup
import shutil
import glob

def makeWebpage(maplayerlist, configs, web_template, shakemap, outfolder=None,
                includeunc=False, cleanup=False):
    """
    :param maplayers: list of maplayer outputs from multiple models
    TODO add in logic to deal with when one of the model types is missing
    """
    # get ShakeMap id
    sm_id = maplayerlist[0]['model']['description']['shakemap']
    if outfolder is None:
        outfolder = os.path.join(os.getcwd(), sm_id)
    fullout = os.path.join(outfolder, 'webpage')
    content = os.path.join(fullout, 'content')
    articles = os.path.join(content, 'articles')
    #hidden_pages = os.path.join(content, 'hidden_pages')
    pages = os.path.join(content, 'pages')
    images = os.path.join(content, 'images')
    #images1 = os.path.join(images, sm_id)
    theme = web_template
    static = os.path.join(theme, 'static')
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
    copy(os.path.join(os.path.dirname(web_template), 'pelicanconf.py'), peliconf)
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
        if 'landslide' in maplayer['model']['description']['parameters']['modeltype'].lower():
            LS.append(maplayer)
            confLS.append(conf)
        elif 'liquefaction' in maplayer['model']['description']['parameters']['modeltype'].lower():
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
        #maskLS = []
        namesLS = []

        for conf, L in zip(confLS, LS):
            #TODO add threshold option for Hagg
            HaggLS.append(computeHagg(L['model']['grid']))
            maxLS.append(np.nanmax(L['model']['grid'].getData()))
            plotorder, logscale, lims, colormaps, maskthreshes = parseConfigLayers(L, conf, keys=['model'])
            logLS.append(logscale[0])
            limLS.append(lims[0])
            colLS.append(colormaps[0])
            #maskLS.append(maskthreshes[0])
            namesLS.append(L['model']['description']['name'])

        mapLS, filenameLS = makemaps.interactiveMap(concM(LS, astitle='model', includeunc=includeunc),
                                                    colormaps=colLS, lims=limLS, clear_zero=False,
                                                    logscale=logLS, separate=False, outfilename='LS_%s' % sm_id,
                                                    mapid='LS', savefiles=True, outputdir=images,
                                                    sepcolorbar=True, floatcb=False)
        write_individual(HaggLS, maxLS, namesLS, articles, 'Landslides',
                         interactivehtml=filenameLS[0], outjsfile=outjsfileLS)

    if len(LQ) > 0:
        HaggLQ = []
        maxLQ = []
        logLQ = []
        limLQ = []
        colLQ = []
        #maskLQ = []
        namesLQ = []

        for conf, L in zip(confLQ, LQ):
            HaggLQ.append(computeHagg(L['model']['grid']))
            maxLQ.append(np.nanmax(L['model']['grid'].getData()))
            plotorder, logscale, lims, colormaps, maskthreshes = parseConfigLayers(L, conf, keys=['model'])
            logLQ.append(logscale[0])
            limLQ.append(lims[0])
            colLQ.append(colormaps[0])
            #maskLQ.append(maskthreshes[0])
            namesLQ.append(L['model']['description']['name'])
        mapLQ, filenameLQ = makemaps.interactiveMap(concM(LQ, astitle='model', includeunc=includeunc),
                                                    colormaps=colLQ, lims=limLQ, clear_zero=False,
                                                    logscale=logLQ, separate=False, outfilename='LQ_%s' % sm_id,
                                                    savefiles=True, mapid='LQ', outputdir=images,
                                                    sepcolorbar=True, floatcb=False)

        write_individual(HaggLQ, maxLQ, namesLQ, articles, 'Liquefaction',
                         interactivehtml=filenameLQ[0], outjsfile=outjsfileLQ)

    #write_scibackground(LS, LQ)
    write_summary(shakemap, pages, images, HaggLS=HaggLS[namesLS=='Nowicki and others (2014)'],
                  HaggLQ=HaggLQ[namesLQ=='Zhu and others (2016)'])

    # run website
    retcode, stdout, stderr = get_command_output(('pelican -s %s -o %s -t %s') %
                                                (peliconf, os.path.join(fullout, 'output'), theme))
    print(stderr)
    #write_static_map(filenameLS, filenameLQ, static)

    if cleanup: # delete everything except what is needed to make website
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

    return os.path.join(fullout, 'output')


def write_individual(Hagg, maxprobs, modelnames, outputdir, modeltype,
                     topimage=None, staticmap=None, interactivehtml=None,
                     outjsfile=None):
    """
    write markdown file for landslides or liquefaction
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
        file1.write('date: %s\n' % datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        file1.write('<center><h2>%s</h2></center>' % modeltype.title())
        if topimage is not None:
            file1.write('  <img src="/images/%s" width="300" />\n' % topimage)
        if interactivehtml is not None:
            # Extract js and move to map.js
            with open(interactivehtml) as f:
                soup = BeautifulSoup(f, 'html.parser')
                soup.prettify(encoding='utf-8')
                scs = soup.find_all('script')
                if scs is not None:
                    for sc in scs:
                        if 'var' in str(sc):
                            temp= str(sc)
                            temp = temp.strip('<script>').strip('</script>')
                            with open(outjsfile, 'a') as f2:
                                f2.write(temp)
            # Embed and link to fullscreen
            fileloc = interactivehtml.split('images')[-1]
            file1.write('<center><a href="images%s">Click here for full interactive map</a></center>\n'
                        % fileloc)

            file1.write('<center><div class="folium-map" id="map_%s"></div></center>\n' % id1)

            cbname = fileloc.split('.html')[0] + '_colorbar' + '.png'
            file1.write('<center><img src="images%s" width="300" href="images%s"/></center>\n' %
                        (cbname, cbname))
        if staticmap is not None:
            #file1.write('<center> <h2>Static Map<h2> </center>\n')
            file1.write('<center><img src="images%s" width="450" href="images%s"/></center>\n' %
                        (staticmap.split('images')[-1], staticmap.split('images')[-1]))

        file1.write('<hr>\n')
        file1.write('<center><h3>Summary</h3></center>')
        file1.write('<table style="width:100%">')
        file1.write('<tr><th>Model</th><th>Aggregate Hazard</th><th>Max. Probability</th></tr>\n')
        for H, m, n in zip(Hagg, maxprobs, modelnames):
            file1.write('<tr><td>%s</td><td>%0.2f km<sup>2</sup></td><td>%0.2f</td></tr>\n' % (n.title(), H, m))
        file1.write('</table>')


def write_scibackground(configLS, configLQ):
    """
    write markdown file describing model background and references
    """
    pass


def write_summary(shakemap, outputdir, imgoutputdir, HaggLS=None, HaggLQ=None):
    edict = ShakeGrid.load(shakemap, adjust='res').getEventDict()
    temp = ShakeGrid.load(shakemap, adjust='res').getShakeDict()
    edict['eventid'] = temp['shakemap_id']
    edict['version'] = temp['shakemap_version']
    alertLS, alertLQ, statement = get_alert(HaggLS, HaggLQ)
    #TODO Make images for alerts
    
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
        file1.write('<hr>')


def get_alert(HaggLS, HaggLQ, binLS=[100., 850., 4000.], binLQ=[70., 120., 1000.]):
    """
    Bin edges (3 values) between Green and Yellow, Yellow and Orange, and Orange and Red
    LS based on Nowicki et al 2014 model results
    Red >4000
    Orange 850-4000
    Yellow 100-850
    Green <100
    LQ
    
    Based on Zhu et al 2017 general model results
    red = ~Loma Prieta >1000
    orange = Christchurch >120 <1000
    yellow = Greece >70 <120
    green = Northern Italy <70
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
    
    statement = 'This earthquake likely returned %s liquefaction and %s landsliding' % (get_word(alertLQ), get_word(alertLS))

    return alertLS, alertLQ, statement


def get_word(color):
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