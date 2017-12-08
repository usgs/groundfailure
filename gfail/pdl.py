#!/usr/bin/env python3
"""
TODO:
    - Add max probability to product properties
    - Add support to specify location of PDL key
    - Potentially add more files (e.g., png, pdf)
"""
import os
import numpy as np
import json
import shutil
from lxml import etree
from impactutils.io.cmd import get_command_output
from mapio.shake import ShakeGrid


def transfer(eventdir, pdl_conf, pdl_bin=None, source="us", dryrun=False):
    """
    This is to transfer the event's 'pdl_directory' to comcat.

    Args:
        eventdir (str): File path to location of results for event
        pdl_conf (str): Path to PDL conf file.
        eventid (str): Event id, if None, assumes that basename of eventdir is
            the eventid
        pdl_bin (str): Path to 'ProductClient.jar'. If None it guesses that it
            is installed in the user's home directory:
            ``~/ProductClient/ProductClient.jar``.
        source (str): PDL 'source'. This is 'us' for products coming from NEIC,
            and is tied to the authentication key but this only matters once
            we start sending to production servers rather than dev servers, at
            which point we'll need to add a configuration option for the key.
        dryrun (bool): If True, the PDL command is constructed and printed but
            not executed.

    Returns:
        dict or str: Dictionary of pdl return code, standard out, and standard
        error for dryrun=False; PDL command for dryrun=True.
    """
    if pdl_bin is None:
        pdl_bin = os.path.join(os.path.expanduser('~'),
                               'ProductClient',
                               'ProductClient.jar')

    pdl_dir = os.path.join(eventdir, 'pdl_directory')

    # Get the shakefile
    shake_file = os.path.join(eventdir, 'shakefile.txt')
    sfile = open(shake_file, "r")
    shakefile = sfile.read()
    sfile.close()
    event_dict = ShakeGrid.load(shakefile, adjust='res').getEventDict()

    # Get some event info for pdl send command
    lat = event_dict['lat']
    lon = event_dict['lon']
    dep = event_dict['depth']
    mag = event_dict['magnitude']
    eventid = event_dict['event_id']
    time_stamp = event_dict['event_timestamp'].strftime('%Y-%m-%dT%H:%M:%SZ')

    # Note that eventsource is like 'catalog' for scenarios, but I think for
    # real events this is just the same as 'source'
    eventsource = event_dict['event_network']
    if eventid.startswith(eventsource):
        code = eventid
        nch = len(eventsource)
        eventsourcecode = eventid[nch:]
    else:
        code = eventsource + eventid
        eventsourcecode = eventid

    pdl_type = 'groundfailure'

    # PDL properties
    title = '"--property-title=Earthquake-Induced Ground Failure"'
    alert_file = os.path.join(eventdir, 'info.json')
    alert_json = json.load(open(alert_file))
    lq_alert = '"--property-alertLQ=%s" ' % alert_json['Liquefaction']['alert']
    ls_alert = '"--property-alertLS=%s" ' % alert_json['Landslides']['alert']
    lq_hagg = '"--property-HaggLQ=%s" ' % np.round(alert_json['Liquefaction']['paramalertLQ'], 0)
    ls_hagg = '"--property-HaggLS=%s" ' % np.round(alert_json['Landslides']['paramalertLS'], 0)

    # Construct PDL command
    pdl_cmd = ('java -jar %s ' % pdl_bin +
               '--send --configFile=%s ' % pdl_conf +
               '--source=%s ' % source +
               '--eventsource=%s ' % eventsource +
               '--code=%s ' % code +
               '--eventsourcecode=%s ' % eventsourcecode +
               '--latitude=%s ' % lat +
               '--longitude=%s ' % lon +
               '--magnitude=%s ' % mag +
               '--depth=%s ' % dep +
               '--eventtime=%s ' % time_stamp +
               '--type=%s ' % pdl_type +
               '--directory=%s ' % pdl_dir +
               title + " " + lq_alert + " " + ls_alert + " " +
               lq_hagg + " " + ls_hagg
               )

    if not dryrun:
        rc, so, se = get_command_output(pdl_cmd)
        return {'rc': rc, 'so': so, 'se': se}
    else:
        print(pdl_cmd)
        return pdl_cmd


def prepare_pdl_directory(eventdir):
    """
    Make director for transferring to comcat.

    Args:
        eventdir (str): Path to event directory
    """

    pdl_dir = os.path.join(eventdir, 'pdl_directory')
    if os.path.exists(pdl_dir):
        shutil.rmtree(pdl_dir)
    os.makedirs(pdl_dir)

    # Copy web page directory
    old_web_dir = os.path.join(eventdir, 'webpage')
    new_web_dir = os.path.join(pdl_dir, 'html')
    shutil.copytree(old_web_dir, new_web_dir)

    # Put geotif files into pdl directory
    all_files = os.listdir(eventdir)
    geotif_files = [a for a in all_files if a.endswith('.tif')]
    for i in range(len(geotif_files)):
        src = geotif_files[i]
        tfile = os.path.basename(src)
        dst = os.path.join(pdl_dir, tfile)
        shutil.move(src, dst)

    # Put json files into pdl directory
    json_files = [os.path.join(eventdir, a) for a in all_files if a.endswith('.json')]
    for i in range(len(json_files)):
        src = json_files[i]
        jfile = os.path.basename(src)
        dst = os.path.join(pdl_dir, jfile)
        shutil.copy(src, dst)
        
    # Put hdf files into pdl directory
    hdf_files = [os.path.join(eventdir, a) for a in all_files if a.endswith('.hdf5')]
    for i in range(len(hdf_files)):
        src = hdf_files[i]
        hfile = os.path.basename(src)
        dst = os.path.join(pdl_dir, hfile)
        shutil.copy(src, dst)
    
    # Make contents.xml
    contents = etree.Element("contents")

    # index.html
    file1 = etree.SubElement(
        contents, "file",
        title="Groundfailure Webpage",
        id='gf_html')
    caption1 = etree.SubElement(file1, "caption")
    caption1.text = 'Groundfailure Webpage'
    etree.SubElement(
        file1, "format",
        href='html/index.html',
        type='text/html')

    # Geotif files
    tif_files = [None] * len(geotif_files)
    file_caps = [None] * len(geotif_files)
    for i in range(len(geotif_files)):
        fname = os.path.basename(geotif_files[i])
        spl = fname.split('_')
        ftitle = spl[1].capitalize() + ' ' + spl[2] + ' Model'
        fid = '_'.join(spl[1:3])+"_gtiff"
        tif_files[i] = etree.SubElement(
            contents, "file",
            title=ftitle,
            id=fid)
        file_caps[i] = etree.SubElement(tif_files[i], "caption")
        file_caps[i].text = ftitle + ' Geotiff file'
        etree.SubElement(
            tif_files[i], "format",
            href=fname,
            type='image/geotiff')
        
    # Json files
    j_files = [None] * len(json_files)
    file_caps = [None] * len(json_files)
    for i in range(len(json_files)):
        fname = os.path.splitext(os.path.basename(json_files[i]))[0]
        if fname in 'info.json':
            ftitle = 'Info'
            fid = 'info_json'
        else:
            spl = fname.split('_')
            ftitle = spl[1].capitalize() + ' ' + spl[2] + ' Model Metadata'
            fid = '_'.join(spl[1:3])+"_json"
        j_files[i] = etree.SubElement(
            contents, "file",
            title=ftitle,
            id=fid)
        file_caps[i] = etree.SubElement(j_files[i], "caption")
        file_caps[i].text = ftitle + ' json file'
        etree.SubElement(
            j_files[i], "format",
            href=fname,
            type='text/json')

    # Hdf5 files
    h_files = [None] * len(hdf_files)
    file_caps = [None] * len(hdf_files)
    for i in range(len(hdf_files)):
        fname = os.path.splitext(os.path.basename(hdf_files[i]))[0]
        spl = fname.split('_')
        ftitle = spl[1].capitalize() + ' ' + spl[2] + ' Model Results (All)'
        fid = '_'.join(spl[1:3])+"_hdf5"
        h_files[i] = etree.SubElement(
            contents, "file",
            title=ftitle,
            id=fid)
        file_caps[i] = etree.SubElement(h_files[i], "caption")
        file_caps[i].text = ftitle + ' hdf5 file'
        etree.SubElement(
            h_files[i], "format",
            href=fname,
            type='text/hdf5')

    # Write result
    out_file = os.path.join(pdl_dir, 'contents.xml')
    etree.ElementTree(contents).write(
        out_file,
        pretty_print=True,
        xml_declaration=True,
        encoding='UTF-8'
    )
