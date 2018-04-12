#!/usr/bin/env python3
"""
TODO:
    - Add max probability to product properties
    - Potentially add more files (e.g., png, pdf)
"""
import os
import json
import shutil
from lxml import etree
from impactutils.io.cmd import get_command_output


def transfer(event_dir, pdl_conf, pdl_bin=None, source="us", dryrun=False):
    """
    This is to transfer the event's 'pdl_directory' to comcat.

    Args:
        event_dir (str): File path to location of results for event
        pdl_conf (str): Path to PDL conf file.
        eventid (str): Event id, if None, assumes that basename of event_dir is
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

    pdl_dir = os.path.join(event_dir, 'pdl_directory')

    # Load info.json
    with open(os.path.join(pdl_dir, 'info.json')) as f:
        info_dict = json.load(f)

    # Get some event info for pdl send command
    lat = info_dict['Summary']['lat']
    lon = info_dict['Summary']['lon']
    dep = info_dict['Summary']['depth']
    mag = info_dict['Summary']['magnitude']
    time_stamp = info_dict['Summary']['time']
    code = info_dict['Summary']['code']
    eventsourcecode = info_dict['Summary']['code']
    eventsource = info_dict['Summary']['net']

    pdl_type = 'groundfailure'

    # PDL properties
    lq_haz_alert = '"--property-lq_haz_alert=%s" ' % \
                   info_dict['Liquefaction'][0]['hazard_alert']
    ls_haz_alert = '"--property-ls_haz_alert=%s" ' % \
                   info_dict['Landslides'][0]['hazard_alert']
    lq_pop_alert = '"--property-lq_pop_alert=%s" ' % \
                   info_dict['Liquefaction'][0]['population_alert']
    ls_pop_alert = '"--property-ls_pop_alert=%s" ' % \
                   info_dict['Landslides'][0]['population_alert']

    lq_haz_alert_level = '"--property-lq_haz_alert_level=%s" ' % \
                         info_dict['Liquefaction'][0]['hazard_alert_value']
    ls_haz_alert_level = '"--property-ls_haz_alert_level=%s" ' % \
                         info_dict['Landslides'][0]['hazard_alert_value']
    lq_pop_alert_level = '"--property-lq_pop_alert_level=%s" ' % \
                         info_dict['Liquefaction'][0]['population_alert_value']
    ls_pop_alert_level = '"--property-ls_pop_alert_level=%s" ' % \
                         info_dict['Landslides'][0]['population_alert_value']

    # Construct PDL command
    pdl_cmd = (
        'java -jar %s ' % pdl_bin +
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
        lq_haz_alert +
        ls_haz_alert +
        lq_pop_alert +
        ls_pop_alert +
        lq_haz_alert_level +
        ls_haz_alert_level +
        lq_pop_alert_level +
        ls_pop_alert_level
    )

    if not dryrun:
        rc, so, se = get_command_output(pdl_cmd)
        return {'rc': rc, 'so': so, 'se': se}
    else:
        print(pdl_cmd)
        return pdl_cmd


def prepare_pdl_directory(event_dir):
    """
    Make director for transferring to comcat.

    Args:
        event_dir (str): Path to event directory
    """

    pdl_dir = os.path.join(event_dir, 'pdl_directory')
    if os.path.exists(pdl_dir):
        shutil.rmtree(pdl_dir)
    os.makedirs(pdl_dir)

    # Put geotif files into pdl directory
    all_files = os.listdir(event_dir)
    geotif_files = [os.path.join(event_dir, a)
                    for a in all_files if a.endswith('.tif')]
    for i in range(len(geotif_files)):
        src = geotif_files[i]
        tfile = os.path.basename(src)
        dst = os.path.join(pdl_dir, tfile)
        shutil.copy(src, dst)

    # Also the png files (default models for website interactive map,
    # not static maps)
    png_files = ['jessee_2017.png', 'zhu_2017.png']
    for i in range(len(png_files)):
        src = os.path.join(event_dir, png_files[i])
        tfile = os.path.basename(src)
        dst = os.path.join(pdl_dir, tfile)
        if os.path.exists(src):
            shutil.copy(src, dst)

    # Put json file into pdl directory
    src = os.path.join(event_dir, 'info.json')
    dst = os.path.join(pdl_dir, 'info.json')
    if os.path.exists(src):
        shutil.copy(src, dst)

    # Put hdf files into pdl directory
    hdf_files = [os.path.join(event_dir, a)
                 for a in all_files if a.endswith('.hdf5')]
    for i in range(len(hdf_files)):
        src = hdf_files[i]
        hfile = os.path.basename(src)
        dst = os.path.join(pdl_dir, hfile)
        shutil.copy(src, dst)

    # Put binary ShakeCast files into pdl directory
    flt_files = [os.path.join(event_dir, a)
                 for a in all_files if a.endswith('.flt')]
    flth_files = [os.path.join(event_dir, a)
                  for a in all_files if a.endswith('.hdr')]
    for f1, f2 in zip(flt_files, flth_files):
        src = f1
        f1file = os.path.basename(src)
        dst = os.path.join(pdl_dir, f1file)
        shutil.copy(src, dst)
        src = f2
        f2file = os.path.basename(src)
        dst = os.path.join(pdl_dir, f2file)
        shutil.copy(src, dst)

    # Make contents.xml
    contents = etree.Element("contents")

    # Geotif files
    tif_tree = [None] * len(geotif_files)
    file_caps = [None] * len(geotif_files)
    for i in range(len(geotif_files)):
        fname = os.path.basename(geotif_files[i])
        spl = fname.split('_')
        ftitle = spl[1].capitalize() + ' ' + spl[2] + ' Model (geotiff)'
        fid = '_'.join(spl[1:3])+"_gtiff"
        tif_tree[i] = etree.SubElement(
            contents, "file",
            title=ftitle,
            id=fid)
        file_caps[i] = etree.SubElement(tif_tree[i], "caption")
        file_caps[i].text = ftitle + ' Geotiff file'
        etree.SubElement(
            tif_tree[i], "format",
            href=fname,
            type='image/geotiff')

    # PNG files
    png_tree = [None] * len(png_files)
    file_caps = [None] * len(png_files)
    for i in range(len(png_files)):
        fname = os.path.basename(png_files[i])
        ftitle = fname
        spl = fname.split('.')
        fid = spl[0] + "_png"
        png_tree[i] = etree.SubElement(
            contents, "file",
            title=ftitle,
            id=fid)
        file_caps[i] = etree.SubElement(png_tree[i], "caption")
        file_caps[i].text = spl[0] + ' png file'
        etree.SubElement(
            png_tree[i], "format",
            href=fname,
            type='image/png')

    # Json info file
    fname = 'info.json'
    ftitle = 'Info'
    fid = 'info_json'
    j_tree = etree.SubElement(
        contents, "file",
        title=ftitle,
        id=fid)
    file_caps = etree.SubElement(j_tree, "caption")
    file_caps.text = ftitle + ' json file'
    etree.SubElement(
        j_tree, "format",
        href=fname,
        type='text/json')

    # Hdf5 files
    h_files = [None] * len(hdf_files)
    file_caps = [None] * len(hdf_files)
    for i in range(len(hdf_files)):
        fname = os.path.basename(hdf_files[i])
        bname = os.path.splitext(fname)[0]
        spl = bname.split('_')
        ftitle = spl[1].capitalize() + ' ' + spl[2] + \
            ' Model Results (hdf5)'
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
            type='application/x-hdf')

    # Flt ShakeCast files
    f_files = [None] * len(flt_files)
    fh_files = [None] * len(flth_files)
    file_caps = [None] * len(flt_files)
    fileh_caps = [None] * len(flth_files)
    for i in range(len(flt_files)):
        fname = os.path.basename(flt_files[i])
        bname = os.path.splitext(fname)[0]
        spl = bname.split('_')
        ftitle = spl[1].capitalize() + ' ' + spl[2] + ' Model (flt)'
        fid = '_'.join(spl[1:3])+"_flt"
        f_files[i] = etree.SubElement(
            contents, "file",
            title=ftitle,
            id=fid)
        file_caps[i] = etree.SubElement(f_files[i], "caption")
        file_caps[i].text = ftitle + ' binary float file'
        etree.SubElement(
            f_files[i], "format",
            href=fname,
            type='application/octet-stream')

        # same for header
        fname = os.path.basename(flth_files[i])
        bname = os.path.splitext(fname)[0]
        spl = bname.split('_')
        ftitle = spl[1].capitalize() + ' ' + spl[2] + ' Model (flt header)'
        fid = '_'.join(spl[1:3])+"_hdr"
        fh_files[i] = etree.SubElement(
            contents, "file",
            title=ftitle,
            id=fid)
        fileh_caps[i] = etree.SubElement(fh_files[i], "caption")
        fileh_caps[i].text = ftitle + ' header for flt file'
        etree.SubElement(
            fh_files[i], "format",
            href=fname,
            type='text/plain')

    # Write result
    out_file = os.path.join(pdl_dir, 'contents.xml')
    etree.ElementTree(contents).write(
        out_file,
        pretty_print=True,
        xml_declaration=True,
        encoding='UTF-8'
    )
