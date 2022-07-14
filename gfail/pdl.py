#!/usr/bin/env python3
"""
TODO:
    - Add max probability to product properties
    - Potentially add more files (e.g., png, pdf)
"""
import os
import glob
import json
import shutil
from lxml import etree
from configobj import ConfigObj

# import zipfile
from impactutils.io.cmd import get_command_output


def transfer(
    event_dir,
    version,
    pdl_conf,
    pdl_bin=None,
    source="us",
    dryrun=False,
    status="UPDATE",
):
    """
    This is to transfer the event's 'pdl_directory' to comcat. PDL must be
    installed separately, see https://usgs.github.io/pdl/ for information.

    Args:
        event_dir (str): File path to location of results for event.
        version (int): Version number of ground-failure run.
        pdl_conf (str): Path to PDL conf file.
        pdl_bin (str): Path to 'ProductClient.jar'. If None it guesses that it
            is installed in the user's home directory:
            ``~/ProductClient/ProductClient.jar``.
        source (str): PDL 'source'. This is 'us' for products coming from NEIC,
            and is tied to the authentication key but this only matters once
            we start sending to production servers rather than dev servers, at
            which point we'll need to add a configuration option for the key.
        dryrun (bool): If True, the PDL command is constructed and printed but
            not executed.
        status (str): Status of ground-failure product being sent to comcat.
            Default is "UPDATE" but can also be "WARNING" so that the product
            page displays the warning banner.

    Returns:
        dict or str: Dictionary of pdl return code, standard out, and standard
        error for dryrun=False; PDL command for dryrun=True.
    """
    if pdl_bin is None:
        pdl_bin = os.path.join(
            os.path.expanduser("~"), "ProductClient", "ProductClient.jar"
        )

    pdl_dir = os.path.join(event_dir, "pdl_directory")

    # Load info.json
    with open(os.path.join(pdl_dir, "info.json")) as f:
        info_dict = json.load(f)

    # Get some event info for pdl send command
    lat = info_dict["Summary"]["lat"]
    lon = info_dict["Summary"]["lon"]
    dep = info_dict["Summary"]["depth"]
    mag = info_dict["Summary"]["magnitude"]
    time_stamp = info_dict["Summary"]["time"]
    code = info_dict["Summary"]["code"]
    eventsourcecode = info_dict["Summary"]["code"]
    eventsource = info_dict["Summary"]["net"]
    shake_version = info_dict["Summary"]["shakemap_version"]
    xmin, xmax, ymin, ymax = info_dict["Summary"]["zoom_extent"]

    pdl_type = "ground-failure"

    # PDL properties
    # Get preferred models to extract PDL properties
    lqs = info_dict["Liquefaction"]
    lq_pref = None
    ls_pref = None
    for lq in lqs:
        if lq["preferred"]:
            lq_pref = lq
    lss = info_dict["Landslides"]
    for ls in lss:
        if ls["preferred"]:
            ls_pref = ls

    ls_alert = '"--property-landslide-alert=%s" ' % ls_pref["alert"]
    lq_alert = '"--property-liquefaction-alert=%s" ' % lq_pref["alert"]

    lq_haz_alert = (
        '"--property-liquefaction-hazard-alert-color=%s" '
        % lq_pref["hazard_alert"]["color"]
    )
    ls_haz_alert = (
        '"--property-landslide-hazard-alert-color=%s" '
        % ls_pref["hazard_alert"]["color"]
    )
    lq_pop_alert = (
        '"--property-liquefaction-population-alert-color=%s" '
        % lq_pref["population_alert"]["color"]
    )
    ls_pop_alert = (
        '"--property-landslide-population-alert-color=%s" '
        % ls_pref["population_alert"]["color"]
    )

    lq_haz_alert_value = (
        '"--property-liquefaction-hazard-alert-value=%s" '
        % lq_pref["hazard_alert"]["value"]
    )
    ls_haz_alert_value = (
        '"--property-landslide-hazard-alert-value=%s" '
        % ls_pref["hazard_alert"]["value"]
    )
    lq_pop_alert_value = (
        '"--property-liquefaction-population-alert-value=%s" '
        % lq_pref["population_alert"]["value"]
    )
    ls_pop_alert_value = (
        '"--property-landslide-population-alert-value=%s" '
        % ls_pref["population_alert"]["value"]
    )

    if "std" in list(lq_pref["hazard_alert"].keys()):
        lq_haz_alert_std = (
            '"--property-liquefaction-hazard-std=%s" ' % lq_pref["hazard_alert"]["std"]
        )
        ls_haz_alert_std = (
            '"--property-landslide-hazard-std=%s" ' % ls_pref["hazard_alert"]["std"]
        )
        lq_pop_alert_std = (
            '"--property-liquefaction-population-std=%s" '
            % lq_pref["population_alert"]["std"]
        )
        ls_pop_alert_std = (
            '"--property-landslide-population-std=%s" '
            % ls_pref["population_alert"]["std"]
        )
        ls_haz_range1s = '"--property-landslide-hazard-1std=%s" ' % str(
            ls_pref["probability"]["hagg_1std"]
        ).replace("[", "").replace("]", "")
        ls_haz_range2s = '"--property-landslide-hazard-2std=%s" ' % str(
            ls_pref["probability"]["hagg_2std"]
        ).replace("[", "").replace("]", "")
        ls_pop_range1s = '"--property-landslide-population-1std=%s" ' % str(
            ls_pref["probability"]["pop_1std"]
        ).replace("[", "").replace("]", "")
        ls_pop_range2s = '"--property-landslide-population-2std=%s" ' % str(
            ls_pref["probability"]["pop_2std"]
        ).replace("[", "").replace("]", "")
        lq_haz_range1s = '"--property-liquefaction-hazard-1std=%s" ' % str(
            lq_pref["probability"]["hagg_1std"]
        ).replace("[", "").replace("]", "")
        lq_haz_range2s = '"--property-liquefaction-hazard-2std=%s" ' % str(
            lq_pref["probability"]["hagg_2std"]
        ).replace("[", "").replace("]", "")
        lq_pop_range1s = '"--property-liquefaction-population-1std=%s" ' % str(
            lq_pref["probability"]["pop_1std"]
        ).replace("[", "").replace("]", "")
        lq_pop_range2s = '"--property-liquefaction-population-2std=%s" ' % str(
            lq_pref["probability"]["pop_2std"]
        ).replace("[", "").replace("]", "")

        stdstr = (
            lq_haz_alert_std
            + ls_haz_alert_std
            + lq_pop_alert_std
            + ls_pop_alert_std
            + ls_haz_range1s
            + ls_haz_range2s
            + ls_pop_range1s
            + ls_pop_range2s
            + lq_haz_range1s
            + lq_haz_range2s
            + lq_pop_range1s
            + lq_pop_range2s
        )
    else:
        stdstr = ""

    lq_haz_alert_par = (
        '"--property-liquefaction-hazard-alert-parameter=%s" '
        % lq_pref["hazard_alert"]["parameter"]
    )
    ls_haz_alert_par = (
        '"--property-landslide-hazard-alert-parameter=%s" '
        % ls_pref["hazard_alert"]["parameter"]
    )
    lq_pop_alert_par = (
        '"--property-liquefaction-population-alert-parameter=%s" '
        % lq_pref["population_alert"]["parameter"]
    )
    ls_pop_alert_par = (
        '"--property-landslide-population-alert-parameter=%s" '
        % ls_pref["population_alert"]["parameter"]
    )

    lq_overlay = '"--property-liquefaction-overlay=%s" ' % lq_pref["overlay"]
    ls_overlay = '"--property-landslide-overlay=%s" ' % ls_pref["overlay"]

    lq_extent = lq_pref["extent"]
    ls_extent = ls_pref["extent"]

    # Liquefaction extent
    lq_xmin = '"--property-liquefaction-minimum-longitude=%s" ' % lq_extent[0]
    lq_xmax = '"--property-liquefaction-maximum-longitude=%s" ' % lq_extent[1]
    lq_ymin = '"--property-liquefaction-minimum-latitude=%s" ' % lq_extent[2]
    lq_ymax = '"--property-liquefaction-maximum-latitude=%s" ' % lq_extent[3]

    # Landslide extent
    ls_xmin = '"--property-landslide-minimum-longitude=%s" ' % ls_extent[0]
    ls_xmax = '"--property-landslide-maximum-longitude=%s" ' % ls_extent[1]
    ls_ymin = '"--property-landslide-minimum-latitude=%s" ' % ls_extent[2]
    ls_ymax = '"--property-landslide-maximum-latitude=%s" ' % ls_extent[3]

    # Product extent -- default zoom extent for interactive map
    prod_xmin = '"--property-minimum-longitude=%s" ' % xmin
    prod_xmax = '"--property-maximum-longitude=%s" ' % xmax
    prod_ymin = '"--property-minimum-latitude=%s" ' % ymin
    prod_ymax = '"--property-maximum-latitude=%s" ' % ymax

    rupt_warn = (
        '"--property-rupture-warning=%s" ' % info_dict["Summary"]["rupture_warning"]
    )

    # Check for PDL key:
    defaults = os.path.join(os.path.expanduser("~"), ".gfail_defaults")
    configs = ConfigObj(defaults)
    if "key" in configs.keys():
        pdl_key = configs["key"]
    else:
        pdl_key = None

    # Construct PDL command
    pdl_cmd = (
        "java -jar %s " % pdl_bin
        + "--send --configFile=%s " % pdl_conf
        + "--source=%s " % source
        + "--eventsource=%s " % eventsource
        + "--code=%s " % code
        + "--status=%s " % status
        + "--eventsourcecode=%s " % eventsourcecode
        + "--version=%s " % version
        + "--latitude=%s " % lat
        + "--longitude=%s " % lon
        + "--magnitude=%s " % mag
        + "--depth=%s " % dep
        + "--eventtime=%s " % time_stamp
        + "--type=%s " % pdl_type
        + "--directory=%s " % pdl_dir
        + ls_alert
        + lq_alert
        + lq_haz_alert
        + ls_haz_alert
        + lq_pop_alert
        + ls_pop_alert
        + lq_haz_alert_value
        + ls_haz_alert_value
        + lq_pop_alert_value
        + ls_pop_alert_value
        + lq_haz_alert_par
        + ls_haz_alert_par
        + lq_pop_alert_par
        + ls_pop_alert_par
        + stdstr
        + lq_overlay
        + ls_overlay
        + lq_xmin
        + lq_xmax
        + lq_ymin
        + lq_ymax
        + ls_xmin
        + ls_xmax
        + ls_ymin
        + ls_ymax
        + prod_xmin
        + prod_xmax
        + prod_ymin
        + prod_ymax
        + rupt_warn
        + '"--property-shakemap-version=%s" ' % shake_version
        + "--signatureVersion=v2"
    )

    if pdl_key is not None:
        pdl_cmd = pdl_cmd + " --privateKey=%s " % pdl_key

    if not dryrun:
        rc, so, se = get_command_output(pdl_cmd)
        print("PDL return code: %s " % rc)
        print("PDL standard output:\n%s " % so)
        print("PDL standard error:\n%s " % se)
        return {"rc": rc, "so": so, "se": se}
    else:
        print(pdl_cmd)
        return {"rc": True, "so": b"", "se": b""}


def prepare_pdl_directory(event_dir):
    """
    Make directory for transferring to comcat.

    Args:
        event_dir (str): Path to event directory

    Returns:
        event_dir containing copies of all of the files that need to
        be sent to comcat
    """

    pdl_dir = os.path.join(event_dir, "pdl_directory")
    if os.path.exists(pdl_dir):
        shutil.rmtree(pdl_dir)
    os.makedirs(pdl_dir)

    # Get the event id prefix that is prepended to strip it off later
    all_files = os.listdir(event_dir)
    an_hdf_file = [f for f in all_files if f.endswith(".hdf5")][0]
    idx = an_hdf_file.find("_201")
    author = an_hdf_file[:idx].split("_")[-1]
    idx2 = an_hdf_file.find(author)
    event_prefix = an_hdf_file[: idx2 - 1]

    # Put geotif files into pdl directory
    geotif_files = [os.path.join(event_dir, a) for a in all_files if a.endswith(".tif")]
    for i in range(len(geotif_files)):
        src = geotif_files[i]
        tfile = os.path.basename(src)
        if tfile.startswith(event_prefix):
            tfile = tfile[len(event_prefix) + 1 :]
        tfile = tfile.replace('_slim', '') # get rid of 'slim' in names if present
        dst = os.path.join(pdl_dir, tfile)
        shutil.copy(src, dst)

        # Zip the files
        # zfile = os.path.join(pdl_dir, os.path.splitext(tfile)[0] + '.zip')
        # with zipfile.ZipFile(zfile, 'w') as zf:
        #     zf.write(os.path.join(pdl_dir, dst), arcname=dst,
        #              compress_type=zipfile.ZIP_DEFLATED)

        # # Remove originals
        # os.remove(os.path.join(pdl_dir, dst))

    # Put kmz files into pdl directory
    kmz_files = [os.path.join(event_dir, a) for a in all_files if a.endswith(".kmz")]
    for i in range(len(kmz_files)):
        src = kmz_files[i]
        tfile = os.path.basename(src)
        if tfile.startswith(event_prefix):
            tfile = tfile[len(event_prefix) + 1 :]
        tfile = tfile.replace('_slim', '') # get rid of 'slim' in names if present
        dst = os.path.join(pdl_dir, tfile)
        shutil.copy(src, dst)

    # Also the png files (default models for website interactive map,
    # not static maps)
    png_files = [os.path.join(event_dir, a) for a in all_files if a.endswith(".png")]
    for i in range(len(png_files)):
        src = os.path.join(event_dir, png_files[i])
        tfile = os.path.basename(src)
        tfile = tfile.replace('_slim', '') # get rid of 'slim' in names if present
        dst = os.path.join(pdl_dir, tfile)
        if os.path.exists(src):
            shutil.copy(src, dst)

    # Put json file into pdl directory
    src = os.path.join(event_dir, "info.json")
    dst = os.path.join(pdl_dir, "info.json")
    if os.path.exists(src):
        shutil.copy(src, dst)

    # Put hdf files into pdl directory
    hdf_files = [os.path.join(event_dir, a) for a in all_files if a.endswith(".hdf5")]
    for i in range(len(hdf_files)):
        src = hdf_files[i]
        hfile = os.path.basename(src)
        if hfile.startswith(event_prefix):
            hfile = hfile[len(event_prefix) + 1 :]
        hfile = hfile.replace('_slim', '') # get rid of 'slim' in names if present
        dst = os.path.join(pdl_dir, hfile)
        # Strip off event id prefix

        shutil.copy(src, dst)

    # Put binary ShakeCast files into pdl directory
    # flt_files = [os.path.join(event_dir, a)
    #              for a in all_files if a.endswith('.flt')]
    # flth_files = [os.path.join(event_dir, a)
    #               for a in all_files if a.endswith('.hdr')]
    # for f1, f2 in zip(flt_files, flth_files):
    #     src = f1
    #     f1file = os.path.basename(src)
    #     if f1file.startswith(event_prefix):
    #         f1file = f1file[len(event_prefix)+1:]
    #     dst = os.path.join(pdl_dir, f1file)
    #     shutil.copy(src, dst)
    #     src = f2
    #     f2file = os.path.basename(src)
    #     if f2file.startswith(event_prefix):
    #         f2file = f2file[len(event_prefix)+1:]
    #     dst = os.path.join(pdl_dir, f2file)
    #     shutil.copy(src, dst)

    # Make contents.xml
    contents = etree.Element("contents")

    json_mime = "text/json"
    hdf_mime = "application/x-hdf"
    gtif_mime = "image/geotiff"
    png_mime = "image/png"
    kmz_mime = "application/vnd.google-earth.kmz"
    # zip_mime = 'application/zip'

    # Json info file
    j_tree = etree.SubElement(contents, "file", title="Info", id="info_json")
    file_caps = etree.SubElement(j_tree, "caption")
    file_caps.text = "Info json file"
    etree.SubElement(j_tree, "format", href="info.json", type=json_mime)

    # Jessee section
    jessee_tree = etree.SubElement(
        contents,
        "file",
        title="Preferred Landslide Model (displayed)",
        id="nowicki_jessee_2018",
    )
    file_caps = etree.SubElement(jessee_tree, "caption")
    file_caps.text = "Nowicki Jessee and others (2018)"
    etree.SubElement(jessee_tree, "format", href="jessee_2018_model.kmz", type=kmz_mime)
    etree.SubElement(
        jessee_tree, "format", href="jessee_2018_model.tif", type=gtif_mime
    )
    etree.SubElement(jessee_tree, "format", href="jessee_2018.hdf5", type=hdf_mime)
    etree.SubElement(jessee_tree, "format", href="jessee_2018.png", type=png_mime)

    # Does an uncertainty file exist?
    file_base = "jessee_2018_beta_sigma"
    if len(glob.glob(os.path.join(pdl_dir, "%s*" % file_base))):
        jessee_uncertainty_tree = etree.SubElement(
            contents,
            "file",
            title="Preferred Landslide Model Uncertainty (not displayed)",
            id="nowicki_jessee_2018uncertainty",
        )
        file_caps = etree.SubElement(jessee_uncertainty_tree, "caption")
        file_caps.text = "Nowicki Jessee and others (2018) uncertainty"
        extensions = [".kmz", ".tif"]
        for extension in extensions:
            mime_type = kmz_mime if "kmz" in extension else gtif_mime
            std_file = file_base + extension
            if os.path.exists(os.path.join(pdl_dir, std_file)):
                etree.SubElement(
                    jessee_uncertainty_tree, "format", href=std_file, type=mime_type
                )

    # Zhu 2017 section
    zhu2017_tree = etree.SubElement(
        contents,
        "file",
        title="Preferred Liquefaction Model (displayed)",
        id="zhu_2017",
    )
    file_caps = etree.SubElement(zhu2017_tree, "caption")
    file_caps.text = "Zhu and others (2017)"
    etree.SubElement(
        zhu2017_tree, "format", href="zhu_2017_general_model.kmz", type=kmz_mime
    )
    etree.SubElement(
        zhu2017_tree, "format", href="zhu_2017_general_model.tif", type=gtif_mime
    )
    etree.SubElement(
        zhu2017_tree, "format", href="zhu_2017_general.hdf5", type=hdf_mime
    )
    etree.SubElement(zhu2017_tree, "format", href="zhu_2017_general.png", type=png_mime)

    # Does an uncertainty file exist?
    file_base = "zhu_2017_general_beta_sigma"
    if len(glob.glob(os.path.join(pdl_dir, "%s*" % file_base))):
        zhu2017_uncertainty_tree = etree.SubElement(
            contents,
            "file",
            title="Preferred Liquefaction Model Uncertainty (not displayed)",
            id="zhu_2017uncertainty",
        )
        file_caps = etree.SubElement(zhu2017_uncertainty_tree, "caption")
        file_caps.text = "Zhu and others (2017) uncertainty"
        extensions = [".kmz", ".tif"]
        for extension in extensions:
            mime_type = kmz_mime if "kmz" in extension else gtif_mime
            std_file = file_base + extension
            if os.path.exists(os.path.join(pdl_dir, std_file)):
                etree.SubElement(
                    zhu2017_uncertainty_tree, "format", href=std_file, type=mime_type
                )

    # Godt section
    godt_tree = etree.SubElement(
        contents, "file", title="Alternative Landslide Model 1 (not displayed)"
    )
    file_caps = etree.SubElement(godt_tree, "caption")
    file_caps.text = "Godt and others (2008)"
    etree.SubElement(godt_tree, "format", href="godt_2008.hdf5", type=hdf_mime)
    etree.SubElement(godt_tree, "format", href="godt_2008_model.tif", type=gtif_mime)

    # Nowicki section
    now_tree = etree.SubElement(
        contents, "file", title="Alternative Landslide Model 2 (not displayed)"
    )
    file_caps = etree.SubElement(now_tree, "caption")
    file_caps.text = "Nowicki and others (2014)"
    etree.SubElement(now_tree, "format", href="nowicki_2014_global.hdf5", type=hdf_mime)
    etree.SubElement(
        now_tree, "format", href="nowicki_2014_global_model.tif", type=gtif_mime
    )

    # zhu 2015 section
    zhu2015_tree = etree.SubElement(
        contents,
        "file",
        title="Alternative Liquefaction Model (not displayed)",
        id="zhu_2015",
    )
    file_caps = etree.SubElement(zhu2015_tree, "caption")
    file_caps.text = "Zhu and others (2015)"
    etree.SubElement(zhu2015_tree, "format", href="zhu_2015.hdf5", type=hdf_mime)
    etree.SubElement(zhu2015_tree, "format", href="zhu_2015_model.tif", type=gtif_mime)

    # Copy over legend files
    #    data_dir = pkg_resources.resource_filename('gfail', 'data')
    #    src = os.path.join(data_dir, 'legend_landslide.png')
    #    dst = os.path.join(pdl_dir, 'legend_landslide.png')
    #    shutil.copy(src, dst)
    #    src = os.path.join(data_dir, 'legend_liquefaction.png')
    #    dst = os.path.join(pdl_dir, 'legend_liquefaction.png')
    #    shutil.copy(src, dst)

    # Write result
    out_file = os.path.join(pdl_dir, "contents.xml")
    etree.ElementTree(contents).write(
        out_file, pretty_print=True, xml_declaration=True, encoding="UTF-8"
    )
