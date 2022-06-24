# stdlib imports
from configobj import ConfigObj
import os
import shutil
import numpy as np
import tempfile
import urllib
import json
from argparse import Namespace
from zipfile import ZipFile
import warnings
import time
import logging

# third party imports

from mapio.shake import getHeaderData
from mapio.gdal import GDALGrid
from mapio.gmt import GMTGrid
from mapio.geodict import GeoDict
from impactutils.io.cmd import get_command_output
from mapio.shake import ShakeGrid
from libcomcat.search import get_event_by_id

# local imports

from gfail.conf import correct_config_filepaths
import gfail.logisticmodel as LM
from gfail.godt import godt2008
from gfail.godt2 import godt2008_2
from gfail.webpage import hazdev, create_kmz
from gfail.utilities import (
    get_event_comcat,
    parseConfigLayers,
    text_to_json,
    savelayers,
    getFileType,
)

from gfail.zhu_2015 import Zhu2015Model
from gfail.nowicki_2014 import Nowicki2014Model
from gfail.zhu_2017 import Zhu2017Model
from gfail.zhu_2017_coastal import Zhu2017ModelCoastal
from gfail.jessee_2018 import Jessee2018Model


MODEL_FACTORY = {
    "zhu_2015": Zhu2015Model,
    "zhu_2017_general": Zhu2017Model,
    "zhu_2017_coastal": Zhu2017ModelCoastal,
    "nowicki_2014_global": Nowicki2014Model,
    "jessee_2018": Jessee2018Model,
}


def run_gfail(args):
    """Runs ground failure.

    Args:
        args: dictionary or argument parser Namespace output by bin/gfail
            program.

    Returns:
        list: Names of created files.

    """
    # if user asked for logging to stdout, set that up here
    if args.debug:
        logging.basicConfig(level=logging.INFO)
    # TODO: ADD CONFIG VALIDATION STEP THAT MAKES SURE ALL THE FILES EXIST
    filenames = []
    # If args is a dictionary, convert to a Namespace
    if isinstance(args, dict):
        args = Namespace(**args)

    if args.set_default_paths:
        set_default_paths(args)
        print("default paths set, continuing...\n")

    if args.list_default_paths:
        list_default_paths()
        return

    if args.reset_default_paths:
        reset_default_paths()
        return

    if args.make_webpage:
        # Turn on GIS and HDF5 flags
        gis = True
        hdf5 = True
        kmz = True
    else:
        gis = args.gis
        hdf5 = args.hdf5
        kmz = args.kmz

    # Figure out what models will be run
    if args.shakefile is not None:  # user intends to actually run some models
        shakefile = args.shakefile

        # make output location for things
        if args.output_filepath is None:
            outdir = os.getcwd()
        else:
            outdir = args.output_filepath

        if hdf5 or gis or kmz:
            if not os.path.exists(outdir):
                os.makedirs(outdir)

        # download if is url
        # cleanup = False
        if not os.path.isfile(shakefile):
            if isURL(shakefile):
                # getGridURL returns a named temporary file object
                shakefile = getGridURL(shakefile)
                # cleanup = True  # Be sure to delete it after
            else:
                raise NameError(
                    'Could not find "%s" as a file or a valid url' % shakefile
                )
        eventid = getHeaderData(shakefile)[0]["event_id"]

        # Get entire path so won't break if running gfail with relative path
        shakefile = os.path.abspath(shakefile)

        if args.extract_contents:
            outfolder = outdir
        else:  # Nest in a folder named by eventid
            outfolder = os.path.join(outdir, eventid)
            if not os.path.exists(outfolder):
                os.makedirs(outfolder)

        # Copy shake grid into output directory
        # --- this is base on advice from Mike that when running in production
        #     the shake grids are not archived and so if we need/want to have
        #     the exact grid used for the calculation later if there's every a
        #     question about how the calculation was done, the safest thing is
        #     to store a copy of it here.
        shake_copy = os.path.join(outfolder, "grid.xml")
        shutil.copyfile(shakefile, shake_copy)

        if args.uncertfile is not None:
            uncertfile = os.path.abspath(args.uncertfile)
            unc_copy = os.path.join(outfolder, "uncertainty.xml")
            shutil.copyfile(uncertfile, unc_copy)
        else:
            uncertfile = None

        # Write shakefile to a file for use later
        shakename = os.path.join(outfolder, "shakefile.txt")
        shake_file = open(shakename, "wt")
        shake_file.write(shake_copy)
        shake_file.close()
        filenames.append(shakename)

        # Check that shakemap bounds do not cross 180/-180 line
        sd = ShakeGrid.getFileGeoDict(shakefile)
        if not args.keep_shakemap_bounds:
            if args.set_bounds is None:
                if sd.xmin > sd.xmax:
                    print(
                        "\nShakeMap crosses 180/-180 line, setting bounds so "
                        "only side with more land area is run"
                    )
                    if sd.xmax + 180.0 > 180 - sd.xmin:
                        set_bounds = "%s, %s, %s, %s" % (
                            sd.ymin,
                            sd.ymax,
                            -180.0,
                            sd.xmax,
                        )
                    else:
                        set_bounds = "%s, %s, %s, %s" % (
                            sd.ymin,
                            sd.ymax,
                            sd.xmin,
                            180.0,
                        )
                    print("Bounds applied: %s" % set_bounds)
                else:
                    set_bounds = args.set_bounds
            else:
                set_bounds = args.set_bounds
        else:  # we're using a logbase model version
            set_bounds = "%s, %s, %s, %s" % (sd.ymin, sd.ymax, sd.xmin, sd.xmax)

        config = args.config

        if args.config_filepath is not None:
            # only add config_filepath if full filepath not given and file
            # ext is .ini
            if not os.path.isabs(config) and os.path.splitext(config)[-1] == ".ini":
                config = os.path.join(args.config_filepath, config)

        if os.path.splitext(config)[-1] == ".ini":
            temp = ConfigObj(config)
            if len(temp) == 0:
                raise Exception("Could not find specified .ini file: %s" % config)
            if args.data_path is not None:
                temp = correct_config_filepaths(args.data_path, temp)
            configs = [temp]
            conffail = []
        else:
            # input is a list of config files
            f = open(config, "r")
            configlist = f.readlines()
            configs = []
            conffail = []
            for conf in configlist:
                conf = conf.strip()
                if not os.path.isabs(conf):
                    # only add config_filepath if full filepath not given
                    conf = os.path.join(args.config_filepath, conf)
                try:
                    temp = ConfigObj(conf)
                    if temp:
                        if args.data_path is not None:
                            temp = correct_config_filepaths(args.data_path, temp)
                        configs.append(temp)
                    else:
                        conffail.append(conf)
                except BaseException:
                    conffail.append(conf)

        print("\nRunning the following models:")

        for conf in configs:
            print("\t%s" % conf.keys()[0])
        if len(conffail) > 0:
            print("Could not find or read in the following config files:\n")
            for conf in conffail:
                print("\t%s" % conf)
            print("\nContinuing...\n")

        if set_bounds is not None:
            if "zoom" in set_bounds:
                temp = set_bounds.split(",")
                print(
                    "Using %s threshold of %1.1f to cut model bounds"
                    % (temp[1].strip(), float(temp[2].strip()))
                )
                bounds = get_bounds(shakefile, temp[1].strip(), float(temp[2].strip()))
            else:
                temp = eval(set_bounds)
                latmin = temp[0]
                latmax = temp[1]
                lonmin = temp[2]
                lonmax = temp[3]
                bounds = {
                    "xmin": lonmin,
                    "xmax": lonmax,
                    "ymin": latmin,
                    "ymax": latmax,
                }
            print(
                "Applying bounds of lonmin %1.2f, lonmax %1.2f, "
                "latmin %1.2f, latmax %1.2f"
                % (bounds["xmin"], bounds["xmax"], bounds["ymin"], bounds["ymax"])
            )
        else:
            bounds = None

        if args.make_webpage:
            results = []

        # pre-read in ocean trimming file polygons so only do this step once
        if args.trimfile is not None:
            if not os.path.exists(args.trimfile):
                print(
                    "trimfile defined does not exist: %s\n"
                    "Ocean will not be trimmed." % args.trimfile
                )
                trimfile = None
            elif os.path.splitext(args.trimfile)[1] != ".shp":
                print("trimfile must be a shapefile, ocean will not be trimmed")
                trimfile = None
            else:
                trimfile = args.trimfile
        else:
            trimfile = None

        # Get finite fault ready, if exists

        ffault = None
        point = True
        if args.finite_fault is not None:
            point = False
            try:
                if os.path.splitext(args.finite_fault)[-1] == ".txt":
                    ffault = text_to_json(args.finite_fault)
                elif os.path.splitext(args.finite_fault)[-1] == ".json":
                    ffault = args.finite_fault
                else:
                    print(
                        "Could not read in finite fault, will "
                        "try to download from comcat"
                    )
                    ffault = None
            except BaseException:
                print(
                    "Could not read in finite fault, will try to "
                    "download from comcat"
                )
                ffault = None

        if ffault is None:
            # Try to get finite fault file, if it exists
            try:
                returned_ev = get_event_comcat(shakefile)
                if returned_ev is not None:
                    testjd, detail, temp = returned_ev
                    evinfo = testjd["input"]["event_information"]
                    if "faultfiles" in evinfo:
                        ffilename = evinfo["faultfiles"]
                        if len(ffilename) > 0:
                            # Download the file
                            with tempfile.NamedTemporaryFile(
                                delete=False, mode="w"
                            ) as f:
                                temp.getContent(ffilename, filename=f.name)
                                ffault = text_to_json(f.name)
                                os.remove(f.name)
                            point = False
                        else:
                            point = True
                else:
                    print(
                        "Unable to determine source type, unknown if finite"
                        " fault or point source"
                    )
                    ffault = None
                    point = False

            except Exception as e:
                print(e)
                print(
                    "Unable to determine source type, unknown if finite"
                    " fault or point source"
                )
                ffault = None
                point = False

        # Loop over config files
        for conf in configs:
            modelname = conf.keys()[0]
            print("\nNow running %s:" % modelname)
            notcov, newbnds = check_input_extents(
                conf, shakefile=shakefile, bounds=bounds
            )
            if len(notcov) > 0:
                print(
                    "\nThe following input layers do not cover"
                    " the area of interest:\n\t%s" % "\n\t".join(notcov)
                )
                if newbnds is None:
                    print("\nCannnot make bounds that work. Skipping to next model\n")
                    continue
                else:
                    pnt = "%s, %s, %s, %s" % (
                        newbnds["xmin"],
                        newbnds["xmax"],
                        newbnds["ymin"],
                        newbnds["ymax"],
                    )
                    print(
                        "Running model for new bounds that are fully covered"
                        " by input layer: %s" % pnt
                    )
                    bounds2 = newbnds
            else:
                bounds2 = bounds

            modelfunc = conf[modelname]["funcname"]
            if modelfunc == "LogisticModel":
                t1 = time.time()
                lm = LM.LogisticModel(
                    shakefile,
                    conf,
                    uncertfile=uncertfile,
                    saveinputs=args.save_inputs,
                    bounds=bounds2,
                    trimfile=trimfile,
                )

                maplayers = lm.calculate()
                t2 = time.time()
                logging.info(f"{modelname} Elapsed: {t2-t1:.1f} seconds")

            elif modelfunc == "LogBase":
                # newer object oriented approach to logistic models
                model_class = MODEL_FACTORY[modelname]
                t1 = time.time()
                samplebounds = bounds2
                if args.keep_shakemap_bounds:
                    samplebounds = None
                logging.info(f"Constructing model {modelname}...")
                # try:
                model = model_class(
                    shakefile,
                    conf[modelname],
                    uncertfile=uncertfile,
                    bounds=samplebounds,
                    trimfile=trimfile,
                    saveinputs=args.save_inputs,
                )
                # except Exception as e:
                #     if str(e).find("Input geodict has no overlap"):
                #         msg = (
                #             "One or more of the input layers "
                #             f"for {modelname} has NO overlap with "
                #             "input ShakeMap. Skipping."
                #         )
                #         logging.info(msg)
                #     else:
                #         msg = f"Unhandled exception for {modelname}: '{str(e)}'."
                #         logging.warning(msg)
                #     continue
                logging.info(f"Calling {modelname} calculate...")
                maplayers = model.calculate()
                t2 = time.time()
                logging.info(f"{modelname} Elapsed: {t2-t1:.1f} seconds")
            elif modelfunc == "godt2008":
                maplayers = godt2008(
                    shakefile,
                    conf,
                    uncertfile=uncertfile,
                    saveinputs=args.save_inputs,
                    bounds=bounds2,
                    trimfile=trimfile,
                )
            elif modelfunc == "godt2008_2":
                maplayers = godt2008_2(
                    shakefile,
                    conf,
                    uncertfile=uncertfile,
                    saveinputs=args.save_inputs,
                    bounds=bounds2,
                    trimfile=trimfile,
                )
            else:
                print(
                    "Unknown model function specified in config for %s "
                    "model, skipping to next config" % modelfunc
                )
                continue

            if args.appendname is not None:
                filename = "%s_%s_%s" % (eventid, modelname, args.appendname)
            else:
                filename = "%s_%s" % (eventid, modelname)
            if hdf5:
                filenameh = filename + ".hdf5"
                if os.path.exists(filenameh):
                    os.remove(filenameh)
                savelayers(maplayers, os.path.join(outfolder, filenameh))
                filenames.append(filenameh)

            if gis or kmz:
                for key in maplayers:
                    # Rename 'std' key to 'beta_sigma'
                    if key == "std":
                        key_label = "beta_sigma"
                    else:
                        key_label = key
                    if gis:
                        filen = os.path.join(
                            outfolder, "%s_%s.bil" % (filename, key_label)
                        )
                        fileh = os.path.join(
                            outfolder, "%s_%s.hdr" % (filename, key_label)
                        )
                        fileg = os.path.join(
                            outfolder, "%s_%s.tif" % (filename, key_label)
                        )

                        GDALGrid.copyFromGrid(maplayers[key]["grid"]).save(filen)
                        cflags = "-co COMPRESS=DEFLATE -co predictor=2"
                        srs = "-a_srs EPSG:4326"
                        cmd = "gdal_translate %s %s -of GTiff %s %s" % (
                            srs,
                            cflags,
                            filen,
                            fileg,
                        )
                        rc, so, se = get_command_output(cmd)
                        # Delete bil file and its header
                        os.remove(filen)
                        os.remove(fileh)
                        filenames.append(fileg)
                    if kmz and (
                        not key.startswith("quantile") and not key.startswith("std")
                    ):
                        (
                            plotorder,
                            logscale,
                            lims,
                            colormaps,
                            maskthresh,
                        ) = parseConfigLayers(maplayers, conf, keys=["model"])
                        maxprob = np.nanmax(maplayers[key]["grid"].getData())
                        if key == "model":
                            qdict = {
                                k: maplayers[k]
                                for k in maplayers.keys()
                                if k.startswith("quantile")
                            }
                        else:
                            qdict = None
                        if maskthresh is None:
                            maskthresh = [0.0]
                        if maxprob >= maskthresh[0]:
                            filen = os.path.join(
                                outfolder, "%s_%s.kmz" % (filename, key_label)
                            )
                            filek = create_kmz(
                                maplayers[key],
                                filen,
                                mask=maskthresh[0],
                                levels=lims[0],
                                qdict=qdict,
                            )
                            filenames.append(filek)
                        else:
                            print(
                                "No unmasked pixels present, skipping kmz "
                                "file creation"
                            )

            if args.make_webpage:
                # Compile into list of results for later
                results.append(maplayers)

                #  # Make binary output for ShakeCast
                #  filef = os.path.join(outfolder, '%s_model.flt'
                #                       % filename)
                #  # And get name of header
                #  filefh = os.path.join(outfolder, '%s_model.hdr'
                #                        % filename)
                #  # Make file
                #  write_floats(filef, maplayers['model']['grid'])
                #  filenames.append(filef)
                #  filenames.append(filefh)

        eventid = getHeaderData(shakefile)[0]["event_id"]
        if not hasattr(args, "eventsource"):
            args.eventsource = "us"
        if not hasattr(args, "eventsourcecode"):
            args.eventsourcecode = eventid

        if args.make_webpage:
            if len(results) == 0:
                raise Exception("No models were run. Cannot make webpages.")
            outputs = hazdev(
                results,
                configs,
                shakefile,
                outfolder=outfolder,
                pop_file=args.popfile,
                pager_alert=args.property_alertlevel,
                eventsource=args.eventsource,
                eventsourcecode=args.eventsourcecode,
                point=point,
                gf_version=args.gf_version,
                pdlcall=args.pdlcall,
            )
            filenames = filenames + outputs

        #        # create transparent png file
        #        outputs = create_png(outdir)
        #        filenames = filenames + outputs
        #
        #        # create info file
        #        infofile = create_info(outdir)
        #        filenames = filenames + infofile

        print("\nFiles created:\n")
        for filen in filenames:
            print("%s" % filen)

        return filenames


def getGridURL(gridurl, fname=None):
    """
    Args:
        gridurl (str): url for Shakemap grid.xml file.
        fname (str): file location name, if None, will create a temporary file

    Returns:
        file object corresponding to the url.
    """

    with urllib.request.urlopen(gridurl) as fh:
        data = fh.read().decode("utf-8")
        if fname is None:
            with tempfile.NamedTemporaryFile(delete=False, mode="w") as f:
                f.write(data)
        else:
            with open(fname, "w") as f:
                f.write(data)

    return f.name


def getShakefiles(event, outdir, uncert=False, version=None, source="preferred"):
    """
    Download the shakemap grid.xml file and the

    Args:
        event event id or URL
    """
    shakefile = os.path.join(outdir, "grid.xml")
    if uncert:
        uncertfile = os.path.join(outdir, "uncertainty.xml")
    else:
        uncertfile = None
    if version is not None:
        includeSuperseded = True
    else:
        includeSuperseded = False
    # If args.event is a url to a shakemap, download from that url
    if isURL(event):
        if version is not None or source != "preferred":
            raise Exception(
                "Cannot set shakemap version or source when URL "
                "of gridfile is provided"
            )
        try:
            shakefile = getGridURL(event, shakefile)
        except Exception as e:
            raise Exception(
                "Could not download shakemap file from provided URL: %s" % e
            )
        # Now get corresponding event detail
        event = getHeaderData(shakefile)[0]["event_id"]
        version = getHeaderData(shakefile)[0]["shakemap_version"]
        source = getHeaderData(shakefile)[0]["shakemap_originator"]
        try:
            detail = get_event_by_id(event, includesuperseded=includeSuperseded)
        except BaseException:
            # Maybe originator is missing from event id, try another way

            try:
                temp = getHeaderData(shakefile)[0]
                temp2 = "%s%s" % (temp["shakemap_originator"], temp["shakemap_id"])
                detail = get_event_by_id(temp2, includesuperseded=includeSuperseded)
                event = temp2
            except Exception as e:
                msg = "Could not get event detail for shakemap at provided URL: %s"
                print(msg % e)

    else:
        detail = get_event_by_id(event, includesuperseded=includeSuperseded)

    # Get most recent version
    if version is None:  # Get current preferred
        shakemap = detail.getProducts("shakemap", source=source)[0]
        shakemap.getContent("grid.xml", shakefile)
    # or get version requested
    else:
        allversions = detail.getProducts("shakemap", version="all", source=source)
        # First try with properties, more reliable
        vers = []
        for allv in allversions:
            if "version" in allv.properties:
                vers.append(int(allv["version"]))
            else:
                vers.append(-999)
        idx = np.where(np.array(vers) == version)[0]
        if len(idx) < 1:
            # Try using libcomcat version, less reliable...
            vers = [allv.version for allv in allversions]
            idx = np.where(np.array(vers) == version)[0]
            if len(idx) == 1:
                # Check info.json to make sure it's right version
                infobytes, url = allversions[idx[0]].getContentBytes("info.json")
                info = json.loads(infobytes.decode("utf-8"))
                if info["processing"]["shakemap_versions"]["map_version"] != version:
                    idx = []
        if len(idx) < 1:
            msg = "Could not find version %d of Shakemap from source %s"
            raise Exception(msg % (version, source))
        if len(idx) > 1:
            msg = "Found more than one ShakeMap with matching source and version. \
Choosing first one."
            warnings.warn(msg)
            print(msg)

        shakemap = allversions[idx[0]]
        shakemap.getContent("grid.xml", shakefile)

    if uncert:
        uncertfile = getUncert(shakemap, uncertfile)

    return detail, shakefile, uncertfile


def getUncert(shakemap, fname=None):
    """
    download and unzip (if needed) the uncertainty grid corresponding to a
    shakemap

    Args:
        shakemap: libcomcat ShakeMap product class for the event and version
        fname (str): file location name, if None, will create a temporary file

    Returns:
        file object corresponding to the url.
    """
    grid_url = shakemap.getContentURL("uncertainty.*")
    if grid_url is None:
        return None

    try:
        ext = grid_url.split(".")[-1]
        if fname is None:
            basedir = tempfile.mkdtemp()
        else:
            basedir = os.path.dirname(fname)
        uncertfile = os.path.join(basedir, "uncertainty.xml")
        if ext == "xml":
            shakemap.getContent("uncertainty.xml", uncertfile)
        if ext == "zip":
            fname = os.path.join(basedir, "uncertainty.xml.zip")
            shakemap.getContent("uncertainty.xml.zip", fname)
            # urllib.request.urlretrieve(grid_url, fname)
            with ZipFile(fname, "r") as zip1:
                # See if it's inside a file structure
                out = zip1.filelist[0]
                # Extract all the contents of zip file in different directory
                zip1.extractall(basedir)
                # move file uncertainty.xml file to base dir if it was in a
                # weird subdir
                if os.path.isdir(os.path.dirname(out.filename)):
                    os.replace(os.path.join(basedir, out.filename), uncertfile)
    except Exception as e:
        uncertfile = None
        print("Unable to download uncertainty.xml: %s" % e)
    return uncertfile


def isURL(gridurl):
    """
    This function determines if the provided string is a valid url

    Args:
        gridurl (str): url to check.

    Returns:
        bool: True if gridurl is a valid url, False otherwise.
    """

    is_url = False
    try:
        urllib.request.urlopen(gridurl)
        is_url = True
    except BaseException:
        pass
    return is_url


def set_default_paths(args):
    """
    Creates a file called .gfail_defaults that contains default path
    information to simplify running gfail. Can be overwritten by any manually
    entered paths. This updates any existing .gfail_defaults file. If
    args.data_path is 'reset' then any existing defaults will be removed.

    Args:
        args (arparser Namespace): Input arguments.

    Returns:
        Updates .gfail_defaults file on users path, or creates new one if
        file does not already exist.
    """
    filename = os.path.join(os.path.expanduser("~"), ".gfail_defaults")
    if os.path.exists(filename):
        D = ConfigObj(filename)
    else:
        D = {}
    if args.data_path is not None:
        if args.data_path == "reset":
            D.pop("data_path")
        else:
            # check that it's a valid path
            if os.path.exists(args.data_path):
                D.update({"data_path": args.data_path})
            else:
                print("Path given for data_path does not exist: %s" % args.data_path)
    if args.output_filepath is not None:
        if args.output_filepath == "reset":
            D.pop("output_filepath")
        else:
            # check that it's a valid path
            if os.path.exists(args.output_filepath):
                D.update({"output_filepath": args.output_filepath})
            else:
                print(
                    "Path given for output_filepath does not exist: %s"
                    % args.output_filepath
                )
    if args.config_filepath is not None:
        if args.config_filepath == "reset":
            D.pop("config_filepath")
        else:
            # check that it's a valid path
            if os.path.exists(args.config_filepath):
                D.update({"config_filepath": args.config_filepath})
            else:
                print(
                    "Path given for config_filepath does not exist: %s"
                    % args.config_filepath
                )
    if args.popfile is not None:
        if args.popfile == "reset":
            D.pop("popfile")
        else:
            # check that it's a valid path
            if os.path.exists(args.popfile):
                D.update({"popfile": args.popfile})
            else:
                print(
                    "Path given for population file does not exist: %s" % args.popfile
                )
    if args.trimfile is not None:
        if args.trimfile == "reset":
            D.pop("trim")
        else:
            # check that it's a valid path and that it's a shapefile
            if os.path.exists(args.trimfile):
                filename4, fileextension = os.path.splitext(args.trimfile)
                if fileextension == ".shp":
                    D.update({"trimfile": args.trimfile})
                else:
                    print("Ocean trimming file is not a shapefile: %s" % args.trimfile)
            else:
                print(
                    "Path given for ocean trimming file does not exist: %s"
                    % args.trimfile
                )
    if args.pdl_config is not None:
        if args.pdl_config == "reset":
            D.pop("pdl_config")
        else:
            # check that it's a valid path
            if os.path.exists(args.pdl_config):
                D.update({"pdl_config": args.pdl_config})
            else:
                print(
                    "Path given for pdl config file does not exist: %s"
                    % args.pdl_config
                )
    if args.log_filepath is not None:
        if args.log_filepath == "reset":
            D.pop("log_filepath")
        else:
            # check that it's a valid path
            if os.path.exists(args.log_filepath):
                D.update({"log_filepath": args.log_filepath})
            else:
                print("Path given for log file does not exist: %s" % args.log_filepath)
    if args.dbfile is not None:
        if args.dbfile == "reset":
            D.pop("dbfile")
        else:
            # check that it's a valid path (file itself doesnt have to exist)
            if os.path.exists(os.path.dirname(args.dbfile)):
                D.update({"dbfile": args.dbfile})
            else:
                print("Path given for database file does not exist: %s" % args.dbfile)

    print("New default paths set.\n")

    if D:
        C = ConfigObj(D)
        C.filename = filename
        C.write()
        list_default_paths()
    else:
        print("no defaults set because no paths were input\n")


def list_default_paths():
    """
    Lists all default paths currently set.
    """
    filename = os.path.join(os.path.expanduser("~"), ".gfail_defaults")
    if os.path.exists(filename):
        D = ConfigObj(filename)
        print("Default paths currently set to:\n")
        for key in D:
            print("\t%s = %s" % (key, D[key]))
    else:
        print("No default paths currently set\n")


def reset_default_paths():
    """
    Clear default path file
    """
    filename = os.path.join(os.path.expanduser("~"), ".gfail_defaults")
    if os.path.exists(filename):
        os.remove(filename)
        print("Default paths cleared\n")
    else:
        print("No default paths currently set\n")


def get_bounds(shakefile, parameter="pga", threshold=2.0):
    """
    Get the boundaries of the shakemap that include all areas with shaking
    above the defined threshold.

    Args:
        shakefile (str): Path to shakemap file.
        parameter (str): Either 'pga' or 'pgv'.
        threshold (float): Minimum value of parameter of interest, in units
            of %g for pga and cm/s for pgv. The default value of 2% g is based
            on minimum pga threshold ever observed to have triggered landslides
            by Jibson and Harp (2016).

    Returns:
        dict: A dictionary with keys 'xmin', 'xmax', 'ymin', and 'ymax' that
        defines the boundaries in geographic coordinates.
    """
    shakemap = ShakeGrid.load(shakefile, adjust="res")
    if parameter == "pga":
        vals = shakemap.getLayer("pga")
    elif parameter == "pgv":
        vals = shakemap.getLayer("pgv")
    else:
        raise Exception("parameter not valid")
    xmin, xmax, ymin, ymax = vals.getBounds()
    lons = np.linspace(xmin, xmax, vals.getGeoDict().nx)
    lats = np.linspace(ymax, ymin, vals.getGeoDict().ny)
    row, col = np.where(vals.getData() > float(threshold))
    lonmin = lons[col].min()
    lonmax = lons[col].max()
    latmin = lats[row].min()
    latmax = lats[row].max()

    # dummy fillers, only really care about bounds
    boundaries1 = {"dx": 100, "dy": 100.0, "nx": 100.0, "ny": 100}

    if xmin < lonmin:
        boundaries1["xmin"] = lonmin
    else:
        boundaries1["xmin"] = xmin
    if xmax > lonmax:
        boundaries1["xmax"] = lonmax
    else:
        boundaries1["xmax"] = xmax
    if ymin < latmin:
        boundaries1["ymin"] = latmin
    else:
        boundaries1["ymin"] = ymin
    if ymax > latmax:
        boundaries1["ymax"] = latmax
    else:
        boundaries1["ymax"] = ymax

    return boundaries1


def check_input_extents(config, shakefile=None, bounds=None):
    """Make sure all input files exist and cover the extent desired

    Args:
        config: configObj of a single model
        shakefile: path to ShakeMap grid.xml file (used for bounds). If not
            provided, bounds must be provided
        bounds: dictionary of bounds with keys: 'xmin', 'xmax', 'ymin', 'ymax'

    Returns:
        tuple containing:
            notcovered: list of files that do not cover the entire area
                defined by bounds or shakefile
            newbounds: new dictionary of bounds of subarea of original
                bounds or shakefile extent that is covered by all input files
    """
    if shakefile is None and bounds is None:
        raise Exception("Must define either a shakemap file or bounds")
    modelname = config.keys()[0]
    # Make dummy geodict to use
    if bounds is None:
        evdict = ShakeGrid.getFileGeoDict(shakefile)
    else:
        evdict = GeoDict.createDictFromBox(
            bounds["xmin"],
            bounds["xmax"],
            bounds["ymin"],
            bounds["ymax"],
            0.00001,
            0.00001,
            inside=False,
        )

    # Check extents of all input layers
    notcovered = []
    notcovgdicts = []
    newbounds = None
    for item, value in config[modelname]["layers"].items():
        if "file" in value.keys():
            filelook = value["file"]
            if getFileType(filelook) == "gmt":
                tmpgd, _ = GMTGrid.getFileGeoDict(filelook)
            else:
                tmpgd, _ = GDALGrid.getFileGeoDict(filelook)
            # See if tempgd contains evdict
            contains = tmpgd.contains(evdict)
            if not contains:
                notcovered.append(filelook)
                notcovgdicts.append(tmpgd)
                # print(filelook)
    if len(notcovered) > 0:
        # Figure out what bounds COULD be run
        xmins = [gd.xmin for gd in notcovgdicts]
        xmaxs = [gd.xmax for gd in notcovgdicts]
        ymins = [gd.ymin for gd in notcovgdicts]
        ymaxs = [gd.ymax for gd in notcovgdicts]

        # Set in by a buffer of 0.05 degrees because mapio doesn't like
        # when bounds are exactly the same for getboundswithin
        newbounds = dict(
            xmin=evdict.xmin + 0.05,
            xmax=evdict.xmax - 0.05,
            ymin=evdict.ymin + 0.05,
            ymax=evdict.ymax - 0.05,
        )
        # Which one is the problem?
        if evdict.xmin < np.max(xmins):
            newbounds["xmin"] = np.max(xmins) + 0.05
        if evdict.xmax > np.min(xmaxs):
            newbounds["xmax"] = np.min(xmaxs) - 0.05
        if evdict.ymin < np.max(ymins):
            newbounds["ymin"] = np.max(ymins) + 0.05
        if evdict.ymax > np.min(ymaxs):
            newbounds["ymax"] = np.min(ymaxs) - 0.05

        # See if this is a possible extent
        try:
            GeoDict.createDictFromBox(
                newbounds["xmin"],
                newbounds["xmax"],
                newbounds["ymin"],
                newbounds["ymax"],
                0.00001,
                0.00001,
                inside=False,
            )
        except BaseException:
            print("Cannot make new bounds that will work")
            newbounds = None

    return notcovered, newbounds
