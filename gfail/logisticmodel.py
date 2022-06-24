#!/usr/bin/env python
"""
This module contains functions and class definitions for running forward
models of models based on logistic regression.
"""

# stdlib imports
import numpy as np
import os.path
import re
import collections
import shutil
import tempfile
from timeit import default_timer as timer
import logging

# third party imports
from mapio.shake import ShakeGrid
from mapio.shake import getHeaderData
from mapio.gmt import GMTGrid
from mapio.gdal import GDALGrid
from mapio.grid2d import Grid2D
from mapio.geodict import GeoDict


# local imports
from gfail.temphdf import TempHdf
from gfail.spatial import quickcut, trim_ocean
from gfail.utilities import getFileType
from gfail.stats import get_rangebeta

# temporary until mapio is updated
import warnings

warnings.filterwarnings("ignore")


PARAM_PATTERN = "b[0-9]+"
LAYER_PATTERN = "_layer"
TERM_PATTERN = "term"

SM_TERMS = ["MW", "YEAR", "MONTH", "DAY", "HOUR", "pga", "pgv", "mmi"]
SM_GRID_TERMS = ["pga", "pgv", "mmi"]
# these will get np. prepended
OPERATORS = ["log", "log10", "arctan", "power", "sqrt", "minimum", "pi"]
FLOATPAT = r"[+-]?(?=\d*[.eE])(?=\.?\d)\d*\.?\d*(?:[eE][+-]?\d+)?"
INTPAT = "[0-9]+"
OPERATORPAT = r"[\+\-\*\/]*"
MONTHS = [
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec",
]


class LogisticModel(object):
    def __init__(
        self,
        shakefile,
        config,
        uncertfile=None,
        saveinputs=False,
        slopefile=None,
        bounds=None,
        slopemod=None,
        trimfile=None,
    ):
        """
        Sets up the logistic model

        Args:
            shakefile (str): Path to shakemap grid.xml file for the event.
            config: configobj object defining the model and its inputs. Only
                one model should be described in each config file.
            uncertfile (str): Path to uncertainty.xml file.
            saveinputs (bool): Save input layers as Grid2D objects in addition
                to the model? If false (the default), it will just output the
                model.
            slopefile (str): Optional path to slopefile that will be resampled
                to the other input files for applying thresholds. OVERWRITES
                VALUE IN CONFIG.
            bounds (dict): Default of None uses ShakeMap boundaries, otherwise
                a dictionary of boundaries to cut to like

                .. code-block:: python

                    bounds = {
                        'xmin': lonmin, 'xmax': lonmax,
                        'ymin': latmin, 'ymax': latmax
                    }
            slopemod (str): How slope input should be modified to be in
                degrees: e.g., ``np.arctan(slope) * 180. / np.pi`` or
                ``slope/100.`` (note that this may be in the config file
                already).
            trimfile (str): shapefile of earth's landmasses to use to cut
                offshore areas.
        """
        mnames = getLogisticModelNames(config)
        if len(mnames) == 0:
            raise Exception("No config file found or problem with config file format")
        if len(mnames) > 1:
            raise Exception(
                "Config file contains more than one model which "
                "is no longer allowed, update your config file "
                "to the newer format"
            )

        self.model = mnames[0]
        self.config = config
        cmodel = config[self.model]
        self.modeltype = cmodel["gfetype"]
        self.coeffs = validateCoefficients(cmodel)
        # key = layer name, value = file name
        self.layers = validateLayers(cmodel)
        self.terms, timeField = validateTerms(cmodel, self.coeffs, self.layers)
        self.interpolations = validateInterpolations(cmodel, self.layers)
        self.units = validateUnits(cmodel)
        self.gmused = [
            value
            for term, value in cmodel["terms"].items()
            if "pga" in value.lower()
            or "pgv" in value.lower()
            or "mmi" in value.lower()
        ]
        self.modelrefs, self.longrefs, self.shortrefs = validateRefs(cmodel)
        # self.numstd = numstd
        self.clips = validateClips(cmodel, self.layers, self.gmused)
        self.notes = ""

        if cmodel["baselayer"] not in list(self.layers.keys()):
            raise Exception(
                "You must specify a base layer corresponding to "
                "one of the files in the layer section."
            )
        self.saveinputs = saveinputs
        if slopefile is None:
            try:
                self.slopefile = cmodel["slopefile"]
            except BaseException:
                # print('Slopefile not specified in config, no slope '
                #      'thresholds will be applied\n')
                self.slopefile = None
        else:
            self.slopefile = slopefile
        if slopemod is None:
            try:
                self.slopemod = cmodel["slopemod"]
            except BaseException:
                self.slopemod = None

        # See if trimfile exists
        if trimfile is not None:
            if not os.path.exists(trimfile):
                print(
                    "trimfile defined does not exist: %s\nOcean will not be "
                    "trimmed" % trimfile
                )
                self.trimfile = None
            elif os.path.splitext(trimfile)[1] != ".shp":
                print("trimfile must be a shapefile, ocean will not be trimmed")
                self.trimfile = None
            else:
                self.trimfile = trimfile
        else:
            self.trimfile = None

        # Get month of event
        griddict, eventdict, specdict, fields, uncertainties = getHeaderData(shakefile)
        MONTH = MONTHS[eventdict["event_timestamp"].month - 1]

        # Figure out how/if need to cut anything
        geodict = ShakeGrid.getFileGeoDict(shakefile, adjust="res")
        if bounds is not None:  # Make sure bounds are within ShakeMap Grid
            if geodict.xmin < geodict.xmax:  # only if signs are not opposite
                if (
                    geodict.xmin > bounds["xmin"]
                    or geodict.xmax < bounds["xmax"]
                    or geodict.ymin > bounds["ymin"]
                    or geodict.ymax < bounds["ymax"]
                ):
                    print(
                        "Specified bounds are outside shakemap area, using "
                        "ShakeMap bounds instead."
                    )
                    bounds = None

        if bounds is not None:
            tempgdict = GeoDict.createDictFromBox(
                bounds["xmin"],
                bounds["xmax"],
                bounds["ymin"],
                bounds["ymax"],
                geodict.dx,
                geodict.dy,
                inside=False,
            )
            # If Shakemap geodict crosses 180/-180 line, fix geodict so
            # things don't break
            if geodict.xmin > geodict.xmax:
                if tempgdict.xmin < 0:
                    geodict._xmin -= 360.0
                else:
                    geodict._xmax += 360.0
            gdict = geodict.getBoundsWithin(tempgdict)
        else:
            gdict = geodict

        # Now find the layer that is our base layer and get the largest bounds
        # we can guarantee not to exceed shakemap bounds
        basefile = self.layers[cmodel["baselayer"]]
        ftype = getFileType(basefile)
        if ftype == "esri":
            basegeodict, firstcol = GDALGrid.getFileGeoDict(basefile)
            if basegeodict == gdict:
                sampledict = gdict
            else:
                sampledict = basegeodict.getBoundsWithin(gdict)
        elif ftype == "gmt":
            basegeodict, firstcol = GMTGrid.getFileGeoDict(basefile)
            if basegeodict == gdict:
                sampledict = gdict
            else:
                sampledict = basegeodict.getBoundsWithin(gdict)
        else:
            raise Exception(
                "All predictor variable grids must be a valid GMT or ESRI file type."
            )

        # Do we need to subdivide baselayer?
        if "divfactor" in self.config[self.model].keys():
            divfactor = float(self.config[self.model]["divfactor"])
            if divfactor != 1.0:
                # adjust sampledict so everything will be resampled
                newxmin = (
                    sampledict.xmin
                    - sampledict.dx / 2.0
                    + sampledict.dx / (2.0 * divfactor)
                )
                newymin = (
                    sampledict.ymin
                    - sampledict.dy / 2.0
                    + sampledict.dy / (2.0 * divfactor)
                )
                newxmax = (
                    sampledict.xmax
                    + sampledict.dx / 2.0
                    - sampledict.dx / (2.0 * divfactor)
                )
                newymax = (
                    sampledict.ymax
                    + sampledict.dy / 2.0
                    - sampledict.dy / (2.0 * divfactor)
                )
                newdx = sampledict.dx / divfactor
                newdy = sampledict.dy / divfactor
                if np.abs(newxmax) > 180.0:
                    newxmax = np.sign(newxmax) * 180.0
                if np.abs(newxmin) > 180.0:
                    newxmin = np.sign(newxmin) * 180.0

                sampledict = GeoDict.createDictFromBox(
                    newxmin, newxmax, newymin, newymax, newdx, newdy, inside=True
                )

        # Find slope thresholds, if applicable
        self.slopemin = "none"
        self.slopemax = "none"
        if self.slopefile is not None:
            try:
                self.slopemin = float(config[self.model]["slopemin"])
                self.slopemax = float(config[self.model]["slopemax"])
            except BaseException:
                print(
                    "Could not find slopemin and/or slopemax in config, "
                    "limits. No slope thresholds will be applied."
                )
                self.slopemin = "none"
                self.slopemax = "none"

        # Make temporary directory for hdf5 pytables file storage
        self.tempdir = tempfile.mkdtemp()

        # now load the shakemap, resampling and padding if necessary
        temp = ShakeGrid.load(shakefile)  # , adjust='res')
        self.shakedict = temp.getShakeDict()
        self.eventdict = temp.getEventDict()
        self.shakemap = {}

        # Read both PGA and PGV in, may need them for thresholds
        for gm in ["pga", "pgv"]:
            junkfile = os.path.join(self.tempdir, "temp.bil")
            GDALGrid.copyFromGrid(temp.getLayer(gm)).save(junkfile)
            if gm in self.interpolations.keys():
                intermeth = self.interpolations[gm]
            else:
                intermeth = "bilinear"
            junkgrid = quickcut(
                junkfile, sampledict, precise=True, method=intermeth, override=True
            )
            if gm in self.clips:
                junkgrid.setData(
                    np.clip(junkgrid.getData(), self.clips[gm][0], self.clips[gm][1])
                )
            self.shakemap[gm] = TempHdf(
                junkgrid, os.path.join(self.tempdir, "%s.hdf5" % gm)
            )
            os.remove(junkfile)
        del temp

        # get updated geodict
        sampledict = junkgrid.getGeoDict()

        # take uncertainties into account, if available
        if uncertfile is not None:
            self.uncert = {}
            # try:
            # Only read in the ones that will be needed
            temp = ShakeGrid.load(uncertfile)
            already = []
            for gm in self.gmused:
                if "pgv" in gm:
                    gmsimp = "pgv"
                elif "pga" in gm:
                    gmsimp = "pga"
                elif "mmi" in gm:
                    gmsimp = "mmi"
                if gmsimp in already:
                    continue
                junkfile = os.path.join(self.tempdir, "temp.bil")
                GDALGrid.copyFromGrid(temp.getLayer("std%s" % gmsimp)).save(junkfile)
                if gmsimp in self.interpolations.keys():
                    intermeth = self.interpolations[gmsimp]
                else:
                    intermeth = "bilinear"
                junkgrid = quickcut(
                    junkfile, sampledict, precise=True, method=intermeth, override=True
                )
                if gmsimp in self.clips:
                    junkgrid.setData(
                        np.clip(
                            junkgrid.getData(),
                            self.clips[gmsimp][0],
                            self.clips[gmsimp][1],
                        )
                    )
                self.uncert["std" + gmsimp] = TempHdf(
                    junkgrid, os.path.join(self.tempdir, "std%s.hdf5" % gmsimp)
                )
                already.append(gmsimp)
                os.remove(junkfile)
            del temp
            # except:
            # print('Could not read uncertainty file, ignoring '
            #       'uncertainties')
            # self.uncert = None
        else:
            self.uncert = None

        # Load the predictor layers, save as hdf5 temporary files, put file
        # locations into a dictionary.

        # Will be replaced in the next section if a slopefile was defined
        self.nonzero = None

        # key = layer name, value = grid object
        self.layerdict = {}

        didslope = False
        for layername, layerfile in self.layers.items():
            start = timer()
            if isinstance(layerfile, list):
                for lfile in layerfile:
                    if timeField == "MONTH":
                        if lfile.find(MONTH) > -1:
                            layerfile = lfile
                            # ftype = getFileType(layerfile)
                            interp = self.interpolations[layername]
                            temp = quickcut(
                                layerfile, sampledict, precise=True, method=interp
                            )
                            if layername in self.clips:
                                temp.setData(
                                    np.clip(
                                        temp.getData(),
                                        self.clips[layername][0],
                                        self.clips[layername][1],
                                    )
                                )
                            self.layerdict[layername] = TempHdf(
                                temp, os.path.join(self.tempdir, "%s.hdf5" % layername)
                            )
                            del temp
            else:
                interp = self.interpolations[layername]
                temp = quickcut(layerfile, sampledict, precise=True, method=interp)
                if layername in self.clips:
                    temp.setData(
                        np.clip(
                            temp.getData(),
                            self.clips[layername][0],
                            self.clips[layername][1],
                        )
                    )
                # Convert unconsolidated sediments to more reasonable coeff
                if layername == "rock":
                    sub1 = temp.getData()
                    # Change to mixed sed rock coeff
                    sub1[sub1 <= -3.21] = -1.36
                    temp.setData(sub1)
                    self.notes += (
                        "unconsolidated sediment coefficient "
                        "changed to -1.36 (weaker) from -3.22 to "
                        "better reflect that this "
                        "unit is not actually strong\n"
                    )
                self.layerdict[layername] = TempHdf(
                    temp, os.path.join(self.tempdir, "%s.hdf5" % layername)
                )
                td = temp.getGeoDict()
                if td != sampledict:
                    raise Exception("Geodictionaries of resampled files do not match")

                if layerfile == self.slopefile:
                    flag = 0
                    if self.slopemin == "none" and self.slopemax == "none":
                        flag = 1
                    if self.slopemod is None:
                        slope1 = temp.getData().astype(float)
                        slope = 0
                    else:
                        try:
                            slope = temp.getData().astype(float)
                            slope1 = eval(self.slopemod)
                        except BaseException:
                            print(
                                "slopemod provided not valid, continuing "
                                "without slope thresholds."
                            )
                            flag = 1
                    if flag == 0:
                        nonzero = np.array(
                            [(slope1 > self.slopemin) & (slope1 <= self.slopemax)]
                        )
                        self.nonzero = nonzero[0, :, :]
                        del slope1
                        del slope
                    else:
                        # Still remove areas where the slope equals exactly
                        # 0.0 to remove offshore liq areas.
                        nonzero = np.array([slope1 != 0.0])
                        self.nonzero = nonzero[0, :, :]
                        del slope1
                    didslope = True
                del temp

            logging.info(f"Loading {layername} layer: {timer() - start:1.2f} sec")

        if didslope is False and self.slopefile is not None:
            # Slope didn't get read in yet
            temp = quickcut(self.slopefile, sampledict, precise=True, method="bilinear")
            flag = 0
            if self.slopemin == "none" and self.slopemax == "none":
                flag = 1
            if self.slopemod is None:
                slope1 = temp.getData().astype(float)
                slope = 0
            else:
                try:
                    slope = temp.getData().astype(float)
                    slope1 = eval(self.slopemod)
                except BaseException:
                    print(
                        "slopemod provided not valid, continuing without "
                        "slope thresholds"
                    )
                    flag = 1
            if flag == 0:
                nonzero = np.array(
                    [(slope1 > self.slopemin) & (slope1 <= self.slopemax)]
                )
                self.nonzero = nonzero[0, :, :]
                del slope1
                del slope
            else:
                # Still remove areas where the slope equals exactly
                # 0.0 to remove offshore liq areas.
                nonzero = np.array([slope1 != 0.0])
                self.nonzero = nonzero[0, :, :]
                del slope1

        self.nuggets = [str(self.coeffs["b0"])]

        ckeys = sorted(self.terms.keys())
        for key in ckeys:
            term = self.terms[key]
            coeff = self.coeffs[key]
            self.nuggets.append("(%g * %s)" % (coeff, term))

        self.equation = " + ".join(self.nuggets)
        self.geodict = sampledict

    def getEquations(self):
        """
        Method for LogisticModel class to extract strings defining the
        equations for the model for median ground motions.

        Returns:
            equation: the equation for median ground motions,

        """
        return self.equation

    def getGeoDict(self):
        """
        Returns the geodictionary of the LogisticModel class defining bounds
        and resolution of model inputs and outputs.

        Returns:
            geodict: mapio geodict object
        """
        return self.geodict

    def calculate(self, cleanup=True, rowmax=2000, colmax=2000):
        """
        Calculate the model.

        Args:
            cleanup (bool): If True, delete temporary hdf5 files
            rowmax (int): Number of rows to compute at once; If None, all rows
                will be computed at once.
            colmax (int): Number of columns to compute at once; If None, all
                columns will be computed at once.
        Returns:
            dict: Dictionary containing the model results (and model inputs if
            saveinputs was set to True). See
            `the description <https://github.com/usgs/groundfailure#api-for-model-output>`_
            of the structure.
        """
        start_calculate = timer()
        tk = list(self.shakemap.keys())[0]
        # Figure out what slices to do
        rowstarts, rowends, colstarts, colends = self.shakemap[tk].getSliceDiv(
            rowmax, colmax
        )

        # Make empty matrix to fill
        X = np.empty([self.geodict.ny, self.geodict.nx])

        ranges = {
            "b1": (-2, 5),
            "b2": (0, 2.5),
            "b3": (-2, -0.6),
            "b4": (0.05, 0.45),
            "b5": (-1.2, 2.2),
            "b6": (-0.6, 1.2),
        }

        # Loop through slices, appending output each time
        for rowstart, rowend, colstart, colend in zip(
            rowstarts, rowends, colstarts, colends
        ):
            X[rowstart:rowend, colstart:colend] = eval(self.equation)
            for key, term in self.terms.items():
                coeff = self.coeffs[key]
                vmin, vmax = ranges[key]
                layer = eval(term) * coeff

        P = 1 / (1 + np.exp(-X))

        if "vs30max" in self.config[self.model].keys():
            vs30 = self.layerdict["vs30"].getSlice(None, None, None, None, name="vs30")
            P[vs30 > float(self.config[self.model]["vs30max"])] = 0.0

        if "minpgv" in self.config[self.model].keys():
            pgv = self.shakemap["pgv"].getSlice(None, None, None, None, name="pgv")
            P[pgv < float(self.config[self.model]["minpgv"])] = 0.0

        if "minpga" in self.config[self.model].keys():
            pga = self.shakemap["pga"].getSlice(None, None, None, None, name="pga")
            P[pga < float(self.config[self.model]["minpga"])] = 0.0

        if self.uncert is not None:  # hard code for now
            if "Zhu and others (2017)" in self.modelrefs["shortref"]:
                if "stddev" in self.layerdict.keys():
                    stdX = self.layerdict["stddev"].getSlice()
                else:
                    stdX = float(self.config[self.model]["default_stddev"])
                varX = stdX**2.0 + (
                    self.coeffs["b1"] ** 2.0 * self.uncert["stdpgv"].getSlice() ** 2.0
                )
                varP = (np.exp(-X) / (np.exp(-X) + 1) ** 2.0) ** 2.0 * varX
                if "coverage" in self.config[self.model].keys():
                    a = 0.4915
                    b = 42.4
                    c = 9.165
                    # ((2*a*b*c*np.exp(2*c*P))/(b+np.exp(c*P))**3.)**2.*varP
                    varL = (
                        (2 * a * b * c * np.exp(-c * P))
                        / ((1 + b * np.exp(-c * P)) ** 3.0)
                    ) ** 2.0 * varP
                    std1 = np.sqrt(varL)
                else:
                    std1 = np.sqrt(varP)
            elif "Jessee" in self.modelrefs["shortref"]:
                if "stddev" in self.layerdict.keys():
                    stdX = self.layerdict["stddev"].getSlice()
                else:
                    stdX = float(self.config[self.model]["default_stddev"])
                cfs = self.coeffs
                slp = self.layerdict["slope"]
                std = self.uncert["stdpgv"]
                varX = stdX**2.0 + (
                    (cfs["b1"] + cfs["b6"] * (np.arctan(slp.getSlice()) * 180 / np.pi))
                    ** 2.0
                    * std.getSlice() ** 2.0
                )
                varP = (np.exp(-X) / (np.exp(-X) + 1) ** 2.0) ** 2.0 * varX
                if "coverage" in self.config[self.model].keys():
                    a = -7.592
                    b = 5.237
                    c = -3.042
                    d = 4.035
                    varL = (
                        np.exp(a + b * P + c * P**2.0 + d * P**3.0)
                        * (b + 2.0 * P * c + 3.0 * d * P**2.0)
                    ) ** 2.0 * varP
                    std1 = np.sqrt(varL)
                else:
                    std1 = np.sqrt(varP)
            else:
                print(
                    "cannot do uncertainty for %s model, skipping"
                    % self.modelrefs["shortref"]
                )
                self.uncert = None
                std1 = None
        else:
            std1 = None

        # P needs to be converted to areal coverage AFTER dealing with uncert
        if "coverage" in self.config[self.model].keys():
            eqn = self.config[self.model]["coverage"]["eqn"]
            P = eval(eqn)

        # Compute quantiles
        compute_quantiles = False
        mconf = self.config[self.model]
        if ("conf_int_probabilities" in mconf) and (std1 is not None):
            compute_quantiles = True
            ci_probabilities = [float(cip) for cip in mconf["conf_int_probabilities"]]

        if compute_quantiles:
            quantile_dict = {}
            pmax = float(mconf["maxprob"])
            beta_p = P / pmax * (((pmax * P - P**2) / std1**2) - 1)
            beta_q = (1 - P / pmax) * (((pmax * P - P**2) / std1**2) - 1)
            for ci_prob in ci_probabilities:
                min_quantile = str(np.round(100 * (1.0 - ci_prob) / 2.0, 1))
                max_quantile = str(np.round(100 * (1 - ((1.0 - ci_prob)) / 2.0), 1))
                min_prob, max_prob = get_rangebeta(
                    beta_p, beta_q, ci_prob, minlim=0, maxlim=pmax
                )
                quantile_dict[min_quantile] = min_prob
                quantile_dict[max_quantile] = max_prob

        if self.slopefile is not None and self.nonzero is not None:
            # Apply slope min/max limits
            logging.info("applying slope thresholds")
            P = P * self.nonzero
            if std1 is not None:
                # No uncert for masked values
                std1[P == 0] = 0.0
                if compute_quantiles:
                    for q in quantile_dict.values():
                        q[P == 0] = 0.0

        # Stuff into Grid2D object
        if "Jessee" in self.modelrefs["shortref"]:
            if "coverage" not in self.config[self.model].keys():
                units5 = "Relative Hazard"
            else:
                units5 = "Proportion of area affected"
        elif "Zhu" in self.modelrefs["shortref"]:
            if (
                "coverage" not in self.config[self.model].keys()
                and "2017" in self.modelrefs["shortref"]
            ):
                units5 = "Relative Hazard"
            else:
                units5 = "Proportion of area affected"
        else:
            units5 = "Probability of any occurrence"

        shakedetail = "%s_ver%s" % (
            self.shakedict["shakemap_id"],
            self.shakedict["shakemap_version"],
        )
        description = {
            "name": self.modelrefs["shortref"],
            "longref": self.modelrefs["longref"],
            "units": units5,
            "shakemap": shakedetail,
            "event_id": self.eventdict["event_id"],
            "parameters": {
                "slopemin": self.slopemin,
                "slopemax": self.slopemax,
                "modeltype": self.modeltype,
                "notes": self.notes,
            },
        }
        if "vs30max" in self.config[self.model].keys():
            description["vs30max"] = float(self.config[self.model]["vs30max"])
        if "minpgv" in self.config[self.model].keys():
            description["minpgv"] = float(self.config[self.model]["minpgv"])

        Pgrid = Grid2D(P, self.geodict)
        if self.trimfile is not None:
            # Turn all offshore cells to nan
            Pgrid = trim_ocean(Pgrid, self.trimfile)
        rdict = collections.OrderedDict()
        rdict["model"] = {
            "grid": Pgrid,
            "label": "%s estimate - %s" % (self.modeltype.capitalize(), units5.title()),
            "type": "output",
            "description": description,
        }
        if self.uncert is not None:
            Stdgrid = Grid2D(std1, self.geodict)
            if self.trimfile is not None:
                Stdgrid = trim_ocean(Stdgrid, self.trimfile)
            rdict["std"] = {
                "grid": Stdgrid,
                "label": (
                    "%s estimate - %s (std)"
                    % (self.modeltype.capitalize(), units5.title())
                ),
                "type": "output",
                "description": description,
            }
            if compute_quantiles:
                logging.info("Computing quantiles")
                start_quantiles = timer()
                for quantile, qgrid in quantile_dict.items():
                    Qgrid = Grid2D(qgrid, self.geodict)
                    qname = "quantile%s" % quantile
                    rdict[qname] = {
                        "grid": Qgrid,
                        "label": (
                            "%s %sth percentile - %s"
                            % (self.modeltype.capitalize(), quantile, units5.title())
                        ),
                        "type": "output",
                        "description": description,
                    }
                logging.info(f"Quantiles elapsed: {timer() - start_quantiles:1.2f}")

        # This step might swamp memory for higher resolution runs
        if self.saveinputs is True:
            logging.info("Saving inputs.")
            start_saveinputs = timer()
            for layername, layergrid in list(self.layerdict.items()):
                units = self.units[layername]
                if units is None:
                    units = ""
                rdict[layername] = {
                    "grid": Grid2D(
                        layergrid.getSlice(None, None, None, None, name=layername),
                        self.geodict,
                    ),
                    "label": "%s (%s)" % (layername, units),
                    "type": "input",
                    "description": {
                        "units": units,
                        "name": self.shortrefs[layername],
                        "longref": self.longrefs[layername],
                    },
                }
            for gmused in self.gmused:
                if "pga" in gmused:
                    units = "%g"
                    getkey = "pga"
                elif "pgv" in gmused:
                    units = "cm/s"
                    getkey = "pgv"
                elif "mmi" in gmused:
                    units = "intensity"
                    getkey = "mmi"
                else:
                    continue
                    # Layer is derived from several input layers, skip
                    # outputting this layer

                if getkey in rdict:
                    continue

                layer = self.shakemap[getkey].getSlice(
                    None, None, None, None, name=getkey
                )
                rdict[getkey] = {
                    "grid": Grid2D(layer, self.geodict),
                    "label": "%s (%s)" % (getkey.upper(), units),
                    "type": "input",
                    "description": {"units": units, "shakemap": shakedetail},
                }
            logging.info(f"Save inputs elapsed: {timer() - start_saveinputs:1.2f}")

        logging.info(
            f"LogisticModel.calculate elapsed: {timer() - start_calculate:1.2f}"
        )
        if cleanup:
            shutil.rmtree(self.tempdir)
        return rdict


def getLogisticModelNames(config):
    """
    Get the names of the models present in the configobj

    Args:
        config: configobj object defining the model and its inputs.

    Returns:
        list: list of model names.
    """
    names = []
    lmodel_space = config
    for key, value in lmodel_space.items():
        if isinstance(value, str):
            continue
        else:  # this is a model
            names.append(key)
    return names


def getAllGridFiles(indir):
    """
    Get list of all gmt or esri (.grd, .bil) files in a directory.

    Args:
        indir (str): Directory to search.
    Returns:
        list: List of file names.
    """
    # TODO MOVE TO MAPIO
    tflist = os.listdir(indir)
    flist = []
    for tf in tflist:
        fullfile = os.path.join(indir, tf)
        ftype = getFileType(fullfile)
        if ftype in ["gmt", "esri"]:
            flist.append(fullfile)
    return flist


def validateCoefficients(cmodel):
    """
    Ensures coefficients provided in model description are valid and outputs
    a dictionary of the coefficients.

    Args:
        cmodel (dict): Sub-dictionary from config for specific model,
            for example:

            .. code-block:: python

                cmodel = config['test_model']

    Returns:
        dict: a dictionary of model coefficients named b0, b1, b2...
    """
    coeffs = {}
    for key, value in cmodel["coefficients"].items():
        if re.search("b[0-9]*", key) is None:
            raise Exception("coefficients must be named b0, b1, ...")
        coeffs[key] = float(value)
    if "b0" not in list(coeffs.keys()):
        raise Exception("coefficients must include an intercept coefficient named b0.")
    return coeffs


def validateClips(cmodel, layers, gmused):
    """
    Ensures coefficients provided in model description are valid and outputs
    a dictionary of the coefficients.

    Args:
        cmodel (dict): Sub-dictionary from config for specific model,
            for example:

            .. code-block:: python

                cmodel = config['test_model']
        layers: dictionary of layer names
        gmused (list): List of ground motion parameters used

    Returns:
        dict: a dictionary of clip values for each layer (if exists)
    """
    clips = {}
    if "clip" in cmodel:
        for key, value in cmodel["clip"].items():
            if key not in layers:
                if key not in gmused:
                    x1 = [par for par in gmused if key in par]
                    if len(x1) == 0:
                        raise Exception(
                            "Clipping key %s does not match any layers" % key
                        )
            clips[key] = (float(value[0]), float(value[1]))
    return clips


def validateLayers(cmodel):
    """
    Ensures all input files required to run the model exist and are valid
    file types. Make sure all layers are available for area of run

    Args:
        cmodel (dict): Sub-dictionary from config for specific model,
            for example,

            .. code-block:: python

                cmodel = config['test_model']

    Returns:
        dict: a dictionary of file names, e.g.

        .. code-block:: python

            {
                'slope': 'slopefile.bil',
                'vs30': 'vs30.grd'
            }

    """
    layers = {}
    longrefs = {}
    shortrefs = {}
    for key in cmodel["layers"].keys():
        for item, value in cmodel["layers"][key].items():
            if item == "file":
                ftype = getFileType(value)
                if ftype == "unknown":
                    raise Exception(
                        "layer file %s is not a valid GMT or ESRI file." % value
                    )
                if ftype == "dir":
                    value = getAllGridFiles(value)
                layers[key] = value
            elif item == "shortref":
                shortrefs[key] = value
            elif item == "longref":
                longrefs[key] = value
    return layers


def validateTerms(cmodel, coeffs, layers):
    """
    Reformats model inputs from config file, replacing functions with numpy
    functions, inserting code for extracting data from each layer (required
    to run eval in the calculate step), addressing any time variables, and
    checks that term names match coefficient names.

    Args:
        cmodel (dict): Sub-dictionary from config for specific model,
            e.g.

            .. code-block:: python

                cmodel = config['test_model']

        coeffs (dict): Dictionary of model coefficients, e.g.

            .. code-block:: python

                {'b0': 3.5, 'b1': -0.01}

        layers (dict): Dictionary of file names for all input layers, e.g.

            .. code-block:: python

                {'slope': 'slopefile.bil', 'vs30': 'vs30.grd'}

    Returns:
        tuple: (terms, timeField), where
            - 'terms' is a dictionary of terms that form the model equation,
              e.g.

            .. code-block:: python

                {
                    'b1': "self.layerdict['friction'].getData()",
                    'b2': "self.layerdict['slope'].getData()/100."
                }

            - 'timeField' indicates the time that is used to know which input
              file to read in, e.g. for monthly average precipitation, 'MONTH'.
    """
    # TODO:
    #    - Return a time field for every term, not just one global one.

    terms = {}
    timeField = None
    for key, value in cmodel["terms"].items():
        if key not in list(coeffs.keys()):
            raise Exception("Term names must match names of coefficients")
        # replace log with np.log, make sure variables are all in layers list,
        # etc.
        term, rem, tTimeField = checkTerm(value, layers)
        if tTimeField is not None:
            timeField = tTimeField
        if len(rem):
            msg = (
                'Term "%s" contains the unknown text fragment "%s". '
                "This may cause the expression to fail."
            )
            tpl = (term, rem)
            raise Exception(msg % tpl)
        terms[key] = term
    return terms, timeField


def validateInterpolations(cmodel, layers):
    """Validate logistic model interpolation.

    Args:
        cmodel (dict): Sub-dictionary from config for specific model.
        layers (dict): Dictionary of file names for all input layers.

    Returns:
        dict: Model interpolation methods.
    """
    interpolations = {}
    for key, value in cmodel["interpolations"].items():
        if key not in list(layers.keys()):
            raise Exception(
                "Interpolation key %s does not match any names of layers" % key
            )
        methods = ["linear", "nearest", "cubic", "bilinear"]
        if value not in methods:
            raise Exception(
                "Interpolation method %s not in approved list of methods: %s"
                % (key, str(methods))
            )
        interpolations[key] = value
    for key in list(layers.keys()):
        if key not in list(interpolations.keys()):
            raise Exception("No interpolation method configured for layer %s" % key)
    return interpolations


def validateUnits(cmodel):
    """Validate model units.

    Args:
        cmodel (dict): Sub-dictionary from config for specific model.

    Returns:
        dict: Model units.
    """
    units = {}
    for key in cmodel["layers"].keys():
        if "units" in cmodel["layers"][key]:
            units[key] = cmodel["layers"][key]["units"]
        else:
            raise Exception("No unit string configured for layer %s" % key)
    return units


def validateRefs(cmodel):
    """Validate references for models and layers.

    Args:
        cmodel (dict): Sub-dictionary from config for specific model.

    Returns:
        tuple: (modelrefs, longrefs, shortrefs) where:
            * modelrefs: dictionary of citation information for model
                keys='longref', 'shortref'
            * shortrefs: dictionary containing short reference for each
                input layer
            * longrefs: dictionary containing full references for each
                input layer

    """
    longrefs = {}
    shortrefs = {}
    modelrefs = {}
    for key in cmodel["layers"].keys():
        if "longref" in cmodel["layers"][key]:
            longrefs[key] = cmodel["layers"][key]["longref"]
        else:
            print("No longref provided for layer %s" % key)
            longrefs[key] = "unknown"
        if "shortref" in cmodel["layers"][key]:
            shortrefs[key] = cmodel["layers"][key]["shortref"]
        else:
            print("No shortref provided for layer %s" % key)
            shortrefs[key] = "unknown"
    try:
        modelrefs["longref"] = cmodel["longref"]
    except BaseException:
        print("No model longref provided")
        modelrefs["longref"] = "unknown"
    try:
        modelrefs["shortref"] = cmodel["shortref"]
    except BaseException:
        print("No model shortref provided")
        modelrefs["shortref"] = "unknown"
    return modelrefs, longrefs, shortrefs


def checkTerm(term, layers):
    """Checks terms of equation and replaces text with machine readable
    operators

    Args:
        term: term from model configuration file
        layers: dictionary of file names for all input layers

    Returns:
        tuple: (term, tterm, timeField) where:
            * term: dictionary of verified terms for equation with keys
                corresponding to each layer name
            * tterm: any unconverted and unverified text that may cause
                expression to fail
            * timeField: if any inputs are time dependent, output is unit of
                time (e.g., 'YEAR'), otherwise, None.
    """
    # startterm = term
    # Strip out everything that isn't: 0-9.() operators, +-/* or layer names.
    # Anything left is an unknown symbol.
    tterm = term
    # remove log, sqrt, etc.
    for op in OPERATORS:
        tterm = tterm.replace(op, "")
    # remove ShakeMap variables
    for sm_term in SM_TERMS:
        tterm = tterm.replace(sm_term, "")
    # remove layer names
    for layer in layers:
        tterm = tterm.replace(layer, "")
    # remove arithmetic operators
    tterm = re.sub(OPERATORPAT, "", tterm)
    # remove floating point numbers
    tterm = re.sub(FLOATPAT, "", tterm)
    # remove integer numbers
    tterm = re.sub(INTPAT, "", tterm)
    # remove parentheses
    tterm = re.sub("[()]*", "", tterm)
    # remove any blank spaces
    tterm = tterm.strip()
    # remove commas
    tterm = tterm.strip(",")
    # anything left *might* cause an error
    for op in OPERATORS:
        if term.find(op) > -1:
            term = term.replace(op, "np." + op)

    for sm_term in SM_GRID_TERMS:
        term = term.replace(
            sm_term,
            "self.shakemap['%s'].getSlice(rowstart, rowend, "
            "colstart, colend, name='%s')" % (sm_term, sm_term),
        )

    # replace the macro MW with the magnitude value from the shakemap
    term = term.replace("MW", "self.eventdict['magnitude']")

    # term.replace('YEAR',"self.shakemap.getEventDict()['event_time'].year")
    # hasTime = False
    timeField = None
    for unit in ["YEAR", "MONTH", "DAY", "HOUR"]:
        if term.find(unit) > -1:
            term = term.replace(unit, "")
            timeField = unit

    for layer in layers:
        if layer == "friction":
            term = term.replace(
                layer,
                "np.nan_to_num(self.layerdict['%s'].getSlice(rowstart, "
                "rowend, colstart, colend, name='%s'))" % (layer, layer),
            )
        else:
            term = term.replace(
                layer,
                "self.layerdict['%s'].getSlice(rowstart, rowend, colstart, "
                "colend, name='%s')" % (layer, layer),
            )
    return term, tterm, timeField
