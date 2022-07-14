# stdlib imports
import pathlib
import tempfile
import shutil
import logging
import collections
from timeit import default_timer as timer

# third party imports
import numpy as np
from mapio.reader import read, get_file_geodict
from mapio.writer import write
from mapio.shake import ShakeGrid
from mapio.geodict import GeoDict
from mapio.grid2d import Grid2D


# local imports
from gfail.spatial import trim_ocean2

COEFFS = {}

TERMS = {}

TERMLAYERS = {}

SHAKE_LAYERS = []


class LogisticModelBase(object):
    """TODO needs documentation

    Args:
        object (_type_): _description_

    Raises:
        Exception: _description_
        Exception: _description_

    Returns:
        _type_: _description_
    """

    COEFFS = {}
    TERMS = {}
    TERMLAYERS = {}
    SHAKELAYERS = []
    do_coverage = False
    prob_units = None
    notes = ""
    slopefile = None
    nonzero = None
    slopemin = "none"
    slopemax = "none"

    def __init__(
        self,
        shakefile,
        config,
        bounds=None,
        uncertfile=None,
        trimfile=None,
        slopefile=None,
        saveinputs=False,
    ):
        """Initialize LogisticModelBase object.

        Args:
            shakefile (str): Path to ShakeMap grid.xml file.
            config (dict): Dict like object as a result of reading in "logbase"
                           format INI files.
            bounds (dict): Fields are 'xmin', 'xmax', 'ymin', 'ymax'.
            uncertfile (str): Path to ShakeMap uncertainty.xml file.
            trimfile (str): Path to shapefile to use for masking ocean pixels.
            slopefile (str): File containing slope data to be used for slope
                             masking.
            saveinputs (bool): Save input layer grids with model output.

        Notes: All input data grids are loaded in one at a time, and saved to
        a temporary folder that is deleted when the object is deleted. The
        file names loaded are stored in the self.layers dictionary.

        """
        logging.info("Initializing model...")
        self.tempdir = tempfile.mkdtemp()
        self.config = config
        self.bounds = bounds
        self.uncertfile = uncertfile
        self.trimfile = trimfile
        self.saveinputs = saveinputs
        self.modeltype = config["gfetype"]
        logging.info("Loading raw ShakeMap...")
        start_shake = timer()
        self.raw_shake_grid = ShakeGrid.load(shakefile, adjust="res")
        self.get_sample_dict()  # sets self.sampledict
        logging.info(f"Loading resampled ShakeMap with geodict: {self.sampledict}")
        self.shake_grid = ShakeGrid.load(
            shakefile,
            samplegeodict=self.sampledict,
            resample=True,
            doPadding=True,
            method="linear",
            adjust="res",
        )
        #logging.info(f"Load shake elapsed: {timer() - start_shake:1.2f}")
        start_error = timer()
        if self.uncertfile is not None:
            self.error_grid = ShakeGrid.load(
                self.uncertfile,
                samplegeodict=self.sampledict,
                resample=True,
                doPadding=True,
                method="linear",
                adjust="res",
            )
        else:
            self.error_grid = None
        #logging.info(f"Load uncertainty elapsed: {timer() - start_error:1.2f}")

        logging.info("Loaded resampled ShakeMap.")
        self.shakedict = self.shake_grid.getShakeDict()
        self.eventdict = self.shake_grid.getEventDict()

        # make sure that if this model requires slope parameters, a slope file
        # has either been passed in or configured as a layer
        logging.info("Checking slope parameters...")
        if slopefile is not None:
            self.slopefile = slopefile
        elif "slopefile" in config:
            self.slopefile = config["slopefile"]
        elif "slope" in config["layers"]:
            self.slopefile = config["layers"]["slope"]["file"]
        else:
            self.slopefile = None

        if self.slopefile is not None:
            self.set_slope_params(config)

        logging.info("Checking references, units, interpolations...")
        (self.modelrefs, self.longrefs, self.shortrefs) = self.validate_refs(config)
        self.units = self.validate_units(config)
        if self.prob_units is None:
            # this normally will be defined in a subclass
            self.prob_units = "Probability of any occurrence"
        self.interpolations = config["interpolations"]

        logging.info("Looping over layers...")
        self.layers = {}
        for key, layer in config["layers"].items():
            layerfile = pathlib.Path(layer["file"])
            logging.info(f'Reading {layer["file"]}...')
            interp = self.interpolations[key]
            grid = read(
                layerfile,
                samplegeodict=self.sampledict,
                method=interp,
                resample=True,
                doPadding=True,
                interp_approach="rasterio",
            )
            self.pre_process(key, grid)
            outfile = pathlib.Path(self.tempdir, f"{key}.cdf")
            logging.info(f"Writing sampled layer {key}...")
            start_write = timer()
            write(grid, outfile, "netcdf")
            #logging.info(f"write {key} elapsed: {timer() - start_write:1.2f}")
            self.layers[key] = outfile

        for key in self.SHAKELAYERS:
            shake_grid = self.shake_grid.getLayer(key)
            # shake_grid = raw_shake_grid.interpolate2(self.sampledict, method="linear")
            outfile = pathlib.Path(self.tempdir, f"{key}.cdf")
            logging.info(f"Writing sampled layer {key}...")
            start_write = timer()
            write(shake_grid, outfile, "netcdf")
            #logging.info(f"write {key} elapsed: {timer() - start_write:1.2f}")

            self.layers[key] = outfile
        if uncertfile is not None:
            start_write = timer()
            self.save_uncertainty_layers()
            #logging.info(f"write uncertainty elapsed: {timer() - start_write:1.2f}")

        if "slope" not in config["layers"] and self.slopefile is not None:
            self.read_slope()
        if self.slopefile is not None:
            # establish the indices of the slope grid
            # where values are > slopemin and <= slopemax
            self.get_slope_non_zero()

    def read_slope(self):
        interp = "linear"
        if "slope" in self.interpolations:
            interp = self.interpolations["slope"]
        fdict = get_file_geodict(self.slopefile)
        if fdict.isAligned(self.sampledict):
            grid = read(self.slopefile, samplegeodict=self.sampledict)
        else:
            grid = read(
                self.slopefile,
                samplegeodict=self.sampledict,
                method=interp,
                resample=True,
                doPadding=True,
                interp_approach="rasterio",
            )
        outfile = pathlib.Path(self.tempdir, "slope.cdf")
        logging.info("Writing sampled layer slope...")
        write(grid, outfile, "netcdf")
        self.layers["slope"] = outfile

    def get_slope_non_zero(self):
        slopemin = self.slopemin
        slopemax = self.slopemax
        if self.slopemin == "none":
            nx = self.sampledict.nx
            ny = self.sampledict.ny
            self.nonzero = np.ones((ny, nx), dtype=bool)
            return
        slopefile = self.layers["slope"]
        slope = read(slopefile)._data
        # do whatever conversions are necessary to get slope in degrees.
        # this method is implemented in the child class.
        slope = self.modify_slope(slope)
        self.nonzero = np.array([(slope > slopemin) & (slope <= slopemax)])
        self.nonzero = self.nonzero[0, :, :]
        del slope

    def modify_slope(self, slope):
        """This method should be implemented by child classes."""
        pass

    def set_slope_params(self, config):
        # Find slope thresholds, if applicable
        self.slopemin = "none"
        self.slopemax = "none"
        if "slopefile" in config:
            try:
                self.slopemin = float(config["slopemin"])
                self.slopemax = float(config["slopemax"])
            except BaseException:
                print(
                    "Could not find slopemin and/or slopemax in config, "
                    "limits. No slope thresholds will be applied."
                )
                self.slopemin = "none"
                self.slopemax = "none"

    def get_units(self, layer):
        """Get units for an input layer.

        Args:
            layer (str): name of layer.

        Returns:
            str: units.
        """
        try:
            # should work for everything except ground motion layers
            units = self.units[layer]
        except BaseException:
            if "pga" in layer.lower():
                units = "%g"
            elif "pgv" in layer.lower():
                units = "cm/s"
            elif "mmi" in layer.lower():
                units = "intensity"
            else:
                units = ""
        return units

    def validate_units(self, cmodel):
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

    def validate_refs(self, cmodel):
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

    def save_uncertainty_layers(self):
        for shakelayer in self.SHAKELAYERS:
            key = f"std{shakelayer}"
            method = "linear"
            if key in self.config["interpolations"].keys():
                method = self.interpolations[key]
            if self.error_grid is not None:
                gm_grid = self.error_grid.getLayer(key)
                gm_grid = gm_grid.interpolate2(self.sampledict, method=method)

                filename = pathlib.Path(self.tempdir) / f"{key}.cdf"
                self.layers[key] = filename
                write(gm_grid, filename, "netcdf")

    def pre_process(self, key, grid):
        """This method should be implemented by child classes."""
        pass

    def get_sample_dict(self):
        baselayer = self.config["baselayer"]
        basefile = self.config["layers"][baselayer]["file"]
        basedict = get_file_geodict(basefile)
        shake_dict = self.raw_shake_grid.getGeoDict()
        intersection = shake_dict.getIntersection(basedict)
        for layername, layer in self.config["layers"].items():
            layerfile = layer["file"]
            layerdict = get_file_geodict(layerfile)
            intersection = layerdict.getIntersection(intersection)
        self.sampledict = intersection
        # self.sampledict = None
        # shake_geodict = self.shake_grid.getGeoDict()
        # if 'baselayer' in self.config:
        #     baselayer = self.config['baselayer']
        #     fname = self.config['layers'][baselayer]['file']
        #     basefile = pathlib.Path(fname)
        #     base_geodict = get_file_geodict(basefile)
        #     if self.bounds is None:
        #         self.sampledict = base_geodict.getBoundsWithin(shake_geodict)
        #     else:
        #         mod_geodict = GeoDict.createDictFromBox(self.bounds['xmin'],
        #                                                 self.bounds['xmax'],
        #                                                 self.bounds['ymin'],
        #                                                 self.bounds['ymax'],
        #                                                 base_geodict.dx,
        #                                                 base_geodict.dy
        #                                                 )
        #         if not shake_geodict.contains(mod_geodict):
        #             msg = ('Desired sample bounds must be within the bounds '
        #                    'of the ShakeMap.')
        #             raise Exception(msg)
        #         # intersection should have the resolution of the base geodict
        #         intersection = shake_geodict.getIntersection(mod_geodict)
        #         self.sampledict = mod_geodict.getBoundsWithin(intersection)
        # else:
        #     self.sampledict = shake_geodict.copy()

        # we may have some config that tells us to resample even the base layer
        # by some factor of resolution.
        self.subdivide()

    def subdivide(self):
        # Do we need to subdivide baselayer?
        if "divfactor" in self.config.keys():
            divfactor = float(self.config["divfactor"])
            if divfactor != 1.0:
                # adjust sampledict so everything will be resampled
                newxmin = (
                    self.sampledict.xmin
                    - self.sampledict.dx / 2.0
                    + self.sampledict.dx / (2.0 * divfactor)
                )
                newymin = (
                    self.sampledict.ymin
                    - self.sampledict.dy / 2.0
                    + self.sampledict.dy / (2.0 * divfactor)
                )
                newxmax = (
                    self.sampledict.xmax
                    + self.sampledict.dx / 2.0
                    - self.sampledict.dx / (2.0 * divfactor)
                )
                newymax = (
                    self.sampledict.ymax
                    + self.sampledict.dy / 2.0
                    - self.sampledict.dy / (2.0 * divfactor)
                )
                newdx = self.sampledict.dx / divfactor
                newdy = self.sampledict.dy / divfactor
                if np.abs(newxmax) > 180.0:
                    newxmax = np.sign(newxmax) * 180.0
                if np.abs(newxmin) > 180.0:
                    newxmin = np.sign(newxmin) * 180.0

                self.sampledict = GeoDict.createDictFromBox(
                    newxmin, newxmax, newymin, newymax, newdx, newdy, inside=True
                )

    def calculate(self):
        """Calculate the probability, and sigma (if possible).

         Returns:
            dict: Must contain at least one sub dictionary with
                  key "model". That dictionary should look like:
                  {'grid': Grid2D object contanining probability values,
                   'label': String used for display label
                   'type': 'output',
                   'description': String description of model.
                  }

        This method uses an "accumulator" array, and loads in data layers
        from the temporary directory one *term* at a time. As some terms can
        be a combination of layers, this may mean that multiple grids are
        loaded into memory simultaneously, in addition to the accumulator grid.
        Layer grids are deleted from memory once they have been added to the
        accumulator array.


        """
        nx = self.sampledict.nx
        ny = self.sampledict.ny

        # dictionary to hold output
        rdict = collections.OrderedDict()

        X = np.ones((ny, nx)) * self.COEFFS["b0"]
        for term, operation in self.TERMS.items():
            logging.info(f"Reading layer {term}")
            start_term = timer()
            coeff = self.COEFFS[term]
            layers = self.TERMLAYERS[term].split(",")
            layers = [layer.strip() for layer in layers]
            for layer in layers:
                layerfile = self.layers[layer]
                loadstr = f'{layer} = read("{layerfile}",apply_nan=True)'
                globaldict = {f"{layer}": layer, "read": read}
                logging.info(f"Reading sampled layer {layer}...")
                exec(loadstr, globaldict)
                # this should be a reference...
                globals()[layer] = globaldict[layer]
                if self.saveinputs:
                    units = self.get_units(layer)
                    rdict[layer] = {
                        "grid": globaldict[layer],
                        "label": f"{layer} ({units})",
                        "type": "input",
                        "description": {"units": units},
                    }
                    if layer in self.shortrefs:
                        rdict[layer]["description"]["name"] = self.shortrefs[layer]
                    if layer in self.longrefs:
                        rdict[layer]["description"]["longref"] = self.longrefs[layer]
            try:
                msg = f"Executing layer operation {operation} * {coeff}..."
                logging.info(msg)
                # replace the macro MW with the magnitude value from
                # the shakemap
                magstr = f'{self.eventdict["magnitude"]:.1f}'
                operation = operation.replace("MW", magstr)
                X += eval(operation) * coeff
            except Exception as e:
                msg = (
                    f"{str(e)}: Unable to calculate term {term}: "
                    f"operation {operation}*{coeff}."
                )
                raise Exception(msg)
            for layer in layers:
                del globals()[layer]
            #logging.info(f"Read term elapsed: {timer() - start_term:1.2f}")

        logging.info("Calculating probability...")
        # save off the X grid for potential use by
        # uncertainty calculations
        outfile = pathlib.Path(self.tempdir, "X.cdf")
        logging.info("Writing intermediate layer X...")
        grid = Grid2D(data=X, geodict=self.sampledict)
        start_writex = timer()
        write(grid, outfile, "netcdf")
        #logging.info(f"write X elapsed: {timer() - start_writex:1.2f}")
        self.layers["X"] = outfile
        P = 1 / (1 + np.exp(-X))
        del X

        start_modify = timer()
        P = self.modify_probability(P)
        #logging.info(f"modify prob elapsed: {timer() - start_modify:1.2f}")
        # save off the intermediate P grid for potential use by
        # uncertainty calculations
        outfile = pathlib.Path(self.tempdir, "P.cdf")
        logging.info("Writing intermediate layer P...")
        grid = Grid2D(data=P, geodict=self.sampledict)
        start_writep = timer()
        write(grid, outfile, "netcdf")
        #logging.info(f"wite P elapsed: {timer() - start_writep:1.2f}")
        self.layers["P"] = outfile

        sigma_grid = None
        if self.uncertfile is not None:
            start_uncertainty = timer()
            sigma_grid = self.calculate_uncertainty()
            #logging.info(f"uncertainty elapsed: {timer() - start_uncertainty:1.2f}")

        # because non-coverage P grid is used for uncertainty,
        # we cannot convert to areal coverage until that is complete.
        if self.do_coverage:
            start_coverage = timer()
            P = self.calculate_coverage(P)
            #logging.info(f"coverage elapsed: {timer() - start_coverage:1.2f}")

        # apply slope cutoffs
        if self.nonzero is not None:
            P = P * self.nonzero

        # create Grid2D object for P
        p_grid = Grid2D(data=P, geodict=self.sampledict)

        description = self.get_description()

        # trim off ocean pixels if user wanted that
        if self.trimfile is not None:
            start_trim = timer()
            # Turn all offshore cells to nan
            p_grid = trim_ocean2(p_grid, self.trimfile)
            if sigma_grid is not None:
                sigma_grid = trim_ocean2(sigma_grid, self.trimfile)
            #logging.info(f"trim elapsed: {timer() - start_trim:1.2f}")

        rdict["model"] = {
            "grid": p_grid,
            "label": "%s estimate - %s"
            % (self.modeltype.capitalize(), self.prob_units.title()),
            "type": "output",
            "description": description,
        }
        if sigma_grid is not None:
            rdict["std"] = {
                "grid": sigma_grid,
                "label": (
                    "%s estimate - %s (std)"
                    % (self.modeltype.capitalize(), self.prob_units.title())
                ),
                "type": "output",
                "description": description,
            }

        return rdict

    def get_description(self):
        shakedetail = "%s_ver%s" % (
            self.shakedict["shakemap_id"],
            self.shakedict["shakemap_version"],
        )
        description = {
            "name": self.modelrefs["shortref"],
            "longref": self.modelrefs["longref"],
            "units": self.units,
            "shakemap": shakedetail,
            "event_id": self.eventdict["event_id"],
            "parameters": {
                "slopemin": self.slopemin,
                "slopemax": self.slopemax,
                "modeltype": self.modeltype,
                "notes": self.notes,
            },
        }
        if "vs30max" in self.config.keys():
            description["vs30max"] = float(self.config["vs30max"])
        if "minpgv" in self.config.keys():
            description["minpgv"] = float(self.config["minpgv"])

        return description

    def calculate_coverage(self, P):
        """This method should be implemented by child classes."""
        pass

    def calculate_uncertainty(self):
        """This method should be implemented by child classes."""
        pass

    def modify_probability(self, P):
        """This method should be implemented by child classes."""
        pass

    def __del__(self):
        if pathlib.Path(self.tempdir).exists():
            shutil.rmtree(self.tempdir)
