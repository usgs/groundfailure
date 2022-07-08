# third party imports
import numpy as np
from mapio.reader import read
from mapio.grid2d import Grid2D

# local imports
from gfail.logbase import LogisticModelBase


SHAKELAYERS = ["pgv", "pga"]

CLIPS = {"pgv": (0.0, 150.0)}  # cm/s

# Coefficients for conversion to coverage
COV_COEFFS = {"a": 0.4915, "b": 42.4, "c": 9.165}


class Zhu2017Model(LogisticModelBase):
    TERMS = {
        "b1": "np.log(pgv._data*(1/(1+np.power(2.71828,-2*(MW-6)))))",
        "b2": "X0._data",
    }

    COEFFS = {
        "b0": 0,
        "b1": 0.334,
        "b2": 1.,
    }

    TERMLAYERS = {
        "b1": "pgv",
        "b2": "X0",
    }

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
        self.COEFFS = self.COEFFS
        self.TERMS = self.TERMS
        self.TERMLAYERS = self.TERMLAYERS
        self.SHAKELAYERS = SHAKELAYERS
        self.do_coverage = True

        self.prob_units = "Proportion of area affected"

        super().__init__(
            shakefile,
            config,
            bounds=bounds,
            uncertfile=uncertfile,
            trimfile=trimfile,
            slopefile=slopefile,
            saveinputs=saveinputs,
        )

    def pre_process(self, key, grid):
        if key in CLIPS:
            clipmin, clipmax = CLIPS[key]
            grid._data = np.clip(grid._data, clipmin, clipmax)
        return

    def modify_slope(self, slope):
        """Perform modifications to slope to convert to degrees."""
        slope = np.arctan(slope) * 180 / np.pi
        return slope

    def calculate_coverage(self, P):
        a = COV_COEFFS["a"]
        b = COV_COEFFS["b"]
        c = COV_COEFFS["c"]        
        P = a / (1 + b * np.exp(-c * P)) ** 2
        return P

    def calculate_uncertainty(self):
        if "stddev" in self.layers:
            varP = read(self.layers["stddev"])._data
        else:
            varP = float(self.config["default_stddev"])

        stdpgv = read(self.layers["stdpgv"])._data
        varP += varP**2 + (self.COEFFS["b1"] ** 2 * stdpgv**2)
        del stdpgv
        X = read(self.layers["X"])._data
        varP = (np.exp(-X) / (np.exp(-X) + 1) ** 2.0) ** 2.0 * varP
        del X
        if self.do_coverage:
            P = read(self.layers["P"])._data
            a = COV_COEFFS["a"]
            b = COV_COEFFS["b"]
            c = COV_COEFFS["c"]
            std1 = (
                (2 * a * b * c * np.exp(-c * P)) / ((1 + b * np.exp(-c * P)) ** 3.0)
            ) ** 2.0 * varP
            std1 = np.sqrt(std1)
            del P
            del varP
        else:
            std1 = np.sqrt(varP)
            del varP

        sigma = Grid2D(data=std1, geodict=self.sampledict)
        return sigma

    def modify_probability(self, P):
        if "vs30max" in self.config.keys():
            vs30max = float(self.config["vs30max"])
            vs30 = read(self.layers["vs30"])._data
            P[vs30 > vs30max] = 0.0

        if "minpgv" in self.config.keys():
            minpgv = float(self.config["minpgv"])
            pgv = read(self.layers["pgv"])._data
            P[pgv < minpgv] = 0.0

        if "minpga" in self.config.keys():
            minpga = float(self.config["minpga"])
            pga = read(self.layers["pga"])._data
            P[pga < minpga] = 0.0
        return P
