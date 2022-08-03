# third party imports
import numpy as np
from mapio.reader import read

# local imports
from gfail.logbase import LogisticModelBase

SHAKELAYERS = ["pgv", "pga"]

CLIPS = {
    "precip": (0.0, 2500.0),  # use mm because clipping is done on original data
    "pgv": (0.0, 150.0),
}  # cm/s

# Coefficients for conversion to coverage
COV_COEFFS = {"a": 0.4208, "b": 62.59, "c": 11.43}


class Zhu2017ModelCoastal(LogisticModelBase):

    TERMS = {
        "b1": "np.log(pgv._data*(1/(1+np.power(2.71828,-2*(MW-6)))))",
        "b2": "np.log(vs30._data)",
        "b3": "precip._data",
        "b4": "np.power(dc._data, 0.5)",
        "b5": "dr._data",
        "b6": "np.power(dc._data, 0.5) * dr._data",
    }

    COEFFS = {
        "b0": 12.435,
        "b1": 0.301,
        "b2": -2.615,
        "b3": 0.0005556,
        "b4": -0.0287,
        "b5": 0.0666,
        "b6": -0.0369,
    }

    TERMLAYERS = {
        "b1": "pgv",
        "b2": "vs30",
        "b3": "precip",
        "b4": "dc",
        "b5": "dr",
        "b6": "dc, dr",
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

    def modify_probability(self, P):
        if "vs30max" in self.config.keys():
            vs30max = float(self.config["vs30max"])
            vs30 = read(self.layers["vs30"])._data
            P[vs30 > vs30max] = np.nan

        if "minpgv" in self.config.keys():
            minpgv = float(self.config["minpgv"])
            pgv = read(self.layers["pgv"])._data
            P[pgv < minpgv] = np.nan

        if "minpga" in self.config.keys():
            minpga = float(self.config["minpga"])
            pga = read(self.layers["pga"])._data
            P[pga < minpga] = np.nan
        return P
