# third party imports
import numpy as np
from mapio.reader import read
from mapio.grid2d import Grid2D

# local imports
from gfail.logbase import LogisticModelBase


COEFFS = {
    "b0": 0.0,  # intercept (already incorporated into X0)
    "b1": 1.65,  # log(pgv)
    "b2": 1.0,  # X0, all inputs without shaking
    "b3": 0.01,  # log(pgv)*arctan(slope)
}

TERMS = {
    "b1": "np.log(pgv._data)",
    "b2": "X0._data",
    "b3": "np.log(pgv._data) * np.arctan(slope._data) * 180/np.pi",
}

TERMLAYERS = {
    "b1": "pgv",
    "b2": "X0",
    "b3": "pgv, slope",
}

SHAKELAYERS = ["pgv"]

CLIPS = {"pgv": (0.0, 211.0)}  # cm/s

# Coefficients for conversion to coverage
COV_COEFFS = {"a": -7.592, "b": 5.237, "c": -3.042, "d": 4.035}


class Jessee2018Model(LogisticModelBase):
    def __init__(
        self,
        shakefile,
        config,
        bounds=None,
        uncertfile=None,
        trimfile=None,
        saveinputs=False,
    ):
        self.COEFFS = COEFFS
        self.TERMS = TERMS
        self.TERMLAYERS = TERMLAYERS
        self.SHAKELAYERS = SHAKELAYERS
        self.do_coverage = True

        self.prob_units = "Proportion of area affected"

        super().__init__(
            shakefile,
            config,
            bounds=bounds,
            uncertfile=uncertfile,
            trimfile=trimfile,
            saveinputs=saveinputs,
        )

    def pre_process(self, key, grid):
        """Correct grids in model specific way."""
        if key in CLIPS:
            clipmin, clipmax = CLIPS[key]
            grid._data = np.clip(grid._data, clipmin, clipmax)
        return

    def calculate_coverage(self, P):
        a = COV_COEFFS["a"]
        b = COV_COEFFS["b"]
        c = COV_COEFFS["c"]
        d = COV_COEFFS["d"]
        P = np.exp(a + b * P + c * P**2 + d * P**3)
        return P

    def modify_slope(self, slope):
        """Perform modifications to slope to convert to degrees."""
        slope = np.arctan(slope) * 180 / np.pi
        return slope

    def calculate_uncertainty(self):
        if "stddev" in self.layers:
            varP = read(self.layers["stddev"])._data
        else:
            varP = float(self.config["default_stddev"])
        # TODO: Sort this out without having to load in all these
        # layers at the same time, if possible...
        slope = read(self.layers["slope"])._data

        varP = np.arctan(slope) * 180 / np.pi
        varP *= self.COEFFS["b6"]
        varP += self.COEFFS["b1"]
        varP **= 2
        del slope
        std_pgv = read(self.layers["stdpgv"])._data
        varP *= std_pgv**2
        del std_pgv
        if "stddev" in self.layers:
            stddev = read(self.layers["stddev"])._data
        else:
            stddev = float(self.config["default_stddev"])
        varP += stddev**2
        del stddev
        X = read(self.layers["X"])._data
        varP *= (np.exp(-X) / (np.exp(-X) + 1) ** 2) ** 2
        del X

        if self.do_coverage:
            P = read(self.layers["P"])._data
            a = COV_COEFFS["a"]
            b = COV_COEFFS["b"]
            c = COV_COEFFS["c"]
            d = COV_COEFFS["d"]

            std1 = (
                np.exp(a + b * P + c * P**2.0 + d * P**3.0)
                * (b + 2.0 * P * c + 3.0 * d * P**2.0)
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
        return P
