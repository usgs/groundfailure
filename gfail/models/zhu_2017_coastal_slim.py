# third party imports
import numpy as np
from mapio.reader import read

# local imports
from gfail.logbase import LogisticModelBase

# Coefficients for conversion to coverage
COV_COEFFS = {"a": 0.4208, "b": 62.59, "c": 11.43}
SHAKELAYERS = ["pgv", "pga"]

class Zhu2017ModelCoastalSlim(LogisticModelBase):

    TERMS = {
        "b1": "np.log(pgv._data*(1/(1+np.power(2.71828,-2*(MW-6)))))",
        "b2": "X0._data",
    }

    COEFFS = {
        "b0": 0.,
        "b1": 0.301,
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

        super().__init__(
            shakefile,
            config,
            bounds=bounds,
            uncertfile=uncertfile,
            trimfile=trimfile,
            slopefile=slopefile,
            saveinputs=saveinputs,
        )


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
        if "minpgv" in self.config.keys():
            minpgv = float(self.config["minpgv"])
            pgv = read(self.layers["pgv"])._data
            P[pgv < minpgv] = np.nan

        if "minpga" in self.config.keys():
            minpga = float(self.config["minpga"])
            pga = read(self.layers["pga"])._data
            P[pga < minpga] = np.nan
        return P