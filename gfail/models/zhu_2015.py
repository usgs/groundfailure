# third party imports
import numpy as np
from mapio.reader import read

# local imports
from gfail.logbase import LogisticModelBase

TERMS = {
    "b1": "np.log((pga._data/100.0)*(np.power(MW,2.56)/np.power(10,2.24)))",
    "b2": "cti._data",
    "b3": "np.log(vs30._data)",
}

COEFFS = {
    "b0": 24.10,
    "b1": 2.067,
    "b2": 0.355,
    "b3": -4.784,
}


TERMLAYERS = {
    "b1": "pga",
    "b2": "cti",
    "b3": "vs30",
}

SHAKELAYERS = ["pga"]

CLIPS = {
    "cti": (0.0, 15.0),
    "pga": (0.0, 270.0),
}


class Zhu2015Model(LogisticModelBase):
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
        P = 0.81 * P  # correction factor to coverage found in Zhu et al 2017 paper
        return P

    def modify_probability(self, P):
        if "vs30max" in self.config.keys():
            vs30max = float(self.config["vs30max"])
            vs30 = read(self.layers["vs30"])._data
            P[vs30 > vs30max] = np.nan
        return P
