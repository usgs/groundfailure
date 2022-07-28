# third party imports
import numpy as np


# local imports
from gfail.logbase import LogisticModelBase


COEFFS = {
    "b0": -3.6490,  # intercept
    "b1": 0.0133,  # log(pgv)
    "b2": 0.0364,  # slope (deg)
    "b3": -0.0635,  # friction
    "b4": -0.0004,  # cti
    "b5": 0.0019,  # interaction term
}

TERMS = {
    "b1": "pga._data",
    "b2": "slope._data / 100", #input file is scaled up, scale back down
    "b3": "friction._data",
    "b4": "cti1._data * 100",
    "b5": "pga._data * slope._data / 100",
}

TERMLAYERS = {
    "b1": "pga",
    "b2": "slope",
    "b3": "friction",
    "b4": "cti1",
    "b5": "pga, slope",
}

SHAKELAYERS = ["pga"]

CLIPS = {
    "pga": (0.0, 170.0),
}  # cm/s


class Nowicki2014Model(LogisticModelBase):
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
        if key in CLIPS:
            clipmin, clipmax = CLIPS[key]
            grid._data = np.clip(grid._data, clipmin, clipmax)
        return

    def calculate_coverage(self, P):
        return P

    def modify_slope(self, slope):
        """Perform modifications to slope to convert to degrees."""
        slope = np.arctan(slope) * 180 / np.pi
        return slope

    def modify_probability(self, P):
        return P
