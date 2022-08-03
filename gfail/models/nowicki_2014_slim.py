# third party imports
import numpy as np


# local imports
from gfail.logbase import LogisticModelBase


COEFFS = {
    "b0": 0.,  # intercept
    "b1": 0.0133,  # log(pgv)
    "b2": 1., #X0
    "b3": 0.0019,  # interaction term
}

TERMS = {
    "b1": "pga._data",
    "b2": "X0._data", 
    "b3": "pga._data * slope._data / 100", #input slope file is scaled up, scale back down
}

TERMLAYERS = {
    "b1": "pga",
    "b2": "X0",
    "b3": "pga, slope",
}

SHAKELAYERS = ["pga"]

CLIPS = {
    "pga": (0.0, 170.0),
}


class Nowicki2014ModelSlim(LogisticModelBase):
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

        self.prob_units = "Probability of any landslide"

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
