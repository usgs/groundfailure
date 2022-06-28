# third party imports
import numpy as np

# local imports
from gfail import Zhu2017Model


class Zhu2017ModelCoastal(Zhu2017Model):

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

        super().__init__(
            shakefile,
            config,
            bounds=bounds,
            uncertfile=uncertfile,
            trimfile=trimfile,
            slopefile=slopefile,
            saveinputs=saveinputs,
        )

    def calculate_coverage(self, P):
        P = 0.4208 / (1 + 62.59 * np.exp(-11.43 * P)) ** 2
        return P
