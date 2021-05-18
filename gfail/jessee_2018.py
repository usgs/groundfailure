# stdlib imports
import logging
import pathlib
import sys
import time

# third party imports
import numpy as np
from configobj import ConfigObj
from mapio.reader import read
from mapio.writer import write
from mapio.grid2d import Grid2D
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# local imports
from gfail.logbase import LogisticModelBase


COEFFS = {'b0': -6.30,   # intercept
          'b1': 1.65,  # log(pgv)
          'b2': 0.06,  # arctan(slope)
          'b3': 1,  # lithology set to 1.0 - coefficients are in glim file
          'b4': 0.03,  # cti
          'b5': 1.0,  # landcover
          'b6': 0.01  # log(pgv)*arctan(slope)
          }

TERMS = {'b1': 'np.log(pgv._data)',
         'b2': 'np.arctan(slope._data) * 180/np.pi',
         'b3': 'rock._data',
         'b4': 'cti._data',
         'b5': 'landcover._data',
         'b6': 'np.log(pgv._data) * np.arctan(slope._data) * 180/np.pi',
         }

TERMLAYERS = {
    'b1': 'pgv',
    'b2': 'slope',
    'b3': 'rock',
    'b4': 'cti',
    'b5': 'landcover',
    'b6': 'pgv, slope'
}

SHAKELAYERS = ['pgv']

# Jessee specific fixes
OLD_UNCONSOLIDATED = -3.21
NEW_UNCONSOLIDATED = -1.36

CLIPS = {'cti': (0., 19.),
         'pgv': (0., 211.)}  # cm/s

ERROR_COEFFS = {'a': -7.592,
                'b': 5.237,
                'c': -3.042,
                'd': 4.035}


class Jessee2018Model(LogisticModelBase):
    def __init__(self, shakefile, config,
                 bounds=None, uncertfile=None, trimfile=None):
        self.COEFFS = COEFFS
        self.TERMS = TERMS
        self.TERMLAYERS = TERMLAYERS
        self.SHAKELAYERS = SHAKELAYERS
        self.do_coverage = True

        self.prob_units = 'Proportion of area affected'

        super().__init__(shakefile, config,
                         bounds=bounds, uncertfile=uncertfile,
                         trimfile=trimfile)

    def pre_process(self, key, grid):
        '''Correct grids in model specific way.

        '''
        if key == 'rock':
            grid._data[grid._data <= OLD_UNCONSOLIDATED] = NEW_UNCONSOLIDATED
            self.notes += ('unconsolidated sediment coefficient '
                           f'changed to {NEW_UNCONSOLIDATED} (weaker) '
                           'from {OLD_UNCONSOLIDATED} to '
                           'better reflect that this '
                           'unit is not actually strong\n')
        if key in CLIPS:
            clipmin, clipmax = CLIPS[key]
            grid._data = np.clip(grid._data, clipmin, clipmax)
        return

    def calculate_coverage(self, P):
        P = np.exp(-7.592 + 5.237 * P - 3.042 * P**2 + 4.035 * P**3)
        return P

    def modify_slope(self, slope):
        '''Perform modifications to slope to convert to degrees.
        '''
        slope = np.arctan(slope) * 180 / np.pi
        return slope

    def calculate_uncertainty(self):
        if 'stddev' in self.layers:
            varP = read(self.layers['stddev'])._data
        else:
            varP = float(self.config['default_stddev'])
        # TODO: Sort this out without having to load in all these
        # layers at the same time, if possible...
        slope = read(self.layers['slope'])._data

        varP = np.arctan(slope) * 180 / np.pi
        varP *= self.COEFFS['b6']
        varP += self.COEFFS['b1']
        varP **= 2
        del slope
        std_pgv = read(self.layers['stdpgv'])._data
        varP *= std_pgv**2
        del std_pgv
        if 'stddev' in self.layers:
            stddev = read(self.layers['stddev'])._data
        else:
            stddev = float(self.config['default_stddev'])
        varP += stddev**2
        del stddev
        X = read(self.layers['X'])._data
        varP *= (np.exp(-X) / (np.exp(-X) + 1)**2)**2
        del X

        if self.do_coverage:
            P = read(self.layers['P'])._data
            a = ERROR_COEFFS['a']
            b = ERROR_COEFFS['b']
            c = ERROR_COEFFS['c']
            d = ERROR_COEFFS['d']

            std1 = (np.exp(a + b * P + c * P**2. + d * P**3.) *
                    (b + 2. * P * c + 3. * d * P**2.))**2. * varP
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
