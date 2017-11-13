#!/usr/bin/env python

import os.path
import os
from configobj import ConfigObj
import numpy as np
from gfail.godt import godt2008
from mapio.geodict import GeoDict
from gfail.conf import correct_config_filepaths

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
datadir = os.path.abspath(os.path.join(homedir, 'data'))

shakefile = os.path.join(datadir, 'test_shakegrid.xml')
shakefileBounds = os.path.join(datadir, 'test_shakegrid_bounds.xml')
uncertfile = os.path.join(datadir, 'test_uncert.xml')
cofile = os.path.join(datadir, 'test_friction.bil')
slopefile = os.path.join(datadir, 'test_slope.bil')
vs30file = os.path.join(datadir, 'test_vs30.bil')
ctifile = os.path.join(datadir, 'test_cti.bil')
precipfolder = os.path.join(datadir, 'test_precip')

fakegeodict = GeoDict({'xmin': 0.5, 'xmax': 1.5,
                       'ymin': 0.5, 'ymax': 1.5,
                       'dx': 1., 'dy': 1.,
                       'ny': 2, 'nx': 2})


def test_godt2008():
    configfile = os.path.join(datadir, 'testconfig_godt.ini')
    config = ConfigObj(configfile)
    bounds = {'xmin': -1.5, 'xmax': 3, 'ymin': -1.5, 'ymax': 3}
    # Test path correction (from conf.py)
    config = correct_config_filepaths(datadir, config)
    maplayers1 = godt2008(shakefile, config, displmodel='J_PGA',
                          uncertfile=uncertfile, bounds=bounds)
    maplayers2 = godt2008(shakefile, config, displmodel='J_PGA_M',
                          saveinputs=True)
    maplayers3 = godt2008(shakefile, config, displmodel='RS_PGA_M')
    maplayers4 = godt2008(shakefile, config, displmodel='RS_PGA_PGV',
                          uncertfile=uncertfile, saveinputs=True)
    maplayers5 = godt2008(shakefile, config)

    # Testing J_PGA with Uncertainty - Model Type 1
    np.testing.assert_allclose(maplayers1['model']['grid'].getData(),
                               np.array([[0., 0.1], [0.9, 0.]]), atol=0.01)
    np.testing.assert_allclose(maplayers1['modelmin']['grid'].getData(),
                               np.array([[0., 0.01], [0.9, 0.]]))
    np.testing.assert_allclose(maplayers1['modelmax']['grid'].getData(),
                               np.array([[0., 0.7], [0.99, 0.]]))

    # Testing J_PGA_M with Uncertainty - Model Type 2
    np.testing.assert_allclose(maplayers2['model']['grid'].getData(),
                               np.array([[0., 0.01], [0.9, 0.]]), atol=0.01)

    # Testing RS_PGA_M with Uncertainty - Model Type 3
    np.testing.assert_allclose(maplayers3['model']['grid'].getData(),
                               np.array([[0., 0.5], [0.99, 0.]]), atol=0.01)

    # Testing RS_PGA_PGV with Uncertainty - Model Type 4
    np.testing.assert_allclose(maplayers4['model']['grid'].getData(),
                               np.array([[0., 0.7], [0.99, 0.]]), atol=0.01)
    np.testing.assert_allclose(maplayers4['modelmin']['grid'].getData(),
                               np.array([[0., 0.1], [0.9, 0.]]))
    np.testing.assert_allclose(maplayers4['modelmax']['grid'].getData(),
                               np.array([[0., 0.99], [0.99, 0.]]))
    np.testing.assert_allclose(maplayers4['pgvmin']['grid'].getData(),
                               np.array([[32.53, 44.72], [23.24, 71.58]]),
                               atol=0.01)

    # Testing default model and no uncertainty
    np.testing.assert_allclose(maplayers5['model']['grid'].getData(),
                               np.array([[0., 0.01], [0.9, 0.]]), atol=0.01)

    # delete reference
    del config['godt_2008']['shortref']
    # reshape bounds
    bounds = {'xmin': 0.5, 'xmax': 1.5, 'ymin': 0.5, 'ymax': 1.5}
    # bad uncertfile
    uncertfile2 = os.path.join(datadir, 'test_uncert_error.xml')

    maplayers6 = godt2008(shakefile, config, bounds=bounds,
                          uncertfile=uncertfile2)

    # Testing with changed config
    np.testing.assert_allclose(maplayers6['model']['grid'].getData(),
                               np.array([[0., 0.01], [0.9, 0.]]), atol=0.01)


if __name__ == "__main__":
    test_godt2008()
    print('godt2008 tests passed')
