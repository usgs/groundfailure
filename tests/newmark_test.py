#!/usr/bin/env python

import os.path
import os
from configobj import ConfigObj
import numpy as np
import groundfailure.newmark as NM
from mapio.geodict import GeoDict
from groundfailure.conf import correct_config_filepaths

#hack the path so that I can debug these functions if I need to
#homedir = '/Users/kallstadt/SecondaryHazards/Codes/groundfailure/tests'
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


def test_hazus():
    configfile = os.path.join(datadir, 'testconfig_hazus.ini')
    config = ConfigObj(configfile)
    # Test path correction (from conf.py)
    config = correct_config_filepaths(datadir, config)
    bounds = {'xmin': -1.5, 'xmax': 3, 'ymin': -1.5, 'ymax': 3}
    maplayers1 = NM.hazus(shakefile, config, modeltype='coverage', saveinputs=True, bounds=bounds, uncertfile=uncertfile)
    maplayers2 = NM.hazus(shakefile, config, modeltype='dn_hazus', uncertfile=uncertfile)
    maplayers3 = NM.hazus(shakefile, config, modeltype='dn_prob', probtype='jibson2000')
    maplayers4 = NM.hazus(shakefile, config, modeltype='ac_classic_dn', displmodel='RS_PGA_PGV', uncertfile=uncertfile,
                          probtype='jibson2000', saveinputs=True)
    maplayers5 = NM.hazus(shakefile, config, modeltype='ac_classic_prob', displmodel='RS_PGA_M',
                          probtype='threshold', uncertfile=uncertfile)
    maplayers6 = NM.hazus(shakefile, config, modeltype='ac_classic_prob', uncertfile=uncertfile, probtype='jibson2000',
                          saveinputs=True)
    maplayers7 = NM.hazus(shakefile, config, modeltype='ac_classic_prob', uncertfile=uncertfile, probtype='errortype',
                          saveinputs=True)

    np.testing.assert_allclose(maplayers1['model']['grid'].getData(),
                               np.array([[0., 0., 0.1],
                                        [0., 0., 0.],
                                        [0.3, 0., 0.3]]))  # NEED TO FIX NAN's ONCE MAPIO IS FIXED
    np.testing.assert_allclose(maplayers1['Ac']['grid'].getData(),
                               np.array([[0.4, 0.35, 0.25],
                                        [0.4, 0.35, 0.25],
                                        [0.05, 0.05, 0.05]]))
    np.testing.assert_allclose(maplayers1['modelmin']['grid'].getData(),
                               np.array([[0., 0., 0.1],
                                        [0., 0., 0.],
                                        [0.3, 0., 0.3]]))  # NEED TO FIX NAN's ONCE MAPIO IS FIXED
    np.testing.assert_allclose(maplayers1['modelmax']['grid'].getData(),
                               np.array([[0., 0., 0.1],
                                        [0., 0., 0.],
                                        [0.3, 0., 0.3]]))  # NEED TO FIX NAN's ONCE MAPIO IS FIXED
    np.testing.assert_allclose(maplayers2['model']['grid'].getData(),
                               np.array([[0., np.nan, 25.9],
                                        [np.nan, np.nan, np.nan],
                                        [97.2, np.nan, 1.65]]), rtol=0.5)
    np.testing.assert_allclose(maplayers3['model']['grid'].getData(),
                               np.array([[0., np.nan, 0.334865],
                                        [np.nan, np.nan, np.nan],
                                        [0.3350, np.nan, 0.0422]]), atol=0.01)
    np.testing.assert_allclose(maplayers4['model']['grid'].getData(),
                               np.array([[4.6e-39, 0., 13.84],
                                        [0., 0., 0.],
                                        [37.97, 0., 25.368]]), atol=0.1)
    np.testing.assert_allclose(maplayers4['modelmin']['grid'].getData(),
                               np.array([[0., 0., 1.78],
                                        [0., 0., 0.],
                                        [21.69, 0., 12.48]]), atol=0.1)
    np.testing.assert_allclose(maplayers4['pgv']['grid'].getData(),
                               np.array([[50, np.nan, 73.],
                                        [np.nan, np.nan, np.nan],
                                        [32., np.nan, 84.]]), atol=0.1)
    np.testing.assert_allclose(maplayers5['model']['grid'].getData(),
                               np.array([[0, 0, 1],
                                        [0, 0, 0],
                                        [1, 0, 0]]), atol=0)
    np.testing.assert_allclose(maplayers6['model']['grid'].getData(),
                               np.array([[0., 0., 0.02],
                                        [0., 0., 0.],
                                        [0.33, 0., 0.01]]), atol=0.1)
    np.testing.assert_allclose(maplayers6['modelmin']['grid'].getData(),
                               np.array([[0., 0., 0.],
                                        [0., 0., 0.],
                                        [0.28, 0., 0.]]), atol=0.1)
    np.testing.assert_allclose(maplayers7['model']['grid'].getData(),
                               np.array([[0., 0., 0.02],
                                        [0., 0., 0.],
                                        [0.33, 0., 0.01]]), atol=0.1)
    np.testing.assert_allclose(maplayers7['modelmin']['grid'].getData(),
                               np.array([[0., 0., 0.],
                                        [0., 0., 0.],
                                        [0.28, 0., 0.]]), atol=0.1)

    # Recut grid for testing
    bounds = {'xmin': 0.5, 'xmax': 1.5, 'ymin': 0.5, 'ymax': 1.5}
    # Null config values
    del config['hazus']['shortref']
    del config['hazus']['parameters']['dnthresh']
    # Set uncertfile to null location
    uncertfile2 = os.path.join(datadir, 'test_uncert_error.xml')

    maplayers8 = NM.hazus(shakefileBounds, config, bounds=bounds, probtype='threshold', uncertfile=uncertfile2)
    np.testing.assert_allclose(maplayers8['model']['grid'].getData(),
                               np.array([[0., 0., 0.1],
                                        [0., 0., 0.],
                                        [0.3, 0., 0.3]]), atol=0.1)


def test_est_disp():
    """This test compares against values pulled from Figure 4.14 in HAZUS manual
    """
    ac = 0.4
    pga = 1.
    ed_low = 3.
    ed_high = 5.
    el, eh = NM.est_disp(ac, pga)
    np.testing.assert_allclose(ed_low, el, atol=0.3)
    np.testing.assert_allclose(ed_high, eh, atol=0.3)
    el, eh = NM.est_disp(0.01, 1.)
    np.testing.assert_allclose(20., el, atol=0.3)
    np.testing.assert_allclose(40., eh, atol=0.3)


def test_classic():
    configfile = os.path.join(datadir, 'testconfig_classic.ini')
    config = ConfigObj(configfile)
    bounds = {'xmin': -1.5, 'xmax': 3, 'ymin': -1.5, 'ymax': 3}
    # Test path correction (from conf.py)
    config = correct_config_filepaths(datadir, config)
    maplayers1 = NM.classic(shakefile, config, bounds=bounds, displmodel='J_PGA', uncertfile=uncertfile)
    maplayers2 = NM.classic(shakefile, config, displmodel='J_PGA_M', saveinputs=True)
    maplayers3 = NM.classic(shakefile, config, displmodel='RS_PGA_M')
    maplayers4 = NM.classic(shakefile, config, displmodel='RS_PGA_PGV', uncertfile=uncertfile, probtype='threshold', saveinputs=True)
    maplayers5 = NM.classic(shakefile, config, probtype='thressshold', uncertfile=uncertfile)

    # Testing J_PGA with Uncertainty - Model Type 1
    np.testing.assert_allclose(maplayers1['model']['grid'].getData(),
                               np.array([[0.004, 0.335], [0.334, 0.]]), atol=0.01)
    np.testing.assert_allclose(maplayers1['modelmin']['grid'].getData(),
                               np.array([[0., 0.334], [0.334, 0.]]), atol=0.01)
    np.testing.assert_allclose(maplayers1['modelmax']['grid'].getData(),
                               np.array([[0.0459, 0.335], [0.335, 0.]]), atol=0.01)

    # Testing J_PGA_M without Uncertainty - Model Type 2
    np.testing.assert_allclose(maplayers2['model']['grid'].getData(),
                               np.array([[0.001, 0.3314], [0.3349, 0.]]), atol=0.01)

    # Testing RS_PGA_M without Uncertainty - Model Type 3
    np.testing.assert_allclose(maplayers3['model']['grid'].getData(),
                               np.array([[0.014, 0.335], [0.335, 0.]]), atol=0.01)

    # Testing RS_PGA_PGV with Uncertainty - Model Type 4
    np.testing.assert_allclose(maplayers4['model']['grid'].getData(),
                               np.array([[0., 1.], [1., 0.]]), atol=0.01)
    np.testing.assert_allclose(maplayers4['modelmin']['grid'].getData(),
                               np.array([[0., 1.], [1., 0.]]), atol=0.01)
    np.testing.assert_allclose(maplayers4['modelmax']['grid'].getData(),
                               np.array([[1., 1.], [1., 0.0]]), atol=0.01)
    np.testing.assert_allclose(maplayers4['pgvmin']['grid'].getData(),
                               np.array([[32.53, 44.72], [23.24, 71.58]]), atol=0.01)

    # Testing default model, probtype breaking with uncertainty
    np.testing.assert_allclose(maplayers5['model']['grid'].getData(),
                               np.array([[0.004, 0.335], [0.334, 0.]]), atol=0.01)
    np.testing.assert_allclose(maplayers5['modelmin']['grid'].getData(),
                               np.array([[0., 0.311], [0.284, 0.]]), atol=0.01)
    np.testing.assert_allclose(maplayers5['modelmax']['grid'].getData(),
                               np.array([[0.013, 0.335], [0.335, 0.]]), atol=0.01)

    # for models with different configs
    del config['classic_newmark']['parameters']['dnthresh']
    del config['classic_newmark']['shortref']
    # add watertablefile
    config['classic_newmark']['parameters']['m'] = 'file'
    #reshape bounds
    bounds = {'xmin': 0.5, 'xmax': 1.5, 'ymin': 0.5, 'ymax': 1.5}
    # bad uncertfile
    uncertfile2 = os.path.join(datadir, 'test_uncert_error.xml')

    maplayers6 = NM.classic(shakefile, config, bounds=bounds, displmodel='J_PGA', probtype='threshold', saveinputs=True)
    maplayers7 = NM.classic(shakefile, config, displmodel='J_PGA', uncertfile=uncertfile2)

    # Testing J_PGA, new config, weird bounds - Model Type 6
    np.testing.assert_allclose(maplayers6['model']['grid'].getData(),
                               np.array([[0., 1.], [1., 0.]]), atol=0.01)
    np.testing.assert_allclose(maplayers6['sat thick prop']['grid'].getData(),
                               np.array([[0., 0.96], [0.88, 0.17]]), atol=0.01)

    # Testing J_PGA, new config, bad uncertainty file - Model Type 7
    np.testing.assert_allclose(maplayers7['model']['grid'].getData(),
                               np.array([[0.004, 0.335], [0.334, 0.]]), atol=0.01)


def test_godt2008():
    configfile = os.path.join(datadir, 'testconfig_godt.ini')
    config = ConfigObj(configfile)
    bounds = {'xmin': -1.5, 'xmax': 3, 'ymin': -1.5, 'ymax': 3}
    # Test path correction (from conf.py)
    config = correct_config_filepaths(datadir, config)
    maplayers1 = NM.godt2008(shakefile, config, displmodel='J_PGA', uncertfile=uncertfile, bounds=bounds)
    maplayers2 = NM.godt2008(shakefile, config, displmodel='J_PGA_M', saveinputs=True)
    maplayers3 = NM.godt2008(shakefile, config, displmodel='RS_PGA_M')
    maplayers4 = NM.godt2008(shakefile, config, displmodel='RS_PGA_PGV', uncertfile=uncertfile, saveinputs=True)
    maplayers5 = NM.godt2008(shakefile, config)

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
                               np.array([[32.53, 44.72], [23.24, 71.58]]), atol=0.01)

    # Testing default model and no uncertainty
    np.testing.assert_allclose(maplayers5['model']['grid'].getData(),
                               np.array([[0., 0.01], [0.9, 0.]]), atol=0.01)

    # delete reference
    del config['godt_2008']['shortref']
    #reshape bounds
    bounds = {'xmin': 0.5, 'xmax': 1.5, 'ymin': 0.5, 'ymax': 1.5}
    # bad uncertfile
    uncertfile2 = os.path.join(datadir, 'test_uncert_error.xml')

    maplayers6 = NM.godt2008(shakefile, config, bounds=bounds, uncertfile=uncertfile2)

    # Testing with changed config
    np.testing.assert_allclose(maplayers6['model']['grid'].getData(),
                               np.array([[0., 0.01], [0.9, 0.]]), atol=0.01)


def test_NMdisp():
    """Tests the implementation of various models to estimate Newmark displacement
    NEED TO ADD CHECK OF UNCERTAINTIES
    """

    #test J_PGA
    #This test compares against values pulled from Figure 2 in Jibson (2007)
    ac = 0.75
    amax = 1.
    logDn = -1.
    val, std, logtype = NM.NMdisp(ac, amax, model='J_PGA')
    np.testing.assert_allclose(val, 10**logDn, atol=0.05)
    # Try another on different part of curve
    ac = 0.2
    logDn = 1.
    val, std, logtype = NM.NMdisp(ac, amax, model='J_PGA')
    np.testing.assert_allclose(np.log10(val), logDn, atol=0.01)

    #test J_PGA_M
    #No data published to compare against so this model is based on solving the equation analytically
    # NEED TO ADD CHECK OF UNCERTAINTIES
    ac = 0.75
    amax = 1.
    M = 7.
    val, std, logtype = NM.NMdisp(ac, amax, model='J_PGA_M', M=M)
    assert(logtype == 'log10'), 'NMdisp did not return logtype expected for J_PGA_M'
    np.testing.assert_allclose(val, 0.1088552, atol=0.001)

    # test_RS_PGA_M():
    # This test compares against values plotted in Figure 7 of Rathje and Saygili 2009
    # NEED TO ADD CHECK OF UNCERTAINTIES
    ac = 0.2
    amax = 0.33
    M = 7.5
    Dn = 2.25
    val, std, logtype = NM.NMdisp(ac, amax, M=M, model='RS_PGA_M')
    np.testing.assert_allclose(val, Dn, atol=0.01)

    #test_RS_PGA_PGV():
    #This test compares against values plotted in Figure 7 of Rathje and Saygili 2009
    ac = 0.15
    amax = 0.33
    pgv = 30.
    Dn = 2.5
    val, std, logtype = NM.NMdisp(ac, amax, PGV=pgv, model='RS_PGA_PGV')
    np.testing.assert_allclose(val, Dn, atol=0.05)
    assert(logtype == 'ln'), 'NMdisp did not return logtype expected for J_PGA_M'

test_hazus()
test_est_disp()
test_classic()
test_godt2008()
test_NMdisp()
print('newmark.py tests passed')
