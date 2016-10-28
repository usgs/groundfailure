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
#logisticmodeldir = os.path.abspath(os.path.join(homedir, 'groundfailure'))
#sys.path.insert(0, logisticmodeldir)  # put this at the front of the system path, ignoring any installed mapio stuff

configfile = os.path.join(datadir, 'testconfig_newmark.ini')
config = ConfigObj(configfile)
# Test path correction (from conf.py)
config = correct_config_filepaths(datadir, config)

cmodel = config['test_model']
layers = []

shakefile = os.path.join(datadir, 'test_shakegrid.xml')
uncertfile = os.path.join(datadir, 'test_uncert.xml')
cofile = os.path.join(datadir, 'test_friction.bil')
slopefile = os.path.join(datadir, 'test_slope.bil')
vs30file = os.path.join(datadir, 'test_vs30.bil')
ctifile = os.path.join(datadir, 'test_cti.bil')
precipfolder = os.path.join(datadir, 'test_precip')

fakegeodict = GeoDict({'xmin': 0.0, 'xmax': 1.0,
                      'ymin': 0.0, 'ymax': 1.0,
                      'dx': 1.0, 'dy': 1.0,
                      'ny': 2, 'nx': 2})


def test_hazus():
    pass


def test_est_disp():
    pass


def test_classic():
    pass


def test_godt2008():
    pass


def test_J_PGA():
    pass


def test_J_PGA_M():
    pass


def test_RS_PGA_M():
    pass


def test_RS_PGA_PGV():
    pass
