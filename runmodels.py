#!/usr/bin/env python

#stdlib imports
from configobj import ConfigObj
import numpy as np

#local imports
from groundfailure.newmark import godt2008
from mapio.shake import ShakeGrid
from groundfailure.makemaps import modelMap

configfile = '/Users/kallstadt/SecondaryHazards/Codes/inputs/confignew.ini'
shakefile = '/Users/kallstadt/SecondaryHazards/haywired.xml'  # URL or filename

edict = ShakeGrid.load(shakefile).getEventDict()
temp = ShakeGrid.load(shakefile).getShakeDict()
edict['eventid'] = temp['shakemap_id']
edict['version'] = temp['shakemap_version']

config = ConfigObj(configfile)

maplayers = godt2008(shakefile, config, saveinputs=True, regressionmodel='J_PGA')  # Note that there are many options for mapping and saving outputs, being run as defaults here by giving no inputs

probonly = {}
probonly['probability'] = maplayers['probability']
modelMap(probonly, edict, config, modelname='Godt et al. 2008_prob', maproads=True, mapcities=True, boundaries='zoom', scaletype='binned', lims=[[0, 0.2, 0.4, 0.6, 0.8, 1.]])


modelMap(maplayers, edict, config, modelname='Godt et al. 2008', maproads=True, mapcities=True, boundaries='zoom', scaletype='binned', lims=[[0, 0.2, 0.4, 0.6, 0.8, 1.], np.linspace(0., 2., 15), None, None, None])
