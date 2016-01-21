#!/usr/bin/env python

#stdlib imports
from configobj import ConfigObj

#local imports
from groundfailure.newmark import godt2008
from mapio.shake import ShakeGrid

configfile = '/Users/kallstadt/SecondaryHazards/Codes/inputs/confignew.ini'
shakefile = '/Users/kallstadt/SecondaryHazards/Northridge.xml'  # URL or filename

edict = ShakeGrid.getEventDict(shakefile)
config = ConfigObj(configfile)

maplayers = godt2008(shakefile, config, saveinputs=True, regressionmodel='J_PGA')  # Note that there are many options for mapping and saving outputs, being run as defaults here by giving no inputs

makeMap(maplayers, edict, configfile, modelname='Godt et al. 2008', maproads=maproads, mapcities=mapcities, colormaps=colormaps, plotorder=plotorder, boundaries=boundaries1)
        #saveMap(maplayers, edict, modelname)

# Need to add these to config file
#plotorder = []
#colormaps = []