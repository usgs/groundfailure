#!/usr/bin/env python

#local imports
from groundfailure.godt2008 import runmodel

configfile = '/Users/kallstadt/SecondaryHazards/Codes/inputs/confignew.ini'
shakefile = '/Users/kallstadt/SecondaryHazards/Northridge.xml'  # URL or filename

maplayers, edict = runmodel(shakefile, configfile, mapshaking=True, mapmaxslope=True, mapcohesion=True, mapfriction=True, bounds='shakemap')  # Note that there are many options for mapping and saving outputs, being run as defaults here by giving no inputs
