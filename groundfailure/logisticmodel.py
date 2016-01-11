#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys
import os.path
from xml.dom import minidom
from neicmap.distance import sdist
import ConfigParser
import re
from neicio.shake import ShakeGrid
from neicio import gmt

PARAM_PATTERN = 'b[0-9]+'
LAYER_PATTERN = '_layer'
TERM_PATTERN = 'term'

SM_TERMS = ['MW','PGA','PGV','MMI']
SM_GRID_TERMS = ['PGA','PGV','MMI']
OPERATORS = ['log','log10','power','sqrt'] #these will get np. prepended
FLOATPAT = '[+-]?(?=\d*[.eE])(?=\.?\d)\d*\.?\d*(?:[eE][+-]?\d+)?'
INTPAT = '[0-9]+'
OPERATORPAT = '[\+\-\*\/]*'

def getModelNames(configfile):
    config = ConfigParser.RawConfigParser()
    config.read(configfile)
    models = []
    for section in config.sections():
        if section.find('MODEL') > -1:
            models.append(section.split('_')[0].lower())
    return models
    

class LogisticModel(object):
    def __init__(self,configfile,shakefile,model):
        if model not in getModelNames(configfile):
            raise Exception,'Could not find a model called "%s" in config file %s.' % (model,configfile)
        #do everything here short of calculations - parse config, assemble eqn strings, load data.
        config = ConfigParser.RawConfigParser()
        config.read(configfile)
        self.model = model
        self.coeffs = {}
        self.layers = {} #key = layer name, value = file name
        self.terms = {}
        self.layerdict = {} #key = layer name, value = grid object

        section = model.upper() + '_MODEL'
        options = config.options(section)
        for option in options:
            isCoeff = re.search(PARAM_PATTERN,option) is not None
            isTerm = re.search(TERM_PATTERN,option) is not None
            isLayer = re.search(LAYER_PATTERN,option) is not None
            if isLayer:
                layername = option.split('_')[0]
                self.layers[layername] = config.get(section,option)
            if isCoeff and not isTerm:
                try:
                    self.coeffs[option] = float(config.get(section,option))
                except:
                    pass
                if not option.endswith('0'):
                    termname = option+'_term'
                    if not config.has_option(section,termname):
                        tpl = (configfile,model,option)
                        msg = 'Input config file %s, model %s is missing the term for coefficient %s.'
                        raise Exception,msg % tpl
                    self.terms[option] = config.get(section,termname)

        #load the layers
        #first temporarily load the shakemap so we can get the bounds
        baselayer = config.get(section,'baselayer')
        geodict = ShakeGrid(shakefile,variable='MMI').geodict
        baselayerfile = self.layers[baselayer]
        #make the bounds for the base grid smaller than the shakemap
        basebounds = (geodict['xmin']+(geodict['xdim']*2),
                      geodict['xmax']-(geodict['xdim']*2),
                      geodict['ymin']+(geodict['ydim']*2),
                      geodict['ymax']-(geodict['ydim']*2))
        #make the bounds for the other grids the same as the shakemap
        layerbounds = (geodict['xmin'],
                       geodict['xmax'],
                       geodict['ymin'],
                       geodict['ymax'])
        basegrid = gmt.GMTGrid(baselayerfile,bounds=basebounds)

        self.layerdict[baselayer] = basegrid
        layernames = self.layers.keys()
        layernames.remove(baselayer)
        for layername in layernames:
            layerfile = self.layers[layername]
            layergrid = gmt.GMTGrid(layerfile,bounds=layerbounds)
            try:
                layergrid.interpolateToGrid(basegrid.geodict)
            except Exception,msg:
                pass
            self.layerdict[layername] = layergrid
                
        if not self.coeffs.has_key('b0'):
            raise Exception,'Input config file %s, model %s is missing the intercept coefficient b0.' % (configfile,model)
        self.nuggets = [str(self.coeffs['b0'])]
        ckeys = self.terms.keys()
        ckeys.sort()
        for key in ckeys:
            term = self.terms[key]
            term,rem = self.checkTerm(term) #replace log with np.log, make sure variables are all in layers list, etc.
            if len(rem):
                msg = 'Term "%s" contains the unknown text fragment "%s".  This may cause the expression to fail.\n'
                tpl = (term,rem)
                sys.stderr.write(msg % tpl)
            coeff = self.coeffs[key]
            self.nuggets.append('(%g * %s)' % (coeff,term))

        self.shakedict = self.getShakeMapParams(config,shakefile,basegrid.geodict)

        self.equation = ' + '.join(self.nuggets)

    def getEquation(self):
        return self.equation

    def calculate(self):
        X = eval(self.equation)
        P = 1/(1 + np.exp(-X))
        return P

    def getShakeMapParams(self,config,shakefile,geodict):
        #get the shakemap params the user wants
        smparams = config.get('SHAKEMAP','variables').split(',')

        #load the shakemap
        intset = set(SM_TERMS).intersection(smparams)
        if not len(intset):
            print 'The allowed ShakeMap variables are: "%s". Your config file has "%s".' % (str(SM_TERMS),str(smparams))
            sys.exit(1)
        shakedict = {}
        shakemap = None
        if 'PGA' in smparams:
            shakemap = ShakeGrid(shakefile,variable='PGA')
            shakemap.interpolateToGrid(geodict)
            tmpgrid = gmt.GMTGrid()
            tmpgrid.loadFromGrid(shakemap)
            shakedict['PGA'] = tmpgrid
        if 'PGV' in smparams:
            shakemap = ShakeGrid(shakefile,variable='PGV')
            shakemap.interpolateToGrid(geodict)
            tmpgrid = gmt.GMTGrid()
            tmpgrid.loadFromGrid(shakemap)
            shakedict['PGV'] = tmpgrid
        if 'MMI' in smparams:
            shakemap = ShakeGrid(shakefile,variable='MMI')
            shakemap.interpolateToGrid(geodict)
            tmpgrid = gmt.GMTGrid()
            tmpgrid.loadFromGrid(shakemap)
            shakedict['MMI'] = tmpgrid
        if 'MW' in smparams:
            if shakemap is None:
                shakemap = ShakeGrid(shakefile,variable='MMI')
            attdict = shakemap.getAttributes()
            shakedict['MW'] = attdict['event']['magnitude']

        return (shakedict)

    def checkTerm(self,term):
        layers = self.layerdict.keys()
        model = self.model
        #Strip out everything that isn't: 0-9.() operators, +-/* or layer names.  Anything left is an unknown symbol.
        tterm = term
        #remove log, sqrt, etc.
        for op in OPERATORS:
            tterm = tterm.replace(op,'')
        #remove ShakeMap variables
        for sm_term in SM_TERMS:
            tterm = tterm.replace(sm_term,'')
        #remove layer names
        for layer in layers:
            tterm = tterm.replace(layer,'')
        #remove arithmetic operators
        tterm = re.sub(OPERATORPAT,'',tterm)
        #remove floating point numbers
        tterm = re.sub(FLOATPAT,'',tterm)
        #remove integer numbers
        tterm = re.sub(INTPAT,'',tterm)
        #remove parentheses
        tterm = re.sub('[()]*','',tterm)
        #remove any blank spaces
        tterm = tterm.strip()
        #remove commas
        tterm = tterm.strip(',')
        #anything left *might* cause an error
        for op in OPERATORS:
            if term.find(op) > -1:
                term = term.replace(op,'np.'+op)
        for sm_term in SM_GRID_TERMS:
            term = term.replace(sm_term,"self.shakedict['%s'].griddata" % sm_term)

        term = term.replace('MW',"self.shakedict['MW']")
            
        for layer in layers:
            term = term.replace(layer,"self.layerdict['%s'].griddata" % layer)
        return (term,tterm)

if __name__ == '__main__':
    shakefile = sys.argv[1]
    configfile = sys.argv[2]
    models = getModelNames(configfile)
    for model in models:
        lm = LogisticModel(configfile,shakefile,model)
        print 'Equation for the %s model:' % model
        print
        print lm.getEquation()
        print
        P = lm.calculate()
    
