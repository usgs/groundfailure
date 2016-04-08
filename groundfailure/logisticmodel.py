#!/usr/bin/env python

#stdlib imports
import matplotlib.pyplot as plt
import numpy as np
import sys
import os.path
from xml.dom import minidom
import ConfigParser
import re
import tempfile
import warnings

#third party imports
from mapio.shake import ShakeGrid,getHeaderData
from mapio.gmt import GMTGrid
from mapio.gdal import GDALGrid
from mapio.grid2d import Grid2D

PARAM_PATTERN = 'b[0-9]+'
LAYER_PATTERN = '_layer'
TERM_PATTERN = 'term'

SM_TERMS = ['MW','YEAR','MONTH','DAY','HOUR','pga','pgv','mmi']
SM_GRID_TERMS = ['pga','pgv','mmi']
OPERATORS = ['log','log10','power','sqrt'] #these will get np. prepended
FLOATPAT = '[+-]?(?=\d*[.eE])(?=\.?\d)\d*\.?\d*(?:[eE][+-]?\d+)?'
INTPAT = '[0-9]+'
OPERATORPAT = '[\+\-\*\/]*'
MONTHS = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

def getLogisticModelNames(config):
    names = []
    lmodel_space = config['logistic_models']
    for key,value in lmodel_space.iteritems():
        if isinstance(value,str) or isinstance(value,unicode):
            continue
        else: #this is a model
            names.append(key)
    return names

def getFileType(filename):
    if os.path.isdir(filename):
        return 'dir'
    ftype = GMTGrid.getFileType(filename)
    if ftype != 'unknown':
        return 'gmt'
    #skip over ESRI header files
    if filename.endswith('.hdr'):
            return 'unknown'
    try:
        GDALGrid.getFileGeoDict(filename)
        return 'esri'
    except:
        pass
    return 'unknown'

def getAllGridFiles(indir):
    tflist = os.listdir(indir)
    flist = []
    for tf in tflist:
        fullfile = os.path.join(indir,tf)
        ftype = getFileType(fullfile)
        if ftype in ['gmt','esri']:
            flist.append(fullfile)
    return flist

def validateCoefficients(cmodel):
    coeffs = {}
    for key,value in cmodel['coefficients'].iteritems():
        if re.search('b[0-9]*',key) is None:
            raise Exception('coefficients must be named b0, b1, ...')
        coeffs[key] = float(value)
    if 'b0' not in coeffs.keys():
        raise Exception('coefficients must include an intercept coefficient named b0.')
    return coeffs

def validateLayers(cmodel):
    layers = {}
    for key,value in cmodel['layers'].iteritems():
        ftype = getFileType(value)
        if ftype == 'unknown':
            raise Exception('layer file %s is not a valid GMT or ESRI file.' % value)
        if ftype == 'dir':
            value = getAllGridFiles(value)
        layers[key] = value
    return layers

def validateTerms(cmodel,coeffs,layers):
    #TODO - return a time field for every term, not just one global one.
    terms = {}
    timeField = None
    for key,value in cmodel['terms'].iteritems():
        if key not in coeffs.keys():
            raise Exception('Term names must match names of coefficients')
        term,rem,tTimeField = checkTerm(value,layers) #replace log with np.log, make sure variables are all in layers list, etc.
        if tTimeField is not None:
            timeField = tTimeField
        if len(rem):
            msg = 'Term "%s" contains the unknown text fragment "%s".  This may cause the expression to fail.\n'
            tpl = (term,rem)
            raise Exception(msg % tpl)
        terms[key] = term
    return (terms,timeField)

def validateInterpolations(cmodel,layers):
    interpolations = {}
    for key,value in cmodel['interpolations'].iteritems():
        if key not in layers.keys():
            raise Exception('Interpolation key %s does not match any names of layers' % key)
        methods = ['linear','nearest','cubic']
        if value not in methods:
            raise Exception('Interpolation method %s not in approved list of methods: %s' % (key,str(methods)))
        interpolations[key] = value
    for key in layers.keys():
        if key not in interpolations.keys():
            raise Exception('No interpolation method configured for layer %s' % key)
    return interpolations

def validateUnits(cmodel,layers):
    units = {}
    for key,value in cmodel['units'].iteritems():
        if key not in layers.keys():
            raise Exception('Interpolation key %s does not match any names of layers' % key)
        
        units[key] = value
    for key in layers.keys():
        if key not in units.keys():
            raise Exception('No unit string configured for layer %s' % key)
    return units

def validateLogisticModels(config):
    mnames = getLogisticModelNames(config)
    for cmodelname in mnames:
        try:
            cmodel = config['logistic_models'][cmodelname]
            coeffs = validateCoefficients(cmodel)
            layers = validateLayers(cmodel)#key = layer name, value = file name
            terms,timeField = validateTerms(cmodel,coeffs,layers)
            if timeField is not None:
                for (layer,layerfile) in layers.items():
                    if isinstance(layerfile,list):
                        for lfile in layerfile:
                            if timeField == 'MONTH':
                                pass
            interpolations = validateInterpolations(cmodel,layers)
            if cmodel['baselayer'] not in layers:
                raise Exception('Model %s missing baselayer parameter.' % cmodelname)
        except Exception as e:
            raise Exception('Validation failed with error: "%s" on model %s' % (str(e),cmodelname))
        
    return True

def checkTerm(term,layers):
    startterm = term
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
        term = term.replace(sm_term,"self.shakemap.getLayer('%s').getData()" % sm_term)

    #replace the macro MW with the magnitude value from the shakemap
    term = term.replace('MW',"self.shakemap.getEventDict()['magnitude']")

    #term.replace('YEAR',"self.shakemap.getEventDict()['event_time'].year")
    hasTime = False
    timeField = None
    for unit in ['YEAR','MONTH','DAY','HOUR']:
        if term.find(unit) > -1:
            term = term.replace(unit,'')
            timeField = unit

    for layer in layers:
        term = term.replace(layer,"self.layerdict['%s'].getData()" % layer)
    return (term,tterm,timeField)

class LogisticModel(object):
    def __init__(self,config,shakefile,model):
        if model not in getLogisticModelNames(config):
            raise Exception,'Could not find a model called "%s" in config %s.' % (model,config)
        #do everything here short of calculations - parse config, assemble eqn strings, load data.
        self.model = model
        cmodel = config['logistic_models'][model]
        self.coeffs = validateCoefficients(cmodel)
        self.layers = validateLayers(cmodel)#key = layer name, value = file name
        self.terms,timeField = validateTerms(cmodel,self.coeffs,self.layers)
        self.interpolations = validateInterpolations(cmodel,self.layers)
        self.units = validateUnits(cmodel,self.layers)

        if not cmodel.has_key('baselayer'):
            raise Exception('You must specify a base layer file in config.')
        if cmodel['baselayer'] not in self.layers.keys():
            raise Exception('You must specify a base layer corresponding to one of the files in the layer section.')
        
        #get the geodict for the shakemap
        geodict = ShakeGrid.getFileGeoDict(shakefile,adjust='res')
        griddict,eventdict,specdict,fields,uncertainties = getHeaderData(shakefile)
        YEAR = eventdict['event_timestamp'].year
        MONTH = MONTHS[(eventdict['event_timestamp'].month)-1]
        DAY = eventdict['event_timestamp'].day
        HOUR = eventdict['event_timestamp'].hour

        #now find the layer that is our base layer and get the largest bounds we can guaranteed not to exceed shakemap bounds
        basefile = self.layers[cmodel['baselayer']]
        ftype = getFileType(basefile)
        if ftype == 'esri':
            basegeodict = GDALGrid.getFileGeoDict(basefile)
            sampledict = basegeodict.getBoundsWithin(geodict)
        elif ftype == 'gmt':
            basegeodict = GMTGrid.getFileGeoDict(basefile)
            sampledict = basegeodict.getBoundsWithin(geodict)
        else:
            raise Exception('All predictor variable grids must be a valid GMT or ESRI file type')

        #now load the shakemap, resampling and padding if necessary
        self.shakemap = ShakeGrid.load(shakefile,samplegeodict=sampledict,resample=True,doPadding=True,adjust='res')
        
        #load the predictor layers into a dictionary
        self.layerdict = {} #key = layer name, value = grid object
        for layername,layerfile in self.layers.iteritems():
            if isinstance(layerfile,list):
                for lfile in layerfile:
                    if timeField == 'MONTH':
                        if lfile.find(MONTH) > -1:
                            layerfile = lfile
                            ftype = getFileType(layerfile)
                            interp = self.interpolations[layername]
                            if ftype == 'gmt':
                                lyr = GMTGrid.load(layerfile,sampledict,resample=True,method=interp,doPadding=True)
                            elif ftype == 'esri':
                                lyr = GDALGrid.load(layerfile,sampledict,resample=True,method=interp,doPadding=True)
                            else:
                                msg = 'Layer %s (file %s) does not appear to be a valid GMT or ESRI file.' % (layername,layerfile)
                                raise Exception(msg)
                            self.layerdict[layername] = lyr
            else:
                #first, figure out what kind of file we have (or is it a directory?)
                ftype = getFileType(layerfile)
                interp = self.interpolations[layername]
                if ftype == 'gmt':
                    lyr = GMTGrid.load(layerfile,sampledict,resample=True,method=interp,doPadding=True)
                elif ftype == 'esri':
                    lyr = GDALGrid.load(layerfile,sampledict,resample=True,method=interp,doPadding=True)
                else:
                    msg = 'Layer %s (file %s) does not appear to be a valid GMT or ESRI file.' % (layername,layerfile)
                    raise Exception(msg)
                self.layerdict[layername] = lyr

        shapes = {}
        for layername,layer in self.layerdict.iteritems():
            shapes[layername] = layer.getData().shape

        x = 1
        self.nuggets = [str(self.coeffs['b0'])]
        ckeys = self.terms.keys()
        ckeys.sort()
        for key in ckeys:
            term = self.terms[key]
            coeff = self.coeffs[key]
            self.nuggets.append('(%g * %s)' % (coeff, term))

        self.equation = ' + '.join(self.nuggets)
        self.geodict = self.shakemap.getGeoDict()

    def getEquation(self):
        return self.equation

    def getGeoDict(self):
        return self.geodict

    def calculate(self):
        X = eval(self.equation)
        P = 1/(1 + np.exp(-X))
        Pgrid = Grid2D(P,self.geodict)
        rdict = {'model':{'grid':Pgrid,
                                'label':'Probability',
                                'type':'output',
                                'description': {'units':'probability'}}}
        for layername,layergrid in self.layerdict.items():
            units = self.units[layername]
            rdict[layername] = {'grid':layergrid,
                                'label':'%s (%s)' % (layername,units),
                                'type':'input',
                                'description': {'units': units}}
        return rdict

    

def _test(shakefile,cofile,slopefile,precipfolder):
    model = {'logistic_models':{'nowicki_2014':{'description':'This is the Nowicki Model of 2014, which uses cohesion and slope max as input.',
                                                'gfetype':'landslide',
                                                'baselayer':'cohesion',
                                                'layers':{'cohesion':'%s' % cofile,
                                                          'slope':'%s' % slfile,
                                                          'precip':'%s' % precipfolder},
                                                'interpolations':{'cohesion':'linear',
                                                                  'slope':'linear',
                                                                  'precip':'nearest'},
                                                'terms':{'b1':'pga',
                                                         'b2':'slope',
                                                         'b3':'precipMONTH',
                                                         'b4':'pga*slope*MW'},
                                                'coefficients':{'b0':-7.15,
                                                                'b1':0.0604,
                                                                'b2':0.000825,
                                                                'b3':0.0201,
                                                                'b4':1.45e-05}}}}

    lm = LogisticModel(model,shakefile,'nowicki_2014')
    print lm.getEquation()
    P = lm.calculate()
    
    
if __name__ == '__main__':
    shakefile = sys.argv[1] #needs to be an event occurring in January
    cofile = sys.argv[2]
    slfile = sys.argv[3]
    precip = sys.argv[4]
    _test(shakefile,cofile,slfile,precip)

