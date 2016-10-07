#!/usr/bin/env python

#stdlib imports
import numpy as np
import sys
import os.path
import re
import collections
import copy

#third party imports
from mapio.shake import ShakeGrid, getHeaderData
from mapio.gmt import GMTGrid
from mapio.gdal import GDALGrid
from mapio.grid2d import Grid2D

PARAM_PATTERN = 'b[0-9]+'
LAYER_PATTERN = '_layer'
TERM_PATTERN = 'term'

SM_TERMS = ['MW', 'YEAR', 'MONTH', 'DAY', 'HOUR', 'pga', 'pgv', 'mmi']
SM_GRID_TERMS = ['pga', 'pgv', 'mmi']
OPERATORS = ['log', 'log10', 'power', 'sqrt', 'minimum']  # these will get np. prepended
FLOATPAT = '[+-]?(?=\d*[.eE])(?=\.?\d)\d*\.?\d*(?:[eE][+-]?\d+)?'
INTPAT = '[0-9]+'
OPERATORPAT = '[\+\-\*\/]*'
MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']


def getLogisticModelNames(config):
    names = []
    lmodel_space = config['logistic_models']
    for key, value in lmodel_space.items():
        if isinstance(value, str):
            continue
        else:  # this is a model
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
        fullfile = os.path.join(indir, tf)
        ftype = getFileType(fullfile)
        if ftype in ['gmt', 'esri']:
            flist.append(fullfile)
    return flist


def validateCoefficients(cmodel):
    coeffs = {}
    for key, value in cmodel['coefficients'].items():
        if re.search('b[0-9]*', key) is None:
            raise Exception('coefficients must be named b0, b1, ...')
        coeffs[key] = float(value)
    if 'b0' not in list(coeffs.keys()):
        raise Exception('coefficients must include an intercept coefficient named b0.')
    return coeffs


def validateLayers(cmodel):
    layers = {}
    for key in cmodel['layers'].keys():
        for item, value in cmodel['layers'][key].items():
            if item == 'file':
                ftype = getFileType(value)
                if ftype == 'unknown':
                    raise Exception('layer file %s is not a valid GMT or ESRI file.' % value)
                if ftype == 'dir':
                    value = getAllGridFiles(value)
                layers[key] = value
    return layers


def validateTerms(cmodel, coeffs, layers):
    #TODO - return a time field for every term, not just one global one.
    terms = {}
    timeField = None
    for key, value in cmodel['terms'].items():
        if key not in list(coeffs.keys()):
            raise Exception('Term names must match names of coefficients')
        # replace log with np.log, make sure variables are all in layers list, etc.
        term, rem, tTimeField = checkTerm(value, layers)
        print(term)
        if tTimeField is not None:
            timeField = tTimeField
        if len(rem):
            msg = 'Term "%s" contains the unknown text fragment "%s".  This may cause the expression to fail.\n'
            tpl = (term, rem)
            raise Exception(msg % tpl)
        terms[key] = term
    return (terms, timeField)


def validateInterpolations(cmodel, layers):
    interpolations = {}
    for key, value in cmodel['interpolations'].items():
        if key not in list(layers.keys()):
            raise Exception('Interpolation key %s does not match any names of layers' % key)
        methods = ['linear', 'nearest', 'cubic']
        if value not in methods:
            raise Exception('Interpolation method %s not in approved list of methods: %s' % (key, str(methods)))
        interpolations[key] = value
    for key in list(layers.keys()):
        if key not in list(interpolations.keys()):
            raise Exception('No interpolation method configured for layer %s' % key)
    return interpolations


def validateUnits(cmodel, layers):
    units = {}
    for key in cmodel['layers'].keys():
        if 'units' in cmodel['layers'][key]:
            units[key] = cmodel['layers'][key]['units']
        else:
            raise Exception('No unit string configured for layer %s' % key)
    return units


def validateLogisticModels(config):
    mnames = getLogisticModelNames(config)
    for cmodelname in mnames:
        try:
            cmodel = config['logistic_models'][cmodelname]
            coeffs = validateCoefficients(cmodel)
            layers = validateLayers(cmodel)  # key = layer name, value = file name
            terms, timeField = validateTerms(cmodel, coeffs, layers)
            if timeField is not None:
                for (layer, layerfile) in list(layers.items()):
                    if isinstance(layerfile, list):
                        for lfile in layerfile:
                            if timeField == 'MONTH':
                                pass
            interpolations = validateInterpolations(cmodel, layers)
            if cmodel['baselayer'] not in layers:
                raise Exception('Model %s missing baselayer parameter.' % cmodelname)
        except Exception as e:
            raise Exception('Validation failed with error: "%s" on model %s' % (str(e), cmodelname))

    return True


def validateRefs(cmodel):
    longrefs = {}
    shortrefs = {}
    modelrefs = {}
    for key in cmodel['layers'].keys():
        if 'longref' in cmodel['layers'][key]:
            longrefs[key] = cmodel['layers'][key]['longref']
        else:
            print('No longref provided for layer %s' % key)
            longrefs[key] = 'unknown'
        if 'shortref' in cmodel['layers'][key]:
            shortrefs[key] = cmodel['layers'][key]['shortref']
        else:
            print('No shortref provided for layer %s' % key)
            shortrefs[key] = 'unknown'
    try:
        modelrefs['longref'] = cmodel['longref']
    except:
        print('No model longref provided')
        modelrefs['longref'] = 'unknown'
    try:
        modelrefs['shortref'] = cmodel['shortref']
    except:
        print('No model shortref provided')
        modelrefs['shortref'] = 'unknown'
    return modelrefs, longrefs, shortrefs


def checkTerm(term, layers):
    #startterm = term
    #Strip out everything that isn't: 0-9.() operators, +-/* or layer names.  Anything left is an unknown symbol.
    tterm = term
    #remove log, sqrt, etc.
    for op in OPERATORS:
        tterm = tterm.replace(op, '')
    #remove ShakeMap variables
    for sm_term in SM_TERMS:
        tterm = tterm.replace(sm_term, '')
    #remove layer names
    for layer in layers:
        tterm = tterm.replace(layer, '')
    #remove arithmetic operators
    tterm = re.sub(OPERATORPAT, '', tterm)
    #remove floating point numbers
    tterm = re.sub(FLOATPAT, '', tterm)
    #remove integer numbers
    tterm = re.sub(INTPAT, '', tterm)
    #remove parentheses
    tterm = re.sub('[()]*', '', tterm)
    #remove any blank spaces
    tterm = tterm.strip()
    #remove commas
    tterm = tterm.strip(',')
    #anything left *might* cause an error
    for op in OPERATORS:
        if term.find(op) > -1:
            term = term.replace(op, 'np.'+op)

    for sm_term in SM_GRID_TERMS:
        term = term.replace(sm_term, "self.shakemap.getLayer('%s').getData()" % sm_term)

    #replace the macro MW with the magnitude value from the shakemap
    term = term.replace('MW', "self.shakemap.getEventDict()['magnitude']")

    #term.replace('YEAR',"self.shakemap.getEventDict()['event_time'].year")
    #hasTime = False
    timeField = None
    for unit in ['YEAR', 'MONTH', 'DAY', 'HOUR']:
        if term.find(unit) > -1:
            term = term.replace(unit, '')
            timeField = unit

    for layer in layers:
        term = term.replace(layer, "self.layerdict['%s'].getData()" % layer)
    return (term, tterm, timeField)


class LogisticModel(object):
    def __init__(self, config, shakefile, model, uncertfile=None):
        if model not in getLogisticModelNames(config):
            raise Exception('Could not find a model called "%s" in config %s.' % (model, config))
        #do everything here short of calculations - parse config, assemble eqn strings, load data.

        self.model = model
        cmodel = config['logistic_models'][model]
        self.modeltype = cmodel['gfetype']
        self.coeffs = validateCoefficients(cmodel)
        self.layers = validateLayers(cmodel)  # key = layer name, value = file name
        self.terms, timeField = validateTerms(cmodel, self.coeffs, self.layers)
        self.interpolations = validateInterpolations(cmodel, self.layers)
        self.units = validateUnits(cmodel, self.layers)
        self.gmused = [value for term, value in cmodel['terms'].items() if 'pga' in value.lower() or 'pgv' in value.lower() or 'mmi' in value.lower()]
        self.modelrefs, self.longrefs, self.shortrefs = validateRefs(cmodel)
        if 'baselayer' not in cmodel:
            raise Exception('You must specify a base layer file in config.')
        if cmodel['baselayer'] not in list(self.layers.keys()):
            raise Exception('You must specify a base layer corresponding to one of the files in the layer section.')

        #get the geodict for the shakemap
        geodict = ShakeGrid.getFileGeoDict(shakefile, adjust='res')
        griddict, eventdict, specdict, fields, uncertainties = getHeaderData(shakefile)
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
        self.shakemap = ShakeGrid.load(shakefile, samplegeodict=sampledict, resample=True, doPadding=True, adjust='res')

        # take uncertainties into account
        if uncertfile is not None:
            try:
                self.uncert = ShakeGrid.load(uncertfile, samplegeodict=sampledict, resample=True, doPadding=True,
                                             method='linear', adjust='res')
            except:
                print('Could not read uncertainty file, ignoring uncertainties')
                self.uncert = None
        else:
            self.uncert = None

        #load the predictor layers into a dictionary
        self.layerdict = {}  # key = layer name, value = grid object
        for layername, layerfile in self.layers.items():
            if isinstance(layerfile, list):
                for lfile in layerfile:
                    if timeField == 'MONTH':
                        if lfile.find(MONTH) > -1:
                            layerfile = lfile
                            ftype = getFileType(layerfile)
                            interp = self.interpolations[layername]
                            if ftype == 'gmt':
                                lyr = GMTGrid.load(layerfile, sampledict, resample=True, method=interp, doPadding=True)
                            elif ftype == 'esri':
                                lyr = GDALGrid.load(layerfile, sampledict, resample=True, method=interp, doPadding=True)
                            else:
                                msg = 'Layer %s (file %s) does not appear to be a valid GMT or ESRI file.' % (layername, layerfile)
                                raise Exception(msg)
                            self.layerdict[layername] = lyr
            else:
                #first, figure out what kind of file we have (or is it a directory?)
                ftype = getFileType(layerfile)
                interp = self.interpolations[layername]
                if ftype == 'gmt':
                    lyr = GMTGrid.load(layerfile, sampledict, resample=True, method=interp, doPadding=True)
                elif ftype == 'esri':
                    lyr = GDALGrid.load(layerfile, sampledict, resample=True, method=interp, doPadding=True)
                else:
                    msg = 'Layer %s (file %s) does not appear to be a valid GMT or ESRI file.' % (layername, layerfile)
                    raise Exception(msg)
                self.layerdict[layername] = lyr

        shapes = {}
        for layername, layer in self.layerdict.items():
            shapes[layername] = layer.getData().shape

        self.nuggets = [str(self.coeffs['b0'])]

        ckeys = list(self.terms.keys())
        ckeys.sort()
        for key in ckeys:
            term = self.terms[key]
            coeff = self.coeffs[key]
            self.nuggets.append('(%g * %s)' % (coeff, term))

        self.equation = ' + '.join(self.nuggets)

        if self.uncert is not None:
            self.nugmin = copy.copy(self.nuggets)
            self.nugmax = copy.copy(self.nuggets)
            # Find the term with the shakemap input and replace for these nuggets
            for k, nug in enumerate(self.nuggets):
                if "self.shakemap.getLayer('pga').getData()" in nug:
                    self.nugmin[k].replace("self.shakemap.getLayer('pga').getData()", "(np.exp(np.log(self.shakemap.getLayer('pga').getData()) - self.uncert.getLayer('stdpga').getData()))")
                    self.nugmax[k].replace("self.shakemap.getLayer('pga').getData()", "(np.exp(np.log('self.shakemap.getLayer('pga').getData()) + self.uncert.getLayer('stdpga').getData()))")
                elif "self.layerdict['pgv'].getData()" in nug:
                    self.nugmin[k].replace("self.shakemap.getLayer('pgv').getData()", "(np.exp(np.log(self.shakemap.getLayer('pgv').getData()) - self.uncert.getLayer('stdpgv').getData()))")
                    self.nugmax[k].replace("self.shakemap.getLayer('pgv').getData()", "(np.exp(np.log('self.shakemap.getLayer('pgv').getData()) + self.uncert.getLayer('stdpgv').getData()))")
                elif "self.layerdict['mmi'].getData()" in nug:
                    self.nugmin[k].replace("self.shakemap.getLayer('mmi').getData()", "(np.exp(np.log(self.shakemap.getLayer('mmi').getData()) - self.uncert.getLayer('stdmmi').getData()))")
                    self.nugmax[k].replace("self.shakemap.getLayer('mmi').getData()", "(np.exp(np.log('self.shakemap.getLayer('mmi').getData()) + self.uncert.getLayer('stdmmi').getData()))")
            self.equationmin = ' + '.join(self.nugmin)
            self.equationmax = ' + '.join(self.nugmin)
        else:
            self.equationmin = None
            self.equationmax = None

        self.geodict = self.shakemap.getGeoDict()

        try:
            self.slopemin = float(config['logistic_models'][model]['slopemin'])
            self.slopemax = float(config['logistic_models'][model]['slopemax'])
        except:
            print('could not find slopemin and/or slopemax in config, no limits will be applied')
            self.slopemin = 0.
            self.slopemax = 90.

    def getEquation(self):
        return self.equation

    def getEquations(self):
        return self.equation, self.equationmin, self.equationmax

    def getGeoDict(self):
        return self.geodict

    def calculate(self, saveinputs=False, slopefile=None, slopediv=1.):
        """
        saveinputs - if True, saves all the input layers as Grid2D objects in addition to the model
          if false, will just output model
        slopefile - optional slopefile that will be resampled to the other input files for applying thresholds
        slopediv - number to divide slope by to get to degrees (usually will be default of 1.)
        """
        X = eval(self.equation)
        P = 1/(1 + np.exp(-X))
        if self.uncert is not None:
            Xmin = eval(self.equationmin)
            Xmax = eval(self.equationmax)
            Pmin = 1/(1 + np.exp(-Xmin))
            Pmax = 1/(1 + np.exp(-Xmax))
        if slopefile is not None:
            ftype = getFileType(slopefile)
            sampledict = self.shakemap.getGeoDict()
            if ftype == 'gmt':
                slope = GMTGrid.load(slopefile, sampledict, resample=True, method='linear', doPadding=True).getData()/slopediv
                # Apply slope min/max limits
                print('applying slope thresholds')
                P[slope > self.slopemax] = 0.
                P[slope < self.slopemin] = 0.
                if self.uncert is not None:
                    Pmin[slope > self.slopemax] = 0.
                    Pmin[slope < self.slopemin] = 0.
                    Pmax[slope > self.slopemax] = 0.
                    Pmax[slope < self.slopemin] = 0.
            elif ftype == 'esri':
                slope = GDALGrid.load(slopefile, sampledict, resample=True, method='linear', doPadding=True).getData()/slopediv
                # Apply slope min/max limits
                print('applying slope thresholds')
                P[slope > self.slopemax] = 0.
                P[slope < self.slopemin] = 0.
                if self.uncert is not None:
                    Pmin[slope > self.slopemax] = 0.
                    Pmin[slope < self.slopemin] = 0.
                    Pmax[slope > self.slopemax] = 0.
                    Pmax[slope < self.slopemin] = 0.
            else:
                print('Slope file does not appear to be a valid GMT or ESRI file, not applying any slope thresholds.' % (slopefile))
        else:
            print('No slope file provided, slope thresholds not applied')
        # Stuff into Grid2D object
        temp = self.shakemap.getShakeDict()
        shakedetail = '%s_ver%s' % (temp['shakemap_id'], temp['shakemap_version'])
        description = {'name': self.modelrefs['shortref'], 'longref': self.modelrefs['longref'], 'units': 'probability',
                       'shakemap': shakedetail, 'parameters': {'slopemin': self.slopemin, 'slopemax': self.slopemax}}
        Pgrid = Grid2D(P, self.geodict)
        rdict = collections.OrderedDict()
        rdict['model'] = {'grid': Pgrid,
                          'label': ('%s Probability') % (self.modeltype.capitalize()),
                          'type': 'output',
                          'description': description}
        rdict['modelmin'] = {'grid': Grid2D(Pmin, self.geodict),
                             'label': ('%s Probability (-1 std ground motion)') % (self.modeltype.capitalize()),
                             'type': 'output',
                             'description': description}
        rdict['modelmax'] = {'grid': Grid2D(Pmax, self.geodict),
                             'label': ('%s Probability (+1 std ground motion') % (self.modeltype.capitalize()),
                             'type': 'output',
                             'description': description}

        if saveinputs is True:
            for layername, layergrid in list(self.layerdict.items()):
                units = self.units[layername]
                rdict[layername] = {'grid': layergrid,
                                    'label': '%s (%s)' % (layername, units),
                                    'type': 'input',
                                    'description': {'units': units, 'shakemap': shakedetail}}
            for gmused in self.gmused:
                if 'pga' in gmused:
                    units = '%g'
                    getkey = 'pga'
                if 'pgv' in gmused:
                    units = 'cm/s'
                    getkey = 'pgv'
                if 'mmi' in gmused:
                    units = 'intensity'
                    getkey = 'mmi'
                layer = self.shakemap.getLayer(getkey)
                rdict[gmused] = {'grid': layer,
                                 'label': '%s (%s)' % (getkey.upper(), units),
                                 'type': 'input',
                                 'description': {'units': units, 'shakemap': shakedetail}}
                if self.uncert is not None:
                    layer1 = np.exp(np.log(layer.getData()) - self.uncert.getLayer('std'+getkey).getData())
                    rdict[gmused + '-1std'] = {'grid': Grid2D(layer1, self.geodict),
                                               'label': '%s (%s)' % (getkey.upper()+' - 1 std', units),
                                               'type': 'input',
                                               'description': {'units': units, 'shakemap': shakedetail}}
                    layer2 = np.exp(np.log(layer.getData()) + self.uncert.getLayer('std'+getkey).getData())
                    rdict[gmused + '+1std'] = {'grid': Grid2D(layer2, self.geodict),
                                               'label': '%s (%s)' % (getkey.upper()+' + 1 std', units),
                                               'type': 'input',
                                               'description': {'units': units, 'shakemap': shakedetail}}

        return rdict


def _test(shakefile, cofile, slopefile, precipfolder):
    model = {'logistic_models': {'nowicki_2014': {'description': 'This is the Nowicki Model of 2014, which uses cohesion and slope max as input.',
                                                  'gfetype': 'landslide',
                                                  'baselayer': 'cohesion',
                                                  'layers': {'cohesion': '%s' % cofile,
                                                             'slope': '%s' % slfile,
                                                             'precip': '%s' % precipfolder},
                                                  'interpolations': {'cohesion': 'linear',
                                                                     'slope': 'linear',
                                                                     'precip': 'nearest'},
                                                  'terms': {'b1': 'pga',
                                                            'b2': 'slope',
                                                            'b3': 'precipMONTH',
                                                            'b4': 'pga*slope*MW'},
                                                  'coefficients': {'b0': -7.15,
                                                                   'b1': 0.0604,
                                                                   'b2': 0.000825,
                                                                   'b3': 0.0201,
                                                                   'b4': 1.45e-05}}}}

    lm = LogisticModel(model, shakefile, 'nowicki_2014_global')
    print(lm.getEquation())
    P = lm.calculate()

if __name__ == '__main__':
    shakefile = sys.argv[1]  # needs to be an event occurring in January
    cofile = sys.argv[2]
    slfile = sys.argv[3]
    precip = sys.argv[4]
    _test(shakefile, cofile, slfile, precip)
