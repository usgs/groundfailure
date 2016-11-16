#!/usr/bin/env python
"""
This module contains functions and class definitions for running forward models of models based on logistic regression.
"""

#stdlib imports
import numpy as np
import os.path
import re
import collections
import copy

#third party imports
from mapio.shake import ShakeGrid
from mapio.shake import getHeaderData
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
    """Get the names of the models present in the configobj

    :param config: configobj (config .ini file read in using configobj) defining the model and its inputs
    :type config: dictionary
    :returns:
        names: list of model names

    """
    names = []
    lmodel_space = config
    for key, value in lmodel_space.items():
        if isinstance(value, str):
            continue
        else:  # this is a model
            names.append(key)
    return names


def getFileType(filename):
    """Determine whether input file is a shapefile or a grid (ESRI or GMT).
    EVENTUALLY WILL BE MOVED TO MAPIO

    :param filename:
      String path to candidate filename.
    :returns:
      String, one of 'shapefile','grid','unknown'.

    """
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
    """Get list of all gmt or esri (.grd, .bil) files in a directory
    EVENTUALLY WILL BE MOVED TO MAPIO

    :param indir: directory to search
    :type indir: string
    :returns:
        flist: list of file names

    """
    tflist = os.listdir(indir)
    flist = []
    for tf in tflist:
        fullfile = os.path.join(indir, tf)
        ftype = getFileType(fullfile)
        if ftype in ['gmt', 'esri']:
            flist.append(fullfile)
    return flist


def validateCoefficients(cmodel):
    """Ensures coefficients provided in model description are valid and outputs a dictionary of the coefficients.

    :param cmodel: subdictionary from config for specific model, e.g. cmodel = config['test_model']
    :type cmodel: dictionary
    :returns:
        coeffs(dictionary): a dictionary of model coefficients named b0, b1, b2...

    """
    coeffs = {}
    for key, value in cmodel['coefficients'].items():
        if re.search('b[0-9]*', key) is None:
            raise Exception('coefficients must be named b0, b1, ...')
        coeffs[key] = float(value)
    if 'b0' not in list(coeffs.keys()):
        raise Exception('coefficients must include an intercept coefficient named b0.')
    return coeffs


def validateLayers(cmodel):
    """Ensures all input files required to run the model exist and are valid file types

    :param cmodel: subdictionary from config for specific model, e.g. cmodel = config['test_model']
    :type cmodel: dictionary
    :returns:
        layers(dictionary): a dictionary of file names (e.g. {'slope': 'slopefile.bil', 'vs30': 'vs30.grd'})

    """
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
    """Reformats model inputs from config file, replacing functions with numpy functions, inserting code for extracting
    data from each layer (required to run eval in the calculate step), addressing any time variables, and checks that
    term names match coefficient names.

    TODO - return a time field for every term, not just one global one.

    :param cmodel: subdictionary from config for specific model, e.g. cmodel = config['test_model']
    :type cmodel: dictionary
    :param coeffs: a dictionary of model coefficients (e.g. {'b0': 3.5, 'b1': -0.01})
    :type coeffs: dictionary
    :param layers: a dictionary of file names for all input layers (e.g. {'slope': 'slopefile.bil', 'vs30': 'vs30.grd'})
    :type layers: dictionary
    :returns:
        tuple (terms, timeField)
        terms: dictionary of terms that form the model equation (e.g. 'b1': "self.layerdict['friction'].getData()",
            'b2': "self.layerdict['slope'].getData()/100.")
        timeField: Field that indicates time that is used to know which input file to read in (e.g. for monthly average
            precipitation, 'MONTH')

    """
    terms = {}
    timeField = None
    for key, value in cmodel['terms'].items():
        if key not in list(coeffs.keys()):
            raise Exception('Term names must match names of coefficients')
        # replace log with np.log, make sure variables are all in layers list, etc.
        term, rem, tTimeField = checkTerm(value, layers)
        #print(term)
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
    if len(mnames) > 1:
        raise Exception('Config file contains more than one model which is no longer allowed,\
                        update your config file to the newer format')
    for cmodelname in mnames:
        try:
            cmodel = config[cmodelname]
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
    def __init__(self, shakefile, config, uncertfile=None, saveinputs=False, slopefile=None, slopediv=1.,
                 bounds=None):
        """Set up the logistic model
        # ADD BOUNDS TO THIS MODEL
        :param config: configobj (config .ini file read in using configobj) defining the model and its inputs. Only one
          model should be described in each config file.
        :type config: dictionary
        :param shakefile: Full file path to shakemap.xml file for the event of interest
        :type shakefile: string
        :param uncertfile: Full file path to xml file of shakemap uncertainties
        :type uncertfile: string
        :param saveinputs: if True, saves all the input layers as Grid2D objects in addition to the model
          if false, it will just output the model
        :type saveinputs: boolean
        :param slopefile: optional file path to slopefile that will be resampled to the other input files for applying
          thresholds OVERWRITES VALUE IN CONFIG
        :type slopefile: string
        :param slopediv: number to divide slope by to get to degrees (usually will be default
          of 1.)
        :type slopediv: float

        """
        mnames = getLogisticModelNames(config)
        if len(mnames) > 1:
            raise Exception('Config file contains more than one model which is no longer allowed,\
                            update your config file to the newer format')
        self.model = mnames[0]
        cmodel = config[self.model]
        self.modeltype = cmodel['gfetype']
        self.coeffs = validateCoefficients(cmodel)
        self.layers = validateLayers(cmodel)  # key = layer name, value = file name
        self.terms, timeField = validateTerms(cmodel, self.coeffs, self.layers)
        self.interpolations = validateInterpolations(cmodel, self.layers)
        self.units = validateUnits(cmodel, self.layers)
        self.gmused = [value for term, value in cmodel['terms'].items() if 'pga' in value.lower() or 'pgv' in
                       value.lower() or 'mmi' in value.lower()]
        self.modelrefs, self.longrefs, self.shortrefs = validateRefs(cmodel)
        if cmodel['baselayer'] not in list(self.layers.keys()):
            raise Exception('You must specify a base layer corresponding to one of the files in the layer section.')
        self.saveinputs = saveinputs
        if slopefile is None:
            try:
                self.slopefile = cmodel['slopefile']
            except:
                print('Could not find slopefile term in config, no slope thresholds will be applied\n')
                self.slopefile = None
        else:
            self.slopefile = slopefile
        self.slopediv = slopediv

        #get the geodict for the shakemap
        geodict = ShakeGrid.getFileGeoDict(shakefile, adjust='res')
        griddict, eventdict, specdict, fields, uncertainties = getHeaderData(shakefile)
        #YEAR = eventdict['event_timestamp'].year
        MONTH = MONTHS[(eventdict['event_timestamp'].month)-1]
        #DAY = eventdict['event_timestamp'].day
        #HOUR = eventdict['event_timestamp'].hour

        #now find the layer that is our base layer and get the largest bounds we can guarantee not to exceed shakemap bounds
        basefile = self.layers[cmodel['baselayer']]
        ftype = getFileType(basefile)
        if ftype == 'esri':
            basegeodict, firstcol = GDALGrid.getFileGeoDict(basefile)
            sampledict = basegeodict.getBoundsWithin(geodict)
        elif ftype == 'gmt':
            basegeodict, firstcol = GMTGrid.getFileGeoDict(basefile)
            sampledict = basegeodict.getBoundsWithin(geodict)
        else:
            raise Exception('All predictor variable grids must be a valid GMT or ESRI file type')

        #now load the shakemap, resampling and padding if necessary
        if ShakeGrid.getFileGeoDict(shakefile, adjust='res') == sampledict:
            self.shakemap = ShakeGrid.load(shakefile, adjust='res')
            flag = 1
        else:
            self.shakemap = ShakeGrid.load(shakefile, samplegeodict=sampledict, resample=True, doPadding=True,
                                           adjust='res')
            flag = 0

        # take uncertainties into account
        if uncertfile is not None:
            try:
                if flag == 1:
                    self.uncert = ShakeGrid.load(uncertfile, adjust='res')
                else:
                    self.uncert = ShakeGrid.load(uncertfile, samplegeodict=sampledict, resample=True, doPadding=True,
                                                 adjust='res')
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
                                if GMTGrid.getFileGeoDict(layerfile)[0] == sampledict:
                                    lyr = GMTGrid.load(layerfile)
                                else:
                                    lyr = GMTGrid.load(layerfile, sampledict, resample=True, method=interp,
                                                       doPadding=True)
                            elif ftype == 'esri':
                                if GDALGrid.getFileGeoDict(layerfile)[0] == sampledict:
                                    lyr = GDALGrid.load(layerfile)
                                else:
                                    lyr = GDALGrid.load(layerfile, sampledict, resample=True, method=interp,
                                                        doPadding=True)
                            else:
                                msg = 'Layer %s (file %s) does not appear to be a valid GMT or ESRI file.' % (layername,
                                                                                                              layerfile)
                                raise Exception(msg)
                            self.layerdict[layername] = lyr
            else:
                #first, figure out what kind of file we have (or is it a directory?)
                ftype = getFileType(layerfile)
                interp = self.interpolations[layername]
                if ftype == 'gmt':
                    if GMTGrid.getFileGeoDict(layerfile)[0] == sampledict:
                        lyr = GMTGrid.load(layerfile)
                    else:
                        lyr = GMTGrid.load(layerfile, sampledict, resample=True, method=interp,
                                           doPadding=True)
                elif ftype == 'esri':
                    if GDALGrid.getFileGeoDict(layerfile)[0] == sampledict:
                        lyr = GDALGrid.load(layerfile)
                    else:
                        lyr = GDALGrid.load(layerfile, sampledict, resample=True, method=interp,
                                            doPadding=True)
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
                    self.nugmin[k] = self.nugmin[k].replace("self.shakemap.getLayer('pga').getData()",
                                                            "(np.exp(np.log(self.shakemap.getLayer('pga').getData())\
                                                             - self.uncert.getLayer('stdpga').getData()))")
                    self.nugmax[k] = self.nugmax[k].replace("self.shakemap.getLayer('pga').getData()",
                                                            "(np.exp(np.log(self.shakemap.getLayer('pga').getData())\
                                                             + self.uncert.getLayer('stdpga').getData()))")
                elif "self.layerdict['pgv'].getData()" in nug:
                    self.nugmin[k] = self.nugmin[k].replace("self.shakemap.getLayer('pgv').getData()",
                                                            "(np.exp(np.log(self.shakemap.getLayer('pgv').getData())\
                                                             - self.uncert.getLayer('stdpgv').getData()))")
                    self.nugmax[k] = self.nugmax[k].replace("self.shakemap.getLayer('pgv').getData()",
                                                            "(np.exp(np.log(self.shakemap.getLayer('pgv').getData())\
                                                             + self.uncert.getLayer('stdpgv').getData()))")
                elif "self.layerdict['mmi'].getData()" in nug:
                    self.nugmin[k] = self.nugmin[k].replace("self.shakemap.getLayer('mmi').getData()",
                                                            "(np.exp(np.log(self.shakemap.getLayer('mmi').getData())\
                                                             - self.uncert.getLayer('stdmmi').getData()))")
                    self.nugmax[k] = self.nugmax[k].replace("self.shakemap.getLayer('mmi').getData()",
                                                            "(np.exp(np.log(self.shakemap.getLayer('mmi').getData())\
                                                             + self.uncert.getLayer('stdmmi').getData()))")
            self.equationmin = ' + '.join(self.nugmin)
            self.equationmax = ' + '.join(self.nugmax)
        else:
            self.equationmin = None
            self.equationmax = None

        self.geodict = self.shakemap.getGeoDict()

        try:
            self.slopemin = float(config[self.model]['slopemin'])
            self.slopemax = float(config[self.model]['slopemax'])
        except:
            print('could not find slopemin and/or slopemax in config, no limits will be applied')
            self.slopemin = 0.
            self.slopemax = 90.

    def getEquation(self):
        """Method for LogisticModel class to extract a string defining the equation for the model which can be run
        using eval()

        :returns:
            equation for model using median ground motions

        """
        return self.equation

    def getEquations(self):
        """Method for LogisticModel class to extract strings defining the equations for the model for median
        ground motions and +/- one standard deviation (3 total)

        :returns:
            tuple of three equations (equation, equationmin, equationmax) where equation is the equation for
            median ground motions, equationmin is the equation for the same model but with median ground motions
            minus 1 standard deviation and equationmax is the same but for plus 1 standard deviation.

        """
        return self.equation, self.equationmin, self.equationmax

    def getGeoDict(self):
        """Returns the geodictionary of the LogisticModel class defining bounds and resolution of model inputs and outputs

        :returns:
            geodict: mapio geodict object

        """
        return self.geodict

    def calculate(self):
        """Calculate the model

        :returns:
            a dictionary containing the model results and model inputs if saveinputs was set to
            True when class was set up, see <https://github.com/usgs/groundfailure#api-for-model-output> for a
            description of the structure of this output

        """
        X = eval(self.equation)
        P = 1/(1 + np.exp(-X))
        if self.uncert is not None:
            Xmin = eval(self.equationmin)
            Xmax = eval(self.equationmax)
            Pmin = 1/(1 + np.exp(-Xmin))
            Pmax = 1/(1 + np.exp(-Xmax))
        if self.slopefile is not None:
            ftype = getFileType(self.slopefile)
            sampledict = self.shakemap.getGeoDict()
            if ftype == 'gmt':
                if GMTGrid.getFileGeoDict(self.slopefile)[0] == sampledict:
                        slope = GMTGrid.load(self.slopefile).getData()/self.slopediv
                else:
                        slope = GMTGrid.load(self.slopefile, sampledict, resample=True, method='linear',
                                             doPadding=True).getData()/self.slopediv
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
                if GDALGrid.getFileGeoDict(self.slopefile)[0] == sampledict:
                        slope = GDALGrid.load(self.slopefile).getData()/self.slopediv
                else:
                        slope = GDALGrid.load(self.slopefile, sampledict, resample=True, method='linear',
                                              doPadding=True).getData()/self.slopediv
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
                print('Slope file does not appear to be a valid GMT or ESRI file, not applying any slope thresholds.'
                      % (self.slopefile))
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
        if self.uncert is not None:
            rdict['modelmin'] = {'grid': Grid2D(Pmin, self.geodict),
                                 'label': ('%s Probability (-1 std ground motion)') % (self.modeltype.capitalize()),
                                 'type': 'output',
                                 'description': description}
            rdict['modelmax'] = {'grid': Grid2D(Pmax, self.geodict),
                                 'label': ('%s Probability (+1 std ground motion)') % (self.modeltype.capitalize()),
                                 'type': 'output',
                                 'description': description}

        if self.saveinputs is True:
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
                                               'label': '%s (%s)' % (getkey.upper()+' -1 std', units),
                                               'type': 'input',
                                               'description': {'units': units, 'shakemap': shakedetail}}
                    layer2 = np.exp(np.log(layer.getData()) + self.uncert.getLayer('std'+getkey).getData())
                    rdict[gmused + '+1std'] = {'grid': Grid2D(layer2, self.geodict),
                                               'label': '%s (%s)' % (getkey.upper()+' +1 std', units),
                                               'type': 'input',
                                               'description': {'units': units, 'shakemap': shakedetail}}

        return rdict
