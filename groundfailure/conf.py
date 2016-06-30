#!/usr/bin/env python

#stdlib imports
import os.path
import tempfile
import textwrap
import sys
import re
import io
import shutil

#third party libraries
from configobj import ConfigObj
from validate import Validator,VdtTypeError,VdtParamError

def __getCustomValidator():
    '''Return a Validator object with the custom types we have defined here.
    :returns:
      Validator object with custom types embedded.
    '''
    fdict = {
        'file_type': __file_type,
        'path_type': __path_type,
        }
    
    validator = Validator(fdict)
    return validator

def __file_type(value):
    '''Describes a file_type from the ShakeMap config spec.
    A file_type object is simply a string that must be a valid file on the system.
    :param value:
      String representing valid path to a file on the local system.
    :return:
      Input string, if a valid file name.
    :raises VdtTypeError:
      When path is invalid.
    '''
    if not os.path.isfile(value):
        raise VdtTypeError(value)
    return value

def __path_type(value):
    '''Describes a path_type from the groundfailure config spec.
    A path_type object is simply a string that must be a valid file OR directory on the system.
    :param value:
      String representing valid path to a file or directory on the local system.
    :return:
      Input string, if a valid file/directory name.
    :raises VdtTypeError:
      When path is invalid.
    '''
    if not os.path.isfile(value) and not os.path.isdir(value):
        raise VdtTypeError(value)
    return value

def filterResults(result):
    #TODO: this function has a problem where some error messages are duplicated...?
    errormsg = ''
    for key,value in result.items():
        if isinstance(value,dict):
            tmpmsg = filterResults(value)
            errormsg += tmpmsg
        else:
            if not isinstance(value,bool):
                errormsg += "Parameter %s failed with error '%s'\n" % (key,value.message)
                
    return errormsg

def validate(configfile):
    '''Return a validated config object.
    :param configspec:
      Path to config spec file, used to define the valid configuration parameters for the system.
    :param configfile:
      Config file to validate.
    :returns:
      A validated ConfigObj object or a dictionary of which section/parameters failed validation.
    '''
    thispath = os.path.dirname(os.path.abspath(__file__))
    configspec = os.path.join(thispath,'configspec.ini')
    config = ConfigObj(configfile,configspec=configspec)
    validator = __getCustomValidator()
    result = config.validate(validator,preserve_errors=True)
    if result == True:
        return config
    else:
        errormsg = filterResults(result)
        raise TypeError(errormsg)

def _test_validate():
    configtext = '''[logistic_models]

  #default_landslide and default_liquefaction parameters below must refer to named models in this file
  default_landslide = nowicki_2014
  default_liquefaction = zhu_2015

  # this can be any string, but it must be a unique and descriptive name of a logistic regression model.
  [[nowicki_2014]]
    #Detailed description of the model, its inputs, etc.
    description = 'This is the Nowicki Model of 2014, which uses vs30 and CTI as input.'
    
    #which type of ground failure model is this? Options are landslide or liquefaction.
    gfetype = 'landslide' 

    #what is the grid to which all other grids in this model will be resampled?
    baselayer = vs30

    #these layer files can be any Grid2D-subclass supported format
    #These include, but may not be limited to:
    #https://github.com/usgs/MapIO/blob/master/mapio/gdal.py
    #https://github.com/usgs/MapIO/blob/master/mapio/gmt.py
    [[[layers]]]
      cohesion = %s
      slope = %s

    [[[terms]]]
      #These terms must be named as b1-bN, where N is the number of coefficients
      #in a logistic regression, which takes the form:
      #1/(1 + e^-eqn)
      #where eqn is a linear equation of the form:
      #b0 + b1*t1 + b2*t2 + ... + bN*tN
      #where t1, t2, ... tN are the right hand side of the parameters below.
      b1 = PGA
      b2 = slope
      b3 = cohesion/10.0
      b4 = PGA*slope

    [[[coefficients]]]
      #These coefficients must be named as b1-bN, where N is the number of coefficients
      #in a logistic regression, which takes the form:
      #1/(1 + e^-eqn)
      #where eqn is a linear equation of the form:
      #b0 + b1*t1 + b2*t2 + ... + bN*tN
      #where t1, t2, ... tN are the right hand side of the parameters below.
      b0 = -7.15
      b1 = 0.0604
      b2 = 0.000825
      b3 = 0.0201
      b4 = 1.45e-05

    #colormaps are optionally specified for each input layer
    #the groundfailure code will make a map of each input layer, and will use the
    #default Matplotlib color map unless otherwise specified here.      
    [[[colormaps]]]
      slope = jet_r
    '''
    configfile = None
    datadir = None
    try:
        #create a temporary data directory and put two files - cohesion.grd and slope.grd in it
        datadir = tempfile.mkdtemp()
        cofile = os.path.join(datadir,'cohesion.grd')
        f = open(cofile,'wt')
        f.write('This is a test\n')
        f.close()
        slfile = os.path.join(datadir,'slope.grd')
        f = open(slfile,'wt')
        f.write('This is a test\n')
        f.close()

        configtext = configtext % (cofile,slfile)
        
        #write the sample config file
        tmp,configfile = tempfile.mkstemp()
        os.close(tmp)
        f = open(configfile,'wt')
        f.write(textwrap.dedent(configtext))
        f.close()
        
        config = validate(configfile)
        print('Passed validation of sample config file.')
    except:
        print('Failed to validate sample config file.')
    finally:
        if os.path.isfile(configfile):
            os.remove(configfile)
        if os.path.isdir(datadir):
            shutil.rmtree(datadir)

def _failValidate():
    configtext = '''[logistic_models]

  #default_landslide and default_liquefaction parameters below must refer to named models in this file
  default_landslide = nowicki_2014
  default_liquefaction = zhu_2015

  # this can be any string, but it must be a unique and descriptive name of a logistic regression model.
  [[nowicki_2014]]
    #Detailed description of the model, its inputs, etc.
    description = 'This is the Nowicki Model of 2014, which uses vs30 and CTI as input.'
    
    #which type of ground failure model is this? Options are landslide or liquefaction.
    gfetype = 'landslide' 

    #what is the grid to which all other grids in this model will be resampled?
    baselayer = vs30 #deliberate error

    #these layer files can be any Grid2D-subclass supported format
    #These include, but may not be limited to:
    #https://github.com/usgs/MapIO/blob/master/mapio/gdal.py
    #https://github.com/usgs/MapIO/blob/master/mapio/gmt.py
    [[[layers]]]
      cohesion = %s
      slope = %s

    [[[terms]]]
      #These terms must be named as b1-bN, where N is the number of coefficients
      #in a logistic regression, which takes the form:
      #1/(1 + e^-eqn)
      #where eqn is a linear equation of the form:
      #b0 + b1*t1 + b2*t2 + ... + bN*tN
      #where t1, t2, ... tN are the right hand side of the parameters below.
      b1 = PGA
      b2 = slope
      b3 = cohesion/10.0
      b4 = PGA*slope

    [[[coefficients]]]
      #These coefficients must be named as b1-bN, where N is the number of coefficients
      #in a logistic regression, which takes the form:
      #1/(1 + e^-eqn)
      #where eqn is a linear equation of the form:
      #b0 + b1*t1 + b2*t2 + ... + bN*tN
      #where t1, t2, ... tN are the right hand side of the parameters below.
      b0 = -7.15
      b1 = fred #deliberate failure
      b2 = 0.000825
      b3 = 0.0201
      b4 = 1.45e-05

    #colormaps are optionally specified for each input layer
    #the groundfailure code will make a map of each input layer, and will use the
    #default Matplotlib color map unless otherwise specified here.      
    [[[colormaps]]]
      slope = jet_r
    ''' % ('foo','bar')
    #the above config text should fail validate in three places: b1, cohesion, and slope.
    configfile = None
    try:
        #write the sample config file
        tmp,configfile = tempfile.mkstemp()
        os.close(tmp)
        f = open(configfile,'wt')
        f.write(textwrap.dedent(configtext))
        f.close()
        try:
            config = validate(configfile)
            print('Validation failed to fail :)')
        except Exception as e:
            print('Validation failed as expected with sample config file:')
            print(str(e))
        
    except:
        print('Expected fail failed to fail.')
    finally:
        if os.path.isfile(configfile):
            os.remove(configfile)

            
if __name__ == '__main__':
    _test_validate()
    _failValidate()
    sys.exit(0)
    
