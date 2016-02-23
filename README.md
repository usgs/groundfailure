# groundfailure

Introduction
------------
groundfailure is a project designed to implement as many methods for calculating landslide and liquefaction probability 
given an input ShakeMap.  

Installation and Dependencies
-----------------------------
This package depends on a number of other libraries, most of which are installed with Anaconda.  We strongly suggest that
users install that software or one of the other Scipy Stack distributions described here:

http://www.scipy.org/install.html

Some of those that may not be installed with a distribution:

 - configobj (installed with Anaconda)
 - fiona

Usually these packages can be installed with either conda (Anaconda installs) or pip.

The final dependency, mapio, can be installed with:

pip install git+git://github.com/usgs/MapIO.git  

To install this package:

pip install git+git://github.com/usgs/groundfailure.git

To upgrade this package:

pip install -U git+git://github.com/usgs/groundfailure.git


Configuration
-------------

There will be a configuration file found in ~/.groundfailure/config.ini, with the format described below.

The config file format is a modified version of the "INI" format.  It is described in detail here:

http://configobj.readthedocs.org/en/latest/configobj.html#config-files

**Note** - The configuration below does not represent any valid landslide or liquefaction model.  The parameters
and layers are shown here for the purpose of explaining how to configure models. 

<pre>
[logistic_models]

  #default_landslide and default_liquefaction parameters below must refer to named models in this file
  default_landslide = nowicki_2014
  default_liquefaction = zhu_2014

  # this can be any string, but it must be a unique and descriptive name of a logistic regression model.
  [[nowicki_2014]]
    #Detailed description of the model, its inputs, etc.
    description = 'This is the Nowicki Model of 2014, which uses cohesion and slope max as input.'
    
    #which type of ground failure model is this? Options are landslide or liquefaction.
    gfetype = landslide

    #what is the grid to which all other grids in this model will be resampled?
    baselayer = cohesion 

    #these layer files can be any Grid2D-subclass supported format
    #These include, but may not be limited to:
    #https://github.com/usgs/MapIO/blob/master/mapio/gdal.py
    #https://github.com/usgs/MapIO/blob/master/mapio/gmt.py
    #OR
    #A layer can point to a directory containing 12 data files (from above format list), one for each month.
    #The input event time will be used to choose the appropriate file from the list of 12.
    #The files MUST be named with the capitalized three-letter abbreviation of the month name, like "precip_Jan.grd", or "slope_May.grd".
    #If some files contain more than one of these three-letter abbreviations, you will get unexpected results. (i.e., "DecimatedSlope_Jan.grd")
    [[[layers]]]
      cohesion = /Users/user/secondary/data/cohesion_10i.grd
      slope = /Users/user/secondary/data/slope_max.grd
      precip = /Users/user/secondary/data/precipdata

    #indicate what kind of interpolation should be used for each of the above layers (nearest, linear, cubic)
    [[[interpolations]]]
      cohesion = nearest
      slope = linear
      precip = nearest
      
    #What are the physical units of the various predictive layers?  These will be displayed on output plots
    #and preserved in the output data files. 
    [[[units]]]
      cohesion = kPa
      slope = degrees
      precip = cm/hr

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

  [[zhu_2014]]
  
    #Detailed description of the model, its inputs, etc.
    description = 'This is the Zhu model of 2014, where we use vs30 and CTI'

    gfetype = liquefaction
  
    baselayer = vs30

    #these layer files can be any Grid2D-subclass supported format
    #These include, but may not be limited to:
    #https://github.com/usgs/MapIO/blob/master/mapio/gdal.py
    #https://github.com/usgs/MapIO/blob/master/mapio/gmt.py
    [[[layers]]]
      vs30 = /Users/user/secondary/data/global_vs30.grd
      cti = /Users/user/secondary/data/globalcti.grd 

    #What are the physical units of the various predictive layers?  These will be displayed on output plots
    #and preserved in the output data files. 
    [[[units]]]
      vs30 = m/s
      cti = cm/m

    [[[interpolations]]]
      vs30 = nearest
      cti = linear

    [[[terms]]]
      b1 = 'log((PGA/100.0)*(power(MW,2.56)/power(10,2.24)))'
      b2 = 'cti/100.0'
      b3 = 'log(vs30)'

    [[[coefficients]]]
      b0 = 24.10
      b1 = 2.067
      b2 = 0.355
      b3 = -4.784

    [[[colormaps]]]
      vs30_colormap = jet_r


[output]
  folder = '/Users/user/failureoutput/
</pre>

API for Model Output
---------------------

Each model should output a single dictionary, which has keys that correspond to the names of the 
input and output layers from the model.

The *description* dictionary should have the following fields:

 - name Short name, suitable for use as a plot title if necessary.
 - longref Full citation, USGS format as described here: http://internal.usgs.gov/publishing/sta/sta28.pdf
 - units Physical units for input data layers, and one of the following for output "probability" layers:

         * index Relative (low to high) index of occurrence in a given cell (not necessarily bounded).
         * probability Probability of event (landslide,liquefaction) of a given size occurring in a given cell (0 to 1).
         * coverage Fractional coverage of groundfailure in a given cell (0 to 1).
         * displacement Distance material will move from or in given cell (unbounded).

The *parameters* dictionary should contain key/value pairs, where the values must be either numbers or strings.



Template Model Function:
<pre>
def failure_model():
    geodict = GeoDict({'xmin':0.5,'xmax':3.5,
                       'ymin':0.5,'ymax':3.5,
                       'dx':1.0,'dy':1.0,
                       'nx':4,'ny':4})
    pgrid = Grid2D(data = np.arange(0,16).reshape(4,4),geodict=geodict)
    cgrid = Grid2D(data = np.arange(1,17).reshape(4,4),geodict=geodict)
    sgrid = Grid2D(data = np.arange(2,18).reshape(4,4),geodict=geodict)

    problayer = {'description':{'name':'Nowicki 2014',
                                'longref':'Nowicki, A., 2014, A logistic regression landslide model: Failure Monthly, v. 2, p. 1-7.',
                                'units':'index',
                                'parameters':{'b0':1.045,
                                              'b1':5.435}},
                 'type':'output',
                 'label':'Relative Index Value',
                 'grid':pgrid,
                 }
    
    layer1 = {'description':{'name':'cohesion',
                             'longref':'Smith J. and Jones, J., 1994, Holding on to things: Journal of Geophysical Sciences, v. 17,  p. 100-105',
                             'units':'kPa'},
              'type':'input',
              'label':'cohesion (kPa)',
              'grid':cgrid}
    
    layer2 = {'description':{'name':'slope',
                             'longref':'Garfunkel, A., and Oates, J., 2001, I'm afraid to look down: Journal of Steepness, v. 8, p. 10-25',
                             'units':'degrees'},
              'type':'input',
              'label':'slope (degrees)',
              'grid':sgrid}

    output = {'probability':problayer,
              'cohesion':layer1,
              'slope':layer2}

    return output
</pre>