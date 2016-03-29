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

**Notes** 
- The configuration below does not represent any valid landslide or liquefaction model.  The parameters
and layers are shown here for the purpose of explaining how to configure models. 
- References and other inputs with commas within them need to be enclosed in quotes or else they will not be read in properly (commas will be used to separate) - for example: 'Verdin, D.W., Godt, J., Funk, C., Pedreros, D., Worstell, B. and Verdin, J., 2007, Development of a global slope dataset for estimation of landslide occurrence resulting from earthquakes: U.S. Geological Survey Open-File Report 2007–1188, 25p.'
- Arrays should be not be enclosed in brackets and should be comma separated, for example: model = 0, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.

<pre>
[output]
  folder = '/Users/user/failureoutput/

[mapdata]
  [[dem]]
    # Optional, don't need if have hillshade, just for making hillshade
    file = /Users/kallstadt/SecondaryHazards/Codes/inputs/md30_gmted_gmt.grd

  [[[hillshade]]]
    file = /Users/kallstadt/SecondaryHazards/Codes/inputs/gmted_global_hillshade.grd

  [[roads]]
    folder = /Users/kallstadt/SecondaryHazards/Codes/inputs/roads
    longref = 'Center for International Earth Science Information Network - CIESIN, 2013, Global Roads Open Access Data Set, Version 1 (gROADSv1): Columbia University, and Information Technology Outreach Services - ITOS - University of Georgia, Palisades, NY, NASA Socioeconomic Data and Applications Center (SEDAC). http://dx.doi.org/10.7927/H4VD6WCT.'
    shortref = 'CIESIN (2013)'

  [[cities]]
    file = /Users/kallstadt/SecondaryHazards/Codes/inputs/cities1000.txt
    longref = GeoNames, http://geonames.org/ Accessed: 2 Sept 2015
    shortref = GeoNames

  [[lims]]
    # Corresponding to different possible layer keys - don't need these, will just use defaults if missing, don't need full name of layer, just something that is part of it
    model = 0, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.
    pga = None
    pgv = None
    FS = 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0
    slope = None
    cohesion = None
    friction = None

  [[colors]]
    # Define basic colors and transparencies
    roadcolor = 808080
    countrycolor = 474747
    watercolor = B8EEFF
    alpha = 0.7
    default = cm.jet
    # Corresponding to different possible layer keys - don't need these, will just use defaults if missing
    model = cm.jet
    pga = cm.jet
    pgv = cm.jet
    FS = cm.jet
    slope = cm.gnuplot2
    cohesion = cm.jet
    friction = cm.jet

  [[logscale]]
    # Corresponding to different possible layer keys - don't need these, will just use defaults if missing, don't need full name of layer, just something that is part of it
    model = False
    pga = True
    pgv = True
    FS = False 
    slope = False
    cohesion = False
    friction = False

  [[maskthresholds]]
    # Corresponding to different possible layer keys - don't need these, will just use defaults if missing, don't need full name of layer, just something that is part of it
    model = 0.
    pga = None
    pgv = None
    FS = None 
    slope = None
    cohesion = None
    friction = None

[mechanistic_models]
  [[godt_2008]]
    #Detailed description of the model, its inputs, etc.
    longref = 'Godt, J.W., Sener, B., Verdin, K.L., Wald, D.J., Earle, P.S., Harp, E.L. and Jibson, R.W., 2008, Rapid Assessment of Earthquake-induced Landsliding: Proceedings of the First World Landslide Forum, United Nations University, Tokyo, Japan, p. 392-395.'
    shortref = 'Godt and others (2008)'

    #which type of ground failure model is this? Options are landslide or liquefaction.
    gfetype = landslide
    
    [[[layers]]]
      [[[[cohesion]]]]
        file = /Users/kallstadt/SecondaryHazards/Datasets/Godt_inputs/cohesion_10i.flt
        units = kPa
        longref = 'Godt, J.W., Sener, B., Verdin, K.L., Wald, D.J., Earle, P.S., Harp, E.L. and Jibson, R.W., 2008, Rapid Assessment of Earthquake-induced Landsliding: Proceedings of the First World Landslide Forum, United Nations University, Tokyo, Japan, p. 392-395.'
        shortref = 'Godt and others (2008)'

      [[[[friction]]]]
        file = /Users/kallstadt/SecondaryHazards/Datasets/Godt_inputs/friction.flt
        units = degrees
        longref = 'Godt, J.W., Sener, B., Verdin, K.L., Wald, D.J., Earle, P.S., Harp, E.L. and Jibson, R.W., 2008, Rapid Assessment of Earthquake-induced Landsliding: Proceedings of the First World Landslide Forum, United Nations University, Tokyo, Japan, p. 392-395.'
        shortref = 'Godt and others (2008)'

      [[[[slope]]]]
        filepath = /Users/kallstadt/SecondaryHazards/Datasets/Verdin_slopes_resampled
        units = degrees*100
        longref = 'Verdin, D.W., Godt, J., Funk, C., Pedreros, D., Worstell, B. and Verdin, J., 2007, Development of a global slope dataset for estimation of landslide occurrence resulting from earthquakes: U.S. Geological Survey Open-File Report 2007–1188, 25p.'
        shortref = 'Verdin et al. (2007)'

    [[[parameters]]]
      #Slope thickness in meters
      thick = 2.4
      #Soil unit weight
      uwt = 15.7 
      #Cohesion value for no_data grid cells
      nodata_cohesion = 1.0
      #Friction angle value for no_data grid cells
      nodata_friction = 26.
      #Newmark displacement threshold in cm
      dnthresh = 5.
      #Minimum Factor of safety allowed (unitless)
      fsthresh = 1.01
      #Minimum critical acceleration allowed (in g's)
      acthresh = 0.05

  [[classic_newmark]]
    longref = 'Jibson, R. W., Harp, E. L. and Michael, J. A., 2000, A method for producing digital probabilistic seismic landslide hazard maps: Engineering Geology, v. 58, p. 271-289.'
    shortref = 'Jibson and others (2000)'

    [[[layers]]]
      [[[[cohesion]]]]
        file = /Users/kallstadt/SecondaryHazards/Datasets/Godt_inputs/cohesion_10i.flt
        interpolation = nearest
        units = kPa
        longref = 'Godt, J.W., Sener, B., Verdin, K.L., Wald, D.J., Earle, P.S., Harp, E.L. and Jibson, R.W., 2008, Rapid Assessment of Earthquake-induced Landsliding: Proceedings of the First World Landslide Forum, United Nations University, Tokyo, Japan, p. 392-395.'
        shortref = 'Godt and others (2008)'

      [[[[friction]]]]
        file = /Users/kallstadt/SecondaryHazards/Datasets/Godt_inputs/friction.flt
        units = degrees
        longref = 'Godt, J.W., Sener, B., Verdin, K.L., Wald, D.J., Earle, P.S., Harp, E.L. and Jibson, R.W., 2008, Rapid Assessment of Earthquake-induced Landsliding: Proceedings of the First World Landslide Forum, United Nations University, Tokyo, Japan, p. 392-395.'
        shortref = 'Godt and others (2008)'

      [[[[slope]]]]
        file = /Users/kallstadt/SecondaryHazards/Datasets/Verdin_slopes_resampled/slope_max.bil
        units = degrees
        longref = 'Verdin, D.W., Godt, J., Funk, C., Pedreros, D., Worstell, B. and Verdin, J., 2007, Development of a global slope dataset for estimation of landslide occurrence resulting from earthquakes: U.S. Geological Survey Open-File Report 2007–1188, 25p.'
        shortref = 'Verdin and others (2007)'

      [[[[watertable]]]]
        file = /Users/kallstadt/SecondaryHazards/Datasets/Fan2013WaterTable/wtd_fan2013_zhu_fil_na.grd
        units = meters
        longref = 'Fan, Y., Li, H., and Miguez-Macho, G., 2013, Global Patterns of Groundwater Table Depth: Science, 339, 940-943.'
        shortref = 'Fan and others (2013)'

     [[[parameters]]]
        #Slope thickness in meters
        thick = 2.4
        #Soil unit weight (dry)
        uwt = 15.7
        # Soil unit weight (wet)
        uwtw = 18.8
        #Cohesion value for no_data grid cells
        nodata_cohesion = 5.0
        #Friction angle value for no_data grid cells
        nodata_friction = 26.
        #Newmark displacement threshold in cm - only need if probtype is threshold
        dnthresh = 5.
        #Minimum Factor of safety allowed (unitless)
        fsthresh = 1.01
        #Minimum critical acceleration allowed (in g's)
        acthresh = 0.05
        # Slope minimum threshold in degrees
        slopethresh = 5.
        # Constant proportion of saturated thickness [0-1] (used if watertable file not specified)
        m = 1.0

  [[hazus]]
    longref = 'Federal Emergency Management Agency, 2013, Hazus - MH MR5 Multi-hazard loss Estimation Methodology Earthquake Model: Dept. of Homeland Security, Federal Emergency Management Agency, Washington D.C., 736p. [available at: www.fema.gov/plan/prevent/hazus]'
    shortref = 'FEMA (2013)'

    [[[layers]]]
      [[[[susceptibility]]]]
        file = /Users/kallstadt/SecondaryHazards/Datasets/Wills_et_al_2011/Resampled/Susceptiblity_dry_WGS84.bil
        units = N/A
        longref = 'Wills, C.J., Perez, F.G., and Gutierrez, C.I., 2011, Susceptibility to Deep-Seated Landslides in California: California Geological Survey Map Sheet 58, 1p.'
        shortref = 'Wills and others (2011)'

    [[[parameters]]]
      #Newmark displacement threshold for failure in cm - only need for some models
      dnthresh = 5

[logistic_models]
  # NEEDS TO BE UPDATED FOR NEW STRUCTURE

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
</pre>

API for Model Output
---------------------

Each model should output a single dictionary, which has keys that correspond to the names of the 
input and output layers from the model.

Each layer in the dictionary is itself a dictionary, with the following fields:
 - *description* A dictionary with the fields:

   * *name* Short name, suitable for use as a plot title if necessary.
   * *longref* Full citation, USGS format as described here: http://internal.usgs.gov/publishing/sta/sta28.pdf
   * *units* Physical units for input data layers, and one of the following for output "probability" layers:

     * *index* Relative (low to high) index of occurrence in a given cell (not necessarily bounded).
     * *probability* Probability of event (landslide,liquefaction) of a given size occurring in a given cell (0 to 1).
     * *coverage* Fractional coverage of groundfailure in a given cell (0 to 1).
     * *displacement* Distance material will move from or in given cell (unbounded).

   * *parameters* (Not required for input layers) A dictionary of key/value pairs, where the values must be either numbers or strings.

 - *type* Indicates whether this grid contains input data or output from a model.

 - *label* What will be written next to the colorbar for the data layer.

 - *grid* Input data or model output, in the form of a Grid2D object. 

A template model function implementation is shown below.
<pre>
def failure_model():
    geodict = GeoDict({'xmin':0.5,'xmax':3.5,
                       'ymin':0.5,'ymax':3.5,
                       'dx':1.0,'dy':1.0,
                       'nx':4,'ny':4})
    pgrid = Grid2D(data = np.arange(0,16).reshape(4,4),geodict=geodict)
    cgrid = Grid2D(data = np.arange(1,17).reshape(4,4),geodict=geodict)
    sgrid = Grid2D(data = np.arange(2,18).reshape(4,4),geodict=geodict)
    mgrid = Grid2D(data = np.arange(3,19).reshape(4,4),geodict=geodict)

    modellayer = {'description':{'name':'Nowicki 2014',
                                'longref':'Nowicki, A., 2014, A logistic regression landslide model: Failure Monthly, v. 2, p. 1-7.',
                                'units':'index',
                                'shakemap': '19940117123055_ver2'
                                'parameters':{'b0':1.045,
                                              'b1':5.435}},
                 'type':'output',
                 'label':'Relative Index Value',
                 'grid':pgrid,
                 }
    
    layer1 = {'description':{'name':'Smith and Jones 1994',
                             'longref':'Smith J. and Jones, J., 1994, Holding on to things: Journal of Geophysical Sciences, v. 17,  p. 100-105',
                             'units':'kPa'},
              'type':'input',
              'label':'cohesion (kPa)',
              'grid':cgrid}
    
    layer2 = {'description':{'name':'Garfunkel and Oates 2001',
                             'longref':'Garfunkel, A., and Oates, J., 2001, I'm afraid to look down: Journal of Steepness, v. 8, p. 10-25',
                             'units':'degrees'},
              'type':'input',
              'label':'slope (degrees)',
              'grid':sgrid}

    layer3 = {'description':{'units':'g'
                             'shakemap': '19940117123055_ver2'},
              'type':'input',
              'label':'PGA (g)',
              'grid':mgrid}

    output = {'model':problayer,
              'cohesion':layer1,
              'slope':layer2,
              'pga':layer3}

    return output
</pre>