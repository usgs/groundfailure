# groundfailure

Introduction
------------
groundfailure is a project designed to implement as many methods for calculating landslide and liquefaction probability 
given an input ShakeMap.  

Code documentation: http://usgs.github.io/groundfailure/

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

For gdalfuncs to work (find proj.4 libraries), you may need to add the following line to your .bash_profile:
export PROJSO=/Library/Frameworks/PROJ.framework/PROJ

Configuration
-------------

Configuration options, such as the output directory, paths to input data files, model coefficients, and map display options
are all set in a configuration file found in ~/.groundfailure/config.ini.
Since many of these options depend on local paths, you need to create and modify this file youself.
The format is described below along with a template. 

The config file format is a modified version of the "INI" format.  It is described in detail here:

http://configobj.readthedocs.org/en/latest/configobj.html#config-files

**Notes** 
- The configuration below does not represent any valid landslide or liquefaction model.  The parameters
and layers are shown here for the purpose of explaining how to configure models. 
- References and other inputs with commas within them need to be enclosed in quotes or else they will not be read in properly (commas will be used to separate) - for example: 'Verdin, D.W., Godt, J., Funk, C., Pedreros, D., Worstell, B. and Verdin, J., 2007, Development of a global slope dataset for estimation of landslide occurrence resulting from earthquakes: U.S. Geological Survey Open-File Report 2007–1188, 25p.'
- Arrays should be not be enclosed in brackets and should be comma separated, for example: model = 0, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.
- Files, filepaths, and folders being used for input layers should all be labeled as 'file' in the config file regardless of actual structure.
- Files for different input and data layers should only include the object's location within the input folder indicated at the top of the file.

<pre>
[output]
  folder = /Users/user/groundfailure/outputs/

[input]
  folder = /Users/user/groundfailure/inputs/

[mapdata]
  [[dem]]
    # Optional, don't need if have hillshade, just for making hillshade
    file = md30_gmted_gmt.grd

  [[roads]]
    file = roads
    longref = 'Center for International Earth Science Information Network - CIESIN, 2013, Global Roads Open Access Data Set, Version 1 (gROADSv1): Columbia University, and Information Technology Outreach Services - ITOS - University of Georgia, Palisades, NY, NASA Socioeconomic Data and Applications Center (SEDAC). http://dx.doi.org/10.7927/H4VD6WCT.'
    shortref = 'CIESIN (2013)'

  [[cities]]
    file = cities1000.txt
    longref = GeoNames, http://geonames.org/ Accessed: 2 Sept 2015
    shortref = GeoNames

  [[oceans]]
    file = Oceans_Natural_Earth/ne_10m_ocean/ne_10m_ocean.shp
    longref = Natural Earth (2016) Ocean polygon http://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-ocean/
    shortref = Natural Earth (acc. 2016)

  [[lims]]
    # Corresponding to different possible layer keys - don't need these, will just use defaults if missing, don't need full name of layer, just something that is part of it
    model = 0, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.
    pga = None
    pgv = None
    FS = 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0
    slope = None
    cohesion = 'np.linspace(0., 40., 11.)'
    friction = 'np.linspace(0., 50., 11.)'
    suscat = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
    Dn = 'np.linspace(0., 10., 6.)'

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
    friction = 
    Dn = 0.

[mechanistic_models]
  [[godt_2008]]
    #Detailed description of the model, its inputs, etc.
    longref = 'Godt, J.W., Sener, B., Verdin, K.L., Wald, D.J., Earle, P.S., Harp, E.L. and Jibson, R.W., 2008, Rapid Assessment of Earthquake-induced Landsliding: Proceedings of the First World Landslide Forum, United Nations University, Tokyo, Japan, p. 392-395.'
    shortref = 'Godt and others (2008)'

    #which type of ground failure model is this? Options are landslide or liquefaction.
    gfetype = landslide
    
    [[[layers]]]
      [[[[cohesion]]]]
        file = Godt_inputs/cohesion_10i.flt
        units = kPa
        longref = 'Godt, J.W., Sener, B., Verdin, K.L., Wald, D.J., Earle, P.S., Harp, E.L. and Jibson, R.W., 2008, Rapid Assessment of Earthquake-induced Landsliding: Proceedings of the First World Landslide Forum, United Nations University, Tokyo, Japan, p. 392-395.'
        shortref = 'Godt and others (2008)'

      [[[[friction]]]]
        file = Godt_inputs/friction.flt
        units = degrees
        longref = 'Godt, J.W., Sener, B., Verdin, K.L., Wald, D.J., Earle, P.S., Harp, E.L. and Jibson, R.W., 2008, Rapid Assessment of Earthquake-induced Landsliding: Proceedings of the First World Landslide Forum, United Nations University, Tokyo, Japan, p. 392-395.'
        shortref = 'Godt and others (2008)'

      [[[[slope]]]]
        file = Verdin_slopes_resampled
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
        file = Godt_inputs/cohesion_10i.flt
        interpolation = nearest
        units = kPa
        longref = 'Godt, J.W., Sener, B., Verdin, K.L., Wald, D.J., Earle, P.S., Harp, E.L. and Jibson, R.W., 2008, Rapid Assessment of Earthquake-induced Landsliding: Proceedings of the First World Landslide Forum, United Nations University, Tokyo, Japan, p. 392-395.'
        shortref = 'Godt and others (2008)'

      [[[[friction]]]]
        file = Godt_inputs/friction.flt
        units = degrees
        longref = 'Godt, J.W., Sener, B., Verdin, K.L., Wald, D.J., Earle, P.S., Harp, E.L. and Jibson, R.W., 2008, Rapid Assessment of Earthquake-induced Landsliding: Proceedings of the First World Landslide Forum, United Nations University, Tokyo, Japan, p. 392-395.'
        shortref = 'Godt and others (2008)'

      [[[[slope]]]]
        file = Verdin_slopes_resampled/slope_max.bil
        units = degrees
        longref = 'Verdin, D.W., Godt, J., Funk, C., Pedreros, D., Worstell, B. and Verdin, J., 2007, Development of a global slope dataset for estimation of landslide occurrence resulting from earthquakes: U.S. Geological Survey Open-File Report 2007–1188, 25p.'
        shortref = 'Verdin and others (2007)'

      [[[[watertable]]]]
        file = Fan2013WaterTable/wtd_fan2013_zhu_fil_na.grd
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
        file = Wills_et_al_2011/Resampled/Susceptiblity_dry_WGS84.bil
        units = N/A
        longref = 'Wills, C.J., Perez, F.G., and Gutierrez, C.I., 2011, Susceptibility to Deep-Seated Landslides in California: California Geological Survey Map Sheet 58, 1p.'
        shortref = 'Wills and others (2011)'

    [[[parameters]]]
      #Newmark displacement threshold for failure in cm - only need for some models
      dnthresh = 5

[logistic_models]

  #default_landslide and default_liquefaction parameters below must refer to named models in this file
  default_landslide = nowicki_2014_global
  default_liquefaction = zhu_2015

  # this can be any string, but it must be a unique and descriptive name of a logistic regression model.
  [[nowicki_2014_global]]
    #Detailed description of the model, its inputs, etc.
    description = 'This is the original landslide model of Nowicki et al 2014 using PGA, Slope, Cohesion, and CTI based on global datasets.'
    longref = 'Nowicki, M.A., Wald, D.J., Hamburger, M.W., Hearne, Michael, and Thompson, E.M., 2014, Development of a globally applicable model for near real-time prediction of seismically induced landslides: Engineering Geology, v. 173, p. 54–65.'
    shortref = 'Nowicki and others (2014)'
    
    #which type of ground failure model is this? Options are landslide or liquefaction.
    gfetype = landslide

    #what is the grid to which all other grids in this model will be resampled?
    baselayer = slope

    slopemin = 5. # in degrees
    slopemax = 90. # in degrees

    #these layer files can be any Grid2D-subclass supported format
    #These include, but may not be limited to:
    #https://github.com/usgs/MapIO/blob/master/mapio/gdal.py
    #https://github.com/usgs/MapIO/blob/master/mapio/gmt.py
    [[[layers]]]
      [[[[slope]]]]
        file = slope_max.bil
        units = degrees*100
        longref = 'Verdin, D.W., Godt, J., Funk, C., Pedreros, D., Worstell, B. and Verdin, J., 2007, Development of a global slope dataset for estimation of landslide occurrence resulting from earthquakes: U.S. Geological Survey Open-File Report 2007–1188, 25p.'
        shortref = 'Verdin et al. (2007)'
      [[[[friction]]]]
        file = friction.flt
        units = degrees
        longref = 'Godt, J.W., Sener, B., Verdin, K.L., Wald, D.J., Earle, P.S., Harp, E.L. and Jibson, R.W., 2008, Rapid Assessment of Earthquake-induced Landsliding: Proceedings of the First World Landslide Forum, United Nations University, Tokyo, Japan, p. 392-395.'
        shortref = 'Godt and others (2008)'
      [[[[cti1]]]]  # Must use cti1 or else the cti in the word friction will get changed
        file = global_cti_fil.grd
        units = index
        longref = 'USGS HYDRO1k geographic database, available at https://lta.cr.usgs.gov/HYDRO1K'
        shortref = 'HYDRO1k'

    [[[interpolations]]]
      slope = linear
      friction = nearest
      cti1 = linear

    [[[terms]]]
      #These terms must be named as b1-bN, where N is the number of coefficients
      #in a logistic regression, which takes the form:
      #1/(1 + e^-eqn)
      #where eqn is a linear equation of the form:
      #b0 + b1*t1 + b2*t2 + ... + bN*tN
      #where t1, t2, ... tN are the right hand side of the parameters below.
      #The terms may include the names of layers and any of the following ShakeMap macros:
      #pga,pgv,mmi,MW
      b1 = pga
      b2 = 'slope / 100.'  # Divide slopes by 100 because Verdin dataset multiplies by 100
      b3 = friction
      b4 = cti1 * 100. # Multiply by 100 because Anna used a layer representing CTI * 100 in her model
      b5 = 'pga * slope / 100.' # Divide slopes by 100 because Verdin dataset multiplies by 100

    [[[coefficients]]]
      #These coefficients must be named as b1-bN, where N is the number of coefficients
      #in a logistic regression, which takes the form:
      #1/(1 + e^-eqn)
      #where eqn is a linear equation of the form:
      #b0 + b1*t1 + b2*t2 + ... + bN*tN
      #where t1, t2, ... tN are the right hand side of the parameters below.
      b0 = -3.6490 #intercept
      b1 = 0.0133 #pga
      b2 = 0.0364 #slope
      b3 = -0.0635 #friction
      b4 = -0.0004 # cti
      b5 = 0.0019 # pga*slope

  [[nowicki_2015]]
    #Detailed description of the model, its inputs, etc.
    description = 'This is the Nowicki Model of 2015, which uses precip, lithology, and land cover.'
    
    #which type of ground failure model is this? Options are landslide or liquefaction.
    gfetype = landslide

    #what is the grid to which all other grids in this model will be resampled?
    baselayer = slope 

    slopemin = 5. # in degrees
    slopemax = 90. # in degrees

    #these layer files can be any Grid2D-subclass supported format
    #These include, but may not be limited to:
    #https://github.com/usgs/MapIO/blob/master/mapio/gdal.py
    #https://github.com/usgs/MapIO/blob/master/mapio/gmt.py
    [[[layers]]]
      [[[[slope]]]]
        file = gted_maxslope_30c.flt
        units = degrees
        longref = 'Global Multi-resolution Terrain Elevation Data 2010 (GMTED2010) available at http://topotools.cr.usgs.gov/gmted_viewer/'
        shortref = 'GMTED2010'
      [[[[rock]]]]
        file = glim_copy.grd
        units = lithology
        longref = 'Hartmann, Jens and Moosdorf, Nils, 2012, The new global lithological map database GLiM: A representation of rock properties at the Earth surface, G3, vol 13, no. 12., 37 p.'
        shortref = 'Hartmann and Moosdorf (2012)'
      [[[[landcover]]]]
        file = modis_30c_copy.grd
        units = none
        longref = 'Moderate resolution imaging spectroradiometer (MODIS) land cover dataset, http://modis.gsfc.nasa.gov/'
        shortref = 'MODIS land cover'
      [[[[precip]]]]
        file = precip
        units = millimeters
        longref = 
        shortref = 
      [[[[cti]]]]
        #file = globalcti.grd
        file = global_cti_fil.grd
        units = index
        longref = 'USGS HYDRO1k geographic database, available at https://lta.cr.usgs.gov/HYDRO1K'
        shortref = 'HYDRO1k'
      [[[[elev]]]]
        file = gted_meanelev_30c.flt
        units = m
        longref = 'Global Multi-resolution Terrain Elevation Data 2010 (GMTED2010) available at http://topotools.cr.usgs.gov/gmted_viewer/'
        shortref = 'GMTED2010'

    [[[interpolations]]]
      slope = linear
      rock = nearest
      landcover = nearest
      precip = nearest
      cti = linear
      elev = linear
      
    [[[terms]]]
      #These terms must be named as b1-bN, where N is the number of coefficients
      #in a logistic regression, which takes the form:
      #1/(1 + e^-eqn)
      #where eqn is a linear equation of the form:
      #b0 + b1*t1 + b2*t2 + ... + bN*tN
      #where t1, t2, ... tN are the right hand side of the parameters below.
      #The terms may include the names of layers and any of the following ShakeMap macros:
      #pga,pgv,mmi,MW
      b1 = log(pgv)
      b2 = slope / 90.
      b3 = rock
      b4 = cti * 100.
      b5 = MW
      b6 = precipMONTH
      b7 = landcover  
      b8 = elev
      b9 = log(pgv) * slope / 90.     

    [[[coefficients]]]
      #These coefficients must be named as b1-bN, where N is the number of coefficients
      #in a logistic regression, which takes the form:
      #1/(1 + e^-eqn)
      #where eqn is a linear equation of the form:
      #b0 + b1*t1 + b2*t2 + ... + bN*tN
      #where t1, t2, ... tN are the right hand side of the parameters below.
      b0 = -8.3453199   # intercept
      b1 = 1.737721 # log(pgv)
      b2 = 2.1773966 #slope
      b3 = 1 #lithology set to 1.0 - coefficients are in glim file
      b4 = 0.0484136 # cti
      b5 = 0.1634385 # moment magnitude
      b6 = 0.000949 # precip
      b7 = 1.0 # landcover
      b8 = 0.0002273 # elev
      b9 = 0.477635 # log(pgv)*slope

  [[zhu_2015]]
  
    #Detailed description of the model, its inputs, etc.

    longref = 'Zhu, Jing; Daley, Davene; Baise, L.G.; Thompson, E.M.; Wald, D.J.; and Knudsen, K.L., 2015b, A geospatial liquefaction model for rapid response and loss estimation: Earthquake Spectra, v. 31, no. 3, p. 1813–1837.'
    shortref = 'Zhu and others (2015)'
    description = 'Zhu coastal model'

    gfetype = liquefaction
  
    baselayer = vs30

    slopemin = 0. # in degrees
    slopemax = 5. # in degrees

    #these layer files can be any Grid2D-subclass supported format
    #These include, but may not be limited to:
    #https://github.com/usgs/MapIO/blob/master/mapio/gdal.py
    #https://github.com/usgs/MapIO/blob/master/mapio/gmt.py

    [[[layers]]]
      [[[[vs30]]]]
        file = global_vs30.grd
        units = m/s
        longref = 'Computed from GMTED2010 using methods of Wald and Allen (2007) based on topographic slope'
        shortref = 'Wald and Allen (2007)'
      [[[[cti]]]]
        #file = globalcti.grd
        file = global_cti_fil.grd
        units = index
        longref = 'USGS HYDRO1k geographic database, available at https://lta.cr.usgs.gov/HYDRO1K'
        shortref = 'HYDRO1k'

    [[[interpolations]]]
      vs30 = nearest
      cti = linear

    [[[terms]]]
      b1 = 'log((pga/100.0)*(power(MW,2.56)/power(10,2.24)))'
      b2 = 'cti'
      b3 = 'log(vs30)'

    [[[coefficients]]]
      b0 = 24.10
      b1 = 2.067
      b2 = 0.355
      b3 = -4.784

  [[zhu_2016_coastal]]
    # Model meant for proximity to coasts/oceans

    longref = 'Zhu, Jing;  Baise, Laurie; Thompson, Eric, 2016, An Updated Geospatial Liquefaction Model for Global Use, in submission.'
    shortref = 'Zhu and others (2016)'
    description = 'Zhu coastal model'

    gfetype = liquefaction
    baselayer = vs30

    slopemin = 0. # in degrees
    slopemax = 5. # in degrees

    # layer files
    [[[layers]]]
      [[[[vs30]]]]
        file = global_vs30.grd
        units = m/s
        longref = 'Computed from GMTED2010 using methods of Wald and Allen (2007) based on topographic slope'
        shortref = 'Wald and Allen (2007)'
      [[[[precip]]]]
        file = global_precip_fil.grd
        units = millimeters
        longref = 'WorldClim database, http://WorldClim.org; last accessed March 2014'
        shortref = WorldClim
      [[[[dc]]]]
        file = global_dc.grd
        units = km
        longref = 'Computed from a global dataset computed by NASA's Ocean Color Group, http://oceancolor.gsfc.nasa.gov/cms/DOCS/DistFromCoast; last accessed January 2014'
        shortref = 'NASA Ocean Color group'
      [[[[dr]]]]
        file = global_dr.grd
        units = km
        longref = 'Computed from HydroSHEDS database, http://hydrosheds.cr.usgs.gov/dataavail.php; last accessed February 2014'
        shortref = 'HydroSHEDS'

    [[[interpolations]]]
      vs30 = linear
      precip = linear
      dc = linear
      dr = linear

    [[[terms]]]
      b1 = log(pgv)
      b2 = log(vs30)
      b3 = precip
      b4 = "power(dc, 0.5)"
      b5 = dr
      b6 = "power(dc, 0.5) * dr"

    [[[coefficients]]]
      b0 = 12.435
      b1 = 0.301
      b2 = -2.615
      b3 = 0.0005556
      b4 = -0.0287
      b5 = 0.0666
      b6 = -0.0369

  [[zhu_2016_general]]
    # Meant as a generalized option to the above model
    longref = 'Zhu, Jing;  Baise, Laurie; Thompson, Eric, 2016, An Updated Geospatial Liquefaction Model for Global Use, in submission.'
    shortref = 'Zhu and others (2016)'
    
    gfetype = liquefaction

    baselayer = vs30

    slopemin = 0. # in degrees
    slopemax = 5. # in degrees

    description = 'This is the general model from Zhu et al. 2015'

    [[[layers]]]
      [[[[vs30]]]]
        file = global_vs30.grd
        units = m/s
        longref = 'Computed from GMTED2010 using methods of Wald and Allen (2007) based on topographic slope'
        shortref = 'Wald and Allen (2007)'
      [[[[precip]]]]
        file = global_precip_fil.grd
        units = millimeters
        longref = 'WorldClim database, http://WorldClim.org; last accessed March 2014'
        shortref = WorldClim
      [[[[dc]]]]
        file = global_dc.grd
        units = km
        longref = 'Computed from a global dataset computed by NASA's Ocean Color Group, http://oceancolor.gsfc.nasa.gov/cms/DOCS/DistFromCoast; last accessed January 2014'
        shortref = 'NASA Ocean Color group'
      [[[[dr]]]]
        file = global_dr.grd
        units = km
        longref = 'Computed from HydroSHEDS database, http://hydrosheds.cr.usgs.gov/dataavail.php; last accessed February 2014'
        shortref = 'HydroSHEDS'
      [[[[wtd]]]]
        file = wtd_fan2013_zhu_fil_na.grd
        units = m
        longref = 'Fan, Y., Li, H., and Miguez-Macho, G., 2013, Global Patterns of Groundwater Table Depth: Science, 339, 940-943.'
        shortref = 'Fan and others (2013)'

    [[[interpolations]]]
      vs30 = linear
      precip = linear
      dc = linear
      dr = linear
      wtd = linear

    [[[terms]]]
      b1 = log(pgv)
      b2 = log(vs30)
      b3 = precip
      b4 = "minimum(dc, dr)"
      b5 = wtd

    [[[coefficients]]]
      b0 = 8.801
      b1 = 0.334
      b2 = -1.918
      b3 = 0.0005408
      b4 = -0.2054
      b5 = -0.0333
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