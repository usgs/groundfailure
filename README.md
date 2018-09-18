[![Build Status](https://travis-ci.org/usgs/groundfailure.svg?branch=master)](https://travis-ci.org/usgs/groundfailure)
[![codecov](https://codecov.io/gh/usgs/groundfailure/branch/master/graph/badge.svg)](https://codecov.io/gh/usgs/groundfailure)

# groundfailure

## Introduction

This software is for calculating earthquake-induced ground failures (i.e.,
landslide and liquefaction). These models are intended for regional
or global scale applications, and are intended to be distributed in
near-real-time, triggered by the Shakemaps. 

## Documentation

Information about the methodology can be found on the 
[Ground Failure Scientific Background webpage](https://earthquake.usgs.gov/data/ground-failure/background.php)
and the corresponding [Ground Failure References webpage](https://earthquake.usgs.gov/data/ground-failure/references.php).

The API docs can be found [here](http://usgs.github.io/groundfailure/). 
Besides the API, there are five command-line programs:

`gfail` - runs ground failure models

`callgf` - automation wrapper for gfail

`gfail_transfer` - transfers model results to USGS comcat

`create_info` - creates info.json required for web rendering

`create_png` - creates transparent png of model results required for web rendering

Documentation for the use of these programs can be seen by calling them
with the `-h` flag. 

## Installation and Dependencies

1. If a current version of conda is not already installed, install Miniconda 
    with Python 3.6 (or greater) following the directions provided on the 
    [conda webpage.](https://conda.io/docs/user-guide/install/index.html). 
    Anaconda will also work, but is a larger installation and is not necessary
    unless you want to use it for other purposes. Take note of the folder name
    where it is installed (e.g., miniconda or miniconda3)

2. The current version of miniconda requires that you manually edit your .bash_profile.
    Make the following changes, updating the path below with whatever folder miniconda was installed in:
    * If the installation added a line that looks like this, delete it:
        export PATH="/Users/YourName/miniconda3/bin:$PATH
    * add this line:
        . $HOME/miniconda3/etc/profile.d/conda.sh
    * Save and exit and either close the terminal and open a new one or source
        the .bash_profile ```source ~/.bash_profile```
    * Type ```which conda``` in terminal to make sure conda is found.

3. Clone the groundfailure repository in the location where you want it installed:
```sh
cd Users/YourName
git clone https://github.com/usgs/groundfailure.git
```
There will now be a folder called groundfailure in Users/YourName that contains
all of the files.

4. Run the install.sh script located in the main repository directory:
```sh
cd groundfailure
bash install.sh
```
This will take a while and will show numerous dependencies being installed.

5. The previous step installs a self-contained virtual environment called gf.
    To ensure the virtual environment was successfully installed,
    type ```conda activate gf```. You will need to activate the gf environment
    every time you want to run groundfailure.

6. With the gf virtual environment active, type ```gfail -h``` to ensure gfail
    was correctly installed. If successful, you will see the help section of
    the gfail code.

### Updating

To ensure all of your dependencies are up to date, reinstall completely starting
at Step 3 above.

To update groundfailure to the current master branch without altering dependencies
(if you have altered the master branch, you will first need to stash your changes):
```sh
cd Users/YourName/groundfailure
git pull
```

### Uninstalling

To uninstall, delete the virtual environment:
```sh
conda remove --name gf --all
```
And remove the groundfailure folder that was cloned in step 3.

### Troubleshooting

* Check step 2 from the installation steps above, make sure paths in .bash_profile are correct
    and point to the actual location of miniconda on your machine.

* Try opening a new terminal in case the updated .bash_profile was not sourced in the current terminal window.

* Uninstall (or move) your current anaconda or conda installation and reinstall from scratch. 
    Due to recent conda updates, older preexisting installations of anaconda or 
    miniconda may not function with our installer.
    
* Ensure that miniconda is in your user directory or somewhere that does not
    require admin permissions.

### Releases

* [v1.0](https://github.com/usgs/groundfailure/releases/tag/1.0) of this software was reviewed and released on 11 May 2018. 
* The Digital Object Identifier of this software is: https://doi.org/10.5066/P91G4NS4

## Dependencies

The install.sh script installs this package and dependencies. It is
regularly tested on OSX and Ubuntu. For a full list of dependencies, refer to
[environment.yml](https://github.com/usgs/groundfailure/blob/master/environment.yml).

Some functions of this program require the use of the USGS Product Distribution
Layer (PDL). This must be installed separately. See the [PDL User Guide](https://usgs.github.io/pdl)
for installation information.

## Configuration

For each model, there is a configuration file that describes the model
parameters, default values, metadata/source details, input file locations,
and display preferences. Default versions with relative file paths are found
in the defaultconfigfiles folder of the
[repository](https://github.com/usgs/groundfailure/tree/master/defaultconfigfiles).
These can be edited but to avoid overwriting your changes each time you update
the groundfailure codes, you should edit copies outside of the repository.

Default options, such as the output directory, paths to input data files, paths
to mapping files etc. can be specified using gfail, see below for details.

### Using gfail to set default paths on system (gfail must be in system path)

#### set default paths

```sh
gfail --set-default-paths \
    -d full/modelinput/data/path \
    -o full/output/location/filepath \
    -c full/filepath/to/model/config/files \
    -m full/filepath/to/mapping/config/file \
    -md full/filepath/to/data/for/mapping
```

#### check default paths that are currently set

```sh
gfail --list-default-paths
```

#### clear all default paths

```sh
gfail --reset-default-paths
```

#### after setting default paths, gfail can be run like this:

```sh
gfail modelconfig.ini shakefile.xml -s -pd -pi
```

* the pd flag outputs static pdfs of model results
* the pi flag creates interactive plots as html files
* type gfail -h to see all options

### Model config file format

The config file format is a modified version of the "INI" format.  It is
described in detail
[here](http://configobj.readthedocs.org/en/latest/configobj.html#config-files).

**Notes** 
* References and other strings with commas within them need to be enclosed in
  quotes. Example:
  * 'Verdin, D.W., Godt, J., Funk, C., Pedreros, D., Worstell, B. and Verdin,
    J., 2007, Development of a global slope dataset for estimation of landslide
    occurrence resulting from earthquakes: U.S. Geological Survey Open-File
    Report 2007–1188, 25p.'
* Arrays should be not be enclosed in brackets and should be comma separated.
  Example:
  * model = 0, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.
* Files, filepaths, and folders being used for input layers should all be labeled
  as 'file' in the config file regardless of actual structure.
* If file paths in the config file are relative to a base folder on the users
  system, and you are running models manually and/or not setting default paths in
  gfail you will need to run correct_config_filepaths() after reading in the
  config file. See notebooks for details.

```ini
[nowicki_2014_global]
  # Detailed description of the model and references.
  description = 'This is the original landslide model of Nowicki et al 2014 using PGA, Slope, Cohesion, and CTI based on global datasets.'
  longref = """Nowicki, M.A., Wald, D.J., Hamburger, M.W., Hearne, Michael, and Thompson, E.M., 2014, Development of a 
            globally applicable model for near real-time prediction of seismically induced landslides: Engineering 
            Geology, v. 173, p. 54–65."""
  shortref = 'Nowicki and others (2014)'
  
  # Ground failure model type, pptions are landslide or liquefaction.
  gfetype = landslide

  # The grid to which all other grids in this model will be resampled
  baselayer = slope

  # Thresholds for filtering results, probabilities will be set to zero below
  # slopemin, above slopemax using the slopefile defined
  slopemin = 5. # in degrees
  slopemax = 90. # in degrees
  slopefile = global_Verdin_slopes_resampled_degx100/slope_max.bil
  slopemod = 'slope/100.'
  # Slopemod indicates how the slopefile should be modified so that it is in
  # degrees (use np.function for any functions)

  # groundfailure function corresponding to this model
  funcname = LogisticModel

  [[layers]]
    [[[slope]]]
      file = global_Verdin_slopes_resampled_degx100/slope_max.bil
      units = degrees*100
      longref = """Verdin, D.W., Godt, J., Funk, C., Pedreros, D., Worstell, B. and Verdin, J., 2007, Development of a 
                global slope dataset for estimation of landslide occurrence resulting from earthquakes: U.S. Geological 
                Survey Open-File Report 2007–1188, 25p."""
      shortref = 'Verdin et al. (2007)'
    [[[friction]]]
      file = global_friction_deg.flt
      units = degrees
      longref = """Godt, J.W., Sener, B., Verdin, K.L., Wald, D.J., Earle, P.S., Harp, E.L. and Jibson, R.W., 2008, 
                   Rapid Assessment of Earthquake-induced Landsliding: Proceedings of the First World Landslide Forum, 
                   United Nations University, Tokyo, Japan, p. 392-395."""
      shortref = 'Godt and others (2008)'
    [[[cti1]]] # Must use cti1 or else the cti in the word friction will get changed
      file = global_cti_fil.grd
      units = index
      longref = 'USGS HYDRO1k geographic database, available at https://lta.cr.usgs.gov/HYDRO1K'
      shortref = 'HYDRO1k'

  [[interpolations]]
    slope = linear
    friction = nearest
    cti1 = linear

  [[terms]]
    # These terms must be named as b1-bN, where N is the number of coefficients
    # in a logistic regression, which takes the form:
    # 1/(1 + e^-eqn)
    # where eqn is a linear equation of the form:
    # b0 + b1*t1 + b2*t2 + ... + bN*tN
    # where t1, t2, ... tN are the right hand side of the parameters below.
    # The terms may include the names of layers and any of the following 
    # ShakeMap macros: pga,pgv,mmi,MW
    b1 = pga
    b2 = 'slope / 100.'
    # Divide slopes by 100 because Verdin dataset multiplies by 100
    b3 = friction
    b4 = cti1 * 100.
    # Multiply by 100 because Anna used a layer representing CTI * 100 in her
    # model.
    b5 = 'pga * slope / 100.'
    # Divide slopes by 100 because Verdin dataset multiplies by 100

  [[coefficients]]
    # These coefficients must be named as b1-bN, where N is the number of coefficients
    # in a logistic regression, which takes the form:
    # 1/(1 + e^-eqn)
    # where eqn is a linear equation of the form:
    # b0 + b1*t1 + b2*t2 + ... + bN*tN
    # where t1, t2, ... tN are the right hand side of the parameters below.
    # intercept
    b0 = -3.6490
    # pga
    b1 = 0.0133
    # slope
    b2 = 0.0364
    # friction
    b3 = -0.0635
    # cti
    b4 = -0.0004
    # pga*slope
    b5 = 0.0019 

  [[display_options]]
    # Display preferences for mapping programs (not required)
    [[[lims]]]  # Optional
      # Corresponding to different possible layer keys - don't need these, will
      # just use defaults if missing, don't need full name of layer, just
      # something that is part of it.
      model = None
      friction = None
      pga = None
      slope = None
      cti1 = None

    [[[colors]]]
      # use matplotlib colormaps, starting with cm.
      default = cm.jet
      alpha = 0.7
      # Corresponding to different possible layer keys - don't need these, will
      # just use defaults if missing.
      model = cm.jet
      friction = cm.jet
      pga = cm.jet
      slope = cm.jet
      cti1 = cm.jet

    [[[logscale]]]
      # Whether to use a log scale in figures for given layer.
      # Corresponding to different possible layer keys - don't need these, will
      # just use defaults if missing, don't need full name of layer, just
      # something that is part of it.
      model = False
      friction = False
      pga = False
      slope = False
      cti1 = False

    [[[maskthresholds]]]
      # Model result cells with values below this threshold will be masked in plots
      # Corresponding to different possible layer keys - don't need these, will just use defaults if missing,
      # don't need full name of layer, just something that is part of it
      model = 0.05
      friction = 0.
      pga = None
      slope = None
      cti1 = None
```

### Mapping config file format
The mapping config file tells the program what files to use (roads, cities etc.)
for creating static maps. It is not required, but static maps may be bland
without these.

```ini
# All map inputs are optional, but desireable. File paths are relative to the
# mapping_inputs folder. Default display options for each input layer are found
# in the model config file. 

[dem]  # Optional, for making hillshade
  file = md30_gmted_gmt.grd

[roads]  # Optional
  file = roads  # Folder in this case
  longref = """Center for International Earth Science Information Network - CIESIN, 2013, Global Roads Open Access Data Set,
   Version 1 (gROADSv1): Columbia University, and Information Technology Outreach Services - ITOS - University of Georgia,
   Palisades, NY, NASA Socioeconomic Data and Applications Center (SEDAC). http://dx.doi.org/10.7927/H4VD6WCT.""
  shortref = 'CIESIN (2013)"""

[cities]  # Optional
  file = cities1000.txt
  longref = GeoNames, http://geonames.org/ Accessed: 2 Sept 2015
  shortref = GeoNames

[ocean]  # Optional
  file = ne_10m_ocean/ne_10m_ocean.shp
  longref = Natural Earth (acc. 2016)
  shortref = Natural Earth (2016) Ocean polygon http://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-ocean/

[colors]
  # Define basic colors
  roadcolor = 808080
  countrycolor = 474747
  watercolor = B8EEFF

```

## API for Model Output

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

```py
def failure_model():
    geodict = GeoDict({
        'xmin':0.5,'xmax':3.5,
        'ymin':0.5,'ymax':3.5,
        'dx':1.0,'dy':1.0,
        'nx':4,'ny':4
    })
    pgrid = Grid2D(data = np.arange(0,16).reshape(4,4),geodict=geodict)
    cgrid = Grid2D(data = np.arange(1,17).reshape(4,4),geodict=geodict)
    sgrid = Grid2D(data = np.arange(2,18).reshape(4,4),geodict=geodict)
    mgrid = Grid2D(data = np.arange(3,19).reshape(4,4),geodict=geodict)

    modellayer = {
        'description':{
	    'name':'Nowicki 2014',
            'longref':'Nowicki, A., 2014, A logistic regression landslide model: Failure Monthly, v. 2, p. 1-7.',
            'units':'index',
            'shakemap': '19940117123055_ver2'
            'parameters':{
	        'b0':1.045,
                'b1':5.435
	    }
        },
        'type':'output',
        'label':'Relative Index Value',
        'grid':pgrid,
    }
    
    layer1 = {
        'description':{
	    'name':'Smith and Jones 1994',
            'longref':'Smith J. and Jones, J., 1994, Holding on to things: Journal of Geophysical Sciences, v. 17,  p. 100-105',
            'units':'kPa'
	},
        'type':'input',
        'label':'cohesion (kPa)',
        'grid':cgrid
    }
    
    layer2 = {
        'description':{
	    'name':'Garfunkel and Oates 2001',
            'longref':'Garfunkel, A., and Oates, J., 2001, I'm afraid to look down: Journal of Steepness, v. 8, p. 10-25',
            'units':'degrees'
	},
        'type':'input',
        'label':'slope (degrees)',
        'grid':sgrid
    }

    layer3 = {
        'description':{
	    'units':'g'
            'shakemap': '19940117123055_ver2'
	},
        'type':'input',
        'label':'PGA (g)',
        'grid':mgrid
    }

    output = {
        'model':problayer,
        'cohesion':layer1,
        'slope':layer2,
        'pga':layer3
    }

    return output
```

## Sources of test datasets

### Datasets for example notebooks

We have extracted the input datasets required to run the models demonstrated in the [example notebooks](https://github.com/usgs/groundfailure/tree/master/notebooks) for the 
1994 Northridge, CA, earthquake, including the [USGS ShakeMap, from the ShakeMap Atlas (v1)](https://earthquake.usgs.gov/earthquakes/eventpage/ci3144585#shakemap).
The sources of the input files are listed in the [default config file](https://github.com/usgs/groundfailure/tree/master/defaultconfigfiles/models) for each model.
The reference for each filename is listed as a  "longref" ([example:](https://github.com/usgs/groundfailure/blob/master/defaultconfigfiles/models/jessee_2017.ini#L27)) in the section of the
config file below the corresponding filename, defined as "file: ([example](https://github.com/usgs/groundfailure/blob/master/defaultconfigfiles/models/jessee_2017.ini#L25)).
A digital terrain model is also provided for mapping purposes. It is extracted from the [GMTED2010 Terrain Elevation model](https://topotools.cr.usgs.gov/gmted_viewer).
The extracted input data files are located with the notebooks in the [data folder](https://github.com/usgs/groundfailure/tree/master/notebooks/data).

### Datasets for testing

Test input datasets are included with the repository in order to run the [tests](https://github.com/usgs/groundfailure/tree/master/tests).
Some tests used artificial datasets, but others use input datasets for a subsection
of the area affected by the 1989 Loma Prieta, CA, earthquake. These extracted
sections of the input datasets are located with the tests in the [data folder](https://github.com/usgs/groundfailure/tree/master/tests/data/loma_prieta).
The input layers for each model can be found in the [default config file](https://github.com/usgs/groundfailure/tree/master/defaultconfigfiles/models) for each model,
as described above. Additional layers used in the tests were extracted from the global input layers defined below:

* ne_10m_ocean: [Natural Earth (2016) Ocean polygon](http://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-ocean) last accessed 17 Nov 2017
* cities1000.txt: Global city information from [GeoNames](http://geonames.org) last accessed 2 Sept 2015.
* md30_gmted_gmt.grd: [GMTED2010 Terrain Elevation model](https://topotools.cr.usgs.gov/gmted_viewer)
* gmted_global_hillshade.grd: Hillshade created from md30_gmted_gmt.grd
* lspop2016_lp.flt: [LandScan (2016)™ High Resolution global Population Data Set](https://landscan.ornl.gov/)
