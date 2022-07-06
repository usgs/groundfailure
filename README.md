
[![codecov](https://codecov.io/gh/usgs/groundfailure/branch/main/graph/badge.svg)](https://codecov.io/gh/usgs/groundfailure)

| OS   | Python | Status |
| :--- | :----   | :--- |
| macOS-latest | 3.8  | [![Build Status](https://dev.azure.com/GHSC-ESI/USGS-groundfailure/_apis/build/status/usgs.groundfailure?branchName=main&jobName=groundfailure&configuration=groundfailure%20MacOS_py38)](https://dev.azure.com/GHSC-ESI/USGS-groundfailure/_build/latest?definitionId=7&branchName=main) |
| macOS-latest | 3.9  | [![Build Status](https://dev.azure.com/GHSC-ESI/USGS-groundfailure/_apis/build/status/usgs.groundfailure?branchName=main&jobName=groundfailure&configuration=groundfailure%20MacOS_py39)](https://dev.azure.com/GHSC-ESI/USGS-groundfailure/_build/latest?definitionId=7&branchName=main) |
| ubuntu-latest | 3.8  | [![Build Status](https://dev.azure.com/GHSC-ESI/USGS-groundfailure/_apis/build/status/usgs.groundfailure?branchName=main&jobName=groundfailure&configuration=groundfailure%20Linux_py38)](https://dev.azure.com/GHSC-ESI/USGS-groundfailure/_build/latest?definitionId=7&branchName=main) |
| ubuntu-latest | 3.9  | [![Build Status](https://dev.azure.com/GHSC-ESI/USGS-groundfailure/_apis/build/status/usgs.groundfailure?branchName=main&jobName=groundfailure&configuration=groundfailure%20Linux_py39)](https://dev.azure.com/GHSC-ESI/USGS-groundfailure/_build/latest?definitionId=7&branchName=main) |

# groundfailure

## Introduction

This software provides spatial estimates of earthquake-induced ground failure hazard (i.e., landslide and liquefaction) along with qualitative hazard and population exposure-based alert levels. These models are intended for regional
or global scale applications, and are intended to be distributed in
near-real-time by the USGS, triggered by the Shakemaps.

## Documentation

The methodology is detailed by [Allstadt and others (2015)](https://doi.org/10.1177/87552930211032685), and summarized on the 
[Ground Failure Scientific Background webpage](https://earthquake.usgs.gov/data/ground-failure/background.php)
and the corresponding [Ground Failure References webpage](https://earthquake.usgs.gov/data/ground-failure/references.php).

The software documentation can be found [here](http://usgs.github.io/groundfailure/). There are seven command-line programs, most non-operational users would only need to use `gfail`:

`gfail` - runs one or more ground failure models

`callgf` - automation wrapper for gfail (A file called "autogf_models" that lists the models to run must be placed in
the data_path directory.)

`gfail_transfer` - transfers model results to USGS comcat

`create_info` - creates info.json required for web rendering

`create_png` - creates transparent png of model results required for web rendering

`cleandb` - cleans up the database file created by callgf if, for example, a run is interrupted or stalls.

`viewdb` - outputs summary information about the database created by callgf.

Documentation for the use of these programs can be seen by calling them
with the `-h` flag. 

## Installation and Dependencies

1. Clone the groundfailure repository in the location where you want it installed:
    ```sh
    cd Users/YourName
    git clone https://github.com/usgs/groundfailure.git
    ```
    There will now be a folder called groundfailure in Users/YourName that contains all of the files.

2. Navigate to the groundfailure folder in terminal. Run the install.sh script located in the main repository directory:
    ```sh
    cd groundfailure
    bash install.sh
    ```
    This will take a while and will show numerous dependencies being installed. Note that this installation script will install miniconda or anaconda first if you do not already have it installed.

3. The previous step installs a self-contained virtual environment called gf.
    To ensure the virtual environment was successfully installed,
    type ```conda activate gf```. You will need to activate the gf environment
    every time you want to run groundfailure.

4. With the gf virtual environment active, type ```gfail -h``` to ensure gfail
    was correctly installed. If successful, you will see the help section of
    the gfail code.

### Updating

To ensure all of your dependencies are up to date, reinstall completely starting at Step 2 above. The installation script will uninstall your current gf environment and then reinstall a new one.

To update groundfailure to the current main branch without altering dependencies (if you have altered the main branch, you will first need to stash your changes):
```sh
cd Users/YourName/groundfailure
git pull
```

### Uninstalling

To uninstall, delete the virtual environment:
```sh
conda remove --name gf --all
```
And remove the groundfailure folder that was cloned in step 1.

### Troubleshooting

* Make sure miniconda or anaconda were properly installed and added to your .bash_profile, this should have been done by the installation script if you didn't already have a conda option installed.

* Be sure you are in the gf conda environment.

* Try opening a new terminal in case the updated .bash_profile was not sourced in the current terminal window.

* Uninstall (or move) your current anaconda or conda installation and reinstall from scratch. Older preexisting installations of anaconda or miniconda may not function with our installer.
    
* Ensure that miniconda is in your user directory or somewhere that does not require admin permissions.

### Releases

* [v1.0](https://github.com/usgs/groundfailure/releases/tag/1.0) of this software was reviewed and released on 11 May 2018. 
* The Digital Object Identifier of this software is: https://doi.org/10.5066/P91G4NS4

## Dependencies

The install.sh script installs this package and dependencies. It is
regularly tested on OSX and Ubuntu. 

Some functions of this program require the use of the USGS Product Distribution Layer (PDL). This must be installed separately. See the [PDL User Guide](https://usgs.github.io/pdl) for installation information.

## Configuration

For each model, there is a configuration file that describes the model
type, input layer names, thresholds, metadata/source details, interpolation choices, and display preferences. Default model configuration files with relative file paths are found in the defaultconfigfiles folder of the
[repository](https://github.com/usgs/groundfailure/tree/main/defaultconfigfiles). These can be edited but to avoid overwriting your changes each time you update the groundfailure codes, you should edit and use copies outside of the repository.

Default options, such as the output directory, paths to input data files, paths to mapping files etc. can be specified using gfail, see below for details.

### Setting default paths on system 

The `gfail` command lets you set paths to all the input files on each run, or you can set some default paths that will be used anytime that the paths are not explicitly specified. Most path options are used only for operational near-real-time use, only the data_path, configfilepath, and output_filepath and optionally the coastlines trimming file (trimfile) are needed for local runs. You can set these defaults either using `gfail`'s `--set-default-paths` flag as demonstrated below: 

```sh
gfail --set-default-paths \
    -d /Users/YourName/model_inputs \
    -o /Users/YourName/outputs \
    -c /Users/YourName/groundfailure/defaultconfigfiles/models \
    -pf /Users/YourName/populationfile.flt \
    -tr /Users/YourName/coastlinefile.shp \
    -pdl /Users/YourName/ProductClient/config.ini \
    -log /Users/YourName/logs \
    -db /Users/YourName/events.db \
```
or alternatively, the user can instead manually create a text file called .gfail_defaults following the format example below and put that in their home directory, the outcome will be the same:

```
data_path = /Users/YourName/model_inputs
config_filepath = /Users/YourName/groundfailure/defaultconfigfiles/models
output_filepath = /Users/YourName/outputs
trimfile = /Users/YourName/coastlinefile.shp
popfile = /Users/YourName/populationfile.flt
log_filepath = /Users/YourName/logs
dbfile = /Users/YourName/events.db
pdl_config = /Users/YourName/ProductClient/config.ini
comcat_config = /Users/YourName/ProductClient/comcat.ini

```

#### Check default paths that are currently set

```sh
gfail --list-default-paths
```

#### Clear all default paths

```sh
gfail --reset-default-paths
```

#### After setting default paths, gfail can be run like this:

```sh
gfail modelconfig.ini shakefile.xml --gis --kmz
```

* the --gis flag outputs geotiff files of the models
* the --kmz flag creates stylized kmz files of the models
* type gfail -h to see all options

### Model config file format

The model config file format is a modified version of the "INI" format.  It is described in detail [here](http://configobj.readthedocs.org/en/latest/configobj.html#config-files).

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
[jessee_2018]
  #Detailed description of the model, its inputs, etc.
  description = 'This is the Nowicki Jessee Model, which uses lithology, and land cover.'
  longref = 'Nowicki Jessee, M.A., Hamburger, H.W., Allstadt, K.E., Wald, D.J., Robeson, S.M., Tanyas, H., Hearne, M., Thompson, E.M., 2018, A Global Empirical Model for Near Real-time Assessment of Seismically Induced Landslides, J. Geophys. Res. Earth Surface, 123, 1835-1859'
  shortref = 'Nowicki Jessee and others (2018)'
  
  #which type of ground failure model is this? Options are landslide or liquefaction.
  gfetype = landslide

  #what is the grid to which all other grids in this model will be resampled?
  baselayer = slope 

  slopemin = 2. # in degrees
  slopemax = 90. # in degrees
  slopefile = global_grad.tif

  # Default standard deviation value to use if no map is available. Units are in logit space. 
  default_stddev = 0.03

  # Model's maximum probability, used for computing beta distribution
  maxprob = 0.256

  # Confidence interval probabilities for computing quantiles. 
  # Note, +/- 1 std is 0.68, 2 std is 0.95. Comment out to turn off.
  # conf_int_probabilities = 0.68, 0.95

  # Location of code corresponding to this model
  funcname = LogBase

  [[layers]]
    [[[slope]]]
      file = global_grad.tif
      units = gradient
      longref = """Global Multi-resolution Terrain Elevation Data 2010 (GMTED2010) available at http://topotools.cr.usgs.gov/gmted_viewer/"""
      shortref = 'GMTED2010'
    [[[rock]]]
      file = GLIM_replace.tif
      units = lithology
      longref = """Hartmann, Jens and Moosdorf, Nils, 2012, The new global lithological map database GLiM: A representation of rock properties at the Earth surface, G3, vol 13, no. 12., 37 p."""
      shortref = 'Hartmann and Moosdorf (2012)'
    [[[landcover]]]
      file = globcover_replace.tif
      units = none
      longref = 'Moderate resolution imaging spectroradiometer (MODIS) land cover dataset, http://modis.gsfc.nasa.gov/'
      shortref = 'MODIS land cover'
    [[[cti]]]
      file = global_cti_fil.grd
      units = index
      longref = 'USGS HYDRO1k geographic database, available at https://lta.cr.usgs.gov/HYDRO1K'
      shortref = 'HYDRO1k'
    [[[stddev]]]
      file = jessee_standard_deviation.tif
      units = none
      longref = ''
      shortref = ''

  [[interpolations]]
    slope = linear
    rock = nearest
    landcover = nearest
    cti = linear
    stddev = linear

  [[display_options]]  # These only get used in mapping programs
    [[[lims]]]  # Optional
      # Corresponding to different possible layer keys - don't need these, will just use defaults if missing,
      # don't need full name of layer, just something that is part of it
      model = 0.002, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5
      pgv = None
      slope = None
      rock = None
      landcover = None
      cti = None

    [[[colors]]]
      default = cm.jet
      alpha = 0.8
      # Corresponding to different possible layer keys - don't need these, will just use defaults if missing
      model = cm.CMRmap_r
      pgv = cm.jet
      slope = cm.gnuplot2
      rock = cm.jet
      landcover = cm.jet
      cti = cm.jet

    [[[logscale]]]
      # Corresponding to different possible layer keys - don't need these, will just use defaults if missing,
      # don't need full name of layer, just something that is part of it
      model = True
      pgv = False
      slope = False
      rock = False
      cti = False
      landcover = False

    [[[maskthresholds]]]
      # Corresponding to different possible layer keys - don't need these, will just use defaults if missing,
      # don't need full name of layer, just something that is part of it
      model = 0.002
      pgv = None
      slope = None
      rock = None
      cti = None
      landcover = None
```

## API for Model Output

Each model should output a single dictionary, which has keys that correspond to the names of the input and output layers from the model.

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

We have extracted the input datasets required to run the models demonstrated in the [example notebooks](https://github.com/usgs/groundfailure/tree/main/notebooks) for the 
1994 Northridge, CA, earthquake, including the [USGS ShakeMap, from the ShakeMap Atlas (v1)](https://earthquake.usgs.gov/earthquakes/eventpage/ci3144585#shakemap).
The sources of the input files are listed in the [default config file](https://github.com/usgs/groundfailure/tree/main/defaultconfigfiles/models) for each model.
The reference for each filename is listed as a  "longref" ([example:](https://github.com/usgs/groundfailure/blob/main/defaultconfigfiles/models/jessee_2017.ini#L27)) in the section of the
config file below the corresponding filename, defined as "file: ([example](https://github.com/usgs/groundfailure/blob/main/defaultconfigfiles/models/jessee_2017.ini#L25)).
A digital terrain model is also provided for mapping purposes. It is extracted from the [GMTED2010 Terrain Elevation model](https://topotools.cr.usgs.gov/gmted_viewer).
The extracted input data files are located with the notebooks in the [data folder](https://github.com/usgs/groundfailure/tree/main/notebooks/data).

### Datasets for testing

Test input datasets are included with the repository in order to run the [tests](https://github.com/usgs/groundfailure/tree/main/tests).
Some tests used artificial datasets, but others use input datasets for a subsection of the area affected by the 1989 Loma Prieta, CA, earthquake. These extracted sections of the input datasets are located with the tests in the [data folder](https://github.com/usgs/groundfailure/tree/main/tests/data/loma_prieta).
The input layers for each model can be found in the [default config file](https://github.com/usgs/groundfailure/tree/main/defaultconfigfiles/models) for each model,
as described above. Additional layers used in the tests were extracted from the global input layers defined below:

* ne_10m_ocean: [Natural Earth (2016) Ocean polygon](http://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-ocean) last accessed 17 Nov 2017
* cities1000.txt: Global city information from [GeoNames](http://geonames.org) last accessed 2 Sept 2015.
* md30_gmted_gmt.grd: [GMTED2010 Terrain Elevation model](https://topotools.cr.usgs.gov/gmted_viewer)
* gmted_global_hillshade.grd: Hillshade created from md30_gmted_gmt.grd
* lspop2016_lp.flt: [LandScan (2016)™ High Resolution global Population Data Set](https://landscan.ornl.gov/)

## References

