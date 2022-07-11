
[![codecov](https://codecov.io/gh/usgs/groundfailure/branch/main/graph/badge.svg)](https://codecov.io/gh/usgs/groundfailure)

| OS   | Python | Status |
| :--- | :----   | :--- |
| macOS-latest | 3.8  | [![Build Status](https://dev.azure.com/GHSC-ESI/USGS-groundfailure/_apis/build/status/usgs.groundfailure?branchName=main&jobName=groundfailure&configuration=groundfailure%20MacOS_py38)](https://dev.azure.com/GHSC-ESI/USGS-groundfailure/_build/latest?definitionId=7&branchName=main) |
| macOS-latest | 3.9  | [![Build Status](https://dev.azure.com/GHSC-ESI/USGS-groundfailure/_apis/build/status/usgs.groundfailure?branchName=main&jobName=groundfailure&configuration=groundfailure%20MacOS_py39)](https://dev.azure.com/GHSC-ESI/USGS-groundfailure/_build/latest?definitionId=7&branchName=main) |
| ubuntu-latest | 3.8  | [![Build Status](https://dev.azure.com/GHSC-ESI/USGS-groundfailure/_apis/build/status/usgs.groundfailure?branchName=main&jobName=groundfailure&configuration=groundfailure%20Linux_py38)](https://dev.azure.com/GHSC-ESI/USGS-groundfailure/_build/latest?definitionId=7&branchName=main) |
| ubuntu-latest | 3.9  | [![Build Status](https://dev.azure.com/GHSC-ESI/USGS-groundfailure/_apis/build/status/usgs.groundfailure?branchName=main&jobName=groundfailure&configuration=groundfailure%20Linux_py39)](https://dev.azure.com/GHSC-ESI/USGS-groundfailure/_build/latest?definitionId=7&branchName=main) |

# groundfailure
Allstadt, K.E., Thompson, E.M., Hearne, M., Biegel, K., 2018, groundfailure v1.0: U.S. Geological Survey Software Release, https://doi.org/10.5066/P91G4NS4.

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

Instead of the command line options, users can also run parts of the codes manually and interactively, for example with IPython or Jupyter notebooks. See the [notebook on manual runs](http://localhost:8889/notebooks/notebooks/Run_groundfailure_manually.ipynb) for examples on how to do so.

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

5. (optional) To excercise the code and run the unit tests, you can run this command
   in the base of the repository:
   ```
   pytest .
   ```

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

The `gfail` command lets you set paths to all the input files on each run, or you can set some default paths that will be used anytime that the paths are not explicitly specified. Most path options are used only for operational near-real-time use, only the data_path, configfilepath, and output_filepath and optionally the coastlines trimming file (trimfile) are needed for local runs. You can set these defaults either using `gfail`'s `--set-default-paths` flag as demonstrated in the [commandline notebook](https://github.com/usgs/groundfailure/blob/main/notebooks/Run_groundfailure_commandline.ipynb). Alternatively, the user can manually create a text file called .gfail_defaults following the format example below and put that in their home directory, the outcome will be the same:

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

#### After setting default paths, gfail can be run from the command line like this:

```sh
gfail modelconfig.ini shakefile.xml --gis --kmz
```

* the --gis flag outputs geotiff files of the models
* the --kmz flag creates stylized kmz files of the models
* type gfail -h to see all options

### Model config file format

The model config file format is a modified version of the "INI" format.  It is described in detail [here](http://configobj.readthedocs.org/en/latest/configobj.html#config-files). See any of the [default configuration files](https://github.com/usgs/groundfailure/tree/main/defaultconfigfiles/models) included with this software release for well-commented examples of the model config file format.

**Useful notes about config file formats** 
* References and other strings with commas within them need to be enclosed in
  quotes. Example:
  * 'Verdin, D.W., Godt, J., Funk, C., Pedreros, D., Worstell, B. and Verdin,
    J., 2007, Development of a global slope dataset for estimation of landslide
    occurrence resulting from earthquakes: U.S. Geological Survey Open-File
    Report 2007–1188, 25p.'
* Arrays should be not be enclosed in brackets and should be comma separated.
  Example:
  * model = 0, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.
* Files, filepaths, and folders being used for input layers should all be labeled as 'file' in the config file regardless of actual structure.
* If file paths in the config file are relative to a base data folder on the user's system, and you are running models manually and/or not setting default paths in gfail, you will need to run correct_config_filepaths() after reading in the config file. See the [notebook on manual runs](https://github.com/usgs/groundfailure/blob/main/notebooks/Run_groundfailure_manually.ipynb) for an example.

## Structure for Model Output

If run manually within python (e.g., [see manual notebook](https://github.com/usgs/groundfailure/blob/main/notebooks/Run_groundfailure_manually.ipynb)), each model outputs a single ordered dictionary, which has keys that correspond to the names of the input and output layers from the model (e.g., 'pga', 'slope', 'friction', 'cti1', 'model'). Keys for the input layers will only be present if saveinputs=True was set when calling the model. The 'model' key is the model output, all other names of input layers come from the names used in the model configuration file.

Each layer in the dictionary is itself a dictionary with some of the following fields:
 - *description* A dictionary with the fields:

   * *name* Short name, suitable for use as a plot title if necessary.
   * *longref* Full citation, USGS format as described here: http://internal.usgs.gov/publishing/sta/sta28.pdf
   * *units* Physical units for input data layers, and one of the following for output "probability" layers:

     * *index* Relative (low to high) index of occurrence in a given cell (not necessarily bounded).
     * *probability* Probability of event (landslide,liquefaction) of a given size occurring in a given cell (0 to 1).
     * *coverage* Fractional coverage of groundfailure in a given cell (0 to 1).
     * *displacement* Distance material will move from or in given cell (unbounded).

   * *parameters* (Model output only) A dictionary of key/value pairs, where the values must be either numbers or strings.

 - *type* Indicates whether this grid contains input data or output from a model.

 - *label* What will be written next to the colorbar for the data layer.

 - *grid* Input data or model output, in the form of a MapIO Grid2D object. 

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

Allstadt, K.E., Thompson, E.M., Jibson, R.W., Wald, D.J., Hearne, M., Hunter, E.J., Fee, J., Schovanec, H., Slosky, D., Haynie, K. L., 2021, The USGS ground failure product: near-real-time estimates of earthquake-triggered landslides and liquefaction, Earthquake Spectra, 38, 5-36, https://doi.org/10.1177%2F87552930211032685.
