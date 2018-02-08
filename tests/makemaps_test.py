#!/usr/bin/env python3
"""
The idea here is to use the very fast (small/fake) test data
to try to hit all (most) of the lines of code in makemaps, which
would take just way to long to do with real data.
"""
import os.path
import os
from configobj import ConfigObj
import gfail.logisticmodel as LM
from mapio.geodict import GeoDict
from gfail.conf import correct_config_filepaths
import gfail.makemaps as makemaps
import shutil

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
datadir = os.path.abspath(os.path.join(homedir, 'data'))
upone = os.path.join(homedir, os.pardir)

configfile = os.path.join(datadir, 'testconfig_logimodel.ini')
config = ConfigObj(configfile)
# Test path correction (from conf.py)
config = correct_config_filepaths(datadir, config)

mapconfigfile = os.path.join(upone, 'defaultconfigfiles', 'mapconfig.ini')
mapconfig = ConfigObj(mapconfigfile)


cmodel = config['test_model']
layers = []

shakefile = os.path.join(datadir, 'test_shakegrid.xml')
uncertfile = os.path.join(datadir, 'test_uncert.xml')
cofile = os.path.join(datadir, 'test_cohesion.bil')
slopefile = os.path.join(datadir, 'test_slope.bil')
vs30file = os.path.join(datadir, 'test_vs30.bil')
ctifile = os.path.join(datadir, 'test_cti1.bil')
precipfolder = os.path.join(datadir, 'test_precip')

fakegeodict = GeoDict({'xmin': 0.5, 'xmax': 1.5,
                       'ymin': 0.5, 'ymax': 1.5,
                       'dx': 1.0, 'dy': 1.0,
                       'ny': 2, 'nx': 2})

modelLQ = {
    'TestModelLQ': {
        'description': 'This is a test liquefaction model',
        'gfetype': 'liquefaction',
        'baselayer': 'vs30',
        'slopemin': 0.,
        'slopemax': 5.,
        'layers': {
            'vs30': {
                'file': vs30file,
                'units': 'm/s',
                'longref': 'more words',
                'shortref': 'words'
            },
            'cti1': {
                'file': ctifile,
                'units': 'unitless',
                'longref': 'more words',
                'shortref': 'words'
            }
        },
        'interpolations': {
            'vs30': 'nearest',
            'cti1': 'linear'
        },
        'terms': {
            'b1': 'log((pga/100.0)*(power(MW,2.)))',
            'b2': 'cti1',
            'b3': 'log(vs30)'
        },
        'coefficients': {
            'b0': 15.,
            'b1': 2.,
            'b2': 0.3,
            'b3': -4.
        }
    }
}


def test_parseMapConfig():
    config = mapconfig
    # fileext is None
    makemaps.parseMapConfig(config)
    # 'ocean' in config, oceanref = config['ocean']['shortref'] fails
    del config['ocean']['shortref']
    makemaps.parseMapConfig(config)
    # 'cities' in config, cityref = config['cities']['shortref'] fails
    del config['cities']['shortref']
    makemaps.parseMapConfig(config)
    # Give an invalid city file
    config = mapconfig
    config['cities']['file'] = os.path.join(datadir,
                                            'loma_prieta/mapping_inputs/gmted_global_hillshade.grd.aux.xml')
    makemaps.parseMapConfig(config)
    # 'alpha' in config['colors']
    config = mapconfig
    config['colors']['alpha'] = 0.5
    makemaps.parseMapConfig(config)


def test_parseConfigLayers():
    lq = LM.LogisticModel(shakefile, modelLQ, saveinputs=True)
    maplayers = lq.calculate()
    makemaps.parseConfigLayers(maplayers, config, keys=None)
    # no lims
    del config['test_model']['display_options']['lims']
    makemaps.parseConfigLayers(maplayers, config, keys=None)
    # no colors
    del config['test_model']['display_options']['colors']
    makemaps.parseConfigLayers(maplayers, config, keys=None)
    # no logscale
    del config['test_model']['display_options']['logscale']
    makemaps.parseConfigLayers(maplayers, config, keys=None)
    # no maskthresholds
    del config['test_model']['display_options']['maskthresholds']
    makemaps.parseConfigLayers(maplayers, config, keys=None)
    # plotorder[0] != 'model' --- causes error
#    tmp = collections.OrderedDict()
#    tmp['vs30'] = maplayers['vs30']
#    tmp['model'] = maplayers['model']
#    makemaps.parseConfigLayers(tmp, config, keys=None)


def test_modelMap(tempdir):
    lq = LM.LogisticModel(shakefile, modelLQ, saveinputs=True)
    maplayers = lq.calculate()
    # suptitle is None
    makemaps.modelMap(maplayers, shakefile, suptitle=None,
                      savepdf=False, savepng=False,
                      outputdir=tempdir)
    # shakefile is None
    makemaps.modelMap(maplayers, suptitle=None,
                      savepdf=False, savepng=False,
                      outputdir=tempdir)
    # scaletype == 'binned'
    makemaps.modelMap(maplayers, scaletype='binned',
                      savepdf=False, savepng=False,
                      outputdir=tempdir)
    # scaletype == 'binned' and logscale=!False
    makemaps.modelMap(maplayers, scaletype='binned',
                      logscale=[False, False, True, True],
                      savepdf=False, savepng=False,
                      outputdir=tempdir)
    # logscale=!False
    makemaps.modelMap(maplayers, logscale=[False, False, True, True],
                      savepdf=False, savepng=False,
                      outputdir=tempdir)


def test_zoom(tempdir):

    # boundaries == 'zoom'
    shakefile = os.path.join(datadir, 'loma_prieta', 'grid.xml')
    conf_file = os.path.join(upone, 'defaultconfigfiles', 'models',
                             'zhu_2015.ini')
    conf = ConfigObj(conf_file)
    data_path = os.path.join(datadir, 'loma_prieta', 'model_inputs')
    conf = correct_config_filepaths(data_path, conf)

    lq = LM.LogisticModel(shakefile, conf, saveinputs=True)
    maplayers = lq.calculate()

    makemaps.modelMap(maplayers, boundaries='zoom', zthresh=0.3,
                      savepdf=False, savepng=False, outputdir=tempdir)

    # bounaries dictionary
    bounds = {'xmin': -122.54, 'xmax': -120.36,
              'ymin': 36.1, 'ymax': 37.0}
    makemaps.modelMap(maplayers, boundaries=bounds,
                      savepdf=False, savepng=False)


if __name__ == "__main__":
    #td1 = tempfile.TemporaryDirectory()
    td1 = os.path.join(datadir, 'temporary1')
    test_parseMapConfig()
    test_parseConfigLayers()
    test_modelMap(td1)
    test_zoom(td1)
    # remove tempdir
    shutil.rmtree(td1)
