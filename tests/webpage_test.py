#!/usr/bin/env python3
"""
The idea here is to use the very fast (small/fake) test data
to try to hit all (most) of the lines of code in makemaps, which
would take just way to long to do with real data.
"""
import os.path
import os
from configobj import ConfigObj
from gfail import Zhu2015Model
from mapio.geodict import GeoDict  # TODO replace with read
import tempfile
from gfail.webpage import create_kmz

# import shutil
import gfail.utilities as utilities


homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
datadir = os.path.abspath(os.path.join(homedir, "data"))
upone = os.path.join(homedir, os.pardir)

configfile = os.path.join(upone, "defaultconfigfiles", "models", "zhu_2015.ini")
config = ConfigObj(configfile)

# # Test path correction (from conf.py)
# config = correct_config_filepaths(datadir, config)

shakefile = os.path.join(datadir, "test_shakegrid.xml")
uncertfile = os.path.join(datadir, "test_uncert.xml")
cofile = os.path.join(datadir, "test_cohesion.bil")
slopefile = os.path.join(datadir, "test_slope.bil")
vs30file = os.path.join(datadir, "test_vs30.bil")
ctifile = os.path.join(datadir, "test_cti1.bil")
precipfolder = os.path.join(datadir, "test_precip")

# Replace files
config["zhu_2015"]["layers"]["vs30"]["file"] = vs30file
config["zhu_2015"]["layers"]["cti"]["file"] = ctifile


fakegeodict = GeoDict(
    {
        "xmin": 0.5,
        "xmax": 1.5,
        "ymin": 0.5,
        "ymax": 1.5,
        "dx": 1.0,
        "dy": 1.0,
        "ny": 2,
        "nx": 2,
    }
)


# modelLS = dict(jessee_2018=modelLQ['zhu_2015'])
# modelLS['jessee_2018']['description'] = 'This is a test landslide model'
# modelLS['jessee_2018']['gfetype'] = 'landslide'


def test_parseConfigLayers():
    lq = Zhu2015Model(
        shakefile,
        config["zhu_2015"],
        uncertfile=uncertfile,
        bounds=None,
        trimfile=None,
        saveinputs=True,
    )
    maplayers = lq.calculate()
    utilities.parseConfigLayers(maplayers, config, keys=None)
    # no lims
    # del config['zhu_2015']['display_options']['lims']
    # utilities.parseConfigLayers(maplayers, config, keys=None)
    # no colors
    del config["zhu_2015"]["display_options"]["colors"]
    utilities.parseConfigLayers(maplayers, config, keys=None)
    # no logscale
    del config["zhu_2015"]["display_options"]["logscale"]
    utilities.parseConfigLayers(maplayers, config, keys=None)
    # no maskthresholds
    del config["zhu_2015"]["display_options"]["maskthresholds"]
    utilities.parseConfigLayers(maplayers, config, keys=None)
    # plotorder[0] != 'model' --- causes error


#    tmp = collections.OrderedDict()
#    tmp['vs30'] = maplayers['vs30']
#    tmp['model'] = maplayers['model']
#    makemaps.parseConfigLayers(tmp, config, keys=None)


def test_maps():
    lq = Zhu2015Model(
        shakefile,
        config["zhu_2015"],
        uncertfile=uncertfile,
        bounds=None,
        trimfile=None,
        saveinputs=True,
    )
    maplayers = lq.calculate()

    # ls = LM.LogisticModel(shakefile, modelLS, saveinputs=False)
    # maplayers2 = ls.calculate()

    # Test create_kmz
    tempdir = tempfile.TemporaryDirectory()
    create_kmz(maplayers["model"], outfile=os.path.join(tempdir.name, "test.kmz"))
    # create_kmz(maplayers2['model'], outfile=os.path.join(tempdir.name,
    #                     'test.kmz'), mask=0.003)


if __name__ == "__main__":
    test_parseConfigLayers()
    test_maps()
