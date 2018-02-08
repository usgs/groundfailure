#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import numpy as np

# third party
import tempfile
from mapio.gdal import GDALGrid
from impactutils.io.cmd import get_command_output
import gfail.pdl as pdl

# where is this script?
homedir = os.path.dirname(os.path.abspath(__file__))
upone = os.path.join(homedir, os.pardir)
datadir = os.path.abspath(os.path.join(homedir, 'data'))


def test_zhu2015(tmpdir):
    shakegrid = os.path.join(datadir, 'loma_prieta', 'grid.xml')
    pathcmd = """
        gfail --set-default-paths \
        -d %s/loma_prieta/model_inputs \
        -o [TMPOUT] \
        -c %s/defaultconfigfiles/models \
        -m %s/defaultconfigfiles/mapconfig.ini \
        -md %s/loma_prieta/mapping_inputs
    """ % (datadir, upone, upone, datadir)

    # Make a copy of current defaults
    default_file = os.path.join(os.path.expanduser("~"), ".gfail_defaults")
    if os.path.exists(default_file):
        shutil.copy(default_file, default_file+'_bak')

    try:
        try:
            p = os.path.join(str(tmpdir.name), "sub")
        except:
            p = os.path.join(str(tmpdir), "sub")
        if not os.path.exists(p):
            os.makedirs(p)

        # Modify paths
        pathcmd = pathcmd.replace('[TMPOUT]', p)
        rc, so, se = get_command_output(pathcmd)

        # Run model
        runcmd = "gfail zhu_2015.ini %s --gis -pn -pi -pd" % (shakegrid)
        rc, so, se = get_command_output(runcmd)

        # Read in the testing data
        test_file = os.path.join(p, '19891018000415',
                                 '19891018000415_zhu_2015_model.tif')
        test_grid = GDALGrid.load(test_file)
        test_data = test_grid.getData()

        # Read in target file
        target_file = os.path.join(datadir, 'loma_prieta', 'targets',
                                   '19891018000415_zhu_2015_model.tif')
    #    # To change target data:
        #test_grid.save(test_file)
        #cmd = 'gdal_translate -a_srs EPSG:4326 -of GTiff %s %s' % (test_file, target_file)
        #rc, so, se = get_command_output(cmd)

        target_grid = GDALGrid.load(target_file)
        target_data = target_grid.getData()

    except Exception as e:  # So that the defaults are always put back if something above fails
        print(e)

    # Put defaults back
    if os.path.exists(default_file+'_bak'):
        shutil.copy(default_file+'_bak', default_file)

    # Convert all nans to zeros in both sets
    target_data[np.isnan(target_data)] = 0.0
    test_data[np.isnan(test_data)] = 0.0

    # Then do test
    np.testing.assert_allclose(target_data, test_data)

    # Remove backup and tempfile
    os.remove(default_file+'_bak')
    shutil.rmtree(p)


def test_zhu2015_web(tmpdir):
    shakegrid = os.path.join(datadir, 'loma_prieta', 'grid.xml')
    pathcmd = """
        gfail --set-default-paths \
        -d %s/loma_prieta/model_inputs \
        -o [TMPOUT] \
        -c %s/defaultconfigfiles/models \
        -m %s/defaultconfigfiles/mapconfig.ini \
        -md %s/loma_prieta/mapping_inputs
    """ % (datadir, upone, upone, datadir)

    # Make a copy of current defaults
    default_file = os.path.join(os.path.expanduser("~"), ".gfail_defaults")
    if os.path.exists(default_file):
        shutil.copy(default_file, default_file+'_bak')

    try:
        try:
            p = os.path.join(str(tmpdir.name), "sub")
        except:
            p = os.path.join(str(tmpdir), "sub")
        if not os.path.exists(p):
            os.makedirs(p)
        else:
            shutil.rmtree(p)
            os.makedirs(p)

        # Modify paths
        pathcmd = pathcmd.replace('[TMPOUT]', p)
        rc, so, se = get_command_output(pathcmd)

        # Run model
        conf = os.path.join(datadir, 'test_conf')
        runcmd = "gfail %s %s -w --alert" % (conf, shakegrid)
        rc, so, se = get_command_output(runcmd)

        # Make PDL directory
        pdldir = os.path.join(p, '19891018000415')
        pdl.prepare_pdl_directory(pdldir)

        # Transfer dry run
        transfer_cmd = pdl.transfer(pdldir, 'None', dryrun=True)
    except Exception as e:  # So that defaults are put back even if something goes wrong
        print(e)

    # Put defaults back
    if os.path.exists(default_file+'_bak'):
        shutil.copy(default_file+'_bak', default_file)

    # Remove backup and tempfile
    os.remove(default_file+'_bak')
    shutil.rmtree(p)

    # Then do test
    assert '--property-alertLQ=green' in transfer_cmd
    assert '--property-alertLS=yellow' in transfer_cmd
    assert '--type=groundfailure' in transfer_cmd
    assert '--property-title=Earthquake-Induced Ground Failure' in transfer_cmd
    assert '--eventsourcecode=19891018000415' in transfer_cmd


if __name__ == "__main__":
    td1 = tempfile.TemporaryDirectory()
    test_zhu2015(td1)
    td2 = tempfile.TemporaryDirectory()
    test_zhu2015_web(td2)
