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
    shutil.copy(default_file, default_file+'_bak')

    p = os.path.join(str(tmpdir), "sub")
    if not os.path.exists(p):
        os.makedirs(p)

    # Modify paths
    pathcmd = pathcmd.replace('[TMPOUT]', p)
    rc, so, se = get_command_output(pathcmd)

    # Run model
    runcmd = "gfail zhu_2015.ini %s --gis -pn -pi -pd" % (shakegrid)
    rc, so, se = get_command_output(runcmd)

    # Read in target file
    target_file = os.path.join(datadir, 'loma_prieta', 'targets',
                               '19891018000415_zhu_2015_model.bil')
    target_grid = GDALGrid.load(target_file)
    target_data = target_grid.getData()

    # Read in the testing data
    test_file = os.path.join(p, '19891018000415',
                             '19891018000415_zhu_2015_model.bil')
    test_grid = GDALGrid.load(test_file)
    test_data = test_grid.getData()

    # Put defaults back
    shutil.copy(default_file+'_bak', default_file)

    # Remove backup and tempfile
    os.remove(default_file+'_bak')
    shutil.rmtree(p)

    # Assert
    np.testing.assert_allclose(target_data, test_data)


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
    shutil.copy(default_file, default_file+'_bak')

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
    runcmd = "gfail %s %s --gis -w --alert" % (conf, shakegrid)
    rc, so, se = get_command_output(runcmd)

    # Make PDL directory
    pdl.prepare_pdl_directory('19891018000415')

    # Transfer dry run
    transfer_cmd = pdl.transfer('19891018000415', 'None', dryrun=True)
    assert '--property-alertLQ=yellow' in transfer_cmd
    assert '--property-alertLS=yellow' in transfer_cmd
    assert '--type=groundfailure' in transfer_cmd
    assert '--property-title=Earthquake-Induced Groundfailure' in transfer_cmd
    assert '--eventsourcecode=19891018000415' in transfer_cmd

    # Put defaults back
    shutil.copy(default_file+'_bak', default_file)

    # Remove backup and tempfile
    os.remove(default_file+'_bak')
    shutil.rmtree(p)


if __name__ == "__main__":
    td1 = tempfile.TemporaryDirectory()
    test_zhu2015(td1)
    td2 = tempfile.TemporaryDirectory()
    test_zhu2015_web(td2)
