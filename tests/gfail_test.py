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
from gfail.gfailrun import getGridURL, isURL

# where is this script?
homedir = os.path.dirname(os.path.abspath(__file__))
upone = os.path.join(homedir, os.pardir)
datadir = os.path.abspath(os.path.join(homedir, 'data'))

changetarget = False  # Change to True if need to recompute target data


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

    trimfile = '%s/loma_prieta/mapping_inputs/ne_10m_ocean/ne_10m_ocean.shp' % datadir

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

        # Clear paths
        rc, so, se = get_command_output('gfail -reset')
        # Modify paths
        pathcmd = pathcmd.replace('[TMPOUT]', p)
        rc1, so1, se1 = get_command_output(pathcmd)

        with open(default_file, "a") as f:
            f.write("popfile = %s"
                    % os.path.join(datadir, 'loma_prieta/lspop2016_lp.flt'))

        # List paths
        rc3, so3, se3 = get_command_output('gfail --list-default-paths')

        # Run model with bounds
        runcmd = "gfail %s/test_conf %s -b 'zoom, pga, 2' --hdf5 -tr  %s -ext" % \
            (datadir, shakegrid, trimfile)
        rc4, so4, se4 = get_command_output(runcmd)

        # Run model
        runcmd = "gfail %s/test_conf %s --gis -pn -pi -pd --hdf5 -ext" % (datadir,
                                                                          shakegrid)
        rc2, so2, se2 = get_command_output(runcmd)

        # Read in the testing data
        test_file = os.path.join(p,
                                 '19891018000415_zhu_2015_model.tif')
        test_grid = GDALGrid.load(test_file)
        test_data = test_grid.getData()

        # Read in target file
        target_file = os.path.join(datadir, 'loma_prieta', 'targets',
                                   '19891018000415_zhu_2015_model.tif')

        if changetarget:
            # To change target data:
            test_grid.save(test_file)
            cmd = 'gdal_translate -a_srs EPSG:4326 -of GTiff %s %s' % (
                test_file, target_file)
            rc, so, se = get_command_output(cmd)

        target_grid = GDALGrid.load(target_file)
        target_data = target_grid.getData()

    except Exception as e:
        print(e)

    # Put defaults back
    if os.path.exists(default_file+'_bak'):
        shutil.copy(default_file+'_bak', default_file)

    # Remove backup and tempfile
    if os.path.exists(default_file+'_bak'):
        os.remove(default_file+'_bak')
    shutil.rmtree(p)

    # Test that everything ran
    np.testing.assert_equal(True, rc, 'gfail reset failed')
    np.testing.assert_equal(True, rc1, 'gfail path modification failed')
    np.testing.assert_equal(True, rc2, se2.decode())
    np.testing.assert_equal(True, rc3, 'gfail list-default-paths failed')
    np.testing.assert_equal(True, rc4, se4.decode())

    # Then do test of values
    np.testing.assert_allclose(target_data, test_data, rtol=1e-3)


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
        shutil.copy(default_file, default_file + '_bak')

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

        # Clear paths
        rc, so, se = get_command_output('gfail -reset')
        np.testing.assert_equal(True, rc, 'gfail reset failed')

        # Modify paths
        pathcmd = pathcmd.replace('[TMPOUT]', p)
        rc, so, se = get_command_output(pathcmd)
        np.testing.assert_equal(True, rc, 'gfail path modification failed')

        with open(default_file, "a") as f:
            f.write("popfile = %s"
                    % os.path.join(datadir, 'loma_prieta/lspop2016_lp.flt'))

        # Run model
        conf = os.path.join(datadir, 'test_conf')
        runcmd = "gfail %s %s -w --hdf5 -ext" % (conf, shakegrid)
        rc, so, se = get_command_output(runcmd)
        np.testing.assert_equal(True, rc, se.decode())

        #event_dir = os.path.join(p, '19891018000415')

        # # Make png
        # cmd = 'create_png -e %s' % event_dir
        # rc1, so1, se1 = get_command_output(cmd)
        # np.testing.assert_equal(True, rc1, se1.decode())

        # # Make info
        # cmd = 'create_info -e %s' % event_dir
        # rc2, so2, se2 = get_command_output(cmd)
        # np.testing.assert_equal(True, rc2, se2.decode())

        # Make PDL directory
        #pdldir = os.path.join(p, '19891018000415')
        pdldir = p
        pdl.prepare_pdl_directory(pdldir)

        # Transfer dry run
        pdl_out = pdl.transfer(pdldir, 1, 'None', dryrun=True)
    except Exception as e:  # So that defaults are put back even if something goes wrong
        print(e)

    # Put defaults back
    if os.path.exists(default_file+'_bak'):
        shutil.copy(default_file+'_bak', default_file)

    # remove backup
    os.remove(default_file+'_bak')

    # Then do test
    assert pdl_out['rc'] is True

    # Remove tempfile
    shutil.rmtree(p)


def misc():
    # Note, the url used below may need to be updated if the location of this external resource changes
    gridurl = 'https://earthquake.usgs.gov/archive/product/shakemap/atlas19891018000415/atlas/1508184913111/download/grid.xml'
    test1 = isURL(gridurl)
    test2 = getGridURL(gridurl)
    np.testing.assert_equal(True, test1, 'isURL did not work')
    np.testing.assert_equal(True, os.path.exists(test2),
                            'URL was not downloaded')


if __name__ == "__main__":
    td1 = tempfile.TemporaryDirectory()
    test_zhu2015(td1)
    td2 = tempfile.TemporaryDirectory()
    test_zhu2015_web(td2)
    misc()
    print('gfail tests passed')
