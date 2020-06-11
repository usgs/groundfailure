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

def test_callgf(tmpdir):
    #url = 'https://raw.githubusercontent.com/usgs/groundfailure/master/tests/data/loma_prieta/grid.xml'
    #shakegrid = os.path.join(datadir, 'loma_prieta', 'grid.xml')
    trimfile = '%s/loma_prieta/mapping_inputs/ne_10m_ocean/ne_10m_ocean.shp' \
               % datadir
    pathcmd = """
        gfail --set-default-paths \
        -d %s/loma_prieta/model_inputs \
        -o [TMPOUT] \
        -c %s/defaultconfigfiles/models \
        -m %s/defaultconfigfiles/mapconfig.ini \
        -tr %s \
        -pdl %s/blank.ini \
        -log [TMPOUT] \
        -db [TMPOUT]/test.db
    """ % (datadir, upone, upone, trimfile, datadir)

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

#        # Run model with url
#        runcmd = "callgf -e %s --dry-run" % url 
#        rc4, so4, se4 = get_command_output(runcmd)
#        np.testing.assert_equal(True, rc4, so4.decode())
        
        # Run event that should fail because over water
        runcmd = "callgf -e us2000hg93 --dry-run" 
        rc3, so3, se3 = get_command_output(runcmd)
        np.testing.assert_equal(True, rc3, se3.decode())
        
        # Run Loma Prieta from url
        runcmd = "callgf -e nc216859 --dry-run" 
        rc2, so2, se2 = get_command_output(runcmd)
        np.testing.assert_equal(True, rc2, se2.decode())

        # Run event that shold fail magnitude check
        runcmd = "callgf -e us6000a8nh --dry-run"
        rc, so, se = get_command_output(runcmd)
        np.testing.assert_equal(True, rc, se.decode())
        
        import pdb; pdb.set_trace()
        #TODO add test to simulate pdl triggering
        #TODO test version and source options

    except Exception as e:
        print(e)

    # Put defaults back
    if os.path.exists(default_file+'_bak'):
        shutil.copy(default_file+'_bak', default_file)

    # Remove backup and tempfile
    if os.path.exists(default_file+'_bak'):
        os.remove(default_file+'_bak')
    shutil.rmtree(p)




if __name__ == "__main__":
    td1 = tempfile.TemporaryDirectory()
    test_callgf(td1)
    #print('callgf tests passed')
