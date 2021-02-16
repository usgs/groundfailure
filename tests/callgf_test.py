#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import numpy as np

# third party
import tempfile
from impactutils.io.cmd import get_command_output

# where is this script?
homedir = os.path.dirname(os.path.abspath(__file__))
upone = os.path.join(homedir, os.pardir, 'defaultconfigfiles')
datadir = os.path.abspath(os.path.join(homedir, 'data'))


def test_callgf(tmpdir):
    modinputs = os.path.join(datadir, 'ci39462536', 'model_inputs')
    pathcmd = """
        gfail --set-default-paths \
        -d %s \
        -o [TMPOUT] \
        -c %s \
        -tr %s \
        -pdl %s \
        -log [TMPOUT] \
        -db %s \
        -pf %s
    """ % (modinputs,
           os.path.join(upone, 'models'),
           os.path.join(modinputs, 'ne_10m_ocean.shp'),
           os.path.join(modinputs, 'blank.ini'),
           os.path.join('[TMPOUT]', 'test.db'),
           os.path.join(modinputs, 'lspop2016.flt'))

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
        np.testing.assert_equal(True, rc1, se1.decode())
        np.testing.assert_equal(True, 'New default paths set' in str(so1), so1.decode())

        # Run event that should fail because of magnitude
        runcmd = "callgf -e ci39473968 --dry-run" 
        rc3, so3, se3 = get_command_output(runcmd)
        np.testing.assert_equal(True, 'Magnitude check failed' in str(so3), se3.decode())

        # Run event that should fail because of magnitude, force to run
        runcmd = "callgf -e ci39473968 --dry-run -f" 
        rc5, so5, se5 = get_command_output(runcmd)
        np.testing.assert_equal(True,
                                'Completed gfail run of ci39473968' in str(so5),
                                se5.decode())

        # Force event to run and use certain version and source
        runcmd = "callgf -e ci39473968 --dry-run -v 4 -s cgs -f" 
        rc2, so2, se2 = get_command_output(runcmd)
        np.testing.assert_equal(True,
                                'Completed gfail run of ci39473968' in str(so5),
                                se5.decode())

        # Add and then remove stopfile
        runcmd = "callgf --stop ci39473968" 
        rc6, so6, se6 = get_command_output(runcmd)
        np.testing.assert_equal(True,
                                'Stopfile added' in str(so6),
                                se6.decode())
        # Add and then remove stopfile
        runcmd = "callgf --unstop ci39473968" 
        rc7, so7, se7 = get_command_output(runcmd)
        np.testing.assert_equal(True,
                                'Stopfile removed' in str(so7),
                                se7.decode())

        # Run model with url
        url = 'https://earthquake.usgs.gov/archive/product/shakemap/ci39473968/ci/1591852561898/download/grid.xml'
        runcmd = "callgf -e %s --dry-run" % url
        rc4, so4, se4 = get_command_output(runcmd)
        np.testing.assert_equal(True, 'Completed gfail run of ci39473968' in str(so4),
                                so4.decode())
        
        #TODO add test to simulate pdl triggering

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
