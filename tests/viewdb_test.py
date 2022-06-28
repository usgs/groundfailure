#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test for viewdatabase and related functions

Created on Thu Mar 19 13:35:28 2020

@author: kallstadt
"""

import os
import shutil
import numpy as np
import tempfile
import matplotlib.pyplot as plt

# third party
from impactutils.io.cmd import get_command_output


# where is this script?
homedir = os.path.dirname(os.path.abspath(__file__))
upone = os.path.join(homedir, os.pardir)
datadir = os.path.abspath(os.path.join(homedir, 'data'))


def test_viewdb(tmpdir):
    # Make a copy of current defaults
    default_file = os.path.join(os.path.expanduser("~"), ".gfail_defaults")
    if os.path.exists(default_file):
        shutil.copy(default_file, default_file+'_bak')

    try:
        # Point to test database file
        pathcmd = 'gfail --set-default-paths -db %s' % \
            os.path.join(datadir, 'testevents.db')
        rc, so, se = get_command_output(pathcmd)

        # Run some examples
        runcmd1 = "viewdb -r -c -n 10 --all -s"
        rc1, so1, se1 = get_command_output(runcmd1)

        runcmd2 = "viewdb -r -c -n 10 -p --LShazmin yellow"
        rc2, so2, se2 = get_command_output(runcmd2)

        runcmd3 = "viewdb -e us1000h3p4 us1000h5el -s" # --timeplots"
        rc3, so3, se3 = get_command_output(runcmd3)

        runcmd4 = "viewdb -p -c --minmag 7.5 --color"
        rc4, so4, se4 = get_command_output(runcmd4)
        breakpoint()
        runcmd5 = "viewdb -s --summaryplotfile %s" % \
            os.path.join(tmpdir.name, 'test.png')
        rc5, so5, se5 = get_command_output(runcmd5)

        runcmd6 = "viewdb -s --csvfile %s" % \
            os.path.join(tmpdir.name, 'test.csv')
        rc6, so6, se6 = get_command_output(runcmd6)

        # Test that everything ran
        np.testing.assert_equal(True, rc1, '%s did not run successfully' %
                                runcmd1)
        np.testing.assert_equal(True, rc2, '%s did not run successfully' %
                                runcmd2)
        np.testing.assert_equal(True, rc3, '%s did not run successfully' %
                                runcmd3)
        np.testing.assert_equal(True, rc4, '%s did not run successfully' %
                                runcmd4)
        np.testing.assert_equal(True, rc5, '%s did not run successfully' %
                                runcmd5)

        # Make sure figure was created
        np.testing.assert_equal(True, os.path.isfile(os.path.join(tmpdir.name,
                                'test_overall.png')))
        # Make sure csv file was created
        np.testing.assert_equal(True, os.path.isfile(os.path.join(tmpdir.name,
                                'test.csv')))

    except Exception as e:
        print(e)

    # Put defaults back
    if os.path.exists(default_file+'_bak'):
        shutil.copy(default_file+'_bak', default_file)

    # Remove backup and tempfile
    if os.path.exists(default_file+'_bak'):
        os.remove(default_file+'_bak')


if __name__ == "__main__":
    td1 = tempfile.TemporaryDirectory()
    test_viewdb(td1)
    # print('viewdb tests passed')
