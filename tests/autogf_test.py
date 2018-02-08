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
upone = os.path.join(homedir, os.pardir)
datadir = os.path.abspath(os.path.join(homedir, 'data'))


def test_autogf(tmpdir):
    """
    Runs with current feed, sometimes will run an event, sometimes not
    """
    config = os.path.join(datadir, 'test_autoconfig.ini')
    gconfig = os.path.join(upone, 'defaultconfigfiles', 'models', 'jessee_2017.ini')

    # Set gfail defaults
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
        #p = os.path.join(str(tmpdir.name), "sub")
        p = os.path.join(str(tmpdir), "sub")
        if not os.path.exists(p):
            os.makedirs(p)

        # Modify gfail defaults
        pathcmd = pathcmd.replace('[TMPOUT]', p)
        rc, so, se = get_command_output(pathcmd)

        # Set command for autogf
        agfcmd = """
            autogf -c %s \
            -d test.db \
            -a %s \
            -t -w
        """ % (gconfig, config)

        rc, so, se = get_command_output(agfcmd)
        temp = so.decode().split('\n')[-2]
    except Exception as e:  # To make sure defaults are replaced
        print(e)

    # Change back gfail defaults
    if os.path.exists(default_file+'_bak'):
        shutil.copy(default_file+'_bak', default_file)
    # Remove backup and tempfile
    os.remove(default_file+'_bak')
    shutil.rmtree(p)

    # Then run tests
    np.testing.assert_equal(True, rc, 'autogf run failed')
    np.testing.assert_string_equal(temp, 'Test successful, cleaning up files')
    np.testing.assert_equal(True, os.path.exists('test.db'))
    os.remove('test.db')


if __name__ == "__main__":
    td1 = tempfile.TemporaryDirectory()
    test_autogf(td1)
    print('autogf tests passed')
    shutil.rmtree(td1)
