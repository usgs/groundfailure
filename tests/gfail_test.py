 #!/usr/bin/env python

import os.path
import os
import shutil

#third party
from impactutils.io.cmd import get_command_output

#hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
datadir = os.path.abspath(os.path.join(homedir, 'data'))

shakefile = os.path.join(datadir, 'test_shakegrid.xml')
configfile = os.path.join(datadir, 'testconfig_classic.ini')
uncertfile = os.path.join(datadir, 'test_uncert.xml')
gfaildefaults = os.path.join(datadir, 'test_gfail_defaults')

# Create gfail defaults file
text = {'config_filepath': os.path.join(homedir, '..', '/defaultconfigfiles/models'),
        'output_filepath': os.path.join(datadir, 'test_outputs'),
        'data_path': datadir,
        'mapconfig': os.path.join(homedir, '..', 'defaultconfigfiles/mapconfig.ini')}
f = open(gfaildefaults, 'w')
for key in text:
    f.write('%s = %s\n' % (key, text[key]))
f.close()

# temporarily switch out .gfail_defaults file
defaultfilepath = os.path.join(os.path.expanduser('~'), '.gfail_defaults')
tempname = os.path.join(os.path.expanduser('~'), 'gfail_defaults_temp')
if os.path.exists(defaultfilepath):
    os.rename(defaultfilepath, tempname)
shutil.copy(gfaildefaults, defaultfilepath)


def test_list_default_paths():
    # test all default path things
    cmd = 'gfail --list-default-paths'
    retcode, stdout, stderr = get_command_output(cmd)
    temp = stdout.decode('utf-8')
    assert('Default paths currently set:' not in temp), 'list_default_paths did not return expected output'


def test_reset_default_paths():
    pass


def test_set_default_paths():
    pass


def test_main():
    cmd = 'gfail -i -pd -pn -pi %s %s' % (configfile, shakefile)
    retcode, stdout, stderr = get_command_output(cmd)
    # SEE IF GIVES EXPECTED OUTPUTS


def test_getGridURL():
    pass


def test_isURL():
    pass

test_list_default_paths()
test_reset_default_paths()
test_set_default_paths()
test_getGridURL()
test_isURL()
test_main()

# Put back defaults if they were there
os.remove(defaultfilepath)
if os.path.exists(tempname):
    os.rename(tempname, defaultfilepath)
# Erase any outputs that were made
os.remove(os.path.join(datadir, 'test_outputs', '*'))
os.rmdir(os.path.join(datadir, 'test_outputs'))

print('gfail tests passed')
