 #!/usr/bin/env python

import os.path
import os
import shutil
import glob

#third party
from impactutils.io.cmd import get_command_output

#hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
upone = os.path.dirname(homedir)
datadir = os.path.abspath(os.path.join(homedir, 'data'))

shakefile = os.path.join(datadir, 'test_shakegrid.xml')
configfile = os.path.join(datadir, 'testconfig_classic.ini')
uncertfile = os.path.join(datadir, 'test_uncert.xml')
gfaildefaults = os.path.join(datadir, 'test_gfail_defaults')

# Create gfail defaults file
text = {'config_filepath': os.path.join(upone, 'defaultconfigfiles/models'),
        'output_filepath': os.path.join(datadir, 'test_outputs'),
        'data_path': datadir,
        'mapconfig': os.path.join(upone, 'defaultconfigfiles/mapconfig.ini')}
if not os.path.exists(os.path.join(datadir, 'test_outputs')):
    os.makedirs(os.path.join(datadir, 'test_outputs'))
with open(gfaildefaults, 'w') as f:
    for key, tx in text.items():
        f.write('%s = %s\n' % (key, tx))
        print('%s = %s\n' % (key, tx))

# temporarily switch out .gfail_defaults file
defaultfilepath = os.path.join(os.path.expanduser('~'), '.gfail_defaults')
tempname = os.path.join(os.path.expanduser('~'), 'gfail_defaults_temp')
if os.path.exists(defaultfilepath) and not os.path.exists(tempname):
    os.rename(defaultfilepath, tempname)
shutil.copy(gfaildefaults, defaultfilepath)


def tes_list_default_paths():
    cmd = 'gfail --list-default-paths'
    retcode, stdout, stderr = get_command_output(cmd)
    temp = stdout.decode('utf-8')
    print(stderr)
    assert('Default paths currently set:' in temp), 'list_default_paths did not return expected output'
    temp2 = temp.split('Default paths currently set:')
    assert('/defaultconfigfiles/models' in temp2[1]), 'list_default_paths did not return expected output'
    assert('test_outputs' in temp2[1]), 'list_default_paths did not return expected output'


def tes_reset_default_paths():
    cmd = 'gfail --reset-default-paths'
    retcode, stdout, stderr = get_command_output(cmd)
    temp = stdout.decode('utf-8')
    assert('Default paths cleared' in temp), 'reset_default_paths did not return expected output'
    retcode, stdout, stderr = get_command_output(cmd)


def tes_set_default_paths():
    cmd = ('gfail --set-default-paths -o %s -d %s -c %s -m %s' % (text['output_filepath'],
           text['data_path'], text['config_filepath'], text['mapconfig']))
    retcode, stdout, stderr = get_command_output(cmd)
    temp = stdout.decode('utf-8')
    assert('Default paths currently set:' in temp), 'set_default_paths did not return expected output'
    temp2 = temp.split('Default paths currently set:')
    assert('defaultconfigfiles/models' in temp2[1]), 'set_default_paths did not return expected output'
    assert('test_outputs' in temp2[1]), 'set_default_paths did not return expected output'
    assert('does not exist' not in temp2[1]), 'set_default_paths could not find one of the paths expected'


def tes_main():
    cmd = 'gfail -i -pd -pn -pi %s %s' % (configfile, shakefile)
    retcode, stdout, stderr = get_command_output(cmd)
    temp = stdout.decode('utf-8')
    print(temp)
    # SEE IF GIVES EXPECTED OUTPUTS


def tes_getGridURL():
    pass


def tes_isURL():
    pass

try:
    tes_list_default_paths()
    tes_reset_default_paths()
    tes_set_default_paths()
    tes_getGridURL()
    tes_isURL()
    tes_main()
    print('gfail tests passed')
finally:
    # Put back defaults if they were there
    os.remove(defaultfilepath)
    if os.path.exists(tempname):
        os.rename(tempname, defaultfilepath)
    if os.path.exists(gfaildefaults):
        os.remove(gfaildefaults)
    # Erase any outputs that were made
    filenames = glob.glob(os.path.join(datadir, 'test_outputs', '*'))
    if filenames:
        for filen in filenames:
            os.remove(filen)
    os.rmdir(os.path.join(datadir, 'test_outputs'))
