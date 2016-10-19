#!/usr/bin/env python

import argparse
import os.path
import shutil
import sys
from distutils.dir_util import copy_tree

from impactutils.io.cmd import get_command_output

DEFAULT_TAG = '0.1'


def main(args):

    #-------------------------------------------------------------
    # where should .rst files, Makefile, _build folder be written?
    #-------------------------------------------------------------
    API_DIR = os.path.join(os.path.expanduser('~'), '__api-doc')
    shutil.rmtree(API_DIR, ignore_errors=True)

    #-------------------------------------------------------------
    # where should the temporary clone of the ground failure gh-pages repo live?
    #-------------------------------------------------------------
    CLONE_DIR = os.path.join(os.path.expanduser('~'), '__gf-doc')
    shutil.rmtree(CLONE_DIR, ignore_errors=True)

    #-------------------------------------------------------------
    # Some additional useful directories
    #-------------------------------------------------------------
    REPO_DIR = os.path.dirname(os.path.abspath(__file__))
    PACKAGE_DIR = os.path.join(REPO_DIR, 'groundfailure')

    #-------------------------------------------------------------
    # what is the package called and who are the authors
    #-------------------------------------------------------------
    PACKAGE = "GroundFailure"
    AUTHORS = 'Kate Allstadt, Eric Thompson, Mike Hearne, Katie Biegel'

    # find the make command on this system
    res, stdout, stderr = get_command_output('which make')
    if not res:
        print('Could not find the "make" command on your system. Exiting.')
        sys.exit(1)
    make_cmd = stdout.decode().strip()

    #-------------------------------------------------------------
    # clone the repository
    #-------------------------------------------------------------
    if args.post:
        sys.stderr.write('Cloning groundfailure gh-pages branch...\n')
        if os.path.isdir(CLONE_DIR):
            shutil.rmtree(CLONE_DIR)
        clonecmd = 'git clone -b gh-pages https://github.com/usgs/'\
                   'groundfailure.git %s' % CLONE_DIR
        res, stdout, stderr = get_command_output(clonecmd)
        if not res:
            raise Exception('Could not clone gh-pages branch.')

        # Delete everything in the repository (except hidden git files)
        cmd = 'rm -fr %s/*' % CLONE_DIR
        res, stdout, stderr = get_command_output(cmd)

    #-------------------------------------------------------------
    # run the api doc command; this creates the .rst files
    #-------------------------------------------------------------
    sys.stderr.write('Building groundfailure API documentation (REST)...\n')
    sphinx_cmd = 'sphinx-apidoc -o %s -f -e -l -d 12 -F -H %s -A "%s"'\
                 ' %s' % (API_DIR, PACKAGE, AUTHORS, PACKAGE_DIR)

    res, stdout, stderr = get_command_output(sphinx_cmd)

    if not res:
        raise Exception('Could not build GroundFailure API documentation'
                        ' - error "%s".' % stderr)

    #-------------------------------------------------------------
    # Edit the conf.py file to include the theme.
    #-------------------------------------------------------------
    fname = os.path.join(API_DIR, 'conf.py')
    f = open(fname, 'at')
    f.write("import os\nimport sys\n")
    f.write("sys.path.insert(0, os.path.abspath('%s'))\n" % (REPO_DIR))
    f.write("sys.path.insert(0, os.path.abspath('.'))\n")
    f.write("temp = sys.executable\n")
    f.write("EXECPATH='/'.join(temp.split('/')[:-2])\n")
    f.write("sys.path.append(os.path.join(os.path.expanduser('~'), EXECPATH, 'lib'))\n")

    #-------------------------------------
    # RTD theme
    #-------------------------------------
    f.write("import sphinx_rtd_theme\n")
    f.write("html_theme = 'sphinx_rtd_theme'\n")
    f.write("html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]\n")
    f.write("html_theme_options = {\n")
    f.write("    'collapse_navigation': False,\n")
    f.write("}\n")
    #-------------------------------------

    # Napolean extension? Supports Goggle and Numpy style docstrings, but it
    # also has some side effects such as restrictions on what sections are
    # allowed and it seems to suppress the [source] link to code; maybe this
    # is a configurable option though.
#    f.write("extensions = ['sphinx.ext.autodoc', 'sphinxcontrib.napoleon']\n")

    # This line is needed to inclue __init__ methods in documentation
    f.write("autoclass_content = 'both'\n")
    f.write("autodoc_member_order = 'bysource'\n")
    f.write("html_show_copyright = False\n")
#    f.write("extensions = extensions + [ 'sphinx.ext.autodoc', "\
#            "'sphinx.ext.napoleon', 'sphinx.ext.todo' ] \n")
    f.write("napoleon_include_special_with_doc = False\n")
    f.write("todo_include_todos = True\n")
    f.close()

    #-------------------------------------------------------------
    # Go to the api directory and build the html
    #-------------------------------------------------------------
    sys.stderr.write('Building groundfailure manual (HTML)...\n')
    os.chdir(API_DIR)
    res, stdout, stderr = get_command_output('%s html' % make_cmd)
    if not res:
        raise Exception('Could not build HTML for API documentation. - error "%s"' % stderr)
    if args.verbose:
        print(stdout.decode('utf-8'))
        print(stderr.decode('utf-8'))

    #-------------------------------------------------------------
    # Copy the generated content to the gh-pages branch we created
    # earlier
    #-------------------------------------------------------------
    htmldir = os.path.join(API_DIR, '_build', 'html')
    if not os.path.isdir(CLONE_DIR):
        os.makedirs(CLONE_DIR)
    copy_tree(htmldir, CLONE_DIR)

    if args.post:
        #-------------------------------------------------------------
        # Post to gh-pages
        #-------------------------------------------------------------

        # cd to directory above where html content was pushed
        os.chdir(CLONE_DIR)
        res, stdout, stderr = get_command_output('touch .nojekyll')
        res1, stdout, stderr1 = get_command_output('git add --all')
        res2, stdout, stderr2 = get_command_output(
            'git commit -am"Pushing to GitHub pages"')
        res3, stdout, stderr3 = get_command_output(
            'git push -u origin +gh-pages')
        if res1 + res2 + res3 < 3:
            stderr = stderr1 + stderr2 + stderr3
            print('Something bad happened when attempting to add, commit, '
                  'or push gh-pages content to GitHub - error "%s". Exiting.'
                  % stderr)
            sys.exit(1)
        print('You can inspect the GroundFailure API docs by looking '
              'here: http://usgs.github.io/groundfailure/index.html')
    else:
        if not args.clean:
            indexpage = os.path.join(CLONE_DIR, 'index.html')
            print('GroundFailure documentation index: %s' % indexpage)


if __name__ == '__main__':
    desc = 'Create API documentation for GroundFailure and optionally post '\
           'HTML output to GitHub.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-p', '--post', action='store_true', default=False,
                        help='Post documentation to the web.')
    parser.add_argument('-c', '--clean', action='store_true', default=False,
                        help='Clean up local directories containing generated '
                             'HTML content.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='Produce more output to the screen. ')

    pargs = parser.parse_args()
    main(pargs)
