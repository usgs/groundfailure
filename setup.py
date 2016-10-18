from distutils.core import setup
import os

setup(name='groundfailure',
      version='0.1dev',
      description='Ground failure code (landslide/liquefaction)',
      author='Kate Allstadt, Mike Hearne, Eric Thompson, Katie Biegel',
      author_email='kallstadt@usgs.gov,mhearne@usgs.gov,emthompson@usgs.gov,kbiegel@usgs.gov',
      url='http://github.com/usgs/groundfailure',
      packages=['groundfailure',
                'tests'],
      package_data={'groundfailure': ['configspec.ini',
                                      os.path.join('tests', 'data', '*'),
                                      ]},
      scripts=['gfail'],
      )
