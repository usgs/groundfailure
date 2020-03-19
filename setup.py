from distutils.core import setup
import os

setup(name='groundfailure',
      version='0.1.dev',
      description='Ground failure code (landslide/liquefaction)',
      author='Kate Allstadt, Mike Hearne, Eric Thompson, Katie Biegel',
      author_email='kallstadt@usgs.gov,mhearne@usgs.gov,emthompson@usgs.gov',
      url='http://github.com/usgs/groundfailure',
      packages=['gfail'],
      package_data={
          'gfail': [
              os.path.join('tests', 'data', '*'),
              os.path.join('gfail', 'data', '*'),
          ]
      },
      scripts=[
          'bin/gfail',
          'bin/callgf',
          'bin/create_info',
          'bin/create_png',
          'bin/gfail_transfer',
          'bin/viewdb'
      ],
      )
