from distutils.core import setup

setup(name='groundfailure',
      version='0.1dev',
      description='Ground failure code (landslide/liquefaction)',
      author='Kate Allstadt, Mike Hearne, Eric Thompson, Katie Biegel',
      author_email='kallstadt@usgs.gov,mhearne@usgs.gov,emthompson@usgs.gov,kbiegel@usgs.gov',
      url='http://github.com/usgs/groundfailure',
      packages=['groundfailure'],
      package_data={'groundfailure': ['configspec.ini']},
      scripts=['gfail'],
      )
