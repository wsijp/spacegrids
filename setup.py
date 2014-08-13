import os

#from distutils.core import setup
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname )).read()

#with open('README.rst') as file:
# long_description = file.read()

setup(name='spacegrids',
      version='1.6.8',
      author='Willem Sijp',
      author_email='w.sijp@unsw.edu.au',
      description='numpy array with grids and associated operations',
      keywords=('climate data','grid data','data on grids','spatial grids', 'Netcdf data analysis', 'climate analysis scripts',"interpreting meta data Netcdf","geophysics tools"),
      packages = find_packages(exclude="tests"),
# package_data = {
# "spacegrids": ['README.rst']
# },
      long_description=read('README'),
      url='https://github.com/willo12/spacegrids',
      license = "BSD",
#      install_requires = ["numpy>=1.6","scipy>=0.10","matplotlib>=1.1"]
      install_requires = []
# extras_require = {
# "ncio": ["netCDF4>=1.0.6"],
# "pandas": ["pandas>=0.8.0"],
# "plotting": ["pandas>=0.8.0","matplotlib>=1.2.1"],
# "interp2d": ["basemap>=1.06"],
# }
      )
