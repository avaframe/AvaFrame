from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'make process faster',
  ext_modules = cythonize('processData.pyx'),
)
# to run:
# python setup.py build_ext --inplace