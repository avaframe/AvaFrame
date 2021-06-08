from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

setup(
    ext_modules=cythonize("avaframe/com1DFA/DFAfunctionsCython.pyx",
                          compiler_directives={'linetrace': True},
                          language_level=3),
    include_dirs=[numpy.get_include()]
)
