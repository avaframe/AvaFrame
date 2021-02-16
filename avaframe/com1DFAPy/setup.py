from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

# ext_modules = [Extension("DFAfunctionsCython",
#                ["avaframe/com1DFAPy/DFAfunctionsCython.pyx"],
#               libraries=["m"],
#               extra_compile_args=["-ffast-math"],
#               include_dirs=[numpy.get_include()])]

setup(
    ext_modules=cythonize("avaframe/com1DFAPy/DFAfunctionsCython.pyx"),
    include_dirs=[numpy.get_include()]
)
