from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

# ext_modules=[Extension("cytest",
#              ["cytest.pyx"],
#               libraries=["m"],
#               extra_compile_args = ["-ffast-math"],
#               include_dirs=[numpy.get_include()])]
#
# setup(
#   name = "cytest",
#   cmdclass = {"build_ext": build_ext},
#   ext_modules = ext_modules)

ext_modules=[Extension("SPHfunctionsCython",
             ["SPHfunctionsCython.pyx"],
              libraries=["m"],
              extra_compile_args = ["-ffast-math"],
              include_dirs=[numpy.get_include()])]

setup(
  name = "SPHfunctionsCython",
  cmdclass = {"build_ext": build_ext},
  ext_modules = ext_modules)


#
# setup(
#     ext_modules=[
#         Extension("my_module", ["my_module.c"],
#                   include_dirs=[numpy.get_include()]),
#     ],
# )
#
# # Or, if you use cythonize() to make the ext_modules list,
# # include_dirs can be passed to setup()
#
# setup(
#     ext_modules=cythonize("my_module.pyx"),
#     include_dirs=[numpy.get_include()]
# )
