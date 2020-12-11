from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules=[Extension("SPHfunctionsCython",
             ["SPHfunctionsCython.pyx"],
              libraries=["m"],
              extra_compile_args = ["-ffast-math"])]

setup(
  name = "SPHfunctionsCython",
  cmdclass = {"build_ext": build_ext},
  ext_modules = ext_modules)
