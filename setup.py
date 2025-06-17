from setuptools import setup, Extension
import numpy

from Cython.Build import cythonize

# Define extension modules conditionally
# ext_modules = []

# List all potential Cython file locations
cython_files = [
    'avaframe/com1DFA/DFAfunctionsCython.pyx',
    'avaframe/com1DFA/DFAToolsCython.pyx',
    'avaframe/com1DFA/damCom1DFA.pyx',
]

extensions = [
    Extension(
        name=file.replace('/', '.').replace('.pyx', ''),
        sources=[file],
        include_dirs=[numpy.get_include()]
    )
    for file in cython_files
]

ext_modules = cythonize(
    extensions,
    compiler_directives={
        'language_level': '3',
        'linetrace': True,
    }
)

setup_options = {"build_ext": {"inplace": True}}

setup(
    options=setup_options,
    ext_modules=ext_modules,
    # install_requires=[
    #     'numpy',
    #     'scipy',
    #     'cython',
    #     'matplotlib',
    #     'pandas'
    # ],
    # python_requires='>=3.8',
)