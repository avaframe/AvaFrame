# -*- coding: utf-8 -*-
"""Setup file for the avaframe package.

important commands:
python setup.py sdist
python setup.py build_ext --inplace
python setup.py bdist_wheel
twine uploade dist/*

To create a release (with github):
- update the version number below
- push to github
- create a release there
- run releaseWithMany.. in the actions tab

"""

# from setuptools import setup, find_packages  # Always prefer setuptools
from setuptools import Extension, setup, find_packages
from pathlib import Path
import sys
import numpy
sys.path.append(str(Path(__file__).parent))

from Cython.Build import cythonize

# Cython part
setup_options = {}
print("Package is built with cythonization.")
setup_options = {"build_ext": {"inplace": True}}

extensions = [
    Extension(
        "avaframe.com1DFA.DFAfunctionsCython",
        ["avaframe/com1DFA/DFAfunctionsCython.pyx"],
        include_dirs=[numpy.get_include()],
    ),
    Extension(
        "avaframe.com1DFA.damCom1DFA",
        ["avaframe/com1DFA/damCom1DFA.pyx"],
        include_dirs=[numpy.get_include()],
    ),
    Extension(
        "avaframe.com1DFA.DFAToolsCython",
        ["avaframe/com1DFA/DFAToolsCython.pyx"],
        include_dirs=[numpy.get_include()],
    ),
]

extensions = cythonize(extensions, compiler_directives={"linetrace": True}, language_level=3)

setup(
    # Find packages automatically
    packages=find_packages(exclude=["docs"]),
    # Include package data
    include_package_data=True,
    # Install dependencies
    # Run build_ext
    options=setup_options,
    # Executable scripts
    entry_points={},
    zip_safe=False,
    ext_modules=extensions,
)
