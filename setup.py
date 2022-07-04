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
from avaframe.version import getVersion

DISTNAME = 'avaframe'
LICENSE = 'EUPL'
AUTHOR = 'AvaFrame Contributors'
AUTHOR_EMAIL = 'felix@avaframe.org'
URL = 'http://avaframe.org'
CLASSIFIERS = [
                # How mature is this project? Common values are
                # 3 - Alpha  4 - Beta  5 - Production/Stable
                'Development Status :: 4 - Beta',
                # Indicate who your project is intended for
                'Intended Audience :: Science/Research',
                'License :: OSI Approved :: European Union Public Licence 1.2 (EUPL 1.2)',
                'Programming Language :: Python :: 3.6',
                'Programming Language :: Python :: 3.7',
                'Programming Language :: Python :: 3.8',
                'Programming Language :: Python :: 3.9',
                'Programming Language :: Python :: 3.10',
               ]

DESCRIPTION = 'The Open Avalanche Framework'

req_packages = ['numpy',
                'matplotlib',
                'pyshp',
                'scipy',
                'cmcrameri',
                'seaborn',
                'cython',
                'pandas',
                'shapely',
                'configUpdater',
                'tabulate'
                ]


# read the contents of your README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()


# Decide whether a cythonization of the pyx-file is required.
# if build_ext is supplied -> use cython
nos = (None, "0", "false")
subcommand = sys.argv[1] if len(sys.argv) > 1 else None
use_cython = ((subcommand == "build_ext")  # This is possibly a bit hacky.
              )

if use_cython:
    print("Package is built with cythonization.")

ext = '.pyx' if use_cython else '.c'

extensions = [Extension("avaframe.com1DFA.DFAfunctionsCython",
                        ["avaframe/com1DFA/DFAfunctionsCython"+ext],
                        include_dirs=[numpy.get_include()]
                        )]

if use_cython:
    from Cython.Build import cythonize
    extensions = cythonize(extensions,
                           compiler_directives={'linetrace': True},
                           language_level=3)

setup(
    # Project info
    name=DISTNAME,
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    version=getVersion(),
    # The project's main homepage.
    url=URL,
    # Author details
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    # License
    license=LICENSE,
    classifiers=CLASSIFIERS,
    # We are a python 3 only shop
    python_requires='>=3.6',
    # Find packages automatically
    packages=find_packages(exclude=['docs']),
    # Include package data
    include_package_data=True,
    # Install dependencies
    install_requires=req_packages,
    # additional groups of dependencies here (e.g. development dependencies).
    extras_require={},
    # Executable scripts
    entry_points={
    },
    zip_safe=False,
    ext_modules=extensions
)
