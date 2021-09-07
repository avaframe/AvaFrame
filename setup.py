# -*- coding: utf-8 -*-
"""Setup file for the avaframe package.
   Adapted from the Python Packaging Authority template."""

from setuptools import setup, find_packages  # Always prefer setuptools
from pathlib import Path
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy


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
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ]

DESCRIPTION = 'The Open Avalanche Framework'

req_packages = ['numpy',
                'matplotlib',
                'pyshp',
                'scipy',
                'cmcrameri',
                'seaborn',
                'cython',
                'pandas'
                ]


# read the contents of your README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    # Project info
    name=DISTNAME,
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    version='0.5.2',
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
    ext_modules=cythonize("avaframe/com1DFA/DFAfunctionsCython.pyx",
                          compiler_directives={'linetrace': True},
                          language_level=3),
    include_dirs=[numpy.get_include()]
)
