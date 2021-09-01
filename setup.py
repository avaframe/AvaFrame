# -*- coding: utf-8 -*-
"""Setup file for the avaframe package.
   Adapted from the Python Packaging Authority template."""

from setuptools import setup, find_packages  # Always prefer setuptools

with open('README.md') as f:
    readme = f.read()

with open('LICENSE.txt') as f:
    license = f.read()

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
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ]

DESCRIPTION = 'The Open Avalanche Framework'

LONG_DESCRIPTION = readme

req_packages = ['numpy',
                ]

setup(
    # Project info
    name=DISTNAME,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    version='0.5',
    # The project's main homepage.
    url=URL,
    # Author details
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    # License
    license=LICENSE,
    classifiers=CLASSIFIERS,
    # We are a python 3 only shop
    python_requires='>=3.5',
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
)
