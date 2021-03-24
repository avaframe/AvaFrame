# -*- coding: utf-8 -*-
"""Setup file for the avaframe package.
   Adapted from the Python Packaging Authority template.
   See MANIFEST.in for non-python files (include an exclude)

"""

from setuptools import setup, find_packages  # Always prefer setuptools

with open('README.md') as f:
    readme = f.read()

with open('LICENSE.txt') as f:
    license = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

DISTNAME = 'avaframe'
LICENSE = 'EUPL'
AUTHOR = 'AvaFrame Contributors'
AUTHOR_EMAIL = 'felix@avaframe.org'
URL = 'http://avaframe.org'
CLASSIFIERS = [
        # How mature is this project? Common values are
        # 3 - Alpha  4 - Beta  5 - Production/Stable
        'Development Status :: 3 - Alpha',
        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        # 'License :: OSI Approved :: EUPL License',
        'Programming Language :: Python :: 3',
    ]

DESCRIPTION = 'The Open Avalanche Framework'

LONG_DESCRIPTION = readme

req_packages = ['pyshp',
                'seaborn',
                'cmocean'
                ]
# with open('requirements.txt') as f:
#     requirements = f.read().splitlines()

def local_scheme(version):
    return ""

setup(
    # Project info
    name=DISTNAME,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    version='0.2.4',
    # The project's main homepage.
    url=URL,
    # Author details
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    use_scm_version={"local_scheme": local_scheme},
    # use_scm_version=True,
    setup_requires=['setuptools_scm'],

    # License
    license=LICENSE,
    classifiers=CLASSIFIERS,
    # We are a python 3 only shop
    python_requires='>=3.5',
    # Find packages automatically
    packages=find_packages(exclude=['docs','benchmarks']),
    # Include package data
    include_package_data=True,
    # Install dependencies
    # install_requires=req_packages,
    install_requires=requirements,
    # additional groups of dependencies here (e.g. development dependencies).
    extras_require={},
    # Executable scripts
    entry_points={
    },
)
