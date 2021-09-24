Release Notes
=============

v0.6 (24. September 2021)
-------------------------

ENHANCEMENTS

- installation via pypi (pip install)
- connection to QGis (via plugin manager) 
- function to interpolate data on mesh of different cellSize using splinesp
- testing via pytest extended
- more pathlib usage 
- ASCII header is read as dict
- documentation contains FAQ page
- reworked installation instructions
- cleaner test reports/inis
- github action to deploy to pypi
- switch to codeclimate
- use consistent thickness attributes (shapefiles etc)
- *com1DFA* any resolution is possible now 
- *com1DFA* split the getWeight function in two: first get cell and then get weights. 
- *com1DFA* avoid possibility of segfault because particles exit too quickly the domain.
- *com1DFA* additional particles info: unique identifier for each particle and parent particles
- *com1DFA* central time step calling
- *com1DFA* additional options to set mass per particle directly or via release thickness
- *com1DFA* interpolation option for initialization of Hpart 
- *com1DFA* read entrainment thickness
- *ana3AIMEC* override option for raster cellsize 
- *ana3AIMEC* mass analysis plot even if more than 2 simulations

FIX

- getTimeIndex problem if dtSave < actual dt
- better way to remove particles
- track particles exiting the computation domain
- fix issue save particles
- read aimec grid info from result files and not from dem
- add reasonString to removal of particles
- fix correct module name in AIMEC 
- com2AB write out to shp 

Contributers:

- **Code: core team**


v0.5 (13. July 2021)
--------------------

ENHANCEMENTS

- filtering functions for com1DFA simulations
- flag to disable print at CFG reading
- new colormaps for ppr, pfd, pfv
- *com1DFA* option to add friction explicitly using the method described in #273 .
- *com1DFA* Resistance force is  added explicitly.
- *com1DFA* New method to get the release area
- *com2AB* function to write results to shapefile
- *ana3AIMEC* warning for empty runout zone
- *ana3AIMEC* enable simulation ordering/filtering

FIX

- beta angle issue i.e. distance below angle
- correct removal of particles 
- AIMEC produces warning on empty runout area
- adapt quickplot to new naming scheme

Contributers:

- **Code: core team**
- **Colormaps: C.Tollinger**

DOI for this release:

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.5094509.svg
   :target: https://doi.org/10.5281/zenodo.5094509


v0.4.1 (9. June 2021)
---------------------

Minor release to fix issue with zenodo

v0.4 (8. June 2021)
-------------------

The switch release

This is a big release: we switched our dense flow module 'com1DFA' to the python
version. This means that you know get to use the python version as default.
However, the original version is still available in the module com1DFAOrig. The
full documentation for the python com1DFA version as well as updated benchmarks
will be released in the next version.

Module com2AB (AlphaBeta) recieved an update allowing for custom parameters.

Simulation naming and identification also recieved a major change, we introduced
unique ID's for each individual configuration.

Contributers:

- **Code: core team**


v0.3 (26. April 2021)
---------------------

The AIMEC and Windows release

This release brings an AIMEC refactor, plenty of improvements related to the
test cases and Windows capabilities. 3 new idealised/generic test case are 
included: flat plane, inclined slope and pyramid.

Com1DFAPy recieved a lot of advancement as well, e.g. parts of it are converted
to cython to speed up computation times.  

Documentation regarding our testing is included, see more at the
`testing <https://docs.avaframe.org/en/latest/testing.html>`_ page. 

Contributers:

- **Code: core team**

DOI for this release:

.. image:: https://zenodo.org/badge/281922740.svg
   :target: https://zenodo.org/badge/latestdoi/281922740


v0.2 (28. Dezember 2020)
------------------------

The testing release

Version 0.2 includes the first real world avalanches. It provides data for 6
avalanches, including topographies, release areas and benchmark results.
To know more about our data sources, head over to
`our data sources documentation
<https://docs.avaframe.org/en/latest/dataSources.html>`_.
The existing test cases also recieved some updates by including multiple release
areas and multiple scenarios per avalanche.  

This release also is the first to include `API documentation
<https://docs.avaframe.org/en/latest/api.html>`_ for our modules and functions.
However not all functions are included yet.

Contributers:

- **Data: M.Granig, C. Tollinger**
- **Data: Land Tirol**
- **Code: core team**


v0.1 (06 November 2020)
-----------------------

Initial release. 

This release is the result of several months of development.

Several people have contributed to this release, either directly or through code
that was used as reference/basis:

- **Peter Sampl**, code base for com1DFA
- **Jan-Thomas Fischer**, code base AIMEC, code related to com1DFA
- **Michael Neuhauser**, code for helper and transformation utilities, com1DFA
- **Andreas Kofler**, code related to AIMEC and com1DFA 

and the core team:

- **Anna Wirbel**
- **Matthias Tonnel**
- **Felix Oesterle**

