######################################
in1Data: Input data utilities
######################################

Get input data
================

:py:mod:`in1Data.getInput` module provides functions to fetch or generate input data for avalanche simulations within AvaFrame.
The main functions are described in more detail in the following sections.
A detailed description of all the small helper functions is provided in :py:mod:`in1Data.getInput`


getInputCom1DFA
----------------

This function fetches the input data from an avalanche directory, that is required to start avalanche simulations with com1DFA.
This data consists of a digital elevation model (DEM), a release scenario shapefiles and optional entrainment and resistance area shapefiles.


Inputs
~~~~~~

* path to avalanche directory
* configuration for the avalanche simulations - provides info on whether entrainment and resistance shall be included


Outputs
~~~~~~~~~

* inputSimFiles - dictionary with:
  - DEM
  - release area scenarios
  - secondary release area scenario
  - entrainment area
  - resistance area
  - entResInfo - dictionary with info if there is an entraiment and/or resistance area


computeFromDistribution
==========================

This collection of functions allows to retrieve a sample of values that are distributed following a uniform or Pert distribution.
This is for example useful if you want to perform avalanche simulations for a range of release thickness values,
that should be distributed following a specific distribution, or another parameter that can be varied.
Detailed information on the functions can be found in:

:py:mod:`in1Data.computeFromDistribution`

Inputs
-------

* configuration for the distribution


Outputs
--------

* sample derived from the chosen distribution
* various plots of the chosen distribution and the sample characteristics


To run
-------

An example of how a parameter set distributed following a uniform or Pert distribution can be
generated, is provided by ``runComputeDist.py``.

* first go to ``AvaFrame/avaframe``
* copy ``computeFromDistributionCfg.py`` to ``local_computeFromDistributionCfg.py`` (if not, the standard settings are used)
* adjust path to the desired ``NameOfAvalanche/`` folder in your local copy of ``avaframeCfg.ini``
* run::

      python3 runComputeDist.py

.. _Theory:

Theory
-----------
