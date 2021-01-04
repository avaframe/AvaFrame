######################################
in1Data: Input data utilities
######################################


getInput
==========================

This function fetches the input data required to start avalanche simulations with com1DFA from the avalanche directory.
The avalanche directory can be generated with ``runInitializeProject.py``.


Inputs
-------

* path to avalanche directory
* configuration for the avalanche simulations
* flagDev - if development release area scenario shall be used (default: False)


Outputs
--------

* DEM
* release area scenarios
* entrainment area
* resistance area
* flag if there is an entraiment and/or resistance area



computeFromDistribution
==========================

This function allows to retrieve a sample of values that are distributed following a uniform or Pert distribution.

Inputs
-------

* configuration for the distribution


Outputs
--------

* sample derived from the chosen distribution
* various plots of the chosen distribution and the sample characteristics


To run
-------

* copy ``computeFromDistributionCfg.py`` to ``local_computeFromDistributionCfg.py`` (if not, the standard settings are used)
* adjust path to the desired ``NameOfAvalanche/`` folder in `AvaFrame/avaframe/avaframeCfg.ini``
* in ``AvaFrame/avaframe/`` run::

      python3 runComputeDist.py

.. _Theory:

Theory
-----------
