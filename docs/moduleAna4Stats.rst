######################################
ana4Stats: Statistical analysis tools
######################################


probAna
==========================

probAna is used to derive simple probability maps for a set of simulations for one avalanche track.
These maps show for each point in space the probability for a chosen parameter to exceed a given threshold.
For example, it is possible to compute the probability map of an avalanche to exceed a
peak pressure of 1kPa.

A set of multiple avalanche simulations is required to generate these maps. The simulations can be generated with com1DFA
using a parameter variation, different release-, entrainment- or resistance scenarios.
``runProbAna.py`` gives a example: avalanche simulations for the hockey topography
are performed with varying release thickness. A probability map based on peak pressure is generated.
The output is a raster file (.asc) with values ranging from 0-1. 0 meaning that no simulation exceeded the threshold
in this point in space. 1 on the contrary means that all simulations exceeded the threshold.


Inputs
-------

* path to avalanche directory with simulation result peak files
* configuration settings for probability analysis
* computational module used to perform avalanche simulations
* optional parameter dictionary used to filter simulations
* optional input directory with simulation result peak files

Outputs
--------

* probability map saved as ascii file


To run
-------
An example on how to generate probability maps for avalanche simulations performed with com1DFA
is given in ``runProbAna.py``, where for avaHockeyChannel simulations are performed
varying release thickness values ranging from 0.75 to 1.75 meters in steps of 0.05 meters.
The resulting simulations are then used to generate the probability map. There is also the option
to filter the simulations further - using the function ``getFilterDict`` which generates a
parameter dictionary for filtering according to the filter criteria set in the
configuration file (``probAnaCfg.ini``) of the ``probAna`` function.

* copy ``probAnaCfg.ini`` to ``local_probAnaCfg.ini``
* uncomment FILTER section in ``local_probAnaCfg.ini`` and insert filter parameters if you want to first filter simulations
* in ``AvaFrame/avaframe/`` run::

      python3 runProbAna.py

.. _Theory:

Theory
-----------
This point-wise probability is expressed by the relative
frequency of avalanche peak flow field exceeding a certain threshold for a set of deterministic avalanche simulations
derived from a range of input parameters (see :cite:`HyBeBu2019`).


getStats
==========================

Here, functions that help to compute statistical properties of simulation results are gathered.
In ``readAimecRunout``, the runout length is read from the ana3Aimec output and saved into a 1D numpy array.
``ExtractMaxValues`` can be used to determine the maximum peak values of the simulation results.


Inputs
-------

readAimecRunout
~~~~~~~~~~~~~~~~

* path to directory of ana3AIMEC results
* name of avalanche
* configuration settings used for ana3AIMEC

ExtractMaxValues
~~~~~~~~~~~~~~~~~~~~

* path to directory with peak files of com1DFA simulation results
* configuration settings used for performing com1DFA simulations
* path to avalanche directory
* parameter that was varied to perform avalanche simulations (e.g. release thickness)
* optional name of simulation scenario (e.g. release Scenario)
* optional parameters dictionary that provides criteria for filtering avalanche simulations


Outputs
--------

readAimecRunout
~~~~~~~~~~~~~~~~~

* aimec analysis including resluts text file and plots, see for further information :ref:`moduleAna3AIMEC:ana3AIMEC: Aimec`

ExtractMaxValues
~~~~~~~~~~~~~~~~~

* dictionary that contains max values for all result parameters for each simulation


To run
-------

An example on how to use these statistical functions is given in ``runStatsExample.py``, where
for avaHockeyChannel simulations are performed for two different release area scenarios and
the release thickness is varied from 0.75 to 1.75 meters in steps of 0.05 meters. The resulting
simulations are then analysed using the extractMaxValues function and plots are generated using the
plotting routines from ``statsPlots.py``.

* copy ``getStats.ini`` to ``local_getStatsCfg.ini``
* uncomment FILTER section in ``local_getStatsCfg.ini`` and insert filter parameters if you want to first filter simulations
* in ``AvaFrame/avaframe/`` run::

      python3 runStatsExample.py
