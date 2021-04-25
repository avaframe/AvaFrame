######################################
ana4Stats: Statistical analysis tools
######################################


probAna
==========================

probAna is used to derive simple probability maps for the point-wise probability of
an avalanche result parameter exceeding a threshold for a given set of
avalanche simulations. For example, a probability map of an avalanche exceeding a
peak pressure of 1kPa can be created.

A set of multiple avalanche simulations is required to generate these maps. This can be generated with com1DFA
using parameter variation, different release-, entrainment- or resistance scenarios.
See runProbAna.py for an example: avalanche simulations for the hockey topography
are performed with varying release thickness and a probability map is generated
ranging from 0-1, which corresponds to a probability of 0-100 %.


Inputs
-------

* path to avalanche directory with simulation result peak files
* optional input directory with simulation result peak files
* configuration settings used for avalanche simulations
* configuration settings for probability analysis
* optional output directory


Outputs
--------

* probability map saved as ascii file


To run
-------

* copy ``probAnaCfg.ini`` to ``local_probAnaCfg.ini`` (if not, the standard settings are used)
* adjust path to the desired ``NameOfAvalanche/`` folder in runProbAna.py line 28; it is also possible to add multiple
* adjust path to the com1DFACfg.ini file that shall be used to perform the avalanche simulations with parameter variation in runProbAna.py line 45
* in ``AvaFrame/avaframe/`` run::

      python3 runProbAna.py

.. _Theory:

Theory
-----------
This point-wise probability is expressed by the relative
frequency of avalanche pressure exceeding 1kPa for a set of deterministic avalanche simulations
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
* optional name of simulation scenario


Outputs
--------

readAimecRunout
~~~~~~~~~~~~~~~~~

* 1D numpy array with runout lengths

ExtractMaxValues
~~~~~~~~~~~~~~~~~

* dictionary that contains max values for all result parameters for each simulation
