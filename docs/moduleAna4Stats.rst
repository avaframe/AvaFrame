######################################
ana4Stats: Statistical analysis tools
######################################


probAna
==========================

probAna is used to derive simple probability maps for a set of simulations for one avalanche track.
These maps show for each point in space the probability for a chosen parameter to exceed a given threshold.
For example, it is possible to compute the probability map of an avalanche to exceed a
peak pressure of 1kPa.

A set of multiple avalanche simulations is required to generate these maps. The simulations can be generated with :py:mod:`com1DFA`
using a parameter variation, different release-, entrainment- or resistance scenarios.
:py:mod:`runScripts.runProbAna` gives an example: avalanche simulations for the hockey topography
are performed with varying release thickness. A probability map based on peak pressure is generated.
The output is a raster file (.asc) with values ranging from 0-1. 0 meaning that no simulation exceeded the threshold
in this point in space. 1 on the contrary means that all simulations exceeded the threshold.
Details on this function, as for example required inputs can be found in: :py:mod:`ana4Stats.probAna`.


To run
-------
An example on how to generate probability maps for avalanche simulations performed with :py:mod:`com1DFA`
is given in :py:mod:`runScripts.runProbAna`, where for *avaHockeyChannel* simulations are performed with
varying release thickness values ranging from 0.75 to 1.75 meters in steps of 0.05 meters.
The resulting simulations are then used to generate the probability map. There is also the option
to filter the simulations further - using the function :py:func:`in3Utils.fileHandlerUtils.getFilterDict` which generates a
parameter dictionary for filtering according to the filter criteria set in the
configuration file (``ana4Stats/probAnaCfg.ini``) of the :py:mod:`ana4Stats.probAna` function.
In order to run this example:

* first go to ``AvaFrame/avaframe``
* copy ``ana4Stats/probAnaCfg.ini`` to ``ana4Stats/local_probAnaCfg.ini``
* uncomment ``'FILTER'`` section in ``local_probAnaCfg.ini`` and insert filter parameters if you want to first filter simulations
* run::

      python3 runScripts/runProbAna.py

.. _Theory:

Theory
-----------
This point-wise probability is expressed by the relative
frequency of avalanche peak flow field exceeding a certain threshold for a set of deterministic avalanche simulations
derived from a range of input parameters (see :cite:`HyBeBu2019`).


getStats
==========================

In :py:mod:`ana4Stats.getStats`, functions that help to compute statistical properties of simulation results are gathered.
In :py:func:`ana4Stats.getStats.readAimecRunout`, the runout length is read from the ana3Aimec output and saved into a 1D numpy array.
:py:func:`ana4Stats.getStats.extractMaxValues` can be used to determine the maximum peak values of the simulation results.
These values can then be plotted using the functions in :py:mod:`out3Plot.statsPlots` in order to visualise the statistics of
a set of avalanche simulations. 
For further details on the specific functions, have a look at: :py:mod:`ana4Stats.getStats`.

To run
-------

An example on how to use these statistical functions is given in :py:mod:`runScripts.runStatsExample`, where
for *avaHockeyChannel* simulations are performed for two different release area scenarios and
the release thickness is varied from 0.75 to 1.75 meters in steps of 0.05 meters. The resulting
simulations are then analysed using the :py:func:`ana4Stats.getStats.extractMaxValues` function and plots are generated using the
plotting routines from :py:mod:`out3Plot.statsPlots` .

* first go to ``AvaFrame/avaframe``
* copy ``ana4Stats/getStats.ini`` to ``ana4Stats/local_getStatsCfg.ini``
* uncomment ``'FILTER'`` section in ``ana4Stats/local_getStatsCfg.ini`` and insert filter parameters if you want to first filter simulations
* run::

      python3 runScripts/runStatsExample.py
