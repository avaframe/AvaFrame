##################################
out3SimpPlot: Module simple plots
##################################



Plot all peak fields
===================

This function generates one plot for each peak field in *Outputs/com1DFA/peakFiles*.
These peak fields represent the peak values of the simulation result parameters (*dynamic peak pressure*, *peak flow depth*, *peak velocity*).


Input
-----

* avalanche directory
* configuration settings of com1DFA
* cfgFlags

where cfgFlags is a configparser object that contains the global flags, whether a plot shall be saved and/or printed to screen

Outputs
-------

* one plot of each peak field, called like the peak file (in *Outputs/out3SimpPlot*)
* dictionary with name of plot and full path to plot

To run
------

* call plotAllPeakFields(avalancheDir, cfg)
