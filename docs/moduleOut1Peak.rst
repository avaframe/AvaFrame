##################################
out1Peak: Peak plots
##################################

:py:mod:`out1Peak.outPlotAllPeak`

Plot all peak fields
====================

This function generates one plot for each peak field in *Outputs/com1DFA/peakFiles*.
These peak fields represent the peak values of the simulation result parameters (*dynamic peak pressure*, *peak flow depth*, *peak velocity*).


Input
-----

* avalanche directory

* configuration settings of com1DFA

* cfgFlags

* optional: name of module that has been used to produce data to be plotted

where cfgFlags is a configparser object that contains the global flags, whether a plot shall be saved and/or printed to screen

Outputs
-------

* one plot of each peak field, called like the peak file (in *Outputs/out1Peak*)
* dictionary with name of plot and full path to plot

To run
------

* call plotAllPeakFields(avalancheDir, cfg)


Plot all fields
====================

This function generates one plot for each simulation result field provided in the specified input directory.


Input
-----

* avalanche directory
* input directory
* output directory where plots shall be saved
* configuration settings


Outputs
-------

* one plot of each field, called like the input field file


To run
------

* call plotAllFields(avalancheDir, inputDir, outDir, cfg)
