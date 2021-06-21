##################################
out1Peak: Peak plots
##################################



Plot all peak fields
====================

This function generates one plot for each peak field in *Outputs/modName/peakFiles*.
These peak fields represent the peak values of the simulation result parameters (*dynamic peak pressure*, *peak flow depth*, *peak velocity*),
and modName corresponds to the name of the computational module that has been used to perform the simualtions.


Input
-----

* avalanche directory
* configuration settings of computational module used to produce data to be plotted
* cfgFlags
* modName name of module that has been used to produce data to be plotted

where cfgFlags is a configparser object that contains the global flags, whether a plot shall be saved and/or printed to screen

Outputs
-------

* one plot of each peak field, called like the peak file (in *Outputs/out1Peak*),
if cfgFlags['ReportDir'] = True, plots are saved to *Outputs/modName/reports*
* dictionary with name of plot and full path to plot

To run
------

* call plotAllPeakFields(avalancheDir, cfg, cfgFlags, modName)


Plot all fields
====================

This function generates one plot for each simulation result field provided in the specified input directory.
This function is designed to work for result fields that follow a certain naming convention in order to provide full functionality:
*releaseAreaName_simulationType_modelType_simulationIdentifier_resultType.asc*

Input
-----

* avalanche directory
* input directory
* output directory where plots shall be saved
* result type as string (e.g. 'ppr') - optional if unit shall be specified in plot


Outputs
-------

* one plot of each field, called like the input field file


To run
------

* call plotAllFields(avalancheDir, inputDir, outDir, resType=resType)
