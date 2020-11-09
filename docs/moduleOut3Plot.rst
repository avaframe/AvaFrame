##################################
out3Plot: Peak plots
##################################



Plot peak field
===================

This function generates a plot for a specific peak field from *Outputs/com1DFA/peakFiles*.
These peak field represents the peak value of the passed result parameter.

Input
-----

* avalanche directory
* simulation name
* result parameter (*ppr*, *pfd*, *pv*)
* unit


Outputs
-------

* plot of peak field, called like the peak file (in *Outputs/out1Peak*)
* dictionary with name of plot and full path to plot

To run
------

* call plotPeakFields(avalancheDir, simName, resType, unit)


Quick plot
===========

This function creates two plots, one plot with three panels, first dataset, second dataset and the absolute difference of the two datasets and
one plot showing a cross- and a longprofile of the two datasets.
According to the characteristics of a simulation, as for example provided by the simulation dictionary created with *fileHandlerUtils.makeSimDict*,
the required input parameters can be defined in order to choose the first dataset.
The second dataset is fetched from the benchmark results using the same filtering characteristics.
The function returns a list with the file paths of the generated plots.


Input
-----

* avalanche directory
* result parameter ('ppr', 'pfd', 'pv')
* parameter (e.g. 'Mu', 'RelTh', 'simType')
* values of the parameter (e.g. '0.155', '1.000', 'entres')
* configuration setting of com1DFA
* configuration setting for outQuickPlot


Outputs
-------

* plot with three panels showing the specified peak field, the reference peak field, and their difference including a text box with absolute Mean, Max and Min of the difference
* plot of a cross- and a longprofile of the two datasets
* list with full file paths of the generated plots

Both plots are saved in *Outputs/out3Plot*.



To run
------

* runQuickPlot.py

In this run script, one can specify the result parameter and a second parameter, like e.g. 'Mu' or 'simType' to specify the datasets to be plotted.


Quick plot simple
=================

This function creates two plots of any two raster datasets (with identical dimensions).
The first plot has three panels, showing the first dataset, the second dataset and the absolute difference of the two datasets and
the second plot shows a cross- and a longprofile of the two datasets.
The function returns a list with the file paths of the generated plots.


Input
-----

* avalanche directory
* input directory (where the two datasets are located)

Outputs
-------

* plot with three panels showing dataset 1, dataset 2, and their difference including a text box with absolute Mean, Max and Min of the difference
* plot of a cross- and a longprofile of the two datasets
* list with full file paths of the generated plots

Both plots are saved in *Outputs/out3Plot*.


To run
------

* specifiy input directory, default is *Work/simplePlot*
* runQuickPlotSimple.py
