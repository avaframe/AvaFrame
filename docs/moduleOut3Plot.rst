##################################
out3Plot: Plots
##################################


quickPlot
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


quickPlotSimple
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



in1DataPlots
=================

in1DataPlots can be used to plot a sample and its characteristics derived with ``computeFromDistribution.py``,
such as: cumulative distribution function (CDF), bar plot of sample values, probability density function (PDF) of the sample,
comparison plot of empirical- and desired CDF and comparison of empirical- and desired PDF.



statsPlots
=================

statsPlots can be used produce scatter plots using a peak dictionary where information on two result parameters of avalanche simulations is saved.
This peak dictionary or a list of those can be created using the function ``extractMaxValues`` of ``ana4Stats/getStats.py``.
This can be used to visualize results of avalanche simulations where parameter variation has been used or for e.g. the case of
multiple release area scenarios. If parameter variation was used to derive the simulation results, the plots indicate the parameter values in color.
If the input data includes information about the 'scenario' that was used, the plots use different colors for each scenario.
There is also the option to add a kde plot for each result paramter as marginal plots.


Input
-----

plotValuesScatter
~~~~~~~~~~~~~~~~~~

* list of peak info dictionaries (one or more)
* result parameter type 1
* result parameter type 2
* parameter that has been varied to perform avalanche simulations
* configuration
* flagShow - if True show plots


plotValuesScatterHist
~~~~~~~~~~~~~~~~~~~~~~

* list of peak info dictionaries (one or more)
* result parameter type 1
* result parameter type 2
* configuration
* path to avalanche directory
* flagShow - if True show plots
* flagHue - if True scenario information is provided and plot uses different colors for each scenario


Outputs
-------

plotValuesScatter
~~~~~~~~~~~~~~~~~~

* scatter plot of result type 1 vs result type 2 with color indicating values of parameter that has been varied


plotValuesScatterHist
~~~~~~~~~~~~~~~~~~~~~~~

* scatter plot with marginal kde plots of result type 1 vs result type 2 with color indicating different scenarios (optional)


To run
------

runStatsExample includes an example run script to perform avalanche simulations for varying release thickness and for two different
release area scenarios. The simulation results are visualized using the plotValuesScatter and plotValuesScatterHist and saved to
``data/avaHockeySmoothChannel``. This script can be used as a template to design your own workflow.
