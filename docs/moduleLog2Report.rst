##################################
log2Reports: Module Create Reports
##################################



Generate Report
===================

The Generate Report function creates a markdown style report of data that is provided in a python dictionary.
The report is structured in blocks according to the keys in the dictionary.
Currently, the report is structured as follows:

* Title
* Simulation name
* List of parameters and respective values
* *Release area* - table with columns for each parameter
* *Entrainment area*  - table with columns for each parameter
* *Resistance area* - table with columns for each parameter
* *images* - images with title
* *text* - text is plotted under Additional info

Input
-----

As input, a path variable and a python inputDictionary is required, optionally you can provide a plotDictionary,
including title and the absolute paths to the image you want to be included.
For the inputDictionary, there are the following options for keys to be turned into blocks:

* 'headerLine' - title of the report
* 'simName' - name of the simulation
* 'simParameters' - parameters of the simulation
* 'Release area' - all data regarding the release area scenario
* 'Entrainment area' - all data regarding the entrainment area scenario
* 'Resistance area' - all data regarding the resistance area scenario
* 'images' - title and paths to figures
* 'text' - topic and text

The plotDictionary has to be structured as follows:

* 'simName' with 'title' : 'path'


Output
-------

The default setting is to produce one report file, *fullSimulationReport.md*, however if one wants to receive on report for each *simName*,
set the flag *reportOneFile = False*.


To run
------

* call generateReport(avalancheDir, inputDictionary, plotDict, reportOneFile=True) in the run script of your computational module
