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
* if key: *Release area* is provided - table with columns for each parameter
* if key: *Entrainment area* is provided - table with columns for each parameter
* if key: *Resistance area* is provided - table with columns for each parameter
* if key: *images* is provided - images are plotted
* if key: *text* is provided - text is plotted under Additional info

Input
-----

As input, a path variable and a python inputDictionary is required, optionally you can provide a plotDictionary,
including with title and absolute paths to the images you want to be included.
For the inputDictionary, the following options for keys are turned into blocks:

* 'headerLine' - title of the report
* 'simName' - name of the simulation
* 'simParameters' - parameters of the simulation
* 'Release area' - all data regarding the release area scenario
* 'Entrainment area' - all data regarding the entrainment area scenario
* 'Resistance area' - all data regarding the resistance area scenario
* 'images' - title and paths to figures
* 'text' - topic and text


Outputs
-------

A report in markdown format, default all info is plotted consecutively into one file, called *fullSimulationReport.md*.
If reportOneFile is set to False when calling generateReport, one report file is created for each simulation, called *simName.md*.



To run
------

* call generateReport(avalancheDir, inputDictionary, plotDict) in the run script of your computational module
