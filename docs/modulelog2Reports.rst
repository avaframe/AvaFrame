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
* if key: *releaseArea* is provided - table with columns for each parameter
* if key: *entrainmentArea* is provided - table with columns for each parameter
* if key: *resistanceArea* is provided - table with columns for each parameter
* if key: *images* is provided - images are plotted

Input
-----

As input, a path variable and a python dictionary is required.
For the dictionary, the following keys are required:

* 'headerLine' - title of the report
* 'simName' - name of the simulation

The following keys are optional:

* 'simParameters' - parameters of the simulation
* 'releaseArea' - all data regarding the release area scenario
* 'entrainmentArea' - all data regarding the entrainment area scenario
* 'resistanceArea' - all data regarding the resistance area scenario
* 'images' - paths to figures


Outputs
-------

A report in markdown format, called:

* simName.md


To run
------

* call generateReport(avalancheDir, inputDictionary) in the run script of your computational module
