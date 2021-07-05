##################################
log2Reports: Create Reports
##################################



Generate Report
===================

:py:mod:`log2Report.generateReport` creates a markdown style report corresponding to the
data provided in the python input dictionary.
The report is structured in blocks according to the types provided in the dictionary.
Currently, the report generation supports the following types:

* title - creates a title header
* avaName - includes a line with the avalanche directory
* simName - includes a line with Simulation name in the header
* time - includes a line with the date
* list - creates a table with two columns where the keys and values are listed and the dictionary name is set as title
* columns - creates a table with two rows and one column for each key and value pair and the dictionary name is set as title
* image - creates a block with the dictionary name as title and prints the key as heading with the prefix *Figure* and plots the image provided as value
* text - creates a text block with the dictionary name as title and prints the key as header with the prefix *Topic* and prints the text provided as value


For the input dictionary, a certain structure is required as illustrated for the following example:

::

  reportDict = {'title Block' : {'type' : 'title', 'title' : 'This is my title'},
                'Simulation Name' : {'type' : 'simName', 'name' : 'This is my desired simulation name'},
                'Simulation Parameters' : {'type' : 'list', 'Parameter 1' : 0.155, 'Parameter 2' : 'red'},
                'Release Area' : {'type' : 'columns', 'scenario' : 'rel1', 'area' : '40000', 'exposition' : 'south'},
                'Special Information' : {'type' : 'text', 'Version' : 'This simulation was performed with a dev version.'}
                }

where the keys and values do not have to follow a predefined form, only the key *type* has to be provided.
This example dictionary results in a report of the form:

  .. figure:: _static/report.png

If a separate plotDictionary shall be included, in addition the key *simName* is required so that it can be appended to the correct simulation dictionary.

The default setting is to produce one report file, *fullSimulationReport.md*, however if one wants to receive on report for each *simName*,
set the flag *reportOneFile = False* in ``avaframeCfg.ini`` or preferabley your local copy ``local_com1DFACfg.ini``.
For further details have a look at :py:mod:`log2Report.generateReport` .


Generate Compare Report
=========================

:py:mod:`log2Report.generateCompareReport` creates a markdown style report where the simulation results are compared to the benchmark results.
The report is structured as follows:

  * name of avalanche
  * simulation name
  * text block with info on test
  * table listing simulation parameter name, value (reference) and value (simulation)
    -> if value of simulation run differs from reference value, it is highlighted in red
  * table listing aimec analysis measure, value (reference) and value (simulation)
    -> if value of simulation run differs from reference value, it is highlighted in red
  * block with plots including a header with plot title
    -> if differences exceed threshold a Warning is printed for each plot

For this purpose, python dictionaries are used to read the required input data. In order to produce the structure described above,
the dictionaries require a certain structure too, for example:

::

    dictionary = {
                  '*simName*': 'name of simulation',
                  '*Simulation Parameters*': {
                  '*type*': '*list*',
                  'name of parameter': 'value of parameter'
                  ...
                                              },
                  '*Test Info*': {
                  '*type*': '*text*',
                  'Title of text box': 'text, text',
                  ...
                                  },
                  '*Simulation Results*': {
                  'desired plot title': '*file path*',
                  ...
                                          }
                                  },

where the required **keys** are indicated by the asterisk sign.
The settings can be found in the respective configuration file ``log2Report/genereateCompareReportCfg.ini``.
For further details visit :py:mod:`log2Report.generateCompareReport`.

An example of creating a comparison report can be found in :py:mod:`runStandardTestsCom1DFA.py`. There, all the standard tests are run sequentially
and a report is created where the com1DFA simulation results
and the benchmark results are compared.
