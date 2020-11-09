Workflow example
================================

The following example should make it easier for you to find your way in AvaFrame
and setup your own AvaFrame workflow. Before starting, make sure you have
installed AvaFrame and the required dependencies as explained in
:ref:`installation:Installation`.

Make sure you change to your ``AvaFrame`` direcotry by::

  cd [YOURDIR]/AvaFrame

Replace ``[YOURDIR]`` with the directory from your installation step.

Initialize project
------------------
To create the folder where the input data lies and where the
output results will be saved, specify the full path to the folder
in the ``local_avaframeCfg.ini`` (which is a copy of
``avaframeCfg.ini`` that you need to create). So::

  cd avaframe
  cp avaframeCfg.ini local_avaframeCfg.ini

and edit ``local_avaframeCfg.ini`` with your favorite text editor and adjust the
variable ``avalancheDir``.

Then run ::

  python runInitializeProject.py

This will create a new directory with the input required by AvaFrame structured as described
in :ref:`moduleIn3Utils:Initialize Project`. 

Input data
------------

Check the input data required by the different modules you want to use and fill
the ``Inputs/`` inside the ``[avalancheDir]`` folder from the initialize step accordingly.

For example the ``com1DFA`` module needs input as described in
:ref:`moduleCom1DFA:Input`. You can also have a look at the default setting for
the module you want to use (for example ``com1DFACfg.ini`` for module com1DFA).
If you want to use different settings, create a ``local_`` copy of the ``.ini``
file and modify the desired parameters.

More information about the configuration can be found here: :ref:`configuration:Configuration`

Building your run script
-------------------------

Create your own workflow by taking the ``runOperational.py`` script as
template.

We suggest you copy it and adjust it to your liking. 

Imports
~~~~~~~~~

The run script starts with importing python modules ::

  # Load python modules
  import os
  import time
  import logging

Then import the avaframe modules you want to use. This list will
grow as you use more avaframe modules. You can refer to the different
computational modules documentation to know which imports are required::

  # Local avaframe imports
  # first import the always useful tools
  from avaframe.in3Utils import cfgUtils
  from avaframe.in3Utils import logUtils
  from avaframe.in3Utils import initializeProject as initProj
  from avaframe.log2Report import generateReport as gR

  # then depending on which computational module you want to use

  # for the Com1DFA module
  from avaframe.com1DFA import com1DFA
  from avaframe.out1Peak import outPlotAllPeak as oP

  # for the Alpha Beta module
  from avaframe.com2AB import com2AB
  from avaframe.out3Plot import outAB

General configuration
~~~~~~~~~~~~~~~~~~~~~

The next step is to get the general settings::

  # Load avalanche directory from general configuration file
  cfgMain = cfgUtils.getGeneralConfig()
  avalancheDir = cfgMain['MAIN']['avalancheDir']

More information about the configuration can be found here: :ref:`configuration:Configuration`


Logs
~~~~

Prepare the logging::

  # log file name; leave empty to use default runLog.log
  logName = 'runOperational'
  # Start logging
  log = logUtils.initiateLogger(avalancheDir, logName)
  log.info('MAIN SCRIPT')
  log.info('Current avalanche: %s', avalancheDir)

If you just created the ``avalancheDir`` this one should be clean but if you
already did some calculations you might want to clean it::

  # Clean input directory(ies) of old work and output files
  initProj.cleanSingleAvaDir(avalancheDir, keep=logName)

More information about logging can be found here: :ref:`logging:Logging`


Computation Modules
~~~~~~~~~~~~~~~~~~~

Now the computation modules are called.

For the Com1DFA (dense flow) routine, use::

  # Load input parameters from configuration file
  cfgCom1DFA = cfgUtils.getModuleConfig(com1DFA)

  # Run dense flow
  reportDictList = com1DFA.com1DFAMain(cfgCom1DFA, avalancheDir)

  # Collect results/plots/report  to a single directory
  # Generate plots for all peakFiles
  plotDict = oP.plotAllPeakFields(avalancheDir, cfgCom1DFA, cfgMain['FLAGS'])

For the Com2AB (alpha beta) routine, do::

  # Load input parameters from configuration file
  cfgAB = cfgUtils.getModuleConfig(com2AB)
  
  # Run Alpha Beta
  resAB = com2AB.com2ABMain(cfgAB, avalancheDir)

  # Collect results/plots/report  to a single directory
  reportDictList, _, _ = outAB.writeABpostOut(resAB, cfgAB, reportDictList)

Timing
~~~~~~

If you want to time the execution, add following at the top of your script::

  # Time the routine
  startTime = time.time()

Followed by, at the end::

  # Print time needed
  endTime = time.time()
  log.info('Took %s seconds to calculate.' % (endTime - startTime))



Report
~~~~~~

It is also possible to generate a report file (:ref:`moduleLog2Report:Generate Report`)
in addition to the ``.log`` file::

  # Set directory for report
  reportDir = os.path.join(avalancheDir, 'Outputs')
  # write report
  gR.writeReport(reportDir, reportDictList, cfgMain['FLAGS'])

In the case of the Com1DFA routine, you could just add plots of the peak fields within the report
by first creating the plots::

  # Generate plots for all peakFiles
  plotDict = oP.plotAllPeakFields(avalancheDir, cfg, cfgMain['FLAGS'])

and inlcude this plots in the report by passing the plotDict::

  # write report
  gR.writeReport(reportDir, reportDictList, cfgMain['FLAGS'], plotDict)



