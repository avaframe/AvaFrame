Configuration
-------------

In order to set the configurations required by all the modules within Avaframe, the python module
`configparser <https://docs.python.org/3/library/configparser.html>`_ is used.

This is done in two steps. The first step fetches the main settings::

  from avaframe.in3Utils import cfgUtils
  # Load avalanche directory from general configuration file
  cfgMain = cfgUtils.getGeneralConfig()
  avalancheDir = cfgMain['MAIN']['avalancheDir']

In the second step the specific settings to a given module are imported::

  from avaframe.tmp1Ex import tmp1Ex
  # Load all input Parameters from config file
  # get the configuration of an already imported module
  # Write config to log file
  cfg = cfgUtils.getModuleConfig(tmp1Ex)

The :py:func:`in3Utils.cfgUtils.getModuleConfig` function reads the settings from a configuration file (``tmpEx.ini``
in our example) and writes these settings to the log file. The default settings can be found in the
configuration file provided within each module.

It is possible to modify these settings, there are two options:

* provide the path to your own configuration file when calling
  ``cfgUtils.getModuleConfig(moduleName, path to config file)``

* create a copy of the module configuration file called ``local_`` followed by
  the name of the original configuration file and set the desired values of the
  individual parameters.

So the order is as follows:

#. if there is a path provided, configuration is read from this file.

#. if there is no path provided, the ``local_...`` configuration file is read if
   it exists.

#. if there is no ``local_...``, the ``getModuleConfig`` function reads the
   settings from the default configuration file with the default settings.


In the configuration file itself, there are multiple options to vary a parameter:

* replace the default parameter value with desired value
* provide a number of parameter values separated by ``|`` (e.g. ``relTh=1.|2.|3.``)
* provide a number of parameter values using ``start:stop:numberOfSteps`` (e.g. ``relTh=1.:3.:3``)) - a
  single value can be added by appending ``&4.0`` for example
