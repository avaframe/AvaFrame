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

The ``getModuleConfig`` function reads the settings from the configuration file (``tmpEx.ini``
in our example) and writes these settings to the log file. The default settings can be found in the
configuration file provided within each module.
It is possible to modify these settings by creating a copy of the configuration
file (``local_`` followed by the name of the original configuration file) and set the
desired values of the individual parameters. If there is a configuration file with the
prefix ``local_``, the ``getModuleConfig`` function reads the settings from this file.
Another option is to read the settings from a different configuration file, therefore the path to this file
has to be passed to the ``getModuleConfig`` function directly. 
