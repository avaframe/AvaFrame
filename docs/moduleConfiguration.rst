
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

The ``getModuleConfig`` function reads the configuration file (``tmpEx.ini``
in our example) and writes the setting just read to the log file.
It is possible to modify those parameters by creating a copy of the configuration
file (``local_`` followed by the name of the original configuration file). In
this case, the ``getModuleConfig`` function overwrites the original settings.
