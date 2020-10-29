
Let's have a look at the simple example in runTmp1Ex.py and tmp1Ex.py on how
this is used within Avaframe.

Basically, in your main script call::

  from avaframe.in3Utils import logUtils
  # log file name; leave empty to use default runLog.log
  logName = 'runTmp1Ex'
  # specify the working directory
  avalancheDir = './'
  # ---------------------------------------------
  # Start logging
  log = logUtils.initiateLogger(avalancheDir, logName)

This will configure the logging (it sets the console output as well as the log file).
In your modules/subscripts add::

  import logging
  log = logging.getLogger(__name__)

So you can use::

  log.debug('Should be here')
  log.info('DEM : %s',variable)

To get output that looks like this in your console::

  tmp1Ex:DEBUG - Should be here
  tmp1Ex:INFO - DEM : /path/to/DEM

And something similar in the ``.log`` file which is saved in ``./runTmp1Ex.log`` in this example.
The logging configuration is set in ``Avaframe/avaframe/in3Utils/logging.conf``.
You can modify this ``logging.conf`` file to modify the levels or format of the messages to display
(`python doc will help you <https://docs.python.org/3/library/logging.config.html>`_).
