"""
    Setup required directory structures by computational modules

    This file is part of Avaframe.
"""

# Load modules
import os
import glob
import subprocess
import shutil
import logging

# Local imports
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.in2Trans import ascUtils as aU

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def initialiseRunDirs(avaDir, modName):
    """ Initialise Simulation Run with input data

        Parameters
        ----------
        avaDir : str
            path to avalanche directoy
        modName : str
            name of module

        Returns
        -------
        workDir : str
            path to Work directory
        outputDir : str
            path to Outputs directory
    """

    # Set directories outputs and current work
    outputDir = os.path.join(avaDir, 'Outputs', modName)
    fU.makeADir(outputDir)
    workDir = os.path.join(avaDir, 'Work', modName)
    # If Work directory already exists - error
    if os.path.isdir(workDir):
        log.error('Work directory %s already exists - delete first!' % (workDir))
    else:
        os.makedirs(workDir)
    log.debug('Directory: %s created' % workDir)

    return workDir, outputDir
