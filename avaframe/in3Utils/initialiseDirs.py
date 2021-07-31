"""
    Setup required directory structures by computational modules
"""

# Load modules
import pathlib
import logging

# Local imports
from avaframe.in3Utils import fileHandlerUtils as fU

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def initialiseRunDirs(avaDir, modName):
    """ Initialise Simulation run with input data

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
    outputDir = pathlib.Path(avaDir, 'Outputs', modName)
    fU.makeADir(outputDir)
    workDir = pathlib.Path(avaDir, 'Work', modName)
    # If Work directory already exists - error
    if workDir.is_dir():
        log.error('Work directory %s already exists - delete first!' % (workDir))
    else:
        workDir.mkdir(parents=True, exist_ok=False)
    log.debug('Directory: %s created' % workDir)

    return workDir, outputDir
