"""
    Setup required directory structures by computational modules
"""

# Load modules
import pathlib
import logging

# Local imports
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.in3Utils.initializeProject as initProj

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def initialiseRunDirs(avaDir, modName, cleanRemeshedRasters):
    """ Initialise Simulation run with input data

        Parameters
        ----------
        avaDir : str
            path to avalanche directory
        modName : str
            name of module
        cleanRemeshedRasters: bool
            if True directory Inputs/DEMremeshed shall be cleaned

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
    configDoneDir = pathlib.Path(avaDir, 'Outputs', modName, 'configurationFiles', 'configurationFilesDone')
    fU.makeADir(configDoneDir)
    # If Work directory already exists - error
    if workDir.is_dir():
        message = 'Work directory %s already exists - delete first!' % (workDir)
        log.error(message)
        raise AssertionError(message)
    else:
        workDir.mkdir(parents=True, exist_ok=False)
    log.debug('Directory: %s created' % workDir)

    if cleanRemeshedRasters is True:
        initProj.cleanRemeshedDir(avaDir)

    # first clean configurationFilesLatest dir and create new one
    initProj.cleanLatestConfigurationsDirAndCreate(avaDir, modName)

    return workDir, outputDir
