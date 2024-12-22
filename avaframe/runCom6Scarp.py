
"""
    Run script for scarp analysis
    In this runscript, the path to the scarp configuartion File has to be defined in the main config file
"""

import pathlib
import argparse

# Local imports
from com6RockAvalanche.scarp import runScarpAnalysis
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from avaframe.in3Utils import initializeProject as initProj


def runScarpAnalysisWorkflow(configFile=''):
    """Run the scarp analysis workflow.

    Parameters
    ----------
    configFile : str
        Path to the configuration file in .ini format.

    Returns
    -------
    None
    """

    logName = 'runScarpAnalysis'

    # Load general configuration file
    cfgMain = cfgUtils.getGeneralConfig()

    # Load configuration file path from general config if not provided
    if configFile != '':
        cfgMain['MAIN']['configFile'] = configFile
    else:
        configFile = cfgMain['MAIN']['configFile']

    configFile = pathlib.Path(configFile)

    # Start logging
    log = logUtils.initiateLogger(configFile.parent, logName)
    log.info('MAIN SCRIPT')
    log.info('Using configuration file: %s', configFile)

    # Clean input directory(ies) of old work files
    initProj.cleanSingleAvaDir(configFile.parent, deleteOutput=False)

    # Run the scarp analysis
    runScarpAnalysis(str(configFile))

    log.info('Scarp analysis completed successfully.')

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run scarp analysis workflow')
    parser.add_argument('configFile', metavar='c', type=str, nargs='?', default='',
                        help='Path to the scarp configuration file')

    args = parser.parse_args()
    runScarpAnalysisWorkflow(str(args.configFile))

