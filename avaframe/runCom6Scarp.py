
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
from com6RockAvalanche import scarp


def runScarpAnalysisWorkflow(avadir=''):
    """Run the scarp analysis workflow.

    Parameters
    ----------
    avadir : str
        Path to the avalanche directory containing input and output folders.

    Returns
    -------
    None
    """

    logName = 'runScarpAnalysis'

    # Load general configuration file
    cfgMain = cfgUtils.getGeneralConfig()
    if avadir:
        cfgMain['MAIN']['avalancheDir'] = avadir
    else:
        avadir = cfgMain['MAIN']['avalancheDir']

    avadir = pathlib.Path(avadir)

    # Start logging
    log = logUtils.initiateLogger(avadir.parent, logName)
    log.info('MAIN SCRIPT')
    log.info('Using configuration file: %s', avadir)

    # Clean input directory(ies) of old work files
    initProj.cleanSingleAvaDir(avadir.parent, deleteOutput=False)

    # Load module-specific configuration for Variable Voellmy
    scarpCfg = cfgUtils.getModuleConfig(scarp)

    # Run the scarp analysis
    runScarpAnalysis(avadir, scarpCfg)

    log.info('Scarp analysis completed successfully.')

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run scarp analysis workflow')
    parser.add_argument('avadir', metavar='c', type=str, nargs='?', default='',
                        help='Path to the avalanche directory')

    args = parser.parse_args()
    runScarpAnalysisWorkflow(str(args.avadir))


