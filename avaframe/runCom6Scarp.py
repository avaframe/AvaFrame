"""
    Run script for scarp analysis
    In this runscript, the path to the scarp configuartion File has to be defined in the main config file
"""

import pathlib
import argparse

# Local imports
from avaframe.com6RockAvalanche import scarp
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from avaframe.in3Utils import initializeProject as initProj


def runScarpAnalysisWorkflow(inputDir=""):
    """Run the scarp analysis workflow.

    Parameters
    ----------
    inputDir: str
        Path to the input directory containing DEM, coordinates and perimeter files.

    Returns
    -------
    None
    """

    logName = 'runScarpAnalysis'

    # Load general configuration file
    cfgMain = cfgUtils.getGeneralConfig()

    # Load configuration file path from general config if not provided
    if inputDir != "":
        cfgMain["MAIN"]["avalancheDir"] = inputDir
    else:
        inputDir = cfgMain["MAIN"]["avalancheDir"]

    # Start logging
    log = logUtils.initiateLogger(inputDir, logName)
    log.info("MAIN SCRIPT")
    log.info("Using inputDir: %s", inputDir)

    # ----------------
    # Clean input directory(ies) of old work files
    initProj.cleanSingleAvaDir(inputDir, deleteOutput=False)

    # load scarp config
    scarpCfg = cfgUtils.getModuleConfig(scarp)

    # Run the scarp analysis
    scarp.scarpAnalysisMain(scarpCfg, str(inputDir))

    log.info('Scarp analysis completed successfully.')

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run scarp analysis workflow')
    parser.add_argument(
        "inputDir", metavar="inputDir", type=str, nargs="?", default="", help="the input directory"
    )

    args = parser.parse_args()
    runScarpAnalysisWorkflow(str(args.inputDir))
