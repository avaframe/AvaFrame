"""
Run the rock avalanche setup of com1DFA
"""

import pathlib
import time
import argparse

# Local imports
# import config and init tools
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.in3Utils.initializeProject as initProj
from avaframe.in3Utils import fileHandlerUtils as fU

# import computation modules
from avaframe.com6RockAvalanche import com6RockAvalanche


def runRockAvalanche(avalancheDir=""):
    """Run com1DFA with rock avalanche parameters with only an avalanche directory as input

    Parameters
    ----------
    avalancheDir: str
        path to avalanche directory (setup e.g. with init scripts)

    Returns
    -------
    peakFilesDF: pandas dataframe
        with info about com1DFA peak file locations
    """
    # Time the whole routine
    startTime = time.time()

    # log file name; leave empty to use default runLog.log
    logName = "runCom6RockAvalanche"

    # Load avalanche directory from general configuration file
    # More information about the configuration can be found here
    # on the Configuration page in the documentation
    cfgMain = cfgUtils.getGeneralConfig()
    if avalancheDir != "":
        cfgMain["MAIN"]["avalancheDir"] = avalancheDir
    else:
        avalancheDir = cfgMain["MAIN"]["avalancheDir"]

    # Start logging
    log = logUtils.initiateLogger(avalancheDir, logName)
    log.info("MAIN SCRIPT")
    log.info("Current avalanche: %s", avalancheDir)

    # ----------------
    # Clean input directory(ies) of old work files
    initProj.cleanSingleAvaDir(avalancheDir, deleteOutput=False)

    # load rock avalanche config
    rockAvalancheCfg = cfgUtils.getModuleConfig(com6RockAvalanche)

    # perform com1DFA simulation with rock avalanche settings
    _, plotDict, reportDictList, _ = com6RockAvalanche.runRockAvalanche(cfgMain, rockAvalancheCfg)

    # Get peakfiles to return to QGIS
    avaDir = pathlib.Path(avalancheDir)
    inputDir = avaDir / "Outputs" / "com1DFA" / "peakFiles"
    peakFilesDF = fU.makeSimDF(inputDir, avaDir=avaDir)

    # Print time needed
    endTime = time.time()
    log.info("Took %6.1f seconds to calculate." % (endTime - startTime))

    return peakFilesDF


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run rock avalanche workflow")
    parser.add_argument(
        "avadir", metavar="avadir", type=str, nargs="?", default="", help="the avalanche directory"
    )

    args = parser.parse_args()
    runRockAvalanche(str(args.avadir))
