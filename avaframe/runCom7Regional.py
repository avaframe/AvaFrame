"""
Run script for running simulations in parallel based on input from a regional folder containing
multiple avaFolders
"""

import time
import argparse

from avaframe.com7Regional import com7Regional as com7
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


def runCom7Regional(avalancheDir=""):
    """Run regional avalanche simulations in parallel.

    Parameters
    ----------
    avalancheDir : str, optional
        Path to the main avalanche directory. If not provided, uses the path from general configuration.

    Returns
    -------
    allPeakFilesDir : pathlib.Path or None
        Path to the directory containing the collected peak files from all sub-avalanche directories, if enabled
    allTimeStepsDir : pathlib.Path or None
        Path to the directory containing consolidated time step files, if enabled
    mergedRastersDir : pathlib.Path or None
        Path to the directory containing merged output rasters, if enabled
    """

    # Time the whole routine
    startTime = time.time()

    # Load the avalanche directory from command input or the general configuration file
    cfgMain = cfgUtils.getGeneralConfig()
    if avalancheDir != "":
        cfgMain["MAIN"]["avalancheDir"] = avalancheDir
    else:
        avalancheDir = cfgMain["MAIN"]["avalancheDir"]

    # Start logging
    log = logUtils.initiateLogger(str(avalancheDir), logName="runCom7Regional")
    log.info("MAIN SCRIPT")

    # Load module configuration
    cfg = cfgUtils.getModuleConfig(com7, fileOverride="", toPrint=False, onlyDefault=False)

    # Call main function
    allPeakFilesDir, mergedRastersDir = com7.com7RegionalMain(cfgMain, cfg)

    # Print time needed
    endTime = time.time()
    log.info("Regional process finished after %.1f seconds." % (endTime - startTime))

    return allPeakFilesDir, mergedRastersDir


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Run regional workflow")
    parser.add_argument(
        "avadir", metavar="avalancheDir", type=str, nargs="?", default="", help="the avalanche directory"
    )

    args = parser.parse_args()
    runCom7Regional(str(args.avadir))
