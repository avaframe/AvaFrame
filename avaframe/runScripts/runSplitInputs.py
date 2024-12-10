"""runScript for splitting avalanche input data into multiple folders based on individual release areas"""

import time
import pathlib
import argparse

from avaframe.in3Utils import cfgUtils, logUtils
from avaframe.in4Region import splitInputs as sI

def runSplitInputs(avalancheDir=''):
    """Main function to split input data for each release area feature with only an avalancheDir as input."""

    # Start logging
    log = logUtils.initiateLogger(avalancheDir, logName='runSplitInputs')
    log.info('MAIN SCRIPT')

    # Time the whole routine
    startTime = time.time()

    # Load the avalanche directory from general configuration file
    cfgMain = cfgUtils.getGeneralConfig()
    if avalancheDir != '':
        cfgMain['MAIN']['avalancheDir'] = avalancheDir
    else:
        avalancheDir = cfgMain['MAIN']['avalancheDir']

    # Load input parameters from splitInputs configuration file
    cfg = cfgUtils.getModuleConfig(sI)

    # Paths to input files
    inputRELDir = pathlib.Path(avalancheDir) / 'Inputs' / 'REL'
    inputShp = next(inputRELDir.glob("*.shp"), None)
    if not inputShp:
        log.error(f"No shapefile found in {inputRELDir}.")
        return
    inputDEM = next((pathlib.Path(avalancheDir) / 'Inputs').glob("*.asc"), None)
    if not inputDEM:
        log.error(f"No DEM file found in {pathlib.Path(avalancheDir) / 'Inputs'}.")
        return

    # Create the split inputs dir
    outputDir = pathlib.Path(avalancheDir) / 'Regional'
    outputDir.mkdir(parents=True, exist_ok=True)

    # Step 1: Create the central list
    folderList = sI.createFolderList(inputShp)

    # Step 2: Set up ava directories
    log.info("Running folder initialization for each entry...")
    sI.createFoldersForReleaseAreas(folderList, outputDir)
    log.info("Finished folder initialization")

    # Step 3: Split and move release areas to each directory
    log.info("Splitting and moving release areas...")
    sI.splitAndMoveReleaseAreas(folderList, inputShp, outputDir)
    log.info("Finished splitting and moving release areas")

    # Step 4: Clip and move DEM
    if cfg['GENERAL'].getboolean('splitDEM'):
        log.info("Clipping and moving DEM...")
        sI.splitDEMByCenterpointAndMove(folderList, inputDEM, outputDir, cfg)
        log.info("Finished clipping and moving of DEM")

    # Print time needed
    endTime = time.time()
    log.info(f"Completed splitting input data into '{len(folderList)}' individual folders after "
             f"{endTime - startTime:.1f} seconds.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run split inputs for avalanche directories")
    parser.add_argument('avalancheDir', type=str, nargs='?', default='',
                        help="Directory containing the main avalanche data")

    args = parser.parse_args()
    runSplitInputs(args.avalancheDir)

