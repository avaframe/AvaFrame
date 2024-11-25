"""runScript for splitting avalanche input data into multiple folders based on individual release areas"""

import time
import pathlib
import argparse

from avaframe.in3Utils import cfgUtils, logUtils
from avaframe.in4Region import splitInputs as si

def runSplitInputs(avalancheDir=''):
    """Main function to split input data for each release area feature with only an avalancheDir as input."""

    # Start logging
    log = logUtils.initiateLogger(avalancheDir, logName='runSplitInputs')
    log.info('MAIN SCRIPT')

    # Time the whole routine
    startTime = time.time()

    # Load the (REGIONAL) avalanche directory from general configuration file
    cfgMain = cfgUtils.getGeneralConfig()
    if avalancheDir != '':
        cfgMain['MAIN']['avalancheDir'] = avalancheDir
    else:
        avalancheDir = cfgMain['MAIN']['avalancheDir']

    # Paths to the input files
    input_rel_dir = pathlib.Path(avalancheDir) / 'Inputs' / 'REL'
    input_shapefile = next(input_rel_dir.glob("*.shp"), None)
    if not input_shapefile:
        log.error(f"No shapefile found in {input_rel_dir}. Exiting.")
        return
    input_dem = next((pathlib.Path(avalancheDir) / 'Inputs').glob("*.asc"), None)
    if not input_dem:
        log.error(f"No DEM file found in {pathlib.Path(avalancheDir) / 'Inputs'}. Skipping DEM processing.")
        return

    # Create the split inputs dir
    output_dir = pathlib.Path(avalancheDir) / 'SplitInputs'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Create the central list
    folder_list = si.createFolderList(input_shapefile)

    # Step 2: Set up ava directories
    log.info("Running folder initialization for each entry...")
    si.createFoldersForReleaseAreas(folder_list, output_dir)
    log.info("Finished folder initialization")

    # Step 3: Split and move release areas to each directory
    log.info("Splitting and moving release areas...")
    si.splitAndMoveReleaseAreas(folder_list, input_shapefile, output_dir)
    log.info("Finished splitting and moving release areas")

    # Step 4: Clip and move DEM
    log.info("Clipping and moving DEM...")
    si.splitDEMByCenterpointAndMove(folder_list, input_dem, output_dir)
    log.info("Finished clipping and moving of DEM")

    # Print time needed
    endTime = time.time()
    log.info(f"Completed splitting input data into '{len(folder_list)}' individual folders after "
             f"{endTime - startTime:.1f} seconds.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run split inputs for avalanche directories")
    parser.add_argument('avalancheDir', type=str, nargs='?', default='',
                        help="Directory containing the main avalanche data")

    args = parser.parse_args()
    runSplitInputs(args.avalancheDir)

