"""runScript for splitting avalanche input data into multiple folders based on individual release areas"""

import time
import pathlib
import argparse
import logging
import shapefile  # pyshp
from avaframe.in3Utils import cfgUtils, logUtils
import avaframe.in3Utils.initializeProject as initProj

# create local logger
log = logging.getLogger(__name__)

def splitReleaseAreas(input_shapefile, output_dir):
    """Split the release area shapefile into individual features, each saved in a separate folder."""
    # Open the shapefile
    with shapefile.Reader(input_shapefile) as src:
        fields = src.fields[1:]  # Skip the deletion flag field
        field_names = [field[0] for field in fields]

        log.info(f"Loaded release areas shapefile with {len(src)} features.")

        # Iterate over each feature (record + shape)
        unnamed_count = 1  # Counter for unnamed folders
        for idx, record in enumerate(src.iterShapeRecords()):
            properties = dict(zip(field_names, record.record))
            folder_name = properties.get('name', '').strip() or f"Scenario {unnamed_count}"
            if not properties.get('name', '').strip():
                unnamed_count += 1

            ava_folder_path = output_dir / folder_name

            # Initialize the folder structure if it doesn't exist
            initProj.initializeFolderStruct(str(ava_folder_path), removeExisting=False)

            # Save each feature to its REL directory
            rel_dir = ava_folder_path / 'Inputs' / 'REL'
            rel_dir.mkdir(parents=True, exist_ok=True)
            feature_output_path = rel_dir / f"{folder_name}.shp"

            # Create a new shapefile for the feature
            with shapefile.Writer(str(feature_output_path)) as dst:
                dst.fields = src.fields[1:]  # Copy fields
                dst.record(*record.record)   # Copy data
                dst.shape(record.shape)      # Copy geometry

            log.info(f"Created folder '{folder_name}' with release area shapefile '{feature_output_path}'.")


def runSplitInputs(avalancheDir=''):
    """Main function to split input data for each release area feature."""
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
    input_shapefile = pathlib.Path(avalancheDir) / 'Inputs' / 'REL' / 'refactored_segmentedPRAs_050_cutOr1.shp'
    output_dir = pathlib.Path(avalancheDir) / 'SplitInputs'

    # Ensure Regional directory exists without creating unnecessary subfolders
    output_dir.mkdir(parents=True, exist_ok=True)

    # Split release area shapefile
    splitReleaseAreas(input_shapefile, output_dir)

    # Print time needed
    endTime = time.time()
    log.info('Took %6.1f seconds to calculate.' % (endTime - startTime))

    log.info("Completed splitting release areas into individual folders.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run split inputs for avalanche directories")
    parser.add_argument('avalancheDir', type=str, help="Directory containing the main avalanche data")

    args = parser.parse_args()
    runSplitInputs(args.avalancheDir)

